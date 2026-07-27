!
! FISOC ocean-model wrapper for LADDIE (the One-Layer Antarctic model for Dynamical
! Downscaling of Ice-ocean Exchanges, the 2D ocean model in the UPSY-models repo).
!
! The wrapper's job is to drive LADDIE through the FISOC/ESMF component interface and to
! translate between LADDIE's internal (UPSY mesh) data structures and ESMF objects.
!
! LADDIE is built as a static library (libLADDIE_library.a, see UPSY-models
! src/LADDIE/CMakeLists.txt, BUILD_LADDIE_LIBRARY=ON) which this wrapper links against,
! together with libUPSY_static_library.a and PETSc/NetCDF/HDF5/MPI libraries.
!
! ==========================================================================================

MODULE FISOC_OM_Wrapper

  USE ESMF

  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  USE FISOC_coupler_MOD,                 ONLY: FISOC_coupler_regridMaskedField

  ! LADDIE / UPSY library modules (mirrors the USE list of LADDIE_program.f90)
  USE petscksp
  USE precisions,                        ONLY: dp
  USE basic_program_info,                ONLY: program_name
  USE mpi_basic,                         ONLY: par, initialise_parallelisation
  USE parameters,                        ONLY: initialise_constants, sec_per_year
  USE call_stack_and_comp_time_tracking, ONLY: initialise_control_and_resource_tracker, &
                                               reset_resource_tracker
  USE model_configuration,               ONLY: C, initialise_model_configuration
  USE mesh_types,                        ONLY: type_mesh
  USE laddie_model_types,                ONLY: type_laddie_model
  USE laddie_forcing_types,              ONLY: type_laddie_forcing
  USE laddie_forcing_main,               ONLY: initialise_forcing
  USE laddie_hydrology,                  ONLY: initialise_transects_SGD
  USE LADDIE_main_model,                 ONLY: initialise_laddie_model, run_laddie_model
  USE netcdf_resource_tracking,          ONLY: create_resource_tracking_file, &
                                               write_to_resource_tracking_file
  USE ice_geometry_model_basic,          ONLY: type_ice_geometry_model

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_OM_Wrapper_Init_Phase1,  FISOC_OM_Wrapper_Init_Phase2,  &
       FISOC_OM_Wrapper_Run, FISOC_OM_Wrapper_Finalize

  ! Persistent LADDIE state, held across the FISOC phases (as the FVCOM/ROMS wrappers hold
  ! their model state in module variables). LADDIE's mesh is built in Init_Phase1 by
  ! initialise_forcing; the model state is allocated in Init_Phase2.
  TYPE(type_mesh),           SAVE, TARGET :: mesh
  TYPE(type_laddie_model),   SAVE :: laddie
  TYPE(type_laddie_forcing), SAVE :: forcing

  ! Local ESMF node position -> LADDIE global vertex index, i.e. the same array as
  ! LADDIE2ESMF_mesh's local "nodeIds", kept around (rather than discarded once the mesh
  ! is built) so later field exchange can index between an ESMF field's per-PET local
  ! array and LADDIE's global-vertex-indexed arrays (laddie%melt, forcing%Hib, ...).
  INTEGER,ALLOCATABLE,SAVE :: OM_localNode2globalVi(:)

  ! is_standalone purely controls LADDIE's own output; no physics change
  LOGICAL,PARAMETER :: is_standalone = .TRUE.
  
CONTAINS

  !--------------------------------------------------------------------------------------
  ! Init phase 1: initialise LADDIE far enough to know its mesh, and hand that mesh back
  ! to FISOC as an ESMF_Mesh. Mirrors the first part of LADDIE_program (parallelisation,
  ! PETSc, constants, resource tracker, configuration, forcing -> builds the mesh).
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase1(FISOC_config,vm,OM_ExpFB,OM_mesh,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)              :: vm          ! ESMF virtual machine (parallel context)
    TYPE(ESMF_mesh),INTENT(OUT)           :: OM_mesh
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ExpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: localPet, perr
    CHARACTER(len=ESMF_MAXSTR)            :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: OM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: OM_configFile, OM_stdoutFile
    LOGICAL                               :: verbose_coupling

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Per-PET LADDIE stdout file (matches the FVCOM wrapper convention). This wrapper's
    ! own minimal status lines are written to unit 0 (not the usual OM_outputUnit) so that
    ! they land in the same file, and interleave correctly, with LADDIE/UPSY's own internal
    ! messages -- those use Fortran's preconnected unit 0 throughout UPSY-models (e.g.
    ! write(0,...)), and reopening it here redirects all of them transparently, with no
    ! changes needed to UPSY-models itself. (Two separate units both targeting the same
    ! file, e.g. unit 0 here and the usual OM_outputUnit, is unsafe -- each keeps its own
    ! buffer/file position, and one can silently clobber the other's writes.
    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_stdoutFile, label='OM_stdoutFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    WRITE (OM_stdoutFile, "(a,I0)") TRIM(OM_stdoutFile), localPet
    OPEN(unit=0, file=OM_stdoutFile, STATUS='REPLACE', ERR=101)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Path to the LADDIE .cfg config file
    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_configFile, label='OM_configFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "OM wrapper (LADDIE).  Init phase 1 method."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    ! --- Drive the LADDIE initialisation sequence (cf. LADDIE_program.f90) ---------------
    program_name = 'LADDIE'

    CALL initialise_parallelisation

    ! PETSc. NOTE: guard with PetscInitialized if another component also initialises PETSc.
    CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)

    CALL initialise_constants
    CALL initialise_control_and_resource_tracker
    CALL initialise_model_configuration( OM_configFile)

    ! Standalone LADDIE_program.f90 writes one resource-tracking record for the whole
    ! run; here, one record gets appended (and the tracker reset) after every
    ! run_laddie_model call instead (see Init_Phase2 and Run below), giving a per-step
    ! timing breakdown across the coupled run rather than a single end-of-run summary.
    CALL create_resource_tracking_file( C%output_dir)

    ! initialise_forcing builds LADDIE's (UPSY) mesh from the reference-geometry file and
    ! sets up the ambient ocean forcing. After this, `mesh` is populated. 'ANT' matches the
    ! region_name convention used everywhere else in this single-region (Antarctica) LADDIE
    ! build (e.g. C%lambda_M_ANT, laddie_forcing_main.f90's own 'ANT' calls).
    CALL initialise_forcing( mesh, forcing, 'ANT')

    ! Translate the LADDIE mesh into an ESMF_Mesh to return to FISOC.
    CALL LADDIE2ESMF_mesh(FISOC_config, mesh, OM_mesh, vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    label = 'FISOC_OM_ReqVars:'
    CALL FISOC_getListFromConfig(FISOC_config, label, OM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_populateFieldBundle(OM_ReqVarList,OM_ExpFB,OM_mesh, &
         init_value=FISOC_missingData, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    RETURN

101 msg = "OM (LADDIE) failed to open stdoutFile "//OM_stdoutFile
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase1


  !--------------------------------------------------------------------------------------
  ! Init phase 2: ISM fields are now available. Stage 1 does no exchange, so here we just
  ! finish LADDIE's initialisation (subglacial-discharge transects + allocate model state),
  ! mirroring the rest of the LADDIE_program init sequence.
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase2(FISOC_config,vm,OM_ImpFB,OM_ExpFB,ISM_ExpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)                :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)           :: OM_ImpFB, OM_ExpFB
    TYPE(ESMF_fieldBundle),INTENT(IN),OPTIONAL     :: ISM_ExpFB
    TYPE(ESMF_VM),INTENT(IN)                       :: vm
    INTEGER,INTENT(OUT),OPTIONAL                   :: rc

    INTEGER   :: localPet
    LOGICAL   :: verbose_coupling
    REAL(dp)  :: LADDIE_time

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "OM wrapper (LADDIE).  Init phase 2 method."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    CALL sendFieldDataToOM(OM_ImpFB, FISOC_config, vm, ISM_ExpFB=ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL initialise_transects_SGD( mesh, forcing)
    CALL initialise_laddie_model( mesh, laddie, forcing, is_standalone)

    ! LADDIE's ocean state was just reset to its configured initial conditions by
    ! initialise_laddie_model, so this is the point to spin it up to quasi-equilibrium
    ! (is_initial=.TRUE. -> the longer C%time_duration_laddie_init cycle), mirroring how
    ! UFEMISM itself does this once, immediately after initialise_laddie_model, before its
    ! own main time-stepping loop begins (see BMB_main.f90/UFEMISM_main_model.f90). Doing
    ! this here means FISOC_OM_Wrapper_Run never needs to reason about is_initial at all --
    ! every Run call is a regular (already-quasi-steady) cycle.
    LADDIE_time = LADDIE_currentTime(rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL run_laddie_model( mesh, laddie, forcing, LADDIE_time, .TRUE., is_standalone)

    ! One resource-tracking record for the spin-up cycle, then reset so the next
    ! record (from the first regular Run call) reflects that step alone.
    CALL write_to_resource_tracking_file( LADDIE_time)
    CALL reset_resource_tracker

    ! Populate the export bundle with a real melt value from the just-completed spin-up,
    ! rather than leaving it at FISOC_missingData for the first coupling step.
    CALL getFieldDataFromOM(OM_ExpFB, FISOC_config, vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase2


  SUBROUTINE FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB,OM_ImpFB,ISM_ExpFB,rc_local)

    TYPE(ESMF_config),INTENT(INOUT)                :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL  :: OM_ExpFB, OM_ImpFB
    TYPE(ESMF_fieldBundle),INTENT(IN),OPTIONAL     :: ISM_ExpFB
    TYPE(ESMF_VM),INTENT(IN)                       :: vm
    INTEGER,INTENT(OUT),OPTIONAL                   :: rc_local

    INTEGER            :: localPet, rc
    LOGICAL            :: verbose_coupling
    REAL(dp)           :: LADDIE_time

    rc_local = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (PRESENT(OM_ImpFB)) THEN
       CALL sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,ISM_ExpFB=ISM_ExpFB,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    IF (localPet.EQ.0) THEN
       WRITE (0,*) 'FISOC is about to call the LADDIE run method.'
    END IF
    CALL ESMF_VMBarrier(vm, rc=rc)

    ! NOTE: run_laddie_model does NOT integrate to `time`; it runs a fixed duration
    ! (C%time_duration_laddie days) to quasi-steady state under the current geometry.
    ! is_initial is always .FALSE. here -- the one-off spin-up (is_initial=.TRUE., the
    ! longer C%time_duration_laddie_init cycle) already happened once in Init_Phase2, right
    ! after LADDIE's state was reset by initialise_laddie_model. `time` is derived from how
    ! far FISOC's own clock has advanced since its start, combined with LADDIE's own
    ! configured start time (C%start_time_of_run) -- see LADDIE_currentTime below. This
    ! matters beyond output timestamps: compute_melt_rate (laddie_physics.f90) branches its
    ! ice-heat-diffusion term on time==C%start_time_of_run, so an uncomputed/fixed `time`
    ! would silently pin every call to that first-step simplification.
    LADDIE_time = LADDIE_currentTime(rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL run_laddie_model( mesh, laddie, forcing, LADDIE_time, .FALSE., is_standalone)

    ! One resource-tracking record per coupling step, then reset for the next one.
    CALL write_to_resource_tracking_file( LADDIE_time)
    CALL reset_resource_tracker

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (localPet.EQ.0) THEN
       WRITE (0,*) 'FISOC has just called the LADDIE run method.'
    END IF

    IF (PRESENT(OM_ExpFB)) THEN
       CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    rc_local = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Finalize(FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)           :: vm
    INTEGER,INTENT(OUT),OPTIONAL       :: rc

    INTEGER                            :: localPet, perr
    LOGICAL                            :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "OM wrapper (LADDIE).  Finalise method."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    ! Finalise PETSc. Do NOT call MPI_FINALIZE here: ESMF/FISOC owns MPI finalisation.
    CALL PetscFinalize( perr)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_Wrapper_Finalize


  !--------------------------------------------------------------------------------------
  ! LADDIE's own "time" argument to run_laddie_model (years, see LADDIE_main_model.f90 and
  ! laddie_physics.f90:compute_melt_rate) is a plain elapsed-time scalar relative to
  ! LADDIE's own configured start (C%start_time_of_run), using LADDIE's own year length
  ! (UPSY's sec_per_year). This combines that configured start with how far FISOC's ESMF
  ! clock has actually advanced since FISOC's own start (module variables FISOC_time,
  ! FISOC_startTime, set up by FISOC_setClocks and kept current by FISOC_parent -- see the
  ! FVCOM wrapper for the same established pattern of reading FISOC_time directly).
  FUNCTION LADDIE_currentTime(rc) RESULT(LADDIE_time)

    REAL(dp)                       :: LADDIE_time
    INTEGER,OPTIONAL,INTENT(OUT)   :: rc

    TYPE(ESMF_TimeInterval)        :: elapsed
    REAL(ESMF_KIND_R8)             :: elapsed_sec
    INTEGER                        :: localrc

    IF (PRESENT(rc)) rc = ESMF_FAILURE

    elapsed = FISOC_time - FISOC_startTime
    CALL ESMF_TimeIntervalGet(elapsed, s_r8=elapsed_sec, rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    LADDIE_time = C%start_time_of_run + elapsed_sec / sec_per_year

    IF (PRESENT(rc)) rc = ESMF_SUCCESS

  END FUNCTION LADDIE_currentTime


  !--------------------------------------------------------------------------------
  ! Build an ESMF_Mesh from LADDIE's (UPSY) mesh.
  !--------------------------------------------------------------------------------
  SUBROUTINE LADDIE2ESMF_mesh(FISOC_config,mesh,OM_mesh,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)  :: FISOC_config
    TYPE(type_mesh),INTENT(IN)       :: mesh
    TYPE(ESMF_mesh),INTENT(OUT)      :: OM_mesh
    TYPE(ESMF_VM),INTENT(IN)         :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    INTEGER                          :: localPet, ti, vi, kk, lid, owner
    INTEGER                          :: nLocalElems, nLocalNodes
    INTEGER,ALLOCATABLE              :: gi2local(:)      ! global vi -> local node index (0 = absent)
    INTEGER,ALLOCATABLE              :: nodeIds(:), nodeOwners(:)
    INTEGER,ALLOCATABLE              :: elemIds(:), elemTypes(:), elemConn(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:)
    LOGICAL                          :: verbose_coupling

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)

    ! 1. Local elements = the triangles this PET owns.
    nLocalElems = mesh%ti2 - mesh%ti1 + 1

    ! 2. Collect the distinct vertices referenced by those triangles -> local node list,
    !    recording a global-to-local index map. ESMF_MeshCreate requires every node a
    !    PET declares to be used by one of that PET's own local elements, so the node
    !    set must be derived from the local triangles rather than from UPSY's
    !    independent vertex-ownership partition (vi1:vi2).
    ALLOCATE(gi2local(mesh%nV))
    gi2local = 0
    nLocalNodes = 0
    DO ti = mesh%ti1, mesh%ti2
       DO kk = 1,3
          vi = mesh%Tri(ti,kk)
          IF (gi2local(vi) == 0) THEN
             nLocalNodes = nLocalNodes + 1
             gi2local(vi) = nLocalNodes
          END IF
       END DO
    END DO

    ! 3. Fill the node arrays. The ESMF owner of a node must itself have a local
    !    element using that node (see point 2), so ownership can't be taken from
    !    UPSY's vertex partition (mesh%V_owning_process) either, since that is an
    !    independent partition from the triangle one. Instead, derive it from the
    !    triangles incident to the vertex: the owner is whichever PET owns the
    !    lowest-numbered such triangle, computed identically on every PET from
    !    UPSY's (globally available) mesh connectivity and triangle ownership data.
    ALLOCATE(nodeIds(nLocalNodes), nodeOwners(nLocalNodes), nodeCoords(2*nLocalNodes))
    DO vi = 1, mesh%nV
       lid = gi2local(vi)
       IF (lid > 0) THEN
          owner = mesh%Tri_owning_process(mesh%iTri(vi,1))
          DO kk = 2, mesh%niTri(vi)
             owner = MIN(owner, mesh%Tri_owning_process(mesh%iTri(vi,kk)))
          END DO
          nodeIds(lid)        = vi                          ! global node id
          nodeOwners(lid)     = owner                       ! owning PET (process rank, 0-based)
          nodeCoords(2*lid-1) = mesh%V(vi,1)                ! x [m]
          nodeCoords(2*lid  ) = mesh%V(vi,2)                ! y [m]
       END IF
    END DO

    ! 4. Fill the element arrays. UPSY stores triangle vertices counter-clockwise, which is
    !    the order ESMF expects. elementConn references LOCAL node positions.
    ALLOCATE(elemIds(nLocalElems), elemTypes(nLocalElems), elemConn(3*nLocalElems))
    elemTypes = ESMF_MESHELEMTYPE_TRI
    kk = 0
    DO ti = mesh%ti1, mesh%ti2
       kk = kk + 1
       elemIds(kk)       = ti
       elemConn(3*kk-2)  = gi2local(mesh%Tri(ti,1))
       elemConn(3*kk-1)  = gi2local(mesh%Tri(ti,2))
       elemConn(3*kk  )  = gi2local(mesh%Tri(ti,3))
    END DO

    ! 5. Create the ESMF mesh in one step.
    OM_mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, &
         coordSys=ESMF_COORDSYS_CART,                        &
         nodeIds=nodeIds, nodeCoords=nodeCoords,             &
         nodeOwners=nodeOwners,                              &
         elementIds=elemIds, elementTypes=elemTypes,         &
         elementConn=elemConn,                               &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       WRITE (0,*) 'LADDIE2ESMF_mesh: created ESMF mesh, ', &
            nLocalNodes,' local nodes, ',nLocalElems,' local elements.'
    END IF

    ! Build OM_localNode2globalVi: local ESMF *field* position -> global LADDIE vertex
    ! index. This is NOT the same as nodeIds -- a Field on this mesh only stores this
    ! PET's *owned* nodes (see the VERIFIED note above), so we need the owned subset of
    ! nodeIds, keeping their relative order, which is exactly the order ESMF uses for a
    ! PET's owned-node field storage.
    ALLOCATE(OM_localNode2globalVi(COUNT(nodeOwners == localPet)))
    lid = 0
    DO kk = 1, nLocalNodes
       IF (nodeOwners(kk) == localPet) THEN
          lid = lid + 1
          OM_localNode2globalVi(lid) = nodeIds(kk)
       END IF
    END DO

    DEALLOCATE(gi2local, nodeIds, nodeOwners, nodeCoords, elemIds, elemTypes, elemConn)

  END SUBROUTINE LADDIE2ESMF_mesh


  ! Rebuilds an ESMF Mesh identical in nodes/elements/ownership to LADDIE2ESMF_mesh's
  ! OM_mesh, but with a nodeMask attached (MASK_DRY at vertices that are grounded ice
  ! or icefree land, 0 elsewhere -- i.e. not a valid subglacial-discharge destination).
  ! Node masking can only be set at ESMF_MeshCreate time (ESMF_MeshSet only supports
  ! updating elementMask in place, not nodeMask -- confirmed against ESMF_refdoc-1.pdf
  ! 34.4.20/34.3.11), and forcing%mask_grounded_ice/mask_icefree_land evolve every
  ! coupling step, so unlike the long-lived OM_mesh, this one is rebuilt fresh every
  ! time it's needed (see FISOC_coupler_regridMaskedField, which likewise rebuilds its
  ! routehandle fresh every call for the same reason). The node construction here is
  ! deterministic and identical to LADDIE2ESMF_mesh's, so the existing module-level
  ! OM_localNode2globalVi mapping still applies to Fields built on this mesh -- no
  ! separate mapping is computed or returned.
  SUBROUTINE LADDIE2ESMF_mesh_masked(mesh,forcing,OM_mesh_masked,rc)

    TYPE(type_mesh),INTENT(IN)               :: mesh
    TYPE(type_laddie_forcing),INTENT(IN)     :: forcing
    TYPE(ESMF_mesh),INTENT(OUT)              :: OM_mesh_masked
    INTEGER,INTENT(OUT),OPTIONAL             :: rc

    INTEGER                          :: ti, vi, kk, lid, owner
    INTEGER                          :: nLocalElems, nLocalNodes
    INTEGER,ALLOCATABLE              :: gi2local(:)
    INTEGER,ALLOCATABLE              :: nodeIds(:), nodeOwners(:), nodeMask(:)
    INTEGER,ALLOCATABLE              :: elemIds(:), elemTypes(:), elemConn(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:)

    nLocalElems = mesh%ti2 - mesh%ti1 + 1

    ALLOCATE(gi2local(mesh%nV))
    gi2local = 0
    nLocalNodes = 0
    DO ti = mesh%ti1, mesh%ti2
       DO kk = 1,3
          vi = mesh%Tri(ti,kk)
          IF (gi2local(vi) == 0) THEN
             nLocalNodes = nLocalNodes + 1
             gi2local(vi) = nLocalNodes
          END IF
       END DO
    END DO

    ALLOCATE(nodeIds(nLocalNodes), nodeOwners(nLocalNodes), nodeCoords(2*nLocalNodes), &
         nodeMask(nLocalNodes))
    DO vi = 1, mesh%nV
       lid = gi2local(vi)
       IF (lid > 0) THEN
          owner = mesh%Tri_owning_process(mesh%iTri(vi,1))
          DO kk = 2, mesh%niTri(vi)
             owner = MIN(owner, mesh%Tri_owning_process(mesh%iTri(vi,kk)))
          END DO
          nodeIds(lid)        = vi
          nodeOwners(lid)     = owner
          nodeCoords(2*lid-1) = mesh%V(vi,1)
          nodeCoords(2*lid  ) = mesh%V(vi,2)
          IF (forcing%mask_grounded_ice(vi) .OR. forcing%mask_icefree_land(vi)) THEN
             nodeMask(lid) = MASK_DRY
          ELSE
             nodeMask(lid) = 0
          END IF
       END IF
    END DO

    ALLOCATE(elemIds(nLocalElems), elemTypes(nLocalElems), elemConn(3*nLocalElems))
    elemTypes = ESMF_MESHELEMTYPE_TRI
    kk = 0
    DO ti = mesh%ti1, mesh%ti2
       kk = kk + 1
       elemIds(kk)       = ti
       elemConn(3*kk-2)  = gi2local(mesh%Tri(ti,1))
       elemConn(3*kk-1)  = gi2local(mesh%Tri(ti,2))
       elemConn(3*kk  )  = gi2local(mesh%Tri(ti,3))
    END DO

    OM_mesh_masked = ESMF_MeshCreate(parametricDim=2, spatialDim=2, &
         coordSys=ESMF_COORDSYS_CART,                        &
         nodeIds=nodeIds, nodeCoords=nodeCoords,             &
         nodeOwners=nodeOwners, nodeMask=nodeMask,           &
         elementIds=elemIds, elementTypes=elemTypes,         &
         elementConn=elemConn,                               &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DEALLOCATE(gi2local, nodeIds, nodeOwners, nodeCoords, nodeMask, elemIds, elemTypes, elemConn)

  END SUBROUTINE LADDIE2ESMF_mesh_masked


  SUBROUTINE getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)     :: OM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)          :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)                 :: vm
    INTEGER,INTENT(OUT),OPTIONAL             :: rc

    INTEGER                               :: fieldCount, ii, nn
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)

    rc = ESMF_FAILURE

    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    fieldLoop: DO nn = 1,fieldCount

       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       SELECT CASE (TRIM(ADJUSTL(fieldName)))

       CASE ('OM_bmb')
          ! ptr holds this PET's *owned* nodes only (see LADDIE2ESMF_mesh), in the same
          ! order as OM_localNode2globalVi.
          DO ii = 1,SIZE(ptr)
             ptr(ii) = laddie%melt( OM_localNode2globalVi(ii))
          END DO

       CASE DEFAULT
          ! Same "trust the config" reasoning as sendFieldDataToOM's CASE DEFAULT: a field
          ! FISOC_OM_ReqVars: names that this wrapper doesn't know how to export is allowed
          ! to just stay at its FISOC_missingData init value rather than aborting the run.
          msg = "WARNING: getFieldDataFromOM: unhandled field: "//TRIM(ADJUSTL(fieldName))
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__, rc=rc)

       END SELECT

    END DO fieldLoop

    DEALLOCATE(fieldList)

    rc = ESMF_SUCCESS

  END SUBROUTINE getFieldDataFromOM


  !--------------------------------------------------------------------------------------
  ! ESMF -> LADDIE import. Only ISM_thick (ice thickness) is handled: LADDIE's forcing
  ! wants a full package of *absolute*, already-derived quantities (Hs, Hib, TAF, several
  ! boolean masks, b-grid draft slopes -- see laddie_forcing_types.f90). update_laddie_forcing
  ! (UFEMISM/basal_mass_balance/BMB_main.f90), the native UFEMISM+LADDIE coupling
  ! equivalent of this routine, gets all of these by copying UFEMISM's own already-computed
  ! ice_model fields; here we instead compute them ourselves, on LADDIE's own mesh, from
  ! Elmer's ice thickness via the same generic UPSY utilities UFEMISM itself uses
  ! (masks_mod::determine_masks; ice_geometry_basics; mesh_disc_apply_operators). Bedrock
  ! and sea level are NOT updated from the ISM -- LADDIE's own reference-geometry starting
  ! values (forcing%Hb; sea level assumed 0, matching typical modern-day setups) are
  ! trusted for the whole run.
  SUBROUTINE sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,ISM_ExpFB,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)          :: OM_ImpFB
    TYPE(ESMF_config),INTENT(INOUT)               :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(IN),OPTIONAL    :: ISM_ExpFB
    TYPE(ESMF_VM),INTENT(IN)                      :: vm
    INTEGER,INTENT(OUT),OPTIONAL                  :: rc

    INTEGER                               :: fieldCount, ii, nn
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName, label, listLabel
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)
    TYPE(type_ice_geometry_model)         :: geom
    TYPE(ESMF_Field)                      :: SGD_srcField, SGD_dstField
    TYPE(ESMF_mesh)                       :: SGD_dstMesh
    REAL(ESMF_KIND_R8),POINTER            :: SGD_ptr(:)
    LOGICAL                               :: gotIceThicknessUpdate

    rc = ESMF_FAILURE
    gotIceThicknessUpdate = .FALSE.

    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    fieldLoop: DO nn = 1,fieldCount

       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       SELECT CASE (TRIM(ADJUSTL(fieldName)))

       CASE ('ISM_thick')
          ! ptr holds this PET's *owned* nodes only (see LADDIE2ESMF_mesh), in the same
          ! order as OM_localNode2globalVi.
          DO ii = 1,SIZE(ptr)
             forcing%Hi( OM_localNode2globalVi(ii)) = ptr(ii)
          END DO
          gotIceThicknessUpdate = .TRUE.

       CASE DEFAULT
          ! Trust the config: a field FISOC hands us that we don't act on is allowed --
          ! e.g. ISM2OM_vars: naming something LADDIE has no use for. Only ISM_thick
          ! actually drives anything here.
          msg = "WARNING: sendFieldDataToOM: unhandled field: "//TRIM(ADJUSTL(fieldName))
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__, rc=rc)

       END SELECT

    END DO fieldLoop

    DEALLOCATE(fieldList)

    ! Trust the config here too: if ISM_thick isn't among whatever FISOC_ISM_ReqVars:/
    ! ISM2OM_vars: put in OM_ImpFB (including if that list is empty -- LADDIE then just
    ! runs standalone-like).
    IF (gotIceThicknessUpdate) THEN

       ! Recompute LADDIE's derived geometry (surface elevation, draft, thickness-above-
       ! floatation, masks, b-grid draft slopes) from the freshly-updated ice thickness,
       ! using the same type-bound ice_geometry_model machinery as native UFEMISM+LADDIE
       ! forcing updates (compare laddie_forcing_main.f90's initialise_forcing and
       ! BMB_main.f90's update_laddie_forcing), just computed here rather than copied from
       ! UFEMISM. geom is a scratch, subroutine-local instance: allocate/fill/consume/let
       ! it finalise (type_ice_geometry_model has a FINAL procedure) each call, rather than
       ! keeping one around persistently, since it's only needed transiently here.
       !
       ! IMPORTANT: geom's fields are allocated at [mesh%vi1:mesh%vi2]/[mesh%ti1:mesh%ti2]
       ! (this process's own local range only -- see allocate_ice_geometry_model), while
       ! forcing's fields are dist_shared at [pai_V%i1_nih:i2_nih]/[pai_Tri%i1_nih:i2_nih]
       ! (this whole shared-memory node's range, including halo -- see
       ! laddie_utilities.f90's allocate_laddie_forcing). These are genuinely different
       ! sizes. A bare whole-array "geom%Hi = forcing%Hi" would silently trigger Fortran's
       ! reallocate-on-assignment for the (allocatable) LHS, resizing just that one field
       ! and leaving geom internally inconsistent; a bare "forcing%Hs = geom%Hs" would hit
       ! a hard runtime bounds-check failure since forcing%Hs is a pointer (no
       ! reallocation). Always copy through the explicit vi1:vi2/ti1:ti2 section instead,
       ! which is valid and matching on both sides.
       CALL geom%allocate( 'ANT', mesh)
       geom%Hi(mesh%vi1:mesh%vi2) = forcing%Hi(mesh%vi1:mesh%vi2)
       geom%Hb(mesh%vi1:mesh%vi2) = forcing%Hb(mesh%vi1:mesh%vi2)
       geom%SL(mesh%vi1:mesh%vi2) = 0.0_dp   ! LADDIE's own sea level is not available, but LADDIE assumes it to be zero, so we do too.

       CALL geom%calc_surface_elevation()
       CALL geom%calc_ice_base_elevation()
       CALL geom%calc_thickness_above_floatation()
       CALL geom%determine_masks()
       CALL geom%calc_ice_base_slopes()

       forcing%Hs(mesh%vi1:mesh%vi2)                 = geom%Hs(mesh%vi1:mesh%vi2)
       forcing%Hib(mesh%vi1:mesh%vi2)                = geom%Hib(mesh%vi1:mesh%vi2)
       forcing%TAF(mesh%vi1:mesh%vi2)                = geom%TAF(mesh%vi1:mesh%vi2)
       forcing%mask(mesh%vi1:mesh%vi2)               = geom%mask(mesh%vi1:mesh%vi2)
       forcing%mask_icefree_land(mesh%vi1:mesh%vi2)  = geom%mask_icefree_land(mesh%vi1:mesh%vi2)
       forcing%mask_icefree_ocean(mesh%vi1:mesh%vi2) = geom%mask_icefree_ocean(mesh%vi1:mesh%vi2)
       forcing%mask_grounded_ice(mesh%vi1:mesh%vi2)  = geom%mask_grounded_ice(mesh%vi1:mesh%vi2)
       forcing%mask_floating_ice(mesh%vi1:mesh%vi2)  = geom%mask_floating_ice(mesh%vi1:mesh%vi2)
       forcing%mask_gl_fl(mesh%vi1:mesh%vi2)         = geom%mask_gl_fl(mesh%vi1:mesh%vi2)
       forcing%dHib_dx_b(mesh%ti1:mesh%ti2)          = geom%dHib_dx_b(mesh%ti1:mesh%ti2)
       forcing%dHib_dy_b(mesh%ti1:mesh%ti2)          = geom%dHib_dy_b(mesh%ti1:mesh%ti2)

       CALL geom%deallocate()

    END IF

    ! Subglacial discharge from the ISM (e.g. Elmer/GlaDS), if configured. Uses
    ! forcing%mask_grounded_ice/mask_icefree_land as currently held -- freshly
    ! updated above if ice thickness came in this step, otherwise whatever they
    ! were last set to (trust the config, same philosophy as ISM_thick above).
    IF (PRESENT(ISM_ExpFB)) THEN
       label = "ISM_SGD_flux"
       listLabel = "ISM2OM_vars"
       IF (FISOC_ConfigStringListContains(FISOC_config,label,listLabel,rc=rc)) THEN

          ! laddie%SGD is only allocated once initialise_laddie_model has run (see
          ! FISOC_OM_Wrapper_Init_Phase2, which calls sendFieldDataToOM -- this routine
          ! -- *before* initialise_laddie_model, since the latter needs this routine's
          ! own ice-thickness/mask refresh above to have already happened). So on that
          ! one very first call, laddie itself doesn't exist yet; skip the SGD update
          ! this one time rather than indexing an unallocated array -- real SGD values
          ! get applied starting from the next call (the first regular Run cycle),
          ! well after initialise_laddie_model has completed.
          IF (.NOT. ASSOCIATED(laddie%SGD)) THEN
             CALL ESMF_LogWrite("sendFieldDataToOM: laddie%SGD not yet allocated "// &
                  "(expected on the first Init Phase 2 call only) -- skipping SGD "// &
                  "update this call.", logmsgFlag=ESMF_LOGMSG_INFO, &
                  line=__LINE__, file=__FILE__, rc=rc)
             rc = ESMF_SUCCESS
             RETURN
          END IF

          CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_SGD_flux", field=SGD_srcField, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

          ! Rebuilt fresh every call: forcing%mask_grounded_ice/mask_icefree_land evolve
          ! every step, and node masking can't be updated on an existing Mesh (see
          ! LADDIE2ESMF_mesh_masked).
          CALL LADDIE2ESMF_mesh_masked(mesh, forcing, SGD_dstMesh, rc=rc)

          SGD_dstField = ESMF_FieldCreate(SGD_dstMesh, typekind=ESMF_TYPEKIND_R8, &
               meshloc=ESMF_MESHLOC_NODE, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

          ! zeroregion=TOTAL: nodes excluded by the mask (grounded/icefree land) are reset
          ! to zero, matching the physical expectation of no subglacial discharge there,
          ! rather than retaining a stale value from a previous step.
          CALL FISOC_coupler_regridMaskedField(SGD_srcField, SGD_dstField, MASK_DRY, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

          CALL ESMF_FieldGet(SGD_dstField, farrayPtr=SGD_ptr, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

          ! ptr holds this PET's *owned* nodes only, in the same order as
          ! OM_localNode2globalVi (see LADDIE2ESMF_mesh_masked - construction is
          ! identical to the long-lived OM_mesh, so the same mapping applies).
          ! Volume flux [m^3 s^-1] -> rate [m s^-1], matching how LADDIE's own
          ! compute_subglacial_discharge normalises by area (laddie_physics.f90).
          DO ii = 1,SIZE(SGD_ptr)
             laddie%SGD( OM_localNode2globalVi(ii)) = SGD_ptr(ii) / mesh%A( OM_localNode2globalVi(ii))
          END DO

          CALL ESMF_FieldDestroy(SGD_dstField, rc=rc)
          CALL ESMF_MeshDestroy(SGD_dstMesh, rc=rc)

       END IF
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE sendFieldDataToOM

END MODULE FISOC_OM_Wrapper
