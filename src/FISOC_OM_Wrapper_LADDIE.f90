!
! FISOC ocean-model wrapper for LADDIE (the One-Layer Antarctic model for Dynamical
! Downscaling of Ice-ocean Exchanges, the 2D ocean model in the UPSY-models repo).
!
! The wrapper's job is to drive LADDIE through the FISOC/ESMF component interface and to
! translate between LADDIE's internal (UPSY mesh) data structures and ESMF objects.
!
! LADDIE is built as a static library (libLADDIE_library.a, see UPSY-models
! src/LADDIE/CMakeLists.txt, BUILD_LADDIE_LIBRARY=ON) which this wrapper links against,
! together with libUPSY_static_library.a and the conda PETSc/NetCDF/HDF5/MPI libraries.
!
! ==========================================================================================
! STAGE 1 (this file): run LADDIE *from* FISOC with NO field exchange, and verify it
! reproduces a standalone LADDIE run. So the four entry points simply drive the same
! initialise/run sequence as the standalone LADDIE_program, and build an ESMF_Mesh from
! LADDIE's mesh. Data exchange (export melt; import ice draft + gradients) is a later stage;
! the export/import helper is stubbed below.
! ==========================================================================================
!
! INTEGRATION POINTS still required to build/run stage 1 (tracked in ../../TODO.md):
!   (1) UPSY side (fisoc branch): mpi_basic::initialise_parallelisation must NOT call
!       MPI_INIT when MPI is already initialised (ESMF initialises MPI). i.e. guard with
!       MPI_INITIALIZED. This is the minimal parallelisation change for stage 1; the full
!       MPI_COMM_WORLD -> par%comm sub-communicator refactor stays parked. For stage 1
!       LADDIE uses the global communicator (MPI_COMM_WORLD = all OM PETs).
!   (2) Build wiring: FISOC Makefile invoked with FISOC_OM=LADDIE and
!       FISOC_OM_GEOM=FISOC_OM_MESH, with FISOC_OM_LIBS/LIBPATH/INCLUDE pointing at the
!       LADDIE + UPSY libraries and their mod_files.
!   (3) PETSc: this wrapper calls PetscInitialize/Finalize (as standalone LADDIE does). If
!       another component (e.g. Elmer) also initialises PETSc, guard with PetscInitialized.
!

MODULE FISOC_OM_Wrapper

  USE ESMF

  USE FISOC_utils_MOD
  USE FISOC_types_MOD

  ! LADDIE / UPSY library modules (mirrors the USE list of LADDIE_program.f90)
  USE petscksp
  USE precisions,                        ONLY: dp
  USE basic_program_info,                ONLY: program_name
  USE mpi_basic,                         ONLY: par, initialise_parallelisation
  USE parameters,                        ONLY: initialise_constants
  USE call_stack_and_comp_time_tracking, ONLY: initialise_control_and_resource_tracker, &
                                               reset_resource_tracker
  USE model_configuration,               ONLY: C, initialise_model_configuration
  USE mesh_types,                        ONLY: type_mesh
  USE laddie_model_types,                ONLY: type_laddie_model
  USE laddie_forcing_types,              ONLY: type_laddie_forcing
  USE laddie_forcing_main,               ONLY: initialise_forcing
  USE laddie_hydrology,                  ONLY: initialise_transects_SGD
  USE LADDIE_main_model,                 ONLY: initialise_laddie_model, run_laddie_model

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_OM_Wrapper_Init_Phase1,  FISOC_OM_Wrapper_Init_Phase2,  &
       FISOC_OM_Wrapper_Run, FISOC_OM_Wrapper_Finalize

  ! Persistent LADDIE state, held across the FISOC phases (as the FVCOM/ROMS wrappers hold
  ! their model state in module variables). LADDIE's mesh is built in Init_Phase1 by
  ! initialise_forcing; the model state is allocated in Init_Phase2.
  TYPE(type_mesh),           SAVE :: mesh
  TYPE(type_laddie_model),   SAVE :: laddie
  TYPE(type_laddie_forcing), SAVE :: forcing

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

    ! Per-PET LADDIE stdout file (matches the FVCOM wrapper convention)
    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_stdoutFile, label='OM_stdoutFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    WRITE (OM_stdoutFile, "(a,I0)") TRIM(OM_stdoutFile), localPet
    OPEN(unit=OM_outputUnit, file=OM_stdoutFile, STATUS='REPLACE', ERR=101)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Path to the LADDIE .cfg config file
    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_configFile, label='OM_configFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"*******************************************************************************"
       PRINT*,"**********       OM wrapper (LADDIE).  Init phase 1 method.     **************"
       PRINT*,"*******************************************************************************"
       PRINT*,""
    END IF

    ! --- Drive the LADDIE initialisation sequence (cf. LADDIE_program.f90) ---------------
    program_name = 'LADDIE'

    ! Parallelisation. NOTE (integration point 1): UPSY's initialise_parallelisation must be
    ! modified to skip MPI_INIT when MPI is already initialised (ESMF owns MPI_INIT). For
    ! stage 1 LADDIE uses the global communicator (MPI_COMM_WORLD = all OM PETs); passing a
    ! FISOC sub-communicator is a later stage (tied to the MPI_COMM_LADDIE work).
    CALL initialise_parallelisation

    ! PETSc. NOTE (integration point 3): guard with PetscInitialized if another component
    ! also initialises PETSc.
    CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)

    CALL initialise_constants
    CALL initialise_control_and_resource_tracker
    CALL initialise_model_configuration( OM_configFile)

    ! initialise_forcing builds LADDIE's (UPSY) mesh from the reference-geometry file and
    ! sets up the ambient ocean forcing. After this, `mesh` is populated.
    CALL initialise_forcing( mesh, forcing)

    ! Translate the LADDIE mesh into an ESMF_Mesh to return to FISOC.
    CALL LADDIE2ESMF_mesh(FISOC_config, mesh, OM_mesh, vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create the (empty) export field bundle on the mesh. The required-variable list comes
    ! from the FISOC config. Stage 1 does no exchange, so the fields are left at missing
    ! data (no getFieldDataFromOM call yet).
    label = 'FISOC_OM_ReqVars:'
    CALL FISOC_getListFromConfig(FISOC_config, label, OM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_populateFieldBundle(OM_ReqVarList,OM_ExpFB,OM_mesh, &
         init_value=FISOC_missingData, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! STAGE >1: CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc=rc)  ! e.g. OM_bmb <- laddie%melt

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
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase2(FISOC_config,vm,OM_ImpFB,OM_ExpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ImpFB, OM_ExpFB
    TYPE(ESMF_VM),INTENT(IN)              :: vm
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER   :: localPet
    LOGICAL   :: verbose_coupling
    ! Stage 1: reproduce a standalone run, so is_standalone = .TRUE. (purely controls
    ! LADDIE's own output; no physics change). Flip to .FALSE. once FISOC drives the output.
    LOGICAL,PARAMETER :: is_standalone = .TRUE.

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********      OM wrapper (LADDIE).  Init phase 2 method.    ****************"
       PRINT*,"******************************************************************************"
       PRINT*,""
    END IF

    ! STAGE >1: receive ISM ice-draft geometry from OM_ImpFB and push it into `forcing`
    ! (forcing%Hib + masks), recomputing b-grid slopes on the LADDIE mesh via the UFEMISM
    ! utility calc_ice_shelf_base_slopes; then halo_exchange. (See TODO.md.)

    CALL initialise_transects_SGD( mesh, forcing)
    CALL initialise_laddie_model( mesh, laddie, forcing, is_standalone)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase2


  !--------------------------------------------------------------------------------------
  ! Run: advance LADDIE. Stage 1 reproduces a standalone run, so we simply call
  ! run_laddie_model (which integrates the ocean layer to quasi-equilibrium over a fixed
  ! duration, C%time_duration_laddie). No import/export exchange yet.
  SUBROUTINE FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB,OM_ImpFB,rc_local)

    TYPE(ESMF_config),INTENT(INOUT)                :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL  :: OM_ExpFB, OM_ImpFB
    TYPE(ESMF_VM),INTENT(IN)                       :: vm
    INTEGER,INTENT(OUT),OPTIONAL                   :: rc_local

    INTEGER            :: localPet, rc
    LOGICAL            :: verbose_coupling
    ! Stage 1 matches the standalone LADDIE_program parameters exactly.
    REAL(dp),PARAMETER :: time          = 0.0_dp   ! [yr] reference time stamp (see run_laddie_model)
    LOGICAL,PARAMETER  :: is_initial    = .FALSE.
    LOGICAL,PARAMETER  :: is_standalone = .TRUE.

    rc_local = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! STAGE >1: IF (PRESENT(OM_ImpFB)) CALL sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc)

    IF (localPet.EQ.0) THEN
       WRITE (OM_outputUnit,*) 'FISOC is about to call the LADDIE run method.'
    END IF
    CALL ESMF_VMBarrier(vm, rc=rc)

    ! NOTE: run_laddie_model does NOT integrate to `time`; it runs a fixed duration
    ! (C%time_duration_laddie days, or _init when is_initial) to quasi-steady state under
    ! the current geometry. For stage-1 reproduce-standalone, configure FISOC to call the
    ! OM once. STAGE >1 will set is_initial=.TRUE. on the first call only and is_standalone
    ! =.FALSE., and derive `time` from the FISOC clock.
    CALL run_laddie_model( mesh, laddie, forcing, time, is_initial, is_standalone)

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (localPet.EQ.0) THEN
       WRITE (OM_outputUnit,*) 'FISOC has just called the LADDIE run method.'
    END IF

    ! STAGE >1: IF (PRESENT(OM_ExpFB)) CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc)

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
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********      OM wrapper (LADDIE).  Finalise method.       *****************"
       PRINT*,"******************************************************************************"
       PRINT*,""
    END IF

    ! Finalise PETSc. Do NOT call MPI_FINALIZE here: ESMF/FISOC owns MPI finalisation.
    CALL PetscFinalize( perr)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_Wrapper_Finalize


  !--------------------------------------------------------------------------------
  ! Build an ESMF_Mesh from LADDIE's (UPSY) mesh.
  !
  ! UPSY uses GLOBAL vertex numbering (vi = 1..nV) with the full mesh geometry replicated
  ! on every process, and each vertex owned by exactly one process: mesh%V_owning_process.
  ! This maps directly onto ESMF's model (each node has a unique global id and one owner),
  ! so - unlike the FVCOM wrapper - NO one-to-many route handle is needed, and no MPI
  ! reduction to discover owners (UPSY already provides V_owning_process globally).
  !
  ! Per PET we present: the triangles this PET owns (ti1:ti2) as ESMF elements, and the
  ! vertices those triangles reference (owned + neighbour/halo) as ESMF nodes, each tagged
  ! with its true owner. elementConn uses LOCAL indices into the per-PET node list.
  !
  ! LADDIE works in projected stereographic x,y (metres), so the mesh is Cartesian.
  !
  ! THINGS TO VERIFY during stage-1 testing (flagged rather than assumed):
  !   * par%i (UPSY process rank) must equal the ESMF localPet for nodeOwners to be correct.
  !     True when LADDIE's communicator is the OM VM communicator (global comm, stage 1).
  !   * UPSY partitions vertices and triangles INDEPENDENTLY (vi1:vi2 vs ti1:ti2 are
  !     separate ranges). Here we list only nodes referenced by locally-owned elements (the
  !     FVCOM-proven convention). If ESMF requires every owned node to also be listed by its
  !     owner, we may need to add owned vertices (vi1:vi2) to the local node list. Verify
  !     against ESMF behaviour / the multi-PET MeshCreate docs.
  !   * Coupling fields live on different LADDIE grids (melt on a-grid/vertices=nodes,
  !     velocity on b-grid/triangles=elements). Stage 1 exchanges nothing, so this is not
  !     handled yet; node-based fields (melt, draft) will use the nodes built here.
  !--------------------------------------------------------------------------------
  SUBROUTINE LADDIE2ESMF_mesh(FISOC_config,mesh,OM_mesh,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)  :: FISOC_config
    TYPE(type_mesh),INTENT(IN)       :: mesh
    TYPE(ESMF_mesh),INTENT(OUT)      :: OM_mesh
    TYPE(ESMF_VM),INTENT(IN)         :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    INTEGER                          :: localPet, ti, vi, kk, lid
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
    !    recording a global-to-local index map.
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

    ! 3. Fill the node arrays.
    ALLOCATE(nodeIds(nLocalNodes), nodeOwners(nLocalNodes), nodeCoords(2*nLocalNodes))
    DO vi = 1, mesh%nV
       lid = gi2local(vi)
       IF (lid > 0) THEN
          nodeIds(lid)        = vi                          ! global node id
          nodeOwners(lid)     = mesh%V_owning_process(vi)   ! owning PET (process rank, 0-based)
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
       WRITE (OM_outputUnit,*) 'LADDIE2ESMF_mesh: created ESMF mesh, ', &
            nLocalNodes,' local nodes, ',nLocalElems,' local elements.'
    END IF

    DEALLOCATE(gi2local, nodeIds, nodeOwners, nodeCoords, elemIds, elemTypes, elemConn)

  END SUBROUTINE LADDIE2ESMF_mesh


  ! ====================================================================================
  ! STAGE >1 STUBS (not called yet; kept here to fix the intended structure).
  ! ====================================================================================
  !
  ! SUBROUTINE getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc)
  !   ! LADDIE -> ESMF export. Loop the export field bundle; for each field copy LADDIE
  !   ! data (on owned vertices vi1:vi2) into the ESMF field pointer. e.g.
  !   !   CASE ('OM_bmb')   ptr(:) = laddie%melt( mesh%vi1:mesh%vi2)
  !   !   CASE ('OM_z_l0')  ptr(:) = forcing%Hib( mesh%vi1:mesh%vi2)   ! ice draft
  ! END SUBROUTINE getFieldDataFromOM
  !
  ! SUBROUTINE sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc)
  !   ! ESMF -> LADDIE import. For each ISM field, write owned vertices then halo-exchange:
  !   !   CASE ('ISM_z_l0') forcing%Hib( mesh%vi1:mesh%vi2) = ptr(:)
  !   ! then UPSY halo_exchange on forcing%Hib, and recompute b-grid draft slopes via the
  !   ! UFEMISM routine calc_ice_shelf_base_slopes. No ESMF one-to-many route handle needed.
  ! END SUBROUTINE sendFieldDataToOM

END MODULE FISOC_OM_Wrapper
