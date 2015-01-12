MODULE FISOC_ISM
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_ISM_register
  
  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_register(FISOC_ISM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

!    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
!         userRoutine=FISOC_ISM_init, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_ISM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_ISM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_ISM_register

  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two stages.  The first stage is for independent 
  ! initialisation of the ISM and the second stage occurrs after the OM has been 
  ! initialised, to allow inter-component consistency checks or use of the OM 
  ! state to complete ISM initalisation.
  SUBROUTINE FISOC_ISM_init_phase1(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_config)      :: config
    TYPE(ESMF_mesh)        :: ISM_mesh
    TYPE(ESMF_field)       :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_field)       :: ISM_z_l0, ISM_z_l1
    TYPE(ESMF_fieldBundle) :: ISM_ImpFB, ISM_ExpFB
    CHARACTER(len=ESMF_MAXSTR) :: ISM_meshFile
    REAL(ESMF_KIND_R8),POINTER :: ISM_temperature_l0_ptr(:),ISM_temperature_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:),ISM_z_l1_ptr(:) 

!    real(ESMF_KIND_R8)        :: ownedNodeCoords(:)
    integer                   :: numOwnedElements
    logical                   :: isMemFreed
    type(ESMF_CoordSys_Flag)  :: coordSys
    integer                   :: parametricDim
    integer                   :: spatialDim
    type(ESMF_DistGrid)       :: nodalDistgrid
    type(ESMF_DistGrid)       :: elementDistgrid
    integer                   :: numOwnedNodes

INTEGER :: fieldcount
    
    rc = ESMF_SUCCESS

    NULLIFY(ISM_temperature_l0_ptr)

    msg = "ISM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! Get mesh file name from the config file
    CALL ESMF_GridCompGet(FISOC_ISM, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    CALL ESMF_ConfigGetAttribute(config, ISM_meshFile, label='ISM_meshFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! create ISM mesh and use it to create zeroed fields for the ISM import and export states
    ISM_mesh = ESMF_MeshCreate(filename=ISM_meshFile, &
            filetypeflag=ESMF_FILEFORMAT_ESMFMESH, &
            rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

! Note weakness: currently regridding in 2d instead of a 2d manifold in 3d space.

!    CALL ESMF_MeshGet(ISM_mesh, parametricDim=parametricDim, spatialDim=spatialDim, &
!                    nodalDistgrid=nodalDistgrid, numOwnedNodes=numOwnedNodes, &
!                    numOwnedElements=numOwnedElements, rc=rc)
!    print *,"check this with elmer mesh header file"
!    print *,"mesh stuff", parametricDim, spatialDim, numOwnedNodes, numOwnedElements !, ownedNodeCoords
    
    ISM_temperature_l0 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_temperature_l0", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l0, localDe=0, farrayPtr=ISM_temperature_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_temperature_l0_ptr(:) = -1.0

    ISM_temperature_l1 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_temperature_l1", rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l1, localDe=0, farrayPtr=ISM_temperature_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_temperature_l1_ptr(:) = -1.1
    
    ISM_z_l0 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_z_l0", rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_z_l0_ptr(:) = -100.0
    
    ISM_z_l1 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_z_l1", rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l1, localDe=0, farrayPtr=ISM_z_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_z_l1_ptr(:) = -75.0
    
    msg = "ISM created mesh and fields"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ISM_ExpFB = ESMF_FieldBundleCreate(name="ISM export fields", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_temperature_l0/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_temperature_l1/), rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_z_l0/), rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_z_l1/), rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=fieldcount, rc=rc)
!    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(ISM_ExpSt, (/ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "ISM bundled fields and added to export state"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

!    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  END SUBROUTINE FISOC_ISM_init_phase1


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_init_phase2(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    msg = "ISM initialise phase 2 allows the ISM access to the OM initial state "// &
         "(does nothing for dummy case)"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_init_phase2
  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_run(FISOC_ISM, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount
    INTEGER, INTENT(OUT) :: rc
    

    rc = ESMF_FAILURE

    msg = "ISM run has started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_run
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_finalise(FISOC_ISM, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    msg = "ISM finalise has started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_finalise
  
END MODULE FISOC_ISM
