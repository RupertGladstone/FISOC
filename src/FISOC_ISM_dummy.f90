MODULE FISOC_ISM
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_ISM_register
  
CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_register(FISOC_ISM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init, rc=rc)
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
  SUBROUTINE FISOC_ISM_init(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc

    TYPE(ESMF_config)    :: config
    TYPE(ESMF_mesh)      :: ISM_mesh
    TYPE(ESMF_field)     :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_field)     :: ISM_z_l0, ISM_z_l1
    INTEGER              :: localrc
    CHARACTER(len=ESMF_MAXSTR) :: ISM_meshFile, msg
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
    
    rc = ESMF_SUCCESS

    NULLIFY(ISM_temperature_l0_ptr)

    msg = "ISM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! Get mesh file name from the config file
    CALL ESMF_GridCompGet(FISOC_ISM, config=config, rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    CALL ESMF_ConfigGetAttribute(config, ISM_meshFile, label='ISM_meshFile:', rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! create ISM mesh and use it to create zeroed fields for the ISM import and export states
    ISM_mesh = ESMF_MeshCreate(filename=ISM_meshFile, &
            filetypeflag=ESMF_FILEFORMAT_ESMFMESH, &
            rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

! Note weakness: currently regridding in 2d instead of a 2d manifold in 3d space.

!    CALL ESMF_MeshGet(ISM_mesh, parametricDim=parametricDim, spatialDim=spatialDim, &
!                    nodalDistgrid=nodalDistgrid, numOwnedNodes=numOwnedNodes, &
!                    numOwnedElements=numOwnedElements, rc=rc)
!    print *,"check this with elmer mesh header file"
!    print *,"mesh stuff", parametricDim, spatialDim, numOwnedNodes, numOwnedElements !, ownedNodeCoords
    
    ISM_temperature_l0 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l0, localDe=0, farrayPtr=ISM_temperature_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_temperature_l0_ptr(:) = -1.0

    ISM_temperature_l1 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l1, localDe=0, farrayPtr=ISM_temperature_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_temperature_l1_ptr(:) = -1.1
    
    ISM_z_l0 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_z_l0_ptr(:) = -100.0
    
    ISM_z_l1 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l1, localDe=0, farrayPtr=ISM_z_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_z_l1_ptr(:) = -75.0
    
    msg = "ISM created mesh and fields"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

print*,"make field bundle"
print*,"add bundle to export state"

!    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!get some mesh info, check it looks sensible
!create a zeroed array on it as a new field
!fields for export:
!level0 temperature
!level0 z coord
!level1 temperature
!level1 z coord
!fields for import:
!level0 melt rate

  END SUBROUTINE FISOC_ISM_init
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_run(FISOC_ISM, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount, localrc
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_SUCCESS

    print *," FISOC_ISM_run needs writing"
    
  END SUBROUTINE FISOC_ISM_run
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_finalise(FISOC_ISM, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount, localrc
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_SUCCESS

    print *," FISOC_ISM_finalise needs writing"
    
  END SUBROUTINE FISOC_ISM_finalise
  
END MODULE FISOC_ISM
