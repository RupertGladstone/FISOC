MODULE FISOC_ISM
  
  USE ESMF
  USE FISOC_ISM_Wrapper
  USE FISOC_utils
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_ISM_register
  
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
    TYPE(ESMF_GridComp)        :: FISOC_ISM
    TYPE(ESMF_State)           :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)           :: FISOC_clock
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_config)          :: config
    TYPE(ESMF_mesh)            :: ISM_mesh
    TYPE(ESMF_fieldBundle)     :: ISM_ExpFB
    CHARACTER(len=ESMF_MAXSTR) :: ISM_meshFile

!    TYPE(ESMF_field)           :: ISM_temperature_l0, ISM_temperature_l1
!    TYPE(ESMF_field)           :: ISM_z_l0, ISM_z_l1
!    TYPE(ESMF_field)           :: ISM_dTdz_l0, ISM_z_l0_previous
!    REAL(ESMF_KIND_R8),POINTER :: ISM_temperature_l0_ptr(:),ISM_temperature_l1_ptr(:) 
!    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:),ISM_z_l1_ptr(:) 
!    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:),ISM_z_l0_previous_ptr(:) 

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: ISM_ReqVarList(:),ISM_DerVarList(:)
    CHARACTER(len=ESMF_MAXSTR) :: label
    INTEGER                    :: ii

    rc = ESMF_SUCCESS

    msg = "ISM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! extract a list of required ISM variables from the FISOC config object
    CALL ESMF_GridCompGet(FISOC_ISM, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    label = 'FISOC_ISM_ReqVars:'
    CALL FISOC_getStringListFromConfig(config, label, ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    label = 'FISOC_ISM_DerVars:' ! also derived ISM variables
    CALL FISOC_getStringListFromConfig(config, label, ISM_DerVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! create empty field bundle
    ISM_ExpFB = ESMF_FieldBundleCreate(name='ISM export fields', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! model-specific initialisation
    CALL FISOC_ISM_Wrapper_Init(ISM_ReqVarList,ISM_ExpFB,ISM_mesh,config,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! The model-specific initialisation added the ISM vars to the field bundle, 
    ! and here we add non-model-specific derived variables. 
    CALL FISOC_populateFieldBundle(ISM_DerVarList,ISM_ExpFB,ISM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Calculate values for derived variables from the model-specific ISM vars.
    CALL FISOC_ISM_calcDerivedFields(ISM_ExpFB,config,rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! we only add the fields to the import state as a way of letting the coupler get hold of the 
    ! grid.  There must be a better way to do this.
    CALL ESMF_StateAdd(ISM_ImpSt, (/ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(ISM_ExpSt, (/ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "ISM initialise phase 1 complete"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

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
  SUBROUTINE FISOC_ISM_run(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER                :: petCount
    INTEGER, INTENT(OUT)   :: rc
    
    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    LOGICAL                :: verbose_coupling
    TYPE(ESMF_config)      :: config

    rc = ESMF_FAILURE

    ! query the FISOC config
    CALL ESMF_GridCompGet(FISOC_ISM, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_ConfigGetAttribute(config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! extract field bundles from import and export states to send to model-specific wrapper

    CALL ESMF_StateGet(ISM_ImpSt, "ISM import fields", ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(ISM_ExpSt, "ISM export fields", ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    CALL FISOC_ISM_Wrapper_Run(ISM_ImpFB,ISM_ExpFB,config,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    msg = "ISM run complete for current timestep"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_run


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_finalise(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount
    INTEGER, INTENT(OUT) :: rc

    TYPE(ESMF_fieldbundle)       :: ISM_ImpFB, ISM_ExpFB
    INTEGER                      :: FieldCount,ii
    TYPE(ESMF_field),ALLOCATABLE :: FieldList(:)
    TYPE(ESMF_mesh)              :: ISM_mesh
    
    rc = ESMF_FAILURE

    msg = "ISM finalise: destroy fields, bundles and mesh"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)


    CALL ESMF_StateGet(ISM_ImpSt, "ISM import fields", ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(FieldList(FieldCount))
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldCount=FieldCount, fieldList=FieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldGet(FieldList(1), mesh=ISM_mesh, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_MeshDestroy(ISM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,FieldCount
       CALL ESMF_FieldDestroy(FieldList(ii), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    DEALLOCATE(FieldList)

    CALL ESMF_FieldBundleDestroy(ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    



    CALL ESMF_StateGet(ISM_ExpSt, "ISM export fields", ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(FieldList(FieldCount))
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=FieldCount, fieldList=FieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,FieldCount
       CALL ESMF_FieldDestroy(FieldList(ii), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    DEALLOCATE(FieldList)

    CALL ESMF_FieldBundleDestroy(ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_finalise
  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_calcDerivedFields(ISM_ExpFB,config,rc)

    TYPE(ESMF_config),INTENT(IN) :: config
    INTEGER, INTENT(OUT)         :: rc
    TYPE(ESMF_fieldBundle)       :: ISM_ExpFB

    TYPE(ESMF_field)       :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_field)       :: ISM_z_l0, ISM_z_l1
    TYPE(ESMF_field)       :: ISM_dTdz_l0, ISM_z_l0_previous
    REAL(ESMF_KIND_R8),POINTER :: ISM_temperature_l0_ptr(:),ISM_temperature_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:),ISM_z_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:),ISM_z_l0_previous_ptr(:) 

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_temperature_l0", field=ISM_temperature_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l0, localDe=0, farrayPtr=ISM_temperature_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_temperature_l1", field=ISM_temperature_l1, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l1, localDe=0, farrayPtr=ISM_temperature_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l1", field=ISM_z_l1, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l1, localDe=0, farrayPtr=ISM_z_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dTdz_l0", field=ISM_dTdz_l0, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dTdz_l0, localDe=0, farrayPtr=ISM_dTdz_l0_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! the simplest approximation for the vertical temperature gradient at the ice base is to 
    ! use the lowest two levels and assume the gradient doesn't change between these two 
    ! levels.
    ISM_dTdz_l0_ptr(:) = (ISM_temperature_l1_ptr-ISM_temperature_l0_ptr) / (ISM_z_l1_ptr-ISM_z_l0_ptr)

!    PRINT*,"temperature gradient: add configurable options for more advanced calculations"

!    print*,"previous B to be calculated only at run time"

  END SUBROUTINE FISOC_ISM_calcDerivedFields

END MODULE FISOC_ISM
