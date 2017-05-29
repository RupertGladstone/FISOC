MODULE FISOC_ISM_MOD
  
  USE ESMF
  USE FISOC_ISM_Wrapper
  USE FISOC_utils_MOD
    
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
         userRoutine=FISOC_ISM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_ISM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_ISM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

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

    TYPE(ESMF_config)          :: FISOC_config
    TYPE(ESMF_mesh)            :: ISM_mesh
    TYPE(ESMF_grid)            :: ISM_grid
    TYPE(ESMF_fieldBundle)     :: ISM_ExpFB
    TYPE(ESMF_VM)              :: vm

    CHARACTER(len=ESMF_MAXSTR) :: label, ISM_gridType
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: ISM_ReqVarList(:),ISM_DerVarList(:)
 
    rc = ESMF_FAILURE

    msg = "ISM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! extract a list of required ISM variables from the FISOC config object
    CALL ESMF_GridCompGet(FISOC_ISM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    label = 'FISOC_ISM_ReqVars:'
    CALL FISOC_getListFromConfig(FISOC_config, label, ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    label = 'FISOC_ISM_DerVars:' ! also derived ISM variables
    CALL FISOC_getListFromConfig(FISOC_config, label, ISM_DerVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, ISM_gridType, label='ISM_gridType:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! create empty field bundle
    ISM_ExpFB = ESMF_FieldBundleCreate(name='ISM export fields', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! The model-specific initialisation adds the ISM vars to the field bundle.
    ! This is followed by adding non-model-specific derived variables. 
    SELECT CASE (ISM_gridType)

    CASE("ESMF_grid","ESMF_Grid")
       CALL FISOC_ISM_Wrapper_Init_Phase1(ISM_ReqVarList,ISM_ExpFB,ISM_grid,&
            FISOC_config,vm,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL FISOC_populateFieldBundle(ISM_DerVarList,ISM_ExpFB,ISM_grid,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CASE("ESMF_mesh","ESMF_Mesh")
       CALL FISOC_ISM_Wrapper_Init_Phase1(ISM_ReqVarList,ISM_ExpFB,ISM_mesh,&
            FISOC_config,vm,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL FISOC_populateFieldBundle(ISM_DerVarList,ISM_ExpFB,ISM_mesh,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CASE DEFAULT
       msg = "ERROR: FISOC does not recognise ISM_gridType: "//ISM_gridType
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT


    ! Calculate values for derived variables from the model-specific ISM vars.
    CALL FISOC_ISM_calcDerivedFields(ISM_ExpFB,FISOC_config,rc)
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

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_init_phase1


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_init_phase2(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)

    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_config)      :: FISOC_config
    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_VM)          :: vm

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_ISM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(ISM_ImpSt, "ISM import fields", ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(ISM_ExpSt, "ISM export fields", ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_ISM_Wrapper_Init_Phase2(ISM_ImpFB,ISM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "ISM initialise phase 2 (allows the ISM access to the OM initial state) "
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_init_phase2
  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_run(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)

    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc
    
    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    LOGICAL                :: verbose_coupling
    TYPE(ESMF_config)      :: FISOC_config
!    INTEGER                :: localPet
    TYPE(ESMF_VM)          :: vm

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_ISM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
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
    
    CALL FISOC_ISM_calcDerivedFields_pre(ISM_ExpFB,FISOC_config,rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_ISM_Wrapper_Run(FISOC_config,vm,ISM_ExpFB,ISM_ImpFB,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    CALL FISOC_ISM_calcDerivedFields_post(ISM_ExpFB,FISOC_config,rc)
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
    INTEGER, INTENT(OUT) :: rc

    INTEGER                      :: localPet
    TYPE(ESMF_VM)                :: vm
    TYPE(ESMF_config)            :: FISOC_config
    TYPE(ESMF_fieldbundle)       :: ISM_ImpFB, ISM_ExpFB
    INTEGER                      :: FieldCount,ii
    TYPE(ESMF_field),ALLOCATABLE :: ImpFieldList(:)
    TYPE(ESMF_field),ALLOCATABLE :: ExpFieldList(:)
    TYPE(ESMF_mesh)              :: ISM_mesh

    CHARACTER(len=ESMF_MAXSTR) :: name
    
    rc = ESMF_FAILURE

    msg = "ISM finalise: destroy fields, bundles and mesh"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompGet(FISOC_ISM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ISM_Wrapper_Finalize(FISOC_config,localPet,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    ! destroy ISM export fields

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
    ALLOCATE(ExpFieldList(FieldCount))
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldList=ExpFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,FieldCount
       
       CALL ESMF_FieldGet(field=ExpFieldList(ii),name=name, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
!       CALL ESMF_VMBarrier(vm, rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
       CALL ESMF_FieldDestroy(ExpFieldList(ii), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO
    
    DEALLOCATE(ExpFieldList)
    
    CALL ESMF_FieldBundleDestroy(ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    
    ! destroy ISM import fields
    
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
    ALLOCATE(ImpFieldList(FieldCount))
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldList=ImpFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
!    CALL ESMF_FieldGet(ImpFieldList(1), mesh=ISM_mesh, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    DO ii = 1,FieldCount
       CALL ESMF_FieldDestroy(ImpFieldList(ii), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO
    
!    CALL ESMF_MeshDestroy(ISM_mesh,rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    DEALLOCATE(ImpFieldList)
    
    CALL ESMF_FieldBundleDestroy(ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_ISM_finalise
  


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_calcDerivedFields(ISM_ExpFB,FISOC_config,rc)

    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL          :: rc
    
    CALL FISOC_ISM_calcDerivedFields_pre(ISM_ExpFB,FISOC_config,rc)
    CALL FISOC_ISM_calcDerivedFields_post(ISM_ExpFB,FISOC_config,rc)

  END SUBROUTINE FISOC_ISM_calcDerivedFields

  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_calcDerivedFields_pre(ISM_ExpFB,FISOC_config,rc)

    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: FISOC_ISM_DerVarList(:),FISOC_ISM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)             :: label
    INTEGER                                :: ii, numDerVars

    rc = ESMF_FAILURE

    ! extract a list of derived ISM variables from the configuration object
    label = 'FISOC_ISM_DerVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_ISM_DerVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    label = 'FISOC_ISM_ReqVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    numDerVars = SIZE(FISOC_ISM_DerVarList)

    DO ii = 1,numDerVars

       SELECT CASE(FISOC_ISM_DerVarList(ii))

       CASE ("ISM_z_l0_previous")
          CALL ISM_derive_z_l0_previous(ISM_ExpFB,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE ("ISM_dddt","ISM_dTdz_l0","ISM_z_l0_linterp")

       CASE DEFAULT
          msg="ERROR: derived variable name not recognised: "//FISOC_ISM_DerVarList(ii)
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END SELECT

    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_calcDerivedFields_pre



  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_calcDerivedFields_post(ISM_ExpFB,FISOC_config,rc)

    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)   :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: FISOC_ISM_DerVarList(:),FISOC_ISM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)             :: label
    INTEGER                                :: ii, numDerVars

    rc = ESMF_FAILURE

    ! extract a list of derived ISM variables from the configuration object
    label = 'FISOC_ISM_DerVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_ISM_DerVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    label = 'FISOC_ISM_ReqVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    numDerVars = SIZE(FISOC_ISM_DerVarList)

    DO ii = 1,numDerVars

       SELECT CASE(FISOC_ISM_DerVarList(ii))

       CASE ("ISM_z_l0_previous","ISM_z_l0_linterp")

       CASE ("ISM_dddt")
          CALL ISM_derive_dddt(ISM_ExpFB,FISOC_config,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE ("ISM_dTdz_l0")
          CALL ISM_derive_dTdz_l0(ISM_ExpFB,rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CASE DEFAULT
          msg="ERROR: derived variable name not recognised: "//FISOC_ISM_DerVarList(ii)
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END SELECT

    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_calcDerivedFields_post



  !------------------------------------------------------------------------------
  SUBROUTINE ISM_derive_z_l0_previous(ISM_ExpFB,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)            :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL                   :: rc

    TYPE(ESMF_field)           :: ISM_z_l0
    TYPE(ESMF_field)           :: ISM_z_l0_previous
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_previous_ptr(:) 

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0_previous", field=ISM_z_l0_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0_previous, localDe=0, farrayPtr=ISM_z_l0_previous_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ISM_z_l0_previous_ptr = ISM_z_l0_ptr

    NULLIFY(ISM_z_l0_previous_ptr)
    NULLIFY(ISM_z_l0_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_z_l0_previous


  !------------------------------------------------------------------------------
  ! dddt, short for dd/dt, is the rate of change of depth with time, for use by 
  ! ocean model.
  !
  SUBROUTINE ISM_derive_dddt(ISM_ExpFB,FISOC_config,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL         :: rc
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config

    TYPE(ESMF_field)           :: ISM_z_l0
    TYPE(ESMF_field)           :: ISM_z_l0_previous
    TYPE(ESMF_field)           :: ISM_dddt
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_previous_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dddt_ptr(:) 
    REAL(ESMF_KIND_R8)         :: ISM_dt
    INTEGER                    :: ISM_dt_int

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0_previous", field=ISM_z_l0_previous, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0_previous, localDe=0, farrayPtr=ISM_z_l0_previous_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dddt", field=ISM_dddt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dddt, localDe=0, farrayPtr=ISM_dddt_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_int, 'ISM_dt_sec',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_dt = REAL(ISM_dt_int,ESMF_KIND_R8)
    ISM_dddt_ptr = (ISM_z_l0_ptr - ISM_z_l0_previous_ptr) / ISM_dt

    ! first time we set it to zero
    WHERE (ISM_z_l0_previous_ptr .EQ. 0.0) ISM_dddt_ptr = 0.0

    NULLIFY(ISM_z_l0_previous_ptr)
    NULLIFY(ISM_z_l0_ptr)
    NULLIFY(ISM_dddt_ptr)

    RC = ESMF_SUCCESS

  END SUBROUTINE ISM_derive_dddt



  !------------------------------------------------------------------------------

  SUBROUTINE ISM_derive_dTdz_l0(ISM_ExpFB,rc)

!    CHARACTER(len=ESMF_MAXSTR),INTENT(IN)   :: FISOC_ISM_ReqVarList(:)
    TYPE(ESMF_fieldBundle),INTENT(INOUT)    :: ISM_ExpFB
    INTEGER, INTENT(OUT),OPTIONAL           :: rc

    TYPE(ESMF_field)           :: ISM_z_l0,ISM_z_l1
    TYPE(ESMF_field)           :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_field)           :: ISM_dTdz_l0
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:), ISM_z_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_temperature_l0_ptr(:),ISM_temperature_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:)

    INTEGER :: localDeCount

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_temperature_l0", field=ISM_temperature_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l0, localDe=0, farrayPtr=ISM_temperature_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_temperature_l1", field=ISM_temperature_l1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l1, localDe=0, farrayPtr=ISM_temperature_l1_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldGet(field=ISM_z_l0, localDeCount=localDeCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l1", field=ISM_z_l1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l1, localDe=0, farrayPtr=ISM_z_l1_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dTdz_l0", field=ISM_dTdz_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dTdz_l0, localDe=0, farrayPtr=ISM_dTdz_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! the simplest approximation for the vertical temperature gradient at the ice base is to 
    ! use the lowest two levels and assume the gradient doesn't change between these two 
    ! levels.
    ISM_dTdz_l0_ptr(:) = (ISM_temperature_l1_ptr-ISM_temperature_l0_ptr) / (ISM_z_l1_ptr-ISM_z_l0_ptr)

    NULLIFY(ISM_dTdz_l0_ptr)
    NULLIFY(ISM_temperature_l1_ptr)
    NULLIFY(ISM_temperature_l0_ptr)
    NULLIFY(ISM_z_l1_ptr)
    NULLIFY(ISM_z_l0_ptr)

    rc = ESMF_SUCCESS
    
  END SUBROUTINE ISM_derive_dTdz_l0
  
END MODULE FISOC_ISM_MOD
