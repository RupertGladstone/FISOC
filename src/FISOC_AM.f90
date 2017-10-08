MODULE FISOC_AM_MOD
  
  USE ESMF
  USE FISOC_utils_MOD
  USE FISOC_AM_Wrapper
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_AM_register
    
CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_AM_register(FISOC_AM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_AM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_AM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_AM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetEntryPoint(FISOC_AM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_AM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_AM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_AM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_AM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_AM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_AM_register


  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two stages.  The first stage is for 
  ! independent initialisation of the copmonents and the second stage allows 
  ! modifications based on initialisations of other components or inter-component 
  ! consistency checks.
  SUBROUTINE FISOC_AM_init_phase1(FISOC_AM, AM_ImpSt, AM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)        :: FISOC_AM
    TYPE(ESMF_State)           :: AM_ImpSt, AM_ExpSt 
    TYPE(ESMF_Clock)           :: FISOC_clock
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_config)          :: FISOC_config
    TYPE(ESMF_grid)            :: AM_grid
    TYPE(ESMF_fieldBundle)     :: AM_ISM_ExpFB, AM_OM_ExpFB!,AM_ExpFBcum
    TYPE(ESMF_VM)              :: vm
    CHARACTER(len=ESMF_MAXSTR) :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: AM_ReqVarList(:)

    rc = ESMF_FAILURE

    msg = "AM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! FISOC_AM is the top level gridded component for the AM (atmosphere model).
    ! vm is the virtual machine, i.e. the parallel context.
    CALL ESMF_GridCompGet(FISOC_AM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! Create empty field bundles to be populated by model-specific code.
    ! There is one export field bundle to send to the ISM and one to send to 
    ! the OM. 
    AM_ISM_ExpFB = ESMF_FieldBundleCreate(name='AM ISM export fields', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    AM_OM_ExpFB = ESMF_FieldBundleCreate(name='AM OM export fields', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

! we may activate the cumulators later...
!    AM_ExpFBcum = ESMF_FieldBundleCreate(name='AM export field cumulator', rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! model-specific initialisation
    CALL FISOC_AM_Wrapper_Init_Phase1(FISOC_config,vm,AM_ISM_ExpFB,AM_OM_ExpFB,AM_grid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!    CALL FISOC_initCumulatorFB(AM_ExpFB,AM_ExpFBcum,rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!    CALL FISOC_zeroBundle(AM_ExpFBcum,rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! We only add the AM export field bundles to the import state as a way of 
    ! letting the coupler get hold of the AM grid.  The coupler will remove them 
    ! after initialisation. In the future it will be possible
    ! to add either a mesh or grid object to a gridComp.
    CALL ESMF_StateAdd(AM_ImpSt, (/AM_ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateAdd(AM_ImpSt, (/AM_OM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(AM_ExpSt, (/AM_ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateAdd(AM_ExpSt, (/AM_OM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!    CALL ESMF_StateAdd(AM_ExpSt, (/AM_ExpFBcum/), rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "AM initialise phase 1 complete"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_AM_init_phase1
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_AM_init_phase2(FISOC_AM, AM_ImpSt, AM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_AM
    TYPE(ESMF_State)       :: AM_ImpSt, AM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_VM)          :: vm
    TYPE(ESMF_config)      :: FISOC_config
    TYPE(ESMF_fieldbundle) :: AM_ISM_ImpFB, AM_OM_ImpFB, AM_ISM_ExpFB, AM_OM_ExpFB

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_AM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

! couplers need updating.  Couplers create import field bundles during init.
!    CALL ESMF_StateGet(AM_ImpSt, "AM ISM import fields", AM_ISM_ImpFB, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
!    CALL ESMF_StateGet(AM_ImpSt, "AM OM import fields", AM_OM_ImpFB, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(AM_ExpSt, "AM ISM export fields", AM_ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    CALL ESMF_StateGet(AM_ExpSt, "AM OM export fields", AM_OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_AM_Wrapper_Init_Phase2(FISOC_config,vm,AM_ISM_ImpFB,AM_OM_ImpFB,AM_ISM_ExpFB,AM_OM_ExpFB,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "AM initialise phase 2 (allows the AM access to ISM and OM initial states) "
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_AM_init_phase2

  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_AM_run(FISOC_AM, AM_ImpSt, AM_ExpSt, FISOC_clock, rc)

    TYPE(ESMF_GridComp)        :: FISOC_AM
    TYPE(ESMF_State)           :: AM_ImpSt, AM_ExpSt 
    TYPE(ESMF_Clock)           :: FISOC_clock
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_VM)              :: vm
    INTEGER                    :: localPet
    TYPE(ESMF_fieldbundle)     :: AM_ISM_ImpFB, AM_OM_ImpFB, AM_ISM_ExpFB, AM_OM_ExpFB
    TYPE(ESMF_config)          :: FISOC_config
    LOGICAL                    :: verbose_coupling

    rc = ESMF_FAILURE
   
    CALL ESMF_GridCompGet(FISOC_AM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(AM_ExpSt, "AM ISM export fields", AM_ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateGet(AM_ExpSt, "AM OM export fields", AM_OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_StateGet(AM_ImpSt, "AM ISM import fields", AM_ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateGet(AM_ImpSt, "AM OM import fields", AM_OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! call the model-specific wrapper
! For now don't pass the import field bundles as the couplers need 
! updating before these will contain usable information.
!    CALL FISOC_AM_Wrapper_Run(FISOC_config,vm,          &
!         AM_ISM_ExpFB=AM_ISM_ExpFB,                     & 
!         AM_OM_ExpFB =AM_OM_ExpFB ,                     & 
!         AM_ISM_ImpFB=AM_ISM_ImpFB,                     &
!         AM_OM_ImpFB =AM_OM_ImpFB ,                     &
!         rc_local=rc)
    CALL FISOC_AM_Wrapper_Run(FISOC_config,vm,          &
         AM_ISM_ExpFB=AM_ISM_ExpFB,                     & 
         AM_OM_ExpFB =AM_OM_ExpFB ,                     & 
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
       CALL FISOC_AM_finalise(FISOC_AM, AM_ImpSt, AM_ExpSt, FISOC_clock, rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_AM_run
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_AM_finalise(FISOC_AM, AM_ImpSt, AM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_AM
    TYPE(ESMF_State)       :: AM_ImpSt, AM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_VM)                :: vm
    TYPE(ESMF_config)            :: FISOC_config
    TYPE(ESMF_fieldbundle)       :: AM_ISM_ImpFB, AM_OM_ImpFB, AM_ISM_ExpFB, AM_OM_ExpFB
    INTEGER                      :: FieldCount,ii, localPet
    TYPE(ESMF_grid)              :: AM_grid
    
    rc = ESMF_FAILURE


    msg = "AM finalise: destroy fields, bundles and grid"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompGet(FISOC_AM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(AM_ImpSt, "AM import fields", AM_ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    CALL ESMF_StateGet(AM_ImpSt, "AM import fields", AM_OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_emptyFB_grid(AM_ISM_ImpFB)
    CALL FISOC_emptyFB_grid(AM_OM_ImpFB)

    CALL ESMF_FieldBundleDestroy(AM_ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldBundleDestroy(AM_OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_StateGet(AM_ExpSt, "AM export fields", AM_ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    CALL ESMF_StateGet(AM_ExpSt, "AM export fields", AM_OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_emptyFB_grid(AM_ISM_ExpFB)
    CALL FISOC_emptyFB_grid(AM_OM_ExpFB)

    CALL ESMF_FieldBundleDestroy(AM_ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldBundleDestroy(AM_OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_AM_Wrapper_Finalize(FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_AM_finalise  



END MODULE FISOC_AM_MOD
