
MODULE  FISOC_parent_mod
  
  USE ESMF

  USE FISOC_ISM, ONLY : FISOC_ISM_register
  USE FISOC_coupler, ONLY : FISOC_coupler_register
  USE FISOC_OM, ONLY : FISOC_OM_register
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_parent_register
  
  TYPE(ESMF_GridComp), SAVE :: FISOC_ISM, FISOC_OM
  TYPE(ESMF_CplComp),  SAVE :: FISOC_coupler
  TYPE(ESMF_State),    SAVE :: ISM_ImpSt, ISM_ExpSt, OM_ImpSt, OM_ExpSt
  
CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_parent_register(FISOC_parent, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_parent
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_init, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_RUN, &
         userRoutine=FISOC_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_finalize, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_parent_register
  
    
!------------------------------------------------------------------------------
  SUBROUTINE FISOC_init(FISOC_parent, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc

    TYPE(ESMF_State)     :: ISM_ImpSt, ISM_ExpSt, OM_ImpSt, OM_ExpSt
    TYPE(ESMF_config)    :: config
    INTEGER              :: petCount, localrc, urc
    CHARACTER(len=ESMF_MAXSTR) :: msg

    rc = ESMF_FAILURE

    msg = "Starting FISOC parent initialisation"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! Get config and petCount from the component object (pet is persistent execution thread)
    CALL ESMF_GridCompGet(FISOC_parent, config=config, petCount=petCount, rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return ! bail out

    ! Create and register child components and routines
    FISOC_ISM = ESMF_GridCompCreate(name="Ice Sheet Model", config=config, rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    FISOC_OM = ESMF_GridCompCreate(name="Ocean Model", config=config, rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    FISOC_coupler = ESMF_CplCompCreate(name="FISOC coupler", config=config, rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetServices(FISOC_ISM, FISOC_ISM_register, &
         userRc=urc, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetServices(FISOC_OM, FISOC_OM_register, &
         userRc=urc, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_CplCompSetServices(FISOC_coupler, FISOC_coupler_register, &
         userRc=urc, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "FISOC child routines created and registered"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ISM_ImpSt = ESMF_StateCreate(name='ISM import state', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_ExpSt = ESMF_StateCreate(name='ISM export state', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    OM_ImpSt = ESMF_StateCreate(name='OM import state', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    OM_ExpSt = ESMF_StateCreate(name='OM export state', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "Empty import and export states created"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompInitialize(FISOC_OM, &
         importState=OM_ImpSt, exportState=OM_ExpSt, &
         clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompInitialize(FISOC_ISM, &
         importState=ISM_ImpSt, exportState=ISM_ExpSt, &
         clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_cplCompInitialize(FISOC_coupler, &
         importState=ISM_ExpSt, exportState=OM_ImpSt, &
         clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompInitialize(FISOC_OM, &
         importState=OM_ImpSt, exportState=OM_ExpSt, &
         clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_cplCompInitialize(FISOC_coupler, &
         importState=OM_ExpSt, exportState=ISM_ImpSt, &
         clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompInitialize(FISOC_ISM, &
         importState=ISM_ImpSt, exportState=ISM_ExpSt, &
         clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "Completed calls to component initialiation routines"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    !***more child comps: OM and coupler
    
    !*** can access and check stuff to do with parallelization here... PET count...
    !*** can set up stuff to do with connectivity - how grids are distributed over pets
    ! (PET = persistent execution thread)
    
    !*** create child grids and their decompositions
    !*** maybe elmer needs to set the grid... can write an elmer solver to create the grid 
    ! in the proper ESMF way, using Elmer's own idea of how the grid should be...
    
    !***ROMS grid: need to see how Ufuk has done it...
    
    !*** create empty import and export states for all grid components, and call their initialize 
    ! routines (which have previously been registered)
    ! note components can have multiple initialisation phases (i.e. multiple initialisation routines 
    ! to be called in order determined by the programmer)

    msg = "FISOC initialise completed"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_init

!------------------------------------------------------------------------------
  SUBROUTINE FISOC_run(FISOC_parent, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc
 
    TYPE(ESMF_Alarm)     :: alarm_ocn, alarm_ice

    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ice", alarm_ice, rc=rc)
    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ocn", alarm_ocn, rc=rc)

    ! main timestepping loop
    DO WHILE (.NOT. ESMF_ClockIsStopTime(FISOC_clock, rc=rc))

       !*** run:
       !*** run the child components
       !*** coupler gets called multiple times with different import and export states.  
!       CALL ESMF_TimeIntervalPrint(ts_ice, rc=rc)
!       CALL ESMF_ClockPrint(FISOC_clock, options="advanceCount string isofrac", rc=rc)
!       CALL ESMF_ClockPrint(FISOC_clock, options="currTime string", rc=rc)
       
!       PRINT*,""
!       PRINT*,"New timestep. We need to call FISOC_proc"
       ! ***the FISOC_proc component accesses and can activate or deactivate alarms.
       ! after alarm adjustments FISOC_proc does accumulation of fluxes.
       ! FISOC_proc must operate on the ocean grid and degeneration.

!       IF (ESMF_AlarmIsRinging(alarm_ocn, rc=rc)) THEN
!          PRINT*, "Ocean alarm, we must call the ocean component"
!          PRINT*, "After the ocean component we must call the coupler"
!       END IF

!       IF (ESMF_AlarmIsRinging(alarm_ice, rc=rc)) THEN
!          PRINT*, "Ice alarm, we must call the ice component"
!          PRINT*, "After the ice component we must call the coupler, phase 2?"
!          CALL ESMF_AlarmRingerOff(alarm_ice, rc=rc)
!       END IF
       
       CALL ESMF_ClockAdvance(FISOC_clock, rc=rc)   
    END DO

  END SUBROUTINE FISOC_run


!------------------------------------------------------------------------------
  SUBROUTINE FISOC_finalize(FISOC_parent, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc
    
!*** finalize:
!*** finalize child components, including coupler
!*** destroy al states and components
  END SUBROUTINE FISOC_finalize

END MODULE FISOC_parent_mod
