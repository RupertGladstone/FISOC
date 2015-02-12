
MODULE  FISOC_parent_mod
  
  USE ESMF

  USE FISOC_ISM, ONLY : FISOC_ISM_register
  USE FISOC_coupler, ONLY : FISOC_coupler_register
  USE FISOC_OM, ONLY : FISOC_OM_register
  USE FISOC_utils
  
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
    LOGICAL              :: verbose_coupling


    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_parent, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    IF (verbose_coupling) THEN
       PRINT*,""
       PRINT*,"************************************************************************************"
       PRINT*,"This is a verbose run, set by the verbose_coupling flag in the FISOC_config.rc file."
       PRINT*,"(this affects only printing to screen, not writing to the log file(s), which occurs "
       PRINT*,"in any case)"
       PRINT*,"************************************************************************************"
       PRINT*,""
    END IF

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

    msg = "FISOC initialise completed"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! put child states into the parent import state just so we can make them available during run phases
    CALL ESMF_StateAdd(importstate, (/OM_ImpSt, OM_ExpSt, ISM_ImpSt, ISM_ExpSt/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_init
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_run(FISOC_parent, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc
 
    TYPE(ESMF_State)     :: ISM_ImpSt, ISMExpSt, OM_ImpSt, OM_ExpSt
    TYPE(ESMF_Alarm)     :: alarm_OM, alarm_ISM
    TYPE(ESMF_Alarm)     :: alarm_ISM_exportAvailable, alarm_OM_output
    LOGICAL              :: verbose_coupling
    TYPE(ESMF_config)    :: config
    INTEGER              :: urc


    CALL ESMF_GridCompGet(FISOC_parent, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_StateGet(importstate, "ISM import state", ISM_ImpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(importstate, "ISM export state", ISM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(importstate, "OM import state", OM_ImpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(importstate, "OM export state", OM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ISM", alarm_ISM, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_OM", alarm_OM, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "FISOC parent run: got alarms and child states, now start timestepping"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)


    ! main timestepping loop
    DO WHILE (.NOT. ESMF_ClockIsStopTime(FISOC_clock, rc=rc))

       IF (verbose_coupling) THEN
          CALL ESMF_ClockPrint(FISOC_clock, options="advanceCount string isofrac", rc=rc)
          CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_OM_output", alarm_OM_output, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ISM_exportAvailable", alarm_ISM_exportAvailable, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          PRINT*,"Alarm status: ",ESMF_AlarmIsRinging(alarm_OM),&
               ESMF_AlarmIsRinging(alarm_OM_output),&
               ESMF_AlarmIsRinging(alarm_ISM),&
               ESMF_AlarmIsRinging(alarm_ISM_exportAvailable)
       END IF
       
       IF (ESMF_AlarmIsRinging(alarm_OM, rc=rc)) THEN
          msg = "FISOC parent run: calling OM"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_GridCompRun(FISOC_OM, &
               importState=OM_ImpSt, exportState=OM_ExpSt, &
               clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

          CALL ESMF_cplCompRun(FISOC_coupler, &
               importState=OM_ExpSt, exportState=ISM_ImpSt, &
               clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
       
       IF (ESMF_AlarmIsRinging(alarm_ISM, rc=rc)) THEN
          msg = "FISOC parent run: calling ISM"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_GridCompRun(FISOC_ISM, &
               importState=ISM_ImpSt, exportState=ISM_ExpSt, &
               clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

          CALL ESMF_cplCompRun(FISOC_coupler, &
               importState=ISM_ExpSt, exportState=OM_ImpSt, &
               clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!          CALL ESMF_AlarmRingerOff(alarm_OM, rc=rc)
!          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!               line=__LINE__, file=__FILE__)) &
!               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       END IF

!       IF (ESMF_AlarmIsRinging(alarm_ISM, rc=rc)) THEN
!          CALL ESMF_AlarmRingerOff(alarm_ISM, rc=rc)
!          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!               line=__LINE__, file=__FILE__)) &
!               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
!       END IF

       CALL ESMF_ClockAdvance(FISOC_clock, rc=rc)   
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    END DO

    msg = "FISOC parent run: complete"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

  END SUBROUTINE FISOC_run


  !------------------------------------------------------------------------------
  ! get all states, finalise child components, destroy states.
  SUBROUTINE FISOC_finalize(FISOC_parent, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc

    INTEGER              :: urc
    TYPE(ESMF_State)     :: ISM_ImpSt, ISMExpSt, OM_ImpSt, OM_ExpSt
    TYPE(ESMF_config)    :: config
    LOGICAL              :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_ISM, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (verbose_coupling) THEN
       PRINT*,"Now some finalising/destroying stuff, not very exciting"
    END IF

    msg = "FISOC parent finalise: getting child states"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_StateGet(importstate, "ISM import state", ISM_ImpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(importstate, "ISM export state", ISM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(importstate, "OM import state", OM_ImpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(importstate, "OM export state", OM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "FISOC parent finalise: calling child finalise routines"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompFinalize(FISOC_OM, &
         importState=OM_ImpSt, exportState=OM_ExpSt, &
         clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompFinalize(FISOC_ISM, &
         importState=ISM_ImpSt, exportState=ISM_ExpSt, &
         clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_cplCompFinalize(FISOC_coupler, &
         importState=OM_ExpSt, exportState=ISM_ImpSt, &
         clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_StateDestroy(ISM_ImpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateDestroy(ISM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateDestroy(OM_ImpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateDestroy(OM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateDestroy(exportstate, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateDestroy(importstate, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_finalize

END MODULE FISOC_parent_mod
