MODULE FISOC_OM_MOD
  
  USE ESMF
  USE FISOC_utils_MOD
  USE FISOC_OM_Wrapper
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_OM_register
    
CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_register(FISOC_OM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_OM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_OM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_OM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_OM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_OM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_OM_register


  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two stages.  The first stage is for independent 
  ! initialisation of the ISM and the second stage occurrs after the OM has been 
  ! initialised, to allow inter-component consistency checks or use of the OM 
  ! state to complete ISM initalisation.
  SUBROUTINE FISOC_OM_init_phase1(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)        :: FISOC_OM
    TYPE(ESMF_State)           :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)           :: FISOC_clock
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_config)          :: FISOC_config
    TYPE(ESMF_grid)            :: OM_grid
    TYPE(ESMF_fieldBundle)     :: OM_ExpFB,OM_ExpFBcum
    TYPE(ESMF_VM)              :: vm

    CHARACTER(len=ESMF_MAXSTR) :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: OM_ReqVarList(:)

    rc = ESMF_FAILURE

    msg = "OM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! get information from the OM gridded component.  vm is virtual machine.
    CALL ESMF_GridCompGet(FISOC_OM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_cavityCheckOptions(FISOC_config, rc)

    ! create empty field bundle to be populated by model-specific code.
    OM_ExpFB = ESMF_FieldBundleCreate(name='OM export fields', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    OM_ExpFBcum = ESMF_FieldBundleCreate(name='OM export field cumulator', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! model-specific initialisation
    CALL FISOC_OM_Wrapper_Init_Phase1(OM_ExpFB,OM_grid,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_initCumulatorFB(OM_ExpFB,OM_ExpFBcum,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_zeroBundle(OM_ExpFBcum,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! we only add the OM export field bundle to the import state as a way of 
    ! letting the coupler get hold of the grid.  There must be a better way 
    ! to do this.
    CALL ESMF_StateAdd(OM_ImpSt, (/OM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(OM_ExpSt, (/OM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(OM_ExpSt, (/OM_ExpFBcum/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "OM initialise phase 1 complete"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_init_phase1
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_init_phase2(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_VM)          :: vm
    TYPE(ESMF_config)      :: FISOC_config
    TYPE(ESMF_fieldbundle) :: OM_ImpFB, OM_ExpFB

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_OM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(OM_ImpSt, "OM import fields", OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(OM_ExpSt, "OM export fields", OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    CALL FISOC_OM_Wrapper_Init_Phase2(OM_ImpFB,OM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "OM initialise phase 2 (allows the OM access to the ISM initial state) "
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_init_phase2

  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_run(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)

    TYPE(ESMF_GridComp)        :: FISOC_OM
    TYPE(ESMF_State)           :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)           :: FISOC_clock
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_FileStatus_Flag) :: NC_status
    CHARACTER(len=ESMF_MAXSTR) :: OutputFileName 
    TYPE(ESMF_VM)              :: vm
    INTEGER(ESMF_KIND_I8)      :: advanceCount
    INTEGER                    :: localPet,advanceCountInt4
    TYPE(ESMF_fieldbundle)     :: OM_ImpFB, OM_ExpFB, OM_ExpFBcum
    TYPE(ESMF_config)          :: FISOC_config
    TYPE(ESMF_Alarm)           :: alarm_OM_output, alarm_ISM, alarm_ISM_exportAvailable
    LOGICAL                    :: verbose_coupling, OM_writeNetcdf

    rc = ESMF_FAILURE
   
    CALL ESMF_GridCompGet(FISOC_OM, config=FISOC_config, vm=vm, rc=rc)
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

    CALL ESMF_ConfigGetAttribute(FISOC_config,  OM_writeNetcdf, label='OM_writeNetcdf:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
       msg="WARNING: OM_writeNetcdf not found in FISOC_config.rc, not writing to NetCDF"
       OM_writeNetcdf = .FALSE.
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    ! Clocks and alarms are used to control the asynchronous timestepping
    ! between the OM and ISM.  Most of the logic is on the OM side since 
    ! it is assumed (actually it is required) that the ISM timestep is 
    ! an exact multiple of the OM timestep.
    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ISM", alarm_ISM, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ISM_exportAvailable", alarm_ISM_exportAvailable, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_OM_output", alarm_OM_output, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! "AdvanceCount" gives the number of (OM) timesteps.  Use it to make NetCDF filename.
    CALL ESMF_ClockGet(FISOC_clock, advanceCount=advanceCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    ! The ocean cavity might need temporal linear interpolation at this point.
    CALL OM_HandleCavity(FISOC_config, FISOC_clock, OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    ! Decide how to call OM run wrapper depending on relevant alarms
    OM_output: IF (ESMF_AlarmIsRinging(alarm_OM_output, rc=rc)) THEN
       
       CALL ESMF_StateGet(OM_ExpSt, "OM export fields", OM_ExpFB, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ISM_input: IF (ESMF_AlarmIsRinging(alarm_ISM_exportAvailable, rc=rc)) THEN
          CALL ESMF_StateGet(OM_ImpSt, "OM import fields", OM_ImpFB, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
          
          CALL FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB=OM_ExpFB,OM_ImpFB=OM_ImpFB,rc_local=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
          
          writeNCimp: IF (OM_writeNetcdf) THEN
             CALL ESMF_VMBarrier(vm, rc=rc)
             IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) &
                  CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
             msg = "Writing NetCDF output from FISOC on ocean grid (OM import)"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                  line=__LINE__, file=__FILE__, rc=rc)
             WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_OM_imp_t", advanceCount, ".nc"
             CALL FISOC_FB2NC(OutputFileName,OM_ImpFB,FISOC_config)
          END IF writeNCimp

          CALL ESMF_VMBarrier(vm, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
          
       ELSE
          CALL FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB=OM_ExpFB,rc_local=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)              
       END IF ISM_input
       
       writeNC: IF (OM_writeNetcdf) THEN
          CALL ESMF_VMBarrier(vm, rc=rc)
          msg = "Writing NetCDF output from FISOC on ocean grid (OM export)"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_OM_exp_t", advanceCount, ".nc"
          CALL FISOC_FB2NC(OutputFileName,OM_ExpFB,FISOC_config)
       END IF writeNC
       
    ELSE
       
       IF (ESMF_AlarmIsRinging(alarm_ISM_exportAvailable, rc=rc)) THEN
          
          CALL ESMF_StateGet(OM_ImpSt, "OM import fields", OM_ImpFB, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
          
          CALL FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ImpFB=OM_ImpFB,rc_local=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
          
          IF (OM_writeNetcdf) THEN
             CALL ESMF_VMBarrier(vm, rc=rc)
             msg = "Writing NetCDF output from FISOC on ocean grid (OM import)"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                  line=__LINE__, file=__FILE__, rc=rc)
             WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_OM_imp_t", advanceCount, ".nc"
             CALL FISOC_FB2NC(OutputFileName,OM_ImpFB,FISOC_config)
          END IF

       ELSE
          CALL FISOC_OM_Wrapper_Run(FISOC_config,vm,rc_local=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)              
          
       END IF
       
    END IF OM_output
    

    ! cumulate the outputs if we have new outputs from the OM this step
    OM_cumulate: IF (ESMF_AlarmIsRinging(alarm_OM_output)) THEN
       CALL ESMF_StateGet(OM_ExpSt, "OM export field cumulator", OM_ExpFBcum, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       
       CALL FISOC_cumulateFB(OM_ExpFB,OM_ExpFBcum,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

       ! ...and process the cumulated outputs if this is an ISM step
       IF (ESMF_AlarmIsRinging(alarm_ISM)) THEN
          CALL FISOC_processCumulator(OM_ExpFB,OM_ExpFBcum,FISOC_config,rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       END IF
       
    END IF OM_cumulate
    

    ! turn off alarms relating to data exchange and cumulating
    IF (ESMF_AlarmIsRinging(alarm_ISM_exportAvailable, rc=rc)) THEN
       CALL ESMF_AlarmRingerOff(alarm_ISM_exportAvailable, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)              
    END IF    
    IF (ESMF_AlarmIsRinging(alarm_OM_output, rc=rc)) THEN
       CALL ESMF_AlarmRingerOff(alarm_OM_output, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_run
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_finalise(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_VM)                :: vm
    TYPE(ESMF_config)            :: FISOC_config
    TYPE(ESMF_fieldbundle)       :: OM_ImpFB, OM_ExpFB
    INTEGER                      :: FieldCount,ii, localPet
    TYPE(ESMF_field),ALLOCATABLE :: FieldList(:)
    TYPE(ESMF_grid)              :: OM_grid
    
    rc = ESMF_FAILURE


    msg = "OM finalise: destroy fields, bundles and grid"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompGet(FISOC_OM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(OM_ImpSt, "OM import fields", OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(FieldList(FieldCount))
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldCount=FieldCount, fieldList=FieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldGet(FieldList(1), grid=OM_grid, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_GridDestroy(OM_grid,rc=rc)
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

    CALL ESMF_FieldBundleDestroy(OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    

    CALL ESMF_StateGet(OM_ExpSt, "OM export fields", OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(FieldList(FieldCount))
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=FieldCount, fieldList=FieldList,rc=rc)
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

    CALL ESMF_FieldBundleDestroy(OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_OM_Wrapper_Finalize(FISOC_config,localPet,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_finalise  



  !------------------------------------------------------------------------------
  SUBROUTINE OM_HandleCavity(FISOC_config, FISOC_clock, OM_ImpFB, rc)

    TYPE(ESMF_config),INTENT(INOUT)          :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)     :: OM_ImpFB 
    INTEGER,INTENT(OUT),OPTIONAL             :: rc
    TYPE(ESMF_Clock),INTENT(IN)              :: FISOC_clock

    CHARACTER(len=ESMF_MAXSTR)               :: OM_cavityUpdate
    TYPE(ESMF_Alarm)                         :: alarm_ISM_exportAvailable
    INTEGER, SAVE                            :: linterpCounter=0
    INTEGER                                  :: dt_ratio

    REAL(ESMF_KIND_R8)          :: linterpFactor
    TYPE(ESMF_field)            :: ISM_z_l0, ISM_z_l0_previous, ISM_z_l0_linterp 
    REAL(ESMF_KIND_R8),POINTER  :: ptr_curr(:,:),ptr_prev(:,:),ptr_linterp(:,:)

    rc = ESMF_FAILURE

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_cavityUpdate,    & 
         label='OM_cavityUpdate:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, dt_ratio,           & 
         label='dt_ratio:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    SELECT CASE (OM_cavityUpdate)

    CASE('Rate','RecentIce')
       ! handled elsewhere, do nothing

    CASE('CorrectedRate')
       msg = "OM_cavityUpdate NYI (corrRate)"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CASE('Linterp')
       msg = "OM_cavityUpdate NYI (linterp)"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       ! If an ISM eport is available this means the ISM cavity has been 
       ! updated.  Each time this happens we can reset a counter.  The 
       ! counter is used to count steps since the last ISM cavity update. 
       ! We rely on dt_ratio governing the total number of steps until 
       ! the next ISM cavity update. 
       CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ISM_exportAvailable", alarm_ISM_exportAvailable, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldBundleGet(OM_ImpFB, fieldname='ISM_z_l0', field=ISM_z_l0, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldGet(field=ISM_z_l0, farrayPtr=ptr_curr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldBundleGet(OM_ImpFB, fieldname='ISM_z_l0_previous', field=ISM_z_l0_previous, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldGet(field=ISM_z_l0_previous, farrayPtr=ptr_prev, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldBundleGet(OM_ImpFB, fieldname='ISM_z_l0_linterp', field=ISM_z_l0_linterp, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(field=ISM_z_l0_linterp, farrayPtr=ptr_linterp, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       IF (ESMF_AlarmIsRinging(alarm_ISM_exportAvailable, rc=rc)) THEN 
          linterpCounter = 0
          linterpFactor  = 1.0
       ELSE
          linterpCounter = linterpCounter + 1
          linterpFactor  = linterpCounter/dt_ratio
       END IF

       ptr_linterp = ptr_prev*(1-linterpFactor) + ptr_curr*(linterpFactor)

       NULLIFY(ptr_prev)
       NULLIFY(ptr_curr)
       NULLIFY(ptr_linterp)

    CASE DEFAULT
       msg = "OM_cavityUpdate not recognised"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE OM_HandleCavity

END MODULE FISOC_OM_MOD
