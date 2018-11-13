MODULE FISOC_OM_MOD
  
  USE ESMF
  USE FISOC_utils_MOD
  USE FISOC_OM_Wrapper
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_OM_register

  CHARACTER(len=ESMF_MAXSTR),SAVE :: OM_impFBname = "OM import fields"
  CHARACTER(len=ESMF_MAXSTR),SAVE :: OM_expFBname = "OM export fields"

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
  ! initialisation of the OM and the second stage occurrs after the ISM has been 
  ! initialised, to allow inter-component consistency checks or use of the ISM 
  ! state to complete OM initalisation.
  SUBROUTINE FISOC_OM_init_phase1(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)        :: FISOC_OM
    TYPE(ESMF_State)           :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)           :: FISOC_clock
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_config)          :: FISOC_config
# if defined(FISOC_OM_GRID)
    TYPE(ESMF_grid)            :: OM_grid
# elif defined(FISOC_OM_MESH)
    TYPE(ESMF_mesh)            :: OM_mesh
# endif
    TYPE(ESMF_fieldBundle)     :: OM_ExpFB,OM_ExpFBcum
    TYPE(ESMF_VM)              :: vm
    LOGICAL                    :: ISM_UseOMGrid, OM_UseISMGrid
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
    
    ! if the ISM and OM are on the same grid then the ISM export field bundle is 
    ! used as the OM import field bundle
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_UseOMGrid, 'ISM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_UseISMGrid, 'OM_UseISMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ISM_UseOMGrid.OR.OM_UseISMGrid) THEN
      OM_impFBname = "ISM export fields"
    END IF

    CALL FISOC_cavityCheckOptions(FISOC_config, rc)

    ! create empty field bundle to be populated by model-specific code.
    OM_ExpFB = ESMF_FieldBundleCreate(name=OM_expFBname, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    OM_ExpFBcum = ESMF_FieldBundleCreate(name='OM export field cumulator', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! model-specific initialisation
# if defined(FISOC_OM_GRID)
    CALL FISOC_OM_Wrapper_Init_Phase1(FISOC_config,vm,OM_ExpFB,OM_grid,rc=rc)
# elif defined(FISOC_OM_MESH)
    CALL FISOC_OM_Wrapper_Init_Phase1(FISOC_config,vm,OM_ExpFB,OM_mesh,rc=rc)
# else
    msg = "ERROR: FISOC does not recognise OM geom type."
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
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

    CALL ESMF_StateGet(OM_ImpSt, OM_impFBname, OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(OM_ExpSt, OM_expFBname, OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    CALL FISOC_OM_Wrapper_Init_Phase2(FISOC_config,vm,OM_ImpFB,OM_ExpFB,rc=rc)
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
    CHARACTER(len=ESMF_MAXSTR) :: OM_NCfreq
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

    IF (OM_writeNetcdf) THEN      
       CALL ESMF_ConfigGetAttribute(FISOC_config,  OM_NCfreq, label='OM_NCfreq:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) THEN
          msg="WARNING: OM_NCfreq not found in FISOC_config.rc, setting to all"
          OM_NCfreq = "all"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__, rc=rc)
       END IF
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

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(OM_ExpSt, OM_expFBname, OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_StateGet(OM_ImpSt, OM_impFBname, OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    ! The ocean cavity might need temporal linear interpolation at this point.
    CALL OM_HandleCavity(FISOC_config, FISOC_clock, OM_ImpFB, OM_ExpFB, localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Decide how to call OM run wrapper depending on relevant alarms
    OM_output: IF (ESMF_AlarmIsRinging(alarm_OM_output, rc=rc)) THEN
       
       ! we have new ISM output available and we need OM output
       ISM_input1: IF (ESMF_AlarmIsRinging(alarm_ISM_exportAvailable, rc=rc)) THEN          
          CALL FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB=OM_ExpFB,OM_ImpFB=OM_ImpFB,rc_local=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) THEN
             CALL FISOC_OM_finalise(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          END IF
          CALL OM_NetcdfWrapper(FISOC_config,OM_writeNetcdf,OM_NCfreq,advanceCount, &
               OM_ExpFB=OM_ExpFB,OM_ImpFB=OM_ImpFB)

       ! no new ISM output is available.  we need OM output
       ELSE
          CALL FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB=OM_ExpFB,rc_local=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) THEN
             CALL FISOC_OM_finalise(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          END IF
          CALL OM_NetcdfWrapper(FISOC_config,OM_writeNetcdf,OM_NCfreq,advanceCount,OM_ExpFB=OM_ExpFB)
       END IF ISM_input1
       
    ELSE       
       ! we have new ISM output available but we do not need OM output
       ISM_input2: IF (ESMF_AlarmIsRinging(alarm_ISM_exportAvailable, rc=rc)) THEN
          CALL FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ImpFB=OM_ImpFB,rc_local=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) THEN
             CALL FISOC_OM_finalise(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          END IF
          CALL OM_NetcdfWrapper(FISOC_config,OM_writeNetcdf,OM_NCfreq,advanceCount,OM_ImpFB=OM_ImpFB)

       ! no new ISM output is available and we do not need OM output
       ELSE
          CALL FISOC_OM_Wrapper_Run(FISOC_config,vm,rc_local=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) THEN
             CALL FISOC_OM_finalise(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          END IF
          
       END IF ISM_input2
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
# if defined(FISOC_OM_GRID)
    TYPE(ESMF_grid)              :: OM_grid
# elif defined(FISOC_OM_MESH)
    TYPE(ESMF_grid)              :: OM_mesh
# endif
    LOGICAL                      :: ISM_UseOMGrid, OM_UseISMGrid
    
    rc = ESMF_FAILURE


    msg = "OM finalise: destroy fields, bundles and grid"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompGet(FISOC_OM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_UseOMGrid, 'ISM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_UseISMGrid, 'OM_UseISMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ( (.NOT.ISM_UseOMGrid).AND.(.NOT.OM_UseISMGrid) ) THEN
      CALL ESMF_StateGet(OM_ImpSt, OM_impFBname, OM_ImpFB, rc=rc)
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
      
# if defined(FISOC_OM_GRID)
      CALL ESMF_FieldGet(FieldList(1), grid=OM_grid, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      CALL ESMF_GridDestroy(OM_grid,rc=rc)
# elif defined(FISOC_OM_MESH)
      CALL ESMF_FieldGet(FieldList(1), mesh=OM_mesh, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      CALL ESMF_MeshDestroy(OM_mesh,rc=rc)
# endif
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
    END IF

    CALL ESMF_StateGet(OM_ExpSt, OM_expFBname, OM_ExpFB, rc=rc)
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
    
    CALL FISOC_OM_Wrapper_Finalize(FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_finalise  


  !------------------------------------------------------------------------------
  SUBROUTINE OM_HandleCavity(FISOC_config, FISOC_clock, OM_ImpFB, OM_ExpFB, localPet, rc)
    
!    USE mod_param, ONLY       : BOUNDS, Ngrids
    
    TYPE(ESMF_config),INTENT(INOUT)          :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)     :: OM_ImpFB, OM_ExpFB
    TYPE(ESMF_Clock),INTENT(IN)              :: FISOC_clock
    INTEGER,INTENT(IN)                       :: localPet
    INTEGER,INTENT(OUT),OPTIONAL             :: rc

    CHARACTER(len=ESMF_MAXSTR)               :: OM_cavityUpdate
    TYPE(ESMF_Alarm)                         :: alarm_ISM_exportAvailable
    INTEGER, SAVE                            :: linterpCounter=0
    INTEGER                                  :: dt_ratio, ISM_dt_int

    REAL(ESMF_KIND_R8)          :: linterpFactor, ISM_dt, OM_WCmin, CavCorr
    TYPE(ESMF_field)            :: ISM_z_l0_previous, ISM_z_l0_linterp 
    TYPE(ESMF_field)            :: ISM_z_l0_field, OM_z_l0_field, dddt_field
    TYPE(ESMF_field)            :: OM_bed_field
    REAL(ESMF_KIND_R8),POINTER  :: ptr_curr(:,:),ptr_prev(:,:),ptr_linterp(:,:)
    REAL(ESMF_KIND_R8),POINTER  :: ISM_z_l0(:,:), OM_z_l0(:,:), dddt(:,:)
    REAL(ESMF_KIND_R8),POINTER  :: OM_bed(:,:)
    INTEGER                     :: ii,jj,arrShape(2)
    INTEGER                     :: IendR, IstrR, JendR, JstrR

    REAL(ESMF_KIND_R8), POINTER :: ptrX(:,:), ptrY(:,:)
# if defined(FISOC_OM_GRID)
    TYPE(ESMF_grid)             :: OM_grid
# elif defined(FISOC_OM_MESH)
    TYPE(ESMF_mesh)             :: OM_mesh
# endif
    INTEGER                     :: exclusiveLBound(2),exclusiveUBound(2),exclusiveCount(2),computationalLBound(2)
    INTEGER                     :: computationalUBound(2),computationalCount(2),totalLBound(2),totalUBound(2)

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
      
       !       msg = "OM_cavityUpdate NYI: "//OM_cavityUpdate
      !       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
      !            line=__LINE__, file=__FILE__, rc=rc)
      !       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
      ! ISM var dddt will be used.  If an ISM export is available, which means 
       ! dddt has just been calculated, we modify dddt here to impose a 
      ! correcting drift. The drift is designed to reduce the ISM-OM 
       ! discrepancy over one ISM timestep.
      CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ISM_exportAvailable", alarm_ISM_exportAvailable, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      IF (ESMF_AlarmIsRinging(alarm_ISM_exportAvailable, rc=rc)) THEN 
        
        ! We need the OM cavity geom... 
        CALL ESMF_FieldBundleGet(OM_ExpFB, fieldName="OM_z_l0", field=OM_z_l0_field, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        CALL ESMF_FieldGet(field=OM_z_l0_field, farrayPtr=OM_z_l0, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
        ! ... the OM bedrock... 
        CALL ESMF_FieldBundleGet(OM_ExpFB, fieldName="OM_bed", field=OM_bed_field, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        CALL ESMF_FieldGet(field=OM_bed_field, farrayPtr=OM_bed, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
        ! ...and the ISM cavity geom.
        CALL ESMF_FieldBundleGet(OM_ImpFB, fieldName="ISM_z_l0", field=ISM_z_l0_field, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        CALL ESMF_FieldGet(field=ISM_z_l0_field, farrayPtr=ISM_z_l0, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
        ! And we need dddt, the cavity rate, in order to update it with the 
        ! drift correction.
        CALL ESMF_FieldBundleGet(OM_ImpFB, fieldName="ISM_dddt", field=dddt_field, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        CALL ESMF_FieldGet(field=dddt_field, farrayPtr=dddt, rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
        ! ISM time step in seconds is used (we are modifying dddt here in 
        ! metres per second)
        CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_int, 'ISM_dt_sec',rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        ISM_dt = REAL(ISM_dt_int,ESMF_KIND_R8)
        
        CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_WCmin, 'OM_WCmin',rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
        CALL FISOC_ConfigDerivedAttribute(FISOC_config, CavCorr, 'OM_CavCorr',rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
        IF ( (SIZE(dddt).NE.SIZE(ISM_z_l0)) .OR.     & 
             (SIZE(dddt).NE.SIZE(OM_z_l0))  .OR.     & 
             (SIZE(dddt).NE.SIZE(OM_bed)) ) THEN
          msg = "ERROR: array size inconsistency in corrected cavity rate calc."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          
        END IF
        
        ! TODO: get and check stagger locs
        ! TODO: move grid/mesh choice to cpp
        ! get the grid or mesh from the OM exp bundle.  
# if defined(FISOC_OM_GRID)
        CALL FISOC_getGridFromFB(OM_expFB,OM_grid,rc=rc)
# elif defined(FISOC_OM_MESH)
        CALL FISOC_getMeshFromFB(OM_expFB,OM_mesh,rc=rc)
# else
        msg="invalid CPP options for OM geom type"
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__, rc=rc)          
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

# if defined(FISOC_OM_GRID)
        CALL ESMF_GridGetCoord(OM_grid,                 &
             coordDim=1,                                &
             exclusiveLBound=exclusiveLBound,           &
             exclusiveUBound=exclusiveUBound,           &
             computationalLBound=computationalLBound,   & 
             computationalUBound=computationalUBound,   & 
             totalLBound=totalLBound,                   & 
             totalUBound=totalUBound,                   & 
             farrayPtr=ptrX,                            &
             rc=rc)
        
        !!          arrShape = SHAPE(dddt)
        !          IstrR=BOUNDS(Ngrids)%IstrR(localPet)
        !          IendR=BOUNDS(Ngrids)%IendR(localPet)
        !          JstrR=BOUNDS(Ngrids)%JstrR(localPet)
        !          JendR=BOUNDS(Ngrids)%JendR(localPet)
        
        DO jj=computationalLBound(2),computationalUBound(2)
          DO ii=computationalLBound(1),computationalUBound(1)
            dddt(ii,jj) = dddt(ii,jj) +                          &
                 CavCorr*( MAX(ISM_z_l0(ii,jj),OM_bed(ii,jj)+    &
                 OM_WCmin)-OM_z_l0(ii,jj) ) / ISM_dt
          END DO
        END DO
        
        !          DO jj=JstrR, JendR
        !            DO ii=IstrR, IendR
        !              dddt(ii,jj) = dddt(ii,jj) +                          &
        !                   CavCorr*( MAX(ISM_z_l0(ii,jj),OM_bed(ii,jj)+    &
        !                   OM_WCmin)-OM_z_l0(ii,jj) ) / ISM_dt
        !            END DO
        !          END DO
# else
        msg="invalid OM geom type for correctedRate in handleCavity"
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__, rc=rc)          
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
        
        NULLIFY(ISM_z_l0)
        NULLIFY(OM_z_l0)
        NULLIFY(OM_bed)
        NULLIFY(dddt)
        NULLIFY(ptrX)
        
      END IF
      
      
    CASE('Linterp')
 !      msg = "OM_cavityUpdate NYI: "//OM_cavityUpdate
 !      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
 !           line=__LINE__, file=__FILE__, rc=rc)
 !      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       ! If an ISM export is available this means the ISM cavity has been 
       ! updated.  Each time this happens we can reset a counter.  The 
       ! counter is used to count steps since the last ISM cavity update. 
       ! We rely on dt_ratio governing the total number of steps until 
       ! the next ISM cavity update. 
       ! Note: we assume the ISM initialised both curr and previous cavity 
       ! geom to the initial cavity geom, otherwise we would pass invalid 
       ! cavity to OM
       CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_ISM_exportAvailable", alarm_ISM_exportAvailable, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldBundleGet(OM_ImpFB, fieldname='ISM_z_l0', field=ISM_z_l0_field, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldGet(field=ISM_z_l0_field, farrayPtr=ptr_curr, rc=rc)
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
          linterpFactor  = 0.0
       ELSE
          linterpCounter = linterpCounter + 1
          linterpFactor  = REAL(linterpCounter,ESMF_KIND_R8)/REAL(dt_ratio,ESMF_KIND_R8)
       END IF
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ptr_linterp = ptr_prev*(1.0-linterpFactor) + ptr_curr*(linterpFactor)

       NULLIFY(ptr_prev)
       NULLIFY(ptr_curr)
       NULLIFY(ptr_linterp)

       ! now turn on ISM export alarm, because the ISM may not have turned 
       ! it on, but the time interpolated cavity is effectively a new ISM 
       ! export, at least from the OM perspective
       CALL ESMF_AlarmRingerOn(alarm_ISM_exportAvailable, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)              
       
    CASE DEFAULT
       msg = "OM_cavityUpdate not recognised"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE OM_HandleCavity


  !------------------------------------------------------------------------------
  SUBROUTINE OM_NetcdfWrapper(FISOC_config,OM_writeNetcdf,OM_NCfreq,advanceCount,OM_ExpFB,OM_ImpFB)

    TYPE(ESMF_config),INTENT(INOUT)               :: FISOC_config
    LOGICAL,INTENT(IN)                            :: OM_writeNetcdf
    CHARACTER(len=ESMF_MAXSTR),INTENT(IN)         :: OM_NCfreq
    INTEGER(ESMF_KIND_I8),INTENT(IN)              :: advanceCount
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL :: OM_ExpFB, OM_ImpFB

    CHARACTER(len=ESMF_MAXSTR) :: OutputFileName
    INTEGER :: rc

    IF(OM_writeNetcdf) THEN
       
       SELECT CASE (OM_NCfreq)

       ! Writing out at all timesteps: write any field bundles we're given
       CASE("all")          
          IF (PRESENT(OM_ImpFB)) THEN
             msg = "Writing NetCDF output from FISOC on ocean grid (OM import)"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                  line=__LINE__, file=__FILE__, rc=rc)
             WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_OM_imp_t", advanceCount, ".nc"
             CALL FISOC_FB2NC(OutputFileName,OM_ImpFB,FISOC_config)
          END IF
          IF (PRESENT(OM_ExpFB)) THEN
             msg = "Writing NetCDF output from FISOC on ocean grid (OM export)"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                  line=__LINE__, file=__FILE__, rc=rc)
             WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_OM_exp_t", advanceCount, ".nc"
             CALL FISOC_FB2NC(OutputFileName,OM_ExpFB,FISOC_config)
          END IF

       ! Writing out at ISM export timesteps: hopefully we're given both field bundles
       CASE("ISM")
          IF (PRESENT(OM_ImpFB)) THEN
             IF (PRESENT(OM_ExpFB)) THEN
                msg = "Writing NetCDF output from FISOC on ocean grid (OM imp and exp)"
                CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                     line=__LINE__, file=__FILE__, rc=rc)
                WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_OM_exp_t", advanceCount, ".nc"
                CALL FISOC_FB2NC(OutputFileName,OM_ExpFB,FISOC_config)
                WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_OM_imp_t", advanceCount, ".nc"
                CALL FISOC_FB2NC(OutputFileName,OM_ImpFB,FISOC_config)
          ELSE
                msg = "OM_NCfreq=ISM but we have only OM_ImpFB and not OM_ExpFB... :("
                CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
                     line=__LINE__, file=__FILE__, rc=rc)
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
             END IF
          END IF

       CASE DEFAULT
          msg = "OM_NCfreq value not recognised: "//OM_NCfreq
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END SELECT
       
    END IF
    
  END SUBROUTINE OM_NetcdfWrapper

END MODULE FISOC_OM_MOD
