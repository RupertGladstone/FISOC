
MODULE  FISOC_parent_mod
  
  USE ESMF
  
  !*** register the components? register all four child components here...
  !  use   InjectorMod, only : Injector_register
  !  use FlowSolverMod, only : FlowSolver_register
  !  use    CouplerMod, only : Coupler_register

  USE FISOC_ISM, ONLY : FISOC_ISM_register
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_parent_register
  
  TYPE(ESMF_GridComp), SAVE :: FISOC_ISM, FISOC_OM, FISOC_proc
  TYPE(ESMF_CplComp),  SAVE :: FISOC_coupler
  TYPE(ESMF_State),    SAVE :: FISOC_ice_import, FISOC_ice_export, &
       FISOC_ocean_import, FISOC_ocean_export, FISOC_proc_import,  &
       FISOC_proc_export ! might remove the proc ones? (time processing can be done within Elmer maybe?)
  
CONTAINS
  
  !------------------------------------------------------------------------------
  ! register subroutines to be called for init, run, finalize for this 
  ! parent component
  SUBROUTINE FISOC_parent_register(FISOC_parent, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_parent
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_init, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_RUN, &
         userRoutine=FISOC_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_finalize, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_parent_register
  
    
!------------------------------------------------------------------------------
  SUBROUTINE FISOC_init(FISOC_parent, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount, localrc, urc
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_SUCCESS

    ! Get petCount from the component (pet is persistent execution thread)
    call ESMF_GridCompGet(FISOC_parent, petCount=petCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return ! bail out

    print *,"petCount ",petCount

    ! Create and register child components and routines
    FISOC_ISM = ESMF_GridCompCreate(name="Ice Sheet Model", rc=localrc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetServices(FISOC_ISM, FISOC_ISM_register, &
         userRc=urc, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


! importState=INimp, exportState=INexp, &
    CALL ESMF_GridCompInitialize(FISOC_ISM, &
         clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    IF (rc /= ESMF_SUCCESS) CALL ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    IF (urc /= ESMF_SUCCESS) CALL ESMF_Finalize(endflag=ESMF_END_ABORT, rc=urc)
    
    PRINT *, "Injection Model Initialize finished, rc = ", rc,"  and urc = ", urc
    

!*** how to set up a mesh

!*** import and export states

!***more child comps

 !*** can access and check stuff to do with parallelization here... PET count...
 !*** can set up stuff to do with connectivity - how grids are distributed over pets
 ! (PET = persistent execution thread)
 !*** can create the child grid components here, and register the relevant subroutines

 !*** create child grids and their decompositions
!*** maybe elmer needs to set the grid... can write an elmer solver to create the grid 
! in the proper ESMF way, using Elmer's own idea of how the grid should be...

!***ROMS grid: need to see how Ufuk has done it...

 !*** create import and export states for all grid components, and call their initialize 
 ! routines (which have previously been registered)
 ! note coupler gets initialised multiple times... seems like one of the child 
 ! components get initialised multiple times???  ah, grid components have multiple 
 ! initialisation stages ... is this set in their register routines?  or somewhere else?

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
       !*** do clock stuff (but not advancing the main clock??? where does that happen???)
       !*** run the child components
       !*** coupler gets called multiple times with different import and export states.  
!       CALL ESMF_TimeIntervalPrint(ts_ice, rc=rc)
!       CALL ESMF_ClockPrint(FISOC_clock, options="advanceCount string isofrac", rc=rc)
!       CALL ESMF_ClockPrint(FISOC_clock, options="currTime string", rc=rc)
       
       PRINT*,""
       PRINT*,"New timestep. We need to call FISOC_proc"
       ! ***the FISOC_proc component accesses and can activate or deactivate alarms.
       ! after alarm adjustments FISOC_proc does accumulation of fluxes.
       ! FISOC_proc must operate on the ocean grid and degeneration.

       IF (ESMF_AlarmIsRinging(alarm_ocn, rc=rc)) THEN
          PRINT*, "Ocean alarm, we must call the ocean component"
          PRINT*, "After the ocean component we must call the coupler"
       END IF

       IF (ESMF_AlarmIsRinging(alarm_ice, rc=rc)) THEN
          PRINT*, "Ice alarm, we must call the ice component"
          PRINT*, "After the ice component we must call the coupler"
          CALL ESMF_AlarmRingerOff(alarm_ice, rc=rc)
       END IF
       
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
