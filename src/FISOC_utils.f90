
MODULE FISOC_utils_MOD

  USE ESMF
  USE FISOC_types_MOD

  IMPLICIT NONE

  PRIVATE

  PUBLIC   FISOC_getListFromConfig, FISOC_populateFieldBundle, &
       FISOC_ConfigDerivedAttribute, FISOC_initCumulatorFB,    &
       FISOC_zeroBundle, FISOC_cumulateFB,                     & 
       FISOC_processCumulator, msg, FISOC_VM_MPI_Comm_dup,     &
       FISOC_FieldRegridStore, FISOC_FB2NC, FISOC_setClocks,   & 
       FISOC_destroyClocks, FISOC_ISM2OM, FISOC_OM2ISM,        &
       FISOC_shrink, FISOC_VMAllGather, Unique1DArray,         &
       FISOC_OneGrid, FISOC_cavityCheckOptions,                & 
       FISOC_getFirstFieldRank, FISOC_makeRHfromFB,            &
       FISOC_State2StateCopyFB, FISOC_regridFB,                &
       FISOC_GridCompRun, FISOC_FieldListGetField,             &
       FISOC_getGridFromFB, FISOC_getMeshFromFB,               &                
       FISOC_ConfigStringListContains, FISOC_locallyOwnedNodes,&
       FISOC_CreateOneToManyRouteHandle,                       & 
       FISOC_ArrayRedistFromField, FISOC_MAPLL,                &
       FISOC_IsDerived, FISOC_MPI_NewComm
!         FISOC_getGridOrMeshFromFB, 

  INTERFACE FISOC_IsDerived
     MODULE PROCEDURE FISOC_IsDerived_Real
  END INTERFACE FISOC_IsDerived

  INTERFACE Unique1DArray
     MODULE PROCEDURE Unique1DArray_I4
     MODULE PROCEDURE Unique1DArray_D
  END INTERFACE Unique1DArray

  INTERFACE FISOC_VMAllGather 
     MODULE PROCEDURE FISOC_VMAllGather_Int   
     MODULE PROCEDURE FISOC_VMAllGather_Real
  END INTERFACE FISOC_VMAllGather
  
  INTERFACE FISOC_shrink
     MODULE PROCEDURE FISOC_shrinkInt
     MODULE PROCEDURE FISOC_shrinkReal
  END INTERFACE FISOC_shrink

  INTERFACE FISOC_getListFromConfig
     MODULE PROCEDURE FISOC_getRealListFromConfig 
     MODULE PROCEDURE FISOC_getIntegerListFromConfig 
     MODULE PROCEDURE FISOC_getStringListFromConfig
  END INTERFACE FISOC_getListFromConfig
  
  INTERFACE FISOC_populateFieldBundle
     MODULE PROCEDURE FISOC_populateFieldBundleOn2dGrid
     MODULE PROCEDURE FISOC_populateFieldBundleOnMesh
  END INTERFACE FISOC_populateFieldBundle
  
  INTERFACE FISOC_ConfigDerivedAttribute
     MODULE PROCEDURE FISOC_ConfigDerivedAttributeInteger
     MODULE PROCEDURE FISOC_ConfigDerivedAttributeStaggerLocArray
     MODULE PROCEDURE FISOC_ConfigDerivedAttributeLogical
     MODULE PROCEDURE FISOC_ConfigDerivedAttributeReal
     MODULE PROCEDURE FISOC_ConfigDerivedAttributeString
     MODULE PROCEDURE FISOC_ConfigDerivedAttributeRegridMethod
     MODULE PROCEDURE FISOC_ConfigDerivedAttributeExtrapMethod
  END INTERFACE

  CHARACTER(len=ESMF_MAXSTR) :: msg
  
CONTAINS
  


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_MPI_NewComm(MPI_COMM_in, SizeNew, MPI_COMM_new)

    USE mpi

    INTEGER,INTENT(IN)  :: MPI_COMM_in, SizeNew
    INTEGER,INTENT(OUT) :: MPI_COMM_new

    INTEGER             :: size, rank, ierror, color, rc
    
    CALL MPI_COMM_SIZE(MPI_COMM_in, size, ierror)
    IF (ierror.NE.0) THEN 
       msg = "FATAL: error returned from MPI_COMM_SIZE()"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
       
    CALL MPI_COMM_RANK(MPI_COMM_in, rank, ierror)
    IF (ierror.NE.0) THEN 
       msg = "FATAL: error returned from MPI_COMM_RANK()"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    IF (SizeNew.GT.size) THEN
       msg = "FATAL: SizeNew GT size"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    IF (rank.LT.SizeNew) THEN
       color = 1
    ELSE
       color = MPI_UNDEFINED
    END IF

    CALL MPI_Comm_split(MPI_COMM_in, color, 1, MPI_COMM_new, ierror)
    IF (ierror.NE.0) THEN 
       msg = "FATAL: error returned from MPI_COMM_SPLIT()"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
      
  END SUBROUTINE FISOC_MPI_NewComm

  
  !------------------------------------------------------------------------------
  ! 
  ! wrapper for the ESMF method allows simple time profiling.
  ! Wraps the failfast error checks too.
  SUBROUTINE FISOC_GridCompRun(FISOC_gc, ImpSt, ExpSt, profiling, &
               vm, FISOC_clock, phase, cumulatedCPUtime)
    TYPE(ESMF_GridComp),INTENT(INOUT)       :: FISOC_gc
    TYPE(ESMF_State),INTENT(INOUT)          :: ImpSt, ExpSt
    LOGICAL,INTENT(IN)                      :: profiling
    TYPE(ESMF_vm),INTENT(IN)                :: vm
    TYPE(ESMF_clock),INTENT(INOUT),OPTIONAL :: FISOC_clock
    INTEGER,INTENT(IN),OPTIONAL             :: phase
    REAL,INTENT(INOUT),OPTIONAL             :: cumulatedCPUtime

    INTEGER                                 :: rc, urc
    REAL                                    :: preCallTime, postCallTime

    IF (profiling) THEN
       IF (.NOT. (PRESENT(cumulatedCPUtime))) THEN
          msg = "ERROR: profiling set to true but time not present"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    END IF
    
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    IF (profiling) CALL cpu_time(preCallTime)

    CALL ESMF_GridCompRun(FISOC_gc, &
         importState=ImpSt, exportState=ExpSt, &
         clock=FISOC_clock, phase=phase, rc=rc, userRc=urc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    IF (profiling) THEN
       CALL cpu_time(postCallTime)
       cumulatedCPUtime = cumulatedCPUtime + (postCallTime - preCallTime)
    END IF

  END SUBROUTINE FISOC_GridCompRun


  !------------------------------------------------------------------------------
  ! check whether a field is 'required' or 'derived'
  LOGICAL FUNCTION FISOC_IsDerived_Real(FieldName,FISOC_config,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: FieldName
    INTEGER,INTENT(OUT)                   :: rc
    
    CHARACTER(len=ESMF_MAXSTR)            :: ListNameReq,ListNameDer

    rc = ESMF_failure

    ListNameReq = 'FISOC_ISM_ReqVars:'
    ListNameDer = 'FISOC_ISM_DerVars:'
    
    IF (FISOC_ConfigStringListContains(FISOC_config,FieldName,ListNameReq,rc)) THEN
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       FISOC_IsDerived_Real = .FALSE.
       rc = ESMF_success
       
    ELSE IF (FISOC_ConfigStringListContains(FISOC_config,FieldName,ListNameDer,rc)) THEN
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       FISOC_IsDerived_Real = .TRUE.
       rc = ESMF_success

    ELSE
       msg= "FATAL: FieldName "//FieldName//" not found in ISM_reqVars or ISM_derVars."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END IF

    RETURN
    
  END FUNCTION FISOC_IsDerived_Real


  !------------------------------------------------------------------------------
  ! get a field from a field list by name
  TYPE(ESMF_field) FUNCTION FISOC_FieldListGetField(FieldList,FieldName,rc)

    TYPE(ESMF_field),DIMENSION(:) :: FieldList
    CHARACTER(len=*)              :: FieldName
    INTEGER,INTENT(OUT)           :: rc

    CHARACTER(len=ESMF_MAXSTR)    :: FieldName_nn
    INTEGER                       :: nn, FieldCount

    rc = ESMF_FAILURE

    FieldCount = SIZE(FieldList)


    DO nn = 1,FieldCount
      CALL ESMF_FieldGet(fieldList(nn), name=FieldName_nn, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      IF (FieldName.EQ.FieldName_nn) THEN
        FISOC_FieldListGetField = FieldList(nn)
        rc = ESMF_SUCCESS
        RETURN
      END IF
    END DO

  END FUNCTION FISOC_FieldListGetField


  !------------------------------------------------------------------------------
  ! 
  ! check whether too many cavity geometry options are being passed to the OM
  SUBROUTINE FISOC_cavityCheckOptions(FISOC_config,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT)                   :: rc

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: ISM2OM_Vars(:)
    CHARACTER(len=ESMF_MAXSTR)            :: label, OM_cavityUpdate
    INTEGER                               :: ii, count
    LOGICAL                               :: cavityUpdateMismatch

    rc = ESMF_FAILURE

    cavityUpdateMismatch = .FALSE.

    label = 'ISM2OM_vars:'
    CALL FISOC_getStringListFromConfig(FISOC_config, label, ISM2OM_Vars,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_cavityUpdate, label='OM_cavityUpdate:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    count = 0
    DO ii=1,SIZE(ISM2OM_Vars)
       
       IF (TRIM('ISM_z_l0').EQ.TRIM(ISM2OM_Vars(ii))) THEN
          IF (OM_cavityUpdate.NE.'RecentIce') cavityUpdateMismatch = .TRUE.
          count = count + 1
       END IF
       
       IF (TRIM('ISM_z_l0_linterp').EQ.TRIM(ISM2OM_Vars(ii))) THEN
          IF (OM_cavityUpdate.NE.'Linterp') cavityUpdateMismatch = .TRUE.
          count = count + 1
       END IF
       
       IF (TRIM('ISM_dddt').EQ.TRIM(ISM2OM_Vars(ii))) THEN
          IF ( (OM_cavityUpdate.NE.'Rate')          &
               .AND.                                &
               (OM_cavityUpdate.NE.'CorrectedRate') &
               ) THEN
             cavityUpdateMismatch = .TRUE.
          END IF
          count = count + 1
       END IF
       
    END DO
    
    IF (count.eq.0) THEN
       msg = 'no ISM cavity variable will be passed to the OM'
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)          
    END IF

    IF (count.gt.1) THEN
       msg = 'only 1 ISM cavity variable should be passed to the OM'
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    IF (cavityUpdateMismatch) THEN
       msg = 'Cavity update mismatch (check ISM2OM vars) for OM_cavityUpdate: '//OM_cavityUpdate
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_cavityCheckOptions


  !------------------------------------------------------------------------------
  ! Copy a field bundle from one state to another state
  SUBROUTINE FISOC_State2StateCopyFB(srcState,destState,FBname,rc)
    TYPE(ESMF_state),INTENT(IN)    :: srcState
    TYPE(ESMF_state),INTENT(INOUT) :: destState
    CHARACTER(len=*),INTENT(IN)    :: FBname
    INTEGER,INTENT(OUT),OPTIONAL   :: rc

    TYPE(ESMF_fieldBundle)         :: fb

    rc = ESMF_FAILURE

    CALL ESMF_StateGet(srcState, FBname, fb, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateAdd(destState, (/fb/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_State2StateCopyFB


  !------------------------------------------------------------------------------
  ! 
  ! check whether user wants this variable (fieldname) to be passed from the ISM 
  ! to the OM.
  LOGICAL FUNCTION FISOC_ISM2OM(fieldName,FISOC_config,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: fieldName
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: ISM2OM_Vars(:)
    CHARACTER(len=ESMF_MAXSTR)            :: label
    INTEGER                               :: ii

    rc = ESMF_FAILURE

    label = 'ISM2OM_vars:'
    CALL FISOC_getStringListFromConfig(FISOC_config, label, ISM2OM_Vars,rc=rc)
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
       msg = "ISM2OM_vars not found in FISOC config file, trying to pass all available vars"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)          
       FISOC_ISM2OM = .TRUE. ! pass all vars if list is not present 
       RETURN
    ELSE IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    FISOC_ISM2OM = FISOC_listContains(fieldName,ISM2OM_Vars,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

    RETURN

  END FUNCTION  FISOC_ISM2OM


  !------------------------------------------------------------------------------
  ! 
  ! check whether user wants this variable (fieldname) to be passed from the OM 
  ! to the ISM.
  ! ***TODO:resolve code duplication between this and ism2omroutine
  LOGICAL FUNCTION FISOC_OM2ISM(fieldName,FISOC_config,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: fieldName
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: OM2ISM_Vars(:)
    CHARACTER(len=ESMF_MAXSTR)            :: label
    INTEGER                               :: ii

    rc = ESMF_FAILURE

    label = 'OM2ISM_vars:'
    CALL FISOC_getStringListFromConfig(FISOC_config, label, OM2ISM_Vars,rc=rc)
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
       msg = "OM2ISM_vars not found in FISOC config file, trying to pass all available vars"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)          
       FISOC_OM2ISM = .TRUE. ! pass all vars if list is not present 
       RETURN
    ELSE IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    FISOC_OM2ISM = FISOC_listContains(fieldName,OM2ISM_Vars,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

    RETURN

  END FUNCTION  FISOC_OM2ISM
  
  
  !------------------------------------------------------------------------------
  LOGICAL FUNCTION FISOC_ConfigStringListContains(FISOC_config,itemName,listName,rc)
    
    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: itemName, listName
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: stringList(:)
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    
    rc = ESMF_FAILURE

    CALL FISOC_getStringListFromConfig(FISOC_config, listName, stringList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    FISOC_ConfigStringListContains = FISOC_listContains(itemName,stringList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (ALLOCATED(stringList)) DEALLOCATE(stringList)

    rc = ESMF_SUCCESS

    RETURN

  END FUNCTION FISOC_ConfigStringListContains
  

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION FISOC_listContains(itemName,list,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: itemName, list(:)
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    INTEGER                               :: ii

    rc = ESMF_FAILURE

    IF (SIZE(list).EQ.0) THEN
       FISOC_listContains = .FALSE. ! return false if list is empty 
    ELSE
       FISOC_listContains = .FALSE.
       DO ii=1,SIZE(list)
          IF (TRIM(itemName).EQ.TRIM(list(ii))) THEN 
             FISOC_listContains = .TRUE.
          END IF
       END DO
    END IF

    rc = ESMF_SUCCESS

    RETURN

  END FUNCTION  FISOC_LISTCONTAINS



  !--------------------------------------------------------------------------------------
  ! use the ESMF VM to access the mpi communicator and return a duplicate
  SUBROUTINE FISOC_VM_MPI_Comm_dup(vm,mpic_dup,rc)

    TYPE(ESMF_VM),INTENT(IN)       :: VM
    INTEGER, INTENT(OUT)           :: mpic_dup
    INTEGER, OPTIONAL, INTENT(OUT) :: rc

    INTEGER                        :: mpic, ierr

    rc = ESMF_FAILURE

    !-------------------------------------------------------------------------------
    ! Get the parallel context, specifically the mpi communicator, for the OM to 
    ! use.
    ! The returned MPI communicator spans the same MPI processes that the VM
    ! is defined on.
    CALL ESMF_VMGet(vm, mpiCommunicator=mpic, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ! Duplicate the MPI communicator not to interfere with ESMF communications.
    ! The duplicate MPI communicator can be used in any MPI call in the user
    ! code. 
    ! The ifdef statements are needed here because the MPI_comm_dup call will cause
    ! compilation to fail if the MPI library is not available.
#ifdef FISOC_MPI
    CALL MPI_Comm_dup(mpic, mpic_dup, ierr)
    IF (ierr.NE.0) THEN
       msg = "ERROR: Failed call to MPI_Comm_dup"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
#else
    mpic_dup = FISOC_mpic_missing
#endif

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_VM_MPI_Comm_dup



  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_destroyClocks(FISOC_clock,rc)
    
    TYPE(ESMF_Clock),INTENT(INOUT)     :: FISOC_clock
    INTEGER,OPTIONAL,INTENT(OUT)       :: rc
    
!TODO: get all alarms for this clock and destroy them first
!
!    CALL ESMF_AlarmDestroy(alarm_ISM, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
!    CALL ESMF_AlarmDestroy(alarm_OM, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_ClockDestroy(FISOC_clock, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_destroyClocks



  !------------------------------------------------------------------------------
  ! setting up clocks and alarms
  !
  SUBROUTINE FISOC_setClocks(FISOC_config, FISOC_clock, rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    TYPE(ESMF_Clock),INTENT(INOUT)     :: FISOC_clock
    INTEGER,OPTIONAL,INTENT(OUT)       :: rc

    INTEGER                 :: ISM_dt_sec, OM_dt_sec, dt_ratio
    INTEGER                 :: start_year, end_year, start_month, end_month
    INTEGER                 :: OM_outputInterval, runLength_ISM_steps, runLength_secs
    TYPE(ESMF_TimeInterval) :: ISM_dt, OM_dt, runLength_timeInterval
    TYPE(ESMF_Time)         :: startTime, endTime
    TYPE(ESMF_Alarm)        :: alarm_OM, alarm_OM_output, alarm_ISM, alarm_ISM_exportAvailable
    LOGICAL                 :: gotRunLength

    CALL ESMF_ConfigGetAttribute(FISOC_config, dt_ratio, label='dt_ratio:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_dt_sec, label='OM_dt_sec:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_sec, 'ISM_dt_sec',rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_ConfigGetAttribute(FISOC_config, start_month, label='start_month:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_ConfigGetAttribute(FISOC_config, start_year, label='start_year:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_ConfigGetAttribute(FISOC_config, runLength_ISM_steps, label='runLength_ISM_steps:', rc=rc)
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
       gotRunLength = .FALSE.
    ELSE IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) THEN
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ELSE
       gotRunLength = .TRUE.
    END IF

    IF (.NOT.gotRunLength) THEN
       CALL ESMF_ConfigGetAttribute(FISOC_config, end_month, label='end_month:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_ConfigGetAttribute(FISOC_config, end_year, label='end_year:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_outputInterval, label='OM_outputInterval', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF (MOD(dt_ratio,OM_outputInterval).NE.0) THEN
       msg = "ERROR: dt_ratio/OM_outputInterval is required to be integer"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    CALL ESMF_TimeIntervalSet(OM_dt, s=OM_dt_sec, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_TimeIntervalSet(ISM_dt, s=ISM_dt_sec, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_TimeSet(startTime, yy=start_year, mm=start_month, dd=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF (gotRunLength) THEN
       runLength_secs = runLength_ISM_steps * ISM_dt_sec
       CALL ESMF_TimeIntervalSet(runLength_timeInterval, s=runLength_secs, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       endTime = startTime + runLength_timeInterval
       msg = "used runLength_ISM_steps to set the end time"  
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)       
    ELSE
       CALL ESMF_TimeSet(endTime, yy=end_year, mm=end_month, dd=1, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    FISOC_clock = ESMF_ClockCreate(OM_dt, startTime, stopTime=endTime, &
         name="FISOC main clock", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    alarm_OM = ESMF_AlarmCreate(clock=FISOC_clock, name="alarm_OM", &
         ringTime=startTime, ringInterval=OM_dt, sticky=.FALSE., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    alarm_OM_output = ESMF_AlarmCreate(clock=FISOC_clock, name="alarm_OM_output", &
         ringTime=(startTime+(OM_dt*(OM_outputInterval-1))), ringInterval=OM_dt*OM_outputInterval, sticky=.FALSE., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    alarm_ISM = ESMF_AlarmCreate(clock=FISOC_clock, name="alarm_ISM", &
         ringTime=((startTime+ISM_dt)-OM_dt), ringInterval=ISM_dt, sticky=.FALSE., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    alarm_ISM_exportAvailable = ESMF_AlarmCreate(clock=FISOC_clock, &
         name="alarm_ISM_exportAvailable", &
         ringTime=(startTime+ISM_dt), ringInterval=ISM_dt, sticky=.FALSE.,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_ClockGet(FISOC_clock, currTime=FISOC_time, startTime=FISOC_startTime, stopTime=FISOC_endTime, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    msg = "created and initialised clocks and alarms"  
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_setClocks


  !--------------------------------------------------------------------------------------
  ! set values of fields in this bundle to zero
  SUBROUTINE FISOC_zeroBundle(fieldBundle,rc)
    
    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: fieldBundle
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
# if defined(FISOC_OM_GRID)
    REAL(ESMF_KIND_R8),POINTER            :: field_ptr(:,:)
# elif defined(FISOC_OM_MESH)
    REAL(ESMF_KIND_R8),POINTER            :: field_ptr(:)
# endif
    INTEGER                               :: fieldCount, ii

    rc = ESMF_FAILURE

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! loop over fields, setting all values to zero
    DO ii=1,fieldCount
      CALL ESMF_FieldGet(field=fieldList(ii), farrayPtr=field_ptr, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      field_ptr = 0.0
      
    END DO
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_zeroBundle

  !--------------------------------------------------------------------------------------
  ! add the values for the field in the bundle to the cumulator.  Later it will be divided 
  ! by the number of cumulation operations to give the average. 
  SUBROUTINE FISOC_cumulateFB(fieldBundle,FBcumulator,rc)

    TYPE(ESMF_fieldbundle),INTENT(IN)     :: fieldBundle
    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: FBcumulator
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    INTEGER                               :: fieldCount, ii, fieldCountCum
    REAL(ESMF_KIND_R8)                    :: init_value
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:),fieldListCum(:)
    TYPE(ESMF_TypeKind_Flag)              :: fieldTypeKind
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName, fieldNameCum
# if defined(FISOC_OM_GRID)
    REAL(ESMF_KIND_R8),POINTER            :: fieldCum_ptr(:,:), field_ptr(:,:) 
# elif defined(FISOC_OM_MESH)
    REAL(ESMF_KIND_R8),POINTER            :: fieldCum_ptr(:), field_ptr(:) 
# endif
    TYPE(ESMF_GRID)                       :: grid

    rc = ESMF_FAILURE

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(FBcumulator, fieldCount=fieldCountCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldListCum(fieldCountCum))
    CALL ESMF_FieldBundleGet(FBcumulator, fieldList=fieldListCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (fieldCountCum .NE. fieldCount) THEN
       msg = 'ERROR: fieldBundle and cumulator fieldBundle contain different number of fields '
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    DO ii=1,fieldCount

       CALL ESMF_FieldGet(fieldList(ii), name=fieldName, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(fieldListCum(ii), name=fieldNameCum, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (TRIM(fieldNameCum) .NE. TRIM(fieldName)//"_cum") THEN
          msg = 'ERROR: cumulator field name does not match field name '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    
       CALL ESMF_FieldGet(field=fieldList(ii), localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(field=fieldListCum(ii), localDe=0, farrayPtr=fieldCum_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (SIZE(field_ptr).NE.SIZE(fieldCum_ptr)) THEN
          msg = 'ERROR: cumulator field array pointer dims do not match '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    
       fieldCum_ptr = fieldCum_ptr + field_ptr

    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_cumulateFB


  !--------------------------------------------------------------------------------------
  ! Calculate the average from the cumulator, and stick it back in the field bundle. 
  ! We may one day wish to do other time processing such as keeping the total intead of 
  ! averaging (cant think why though...).
  ! Note: there is a lot of code duplication here of FISOC_cumulateFB.  Is there potential 
  ! for some kind of operatore overloading for field bundle operations?  Perhaps this 
  ! already exists in ESMF?
  SUBROUTINE FISOC_processCumulator(fieldBundle,FBcumulator,FISOC_config,rc)

    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: FBcumulator, fieldBundle
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

!    TYPE(ESMF_field)                      :: cumulatorField
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName, fieldNameCum
    INTEGER                               :: OM_cum_steps, ii, fieldCount, fieldCountCum 
    TYPE(ESMF_field)                      :: fieldCum, field
# if defined(FISOC_OM_GRID)
    REAL(ESMF_KIND_R8),POINTER            :: fieldCum_ptr(:,:), field_ptr(:,:) 
# elif defined(FISOC_OM_MESH)
    REAL(ESMF_KIND_R8),POINTER            :: fieldCum_ptr(:), field_ptr(:) 
# endif
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:),fieldListCum(:)
    TYPE(ESMF_TypeKind_Flag)              :: fieldTypeKind

    rc = ESMF_FAILURE

    ! get number of cumulator steps: ts_ratio / OM_outputInterval
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_cum_steps, 'OM_cum_steps', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(FBcumulator, fieldCount=fieldCountCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldListCum(fieldCountCum))
    CALL ESMF_FieldBundleGet(FBcumulator, fieldList=fieldListCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (fieldCountCum .NE. fieldCount) THEN
       msg = 'ERROR: fieldBundle and cumulator fieldBundle contain different number of fields '
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    DO ii=1,fieldCount

       CALL ESMF_FieldGet(fieldList(ii), name=fieldName, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(fieldListCum(ii), name=fieldNameCum, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (TRIM(fieldNameCum) .NE. TRIM(fieldName)//"_cum") THEN
          msg = 'ERROR: cumulator field name does not match field name '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    
       CALL ESMF_FieldGet(field=fieldList(ii), localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(field=fieldListCum(ii), localDe=0, farrayPtr=fieldCum_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (SIZE(field_ptr).NE.SIZE(fieldCum_ptr)) THEN
          msg = 'ERROR: cumulator field array pointer dims do not match '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    
       field_ptr = fieldCum_ptr / OM_cum_steps

    END DO

    CALL FISOC_zeroBundle(FBcumulator,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_processCumulator


  !--------------------------------------------------------------------------------------
  ! populate a field bundle for cumulating outputs using the field info from an 
  ! existing field bundle.
  SUBROUTINE FISOC_initCumulatorFB(fieldBundle,FBcumulator,rc)

    TYPE(ESMF_fieldbundle),INTENT(IN)     :: fieldBundle
    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: FBcumulator
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    INTEGER                               :: fieldCount, ii
    REAL(ESMF_KIND_R8)                    :: init_value
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:),fieldListCum(:)
    TYPE(ESMF_TypeKind_Flag)              :: fieldTypeKind
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName, fieldNameCum
# if defined(FISOC_OM_GRID)
    TYPE(ESMF_GRID)                       :: grid
# elif defined(FISOC_OM_MESH)
    TYPE(ESMF_MESH)                       :: mesh
# endif

    rc = ESMF_FAILURE

    init_value = 0.0

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! loop over fields, creating a cumulator for each to add to the cumulator field bundle
    ! (cumulator is just another field in which we will sum the field over time)
    ALLOCATE(fieldListCum(fieldCount))
    DO ii=1,fieldCount

       CALL ESMF_FieldGet(fieldList(ii), name=fieldName, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       fieldNameCum = TRIM(fieldName)//"_cum"

# if defined(FISOC_OM_GRID)
       CALL ESMF_FieldGet(fieldList(ii), grid=grid, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       fieldListCum(ii) = ESMF_FieldCreate(grid, typekind=fieldTypeKind, name=fieldNameCum, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

# elif defined(FISOC_OM_MESH)
       CALL ESMF_FieldGet(fieldList(ii), mesh=mesh, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       fieldListCum(ii) = ESMF_FieldCreate(mesh, typekind=fieldTypeKind, name=fieldNameCum, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# else
       msg = "ERROR: FISOC does not recognise OM geom type."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif

    END DO

    CALL ESMF_FieldBundleAdd(FBcumulator, fieldListCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_initCumulatorFB


  !--------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE FISOC_ConfigDerivedAttributeInteger(FISOC_config, derivedAttribute, label,rc)
    
    CHARACTER(len=*),INTENT(IN)           :: label
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT)                   :: derivedAttribute
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc
    
    INTEGER                               :: OM_dt_sec, dt_ratio, OM_outputInterval, OM_AFF
    LOGICAL                               :: APPLY_OM_AFF    

    rc = ESMF_FAILURE

    SELECT CASE(label)
       
    CASE('OM_OPEN_OCEAN') ! used in masking
      CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_OPEN_OCEAN:', rc=rc)
      IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
        derivedAttribute = 0
        msg = "WARNING: OM_OPEN_OCEAN not found, setting to 0"
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__)
      ELSE
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__, rcToReturn=rc)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      END IF

    CASE('ISM_dt_sec')
      CALL ESMF_ConfigGetAttribute(FISOC_config, dt_ratio, label='dt_ratio:', rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__, rcToReturn=rc)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      CALL ESMF_ConfigGetAttribute(FISOC_config, OM_dt_sec, label='OM_dt_sec:', rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__, rcToReturn=rc)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      derivedAttribute = OM_dt_sec * dt_ratio

    CASE('OM_outputInterval')
      CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_outputInterval:', rc=rc)
      IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
        derivedAttribute = 1
        msg = "WARNING: OM_outputInterval not found, setting to 1."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__)
      ELSE
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__, rcToReturn=rc)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      END IF

    CASE('EntriesPerFile')
      CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='EntriesPerFile:', rc=rc)
      IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
        derivedAttribute = 1
        msg = "WARNING: EntriesPerFile not found, setting to 1."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__)
      ELSE
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__, rcToReturn=rc)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      END IF

    CASE('OM_cum_steps')
       CALL ESMF_ConfigGetAttribute(FISOC_config, dt_ratio, label='dt_ratio:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_outputInterval, label='OM_outputInterval', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       derivedAttribute = dt_ratio / OM_outputInterval
       
    CASE('OM_dt_sec','OM_dt_sec:') ! OM run interval (this is "derived" when using accelerated forcing).
      CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_dt_sec:', rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      CALL FISOC_ConfigDerivedAttribute(FISOC_config, APPLY_OM_AFF, 'APPLY_OM_AFF:',rc=rc) 
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      IF (APPLY_OM_AFF) THEN
        CALL ESMF_ConfigGetAttribute(FISOC_config, OM_AFF, label='OM_AFF:', rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        IF ( MOD(derivedAttribute,OM_AFF) .EQ. 0) THEN
          derivedAttribute = derivedAttribute / OM_AFF
        ELSE
          msg = "OM_dt_sec must be divisible by OM_AFF"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        END IF
      END IF
      
    CASE DEFAULT
      msg = 'ERROR: unrecognised derived config attribute label '//label
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
           line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_ConfigDerivedAttributeInteger
  

  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ConfigDerivedAttributeReal(FISOC_config, derivedAttribute, label,rc)
    
    CHARACTER(len=*),INTENT(IN)           :: label
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    REAL(ESMF_KIND_R8),INTENT(OUT)        :: derivedAttribute
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    
    rc = ESMF_FAILURE

    SELECT CASE(label)
      
    CASE('OM_WCmin','OM_WCmin:') ! OM water column minimum thickness
      CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_WCmin:', rc=rc)
      IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
        derivedAttribute = 0.0
        msg = "WARNING: OM_WCmin not found, setting to zero."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__)
      ELSE
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__, rcToReturn=rc)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      END IF
      
    CASE('OM_CavCorr','OM_CavCorr:') ! cavity drift corrrection factor
      CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_CavCorr:', rc=rc)
      IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
        derivedAttribute = 0.2
        msg = "WARNING: OM_CavCorr not found, setting to 0.2."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__)
      ELSE
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__, rcToReturn=rc)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      END IF
      
    CASE DEFAULT
      msg = 'ERROR: unrecognised derived config real attribute label '//label
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
           line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
    END SELECT
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_ConfigDerivedAttributeReal


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ConfigDerivedAttributeRegridMethod(FISOC_config, derivedAttribute, label,rc)
    
    CHARACTER(len=*),INTENT(IN)             :: label
    TYPE(ESMF_config),INTENT(INOUT)         :: FISOC_config
    TYPE(ESMF_RegridMethod_Flag),INTENT(OUT):: derivedAttribute
    INTEGER,OPTIONAL,INTENT(OUT)            :: rc
    
    CHARACTER(len=ESMF_MAXSTR)              :: regridMethodChar
    INTEGER                                 :: rc_local

    rc = ESMF_FAILURE

    SELECT CASE(label)

    CASE('ISM2OM_regrid:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, regridMethodChar, label='ISM2OM_regrid:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          regridMethodChar = "ESMF_REGRIDMETHOD_BILINEAR"
          msg = "WARNING: ISM2OM_regrid not found, setting to bilinear."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE('OM2ISM_regrid:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, regridMethodChar, label='OM2ISM_regrid:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          regridMethodChar = "ESMF_REGRIDMETHOD_BILINEAR"
          msg = "WARNING: OM2ISM_regrid not found, setting to bilinear."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE('Regrid_method')
       CALL ESMF_ConfigGetAttribute(FISOC_config, regridMethodChar, label='Regrid_method:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          regridMethodChar = "ESMF_REGRIDMETHOD_BILINEAR"
          msg = "WARNING: Regrid_method not found, setting to bilinear."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE DEFAULT
       msg = 'ERROR: unrecognised derived config attribute label: '//label
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END SELECT
    
    SELECT CASE (regridMethodChar)
    CASE("ESMF_REGRIDMETHOD_BILINEAR")
       derivedAttribute = ESMF_REGRIDMETHOD_BILINEAR
    CASE("ESMF_REGRIDMETHOD_NEAREST_DTOS")
       derivedAttribute = ESMF_REGRIDMETHOD_NEAREST_DTOS
    CASE("ESMF_REGRIDMETHOD_NEAREST_STOD")
       derivedAttribute = ESMF_REGRIDMETHOD_NEAREST_STOD
    CASE("ESMF_REGRIDMETHOD_CONSERVE")
       derivedAttribute = ESMF_REGRIDMETHOD_CONSERVE
!    CASE("ESMF_REGRIDMETHOD_CONSERVE_2ND")
!       derivedAttribute = ESMF_REGRIDMETHOD_CONSERVE_2ND
    CASE("ESMF_REGRIDMETHOD_PATCH")
       derivedAttribute = ESMF_REGRIDMETHOD_PATCH
    CASE DEFAULT
       msg = 'ERROR: regrid method NYI in FISOC: '//regridMethodChar
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)          
    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ConfigDerivedAttributeRegridMethod


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ConfigDerivedAttributeExtrapMethod(FISOC_config, derivedAttribute, label,rc)
    
    CHARACTER(len=*),INTENT(IN)             :: label
    TYPE(ESMF_config),INTENT(INOUT)         :: FISOC_config
    TYPE(ESMF_ExtrapMethod_Flag),INTENT(OUT):: derivedAttribute
    INTEGER,OPTIONAL,INTENT(OUT)            :: rc
    
    CHARACTER(len=ESMF_MAXSTR)              :: extrapMethodChar
    INTEGER                                 :: rc_local

    rc = ESMF_FAILURE

    SELECT CASE(label)

    CASE('ISM2OM_extrap:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, extrapMethodChar, label='ISM2OM_extrap:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          extrapMethodChar = "ESMF_EXTRAPMETHOD_NONE"
          msg = "WARNING: ISM2OM_extrap not found, setting to ESMF_EXTRAPMETHOD_NONE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE('OM2ISM_extrap:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, extrapMethodChar, label='OM2ISM_extrap:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          extrapMethodChar = "ESMF_EXTRAPMETHOD_NONE"
          msg = "WARNING: OM2ISM_extrap not found, setting to ESMF_EXTRAPMETHOD_NONE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE DEFAULT
       msg = 'ERROR: unrecognised derived config attribute label: '//label
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END SELECT

    
    SELECT CASE (extrapMethodChar)
      
    CASE("ESMF_EXTRAPMETHOD_NONE")
      derivedAttribute = ESMF_EXTRAPMETHOD_NONE
    CASE("ESMF_EXTRAPMETHOD_NEAREST_IDAVG")
      derivedAttribute = ESMF_EXTRAPMETHOD_NEAREST_IDAVG
    CASE("ESMF_EXTRAPMETHOD_NEAREST_STOD")
      derivedAttribute = ESMF_EXTRAPMETHOD_NEAREST_STOD

    CASE DEFAULT
      msg = 'ERROR: extrap method NYI in FISOC: '//extrapMethodChar
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
           line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)          
    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ConfigDerivedAttributeExtrapMethod

  
  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ConfigDerivedAttributeLogical(FISOC_config, derivedAttribute, label,rc)
    
    CHARACTER(len=*),INTENT(IN)            :: label
    TYPE(ESMF_config),INTENT(INOUT)        :: FISOC_config
    LOGICAL,INTENT(OUT)                    :: derivedAttribute
    INTEGER,OPTIONAL,INTENT(OUT)           :: rc

    INTEGER                                :: rc_local, length
    REAL(ESMF_KIND_R8)                     :: realAttribute
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: stringList(:)
    CHARACTER(len=ESMF_MAXSTR)             :: listName
    
    rc = ESMF_FAILURE
    
    SELECT CASE(label)
      
      ! Apply accelerated ocean forcing only if OM_AFF is given and is a positive number other than 1.0
    CASE('APPLY_OM_AFF','APPLY_OM_AFF:')
      CALL ESMF_ConfigGetAttribute(FISOC_config, realAttribute, label='OM_AFF:', rc=rc_local)
      IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
        derivedAttribute = .FALSE.
        msg = "WARNING: OM_AFF not found, not applying accelerated forcing"
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__)
      ELSE
        IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        IF ( realAttribute.GT.0.0) THEN
          IF ( realAttribute.EQ.1.0 ) THEN
            derivedAttribute = .FALSE.
            msg = "WARNING: OM_AFF set to 1.0, not applying accelerated forcing"
            CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
                 line=__LINE__, file=__FILE__)
          ELSE
            derivedAttribute = .TRUE.
          END IF
        ELSE
          msg = "ERROR: OM_AFF must be positive or omitted"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        END IF
      END IF
      
    CASE('WET2DRY','WET2DRY:')
      listName = "WET2DRY_vars:"
      CALL FISOC_getStringListFromConfig(FISOC_config,listName,stringList,rc,returnCount=length)
      IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
        msg = "WET2DRY_vars not found in FISOC config file, setting WET2DRY to false."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__, rc=rc)          
        derivedAttribute = .FALSE.
      ELSE IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) THEN
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      ELSE IF (length.EQ.0) THEN
        msg = "WET2DRY_vars has zero length, setting WET2DRY to false."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
             line=__LINE__, file=__FILE__, rc=rc)          
        derivedAttribute = .FALSE.
      ELSE
        derivedAttribute = .TRUE.
      END IF
      
    CASE('profiling','profiling:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='profiling:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = .FALSE.
          msg = "WARNING: profiling not found, setting to .FALSE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE('OM_initCavityFromISM','OM_initCavityFromISM:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_initCavityFromISM:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = .FALSE.
          msg = "WARNING: OM_initCavityFromISM not found, setting to .FALSE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE('OM_initFrontFromISM','OM_initFrontFromISM:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_initFrontFromISM:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = .FALSE.
          msg = "WARNING: OM_initFrontFromISM not found, setting to .FALSE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
       
    CASE('ISM2OM_init_vars','ISM2OM_init_vars:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='ISM2OM_init_vars:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = .TRUE.
          msg = "WARNING: ISM2OM_init_vars not found, setting to .TRUE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE('ISM_maskOMvars','ISM_maskOMvars:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='ISM_maskOMvars:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = .FALSE.
          msg = "WARNING: ISM_maskOMvars not found, setting to .FALSE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE('ISM_UseOMGrid','ISM_UseOMGrid:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='ISM_UseOMGrid:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = .FALSE.
          msg = "ISM_UseOMGrid not found, setting to .FALSE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE('OM_UseISMGrid','OM_UseISMGrid:')
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_UseISMGrid:', rc=rc_local)
       IF (rc_local.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = .FALSE.
          msg = "OM_UseISMGrid not found, setting to .FALSE."
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc_local, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE DEFAULT
       msg = 'ERROR: unrecognised derived config attribute label: '//label
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ConfigDerivedAttributeLogical


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ConfigDerivedAttributeString(FISOC_config, derivedAttribute, label,rc)
    
    CHARACTER(len=*),INTENT(IN)           :: label
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    CHARACTER(len=ESMF_MAXSTR),INTENT(OUT):: derivedAttribute
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc
    
    CHARACTER(len=ESMF_MAXSTR)            :: OM_cavityUpdate

    rc = ESMF_FAILURE

    SELECT CASE(label)

    CASE("OM_coords")
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='OM_coords:', rc=rc)
       IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = "Cartesian"
          msg = "OM_coords not found, setting to Cartesian"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

    CASE("FileStyle")
       CALL ESMF_ConfigGetAttribute(FISOC_config, derivedAttribute, label='FileStyle:', rc=rc)
       IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
          derivedAttribute = "None"
          msg = "FileStyle not found, setting to None"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__)
       ELSE
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF

       
    CASE("IceDraft")
       CALL ESMF_ConfigGetAttribute(FISOC_config, OM_cavityUpdate, label='OM_cavityUpdate:', rc=rc)
print*,'catch error and set default if missing att'
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!          derivedAttribute = 'draft'
       SELECT CASE(OM_cavityUpdate)
       CASE('RecentIce','Linterp')
          derivedAttribute = 'actual'
       CASE('Rate', 'CorrectedRate')
          derivedAttribute = 'rate'
       CASE DEFAULT
       END SELECT
 
    CASE DEFAULT
       msg = 'ERROR: unrecognised derived config attribute label '
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ConfigDerivedAttributeString


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ConfigDerivedAttributeStaggerLocArray(FISOC_config, derivedAttribute, label, rc)
    
    CHARACTER(len=*),INTENT(IN)           :: label
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_staggerLoc),INTENT(INOUT)   :: derivedAttribute(:)
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc
    
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: attribute_stringList(:)
    INTEGER                               :: OM_dt_sec, dt_ratio, OM_outputInterval

    rc = ESMF_FAILURE

    SELECT CASE(label)

    CASE('OM_ReqVars_stagger:')
       CALL FISOC_getStringListFromConfig(FISOC_config, label, attribute_stringList,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL FISOC_OM_staggerCodes(derivedAttribute,attribute_stringList,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CASE DEFAULT
       msg = 'ERROR: unrecognised derived config attribute label '
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ConfigDerivedAttributeStaggerLocArray


  !--------------------------------------------------------------------------------------
  ! convert a list of strings describning stagger location to ESMF stagger location 
  ! integer codes (for use in ESMF operations such as creating fields on a grid).
  SUBROUTINE FISOC_OM_staggerCodes(staggerLoc,staggerChar,rc)
    
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE,INTENT(IN) :: staggerChar(:)
    TYPE(ESMF_staggerLoc),INTENT(OUT)                 :: staggerLoc(:)
    INTEGER,INTENT(OUT),OPTIONAL                      :: rc

    INTEGER                                           :: ii

    rc = ESMF_FAILURE
    
    IF (size(staggerLoc) .ne. size(staggerChar)) THEN
       msg = 'ERROR: stagger lists must be same length'
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
       
    DO ii = 1,size(staggerLoc)
       SELECT CASE(staggerChar(ii))
       CASE("EDGE1")
          staggerLoc(ii) = ESMF_STAGGERLOC_EDGE1
       CASE("EDGE2")
          staggerLoc(ii) = ESMF_STAGGERLOC_EDGE2
       CASE("CENTER")
          staggerLoc(ii) = ESMF_STAGGERLOC_CENTER
       CASE("CORNER")
          staggerLoc(ii) = ESMF_STAGGERLOC_CORNER
       CASE DEFAULT
          msg = 'ERROR: unrecognised staggerLoc string '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END SELECT
    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_staggerCodes


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_populateFieldBundleOn2dGrid(fieldNames,fieldBundle,grid,init_value,&
       fieldStagger,TLW,TUW,RouteHandle,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN)     :: fieldNames(:)
    TYPE(ESMF_grid),INTENT(IN)                :: grid
    REAL(ESMF_KIND_R8),INTENT(IN),OPTIONAL    :: init_value
    TYPE(ESMF_staggerLoc),INTENT(IN),OPTIONAL :: fieldStagger(:)
    INTEGER,INTENT(IN),OPTIONAL               :: TLW(2), TUW(2)
    TYPE(ESMF_RouteHandle),INTENT(IN),OPTIONAL:: RouteHandle

    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: fieldBundle
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: ii, jj, localDECount
    REAL(ESMF_KIND_R8)                    :: initial_value
    TYPE(ESMF_field)                      :: field
    REAL(ESMF_KIND_R8),POINTER            :: field_ptr(:,:) 

    rc = ESMF_FAILURE

    NULLIFY(field_ptr)

    IF (PRESENT(init_value)) THEN
       initial_value = init_value
    ELSE
       initial_value = 0.0
    END IF

    IF ( ( (PRESENT(TUW)) .AND. (.NOT.PRESENT(TLW)) )   &
         .OR.                                           &
         ( (PRESENT(TLW)) .AND. (.NOT.PRESENT(TUW)) ) ) &
         THEN
       msg = "Expecting neither or both of TUW and TLW"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    IF (PRESENT(RouteHandle).AND.(.NOT.PRESENT(TUW))) THEN 
       msg = "Expecting TUW and TLW if computing halo route handle"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    CALL ESMF_GridGet(grid, localDECount=localDECount, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii=1,SIZE(fieldNames)
       IF (PRESENT(fieldStagger).AND.PRESENT(TUW)) THEN
          field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8,      &
               name=TRIM(fieldNames(ii)),                                &
               totalLWidth=TLW,totalUWidth=TUW,                          &
               staggerloc=fieldStagger(ii),                              &
               rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ELSEIF (PRESENT(fieldStagger)) THEN
          field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8,      &
               name=TRIM(fieldNames(ii)),                                &
               staggerloc=fieldStagger(ii),                              &
               rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ELSEIF (PRESENT(TUW)) THEN
          field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8,      &
               name=TRIM(fieldNames(ii)),                                &
               totalLWidth=TLW,totalUWidth=TUW,                          &
               indexflag=ESMF_INDEX_GLOBAL,             &
               rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ELSE 
          field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8,      &
               name=TRIM(fieldNames(ii)),                                &
               rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END IF

       IF ( (ii.EQ.1) .AND. (PRESENT(RouteHandle)) ) THEN
         print*,"*** need RH ***"
         print*,"*** need RH ***"
       END IF

       CALL ESMF_FieldGet(field=field, localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       field_ptr      = initial_value
       NULLIFY(field_ptr)

       CALL ESMF_FieldBundleAdd(fieldBundle, (/field/), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END DO

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_populateFieldBundleOn2dGrid
  

  !--------------------------------------------------------------------------------------
  ! Note: this subroutine is almost identical to the grid version, and should probably 
  ! be auto generated in a precompile step rather than the current hard coded duplication.
  SUBROUTINE FISOC_populateFieldBundleOnMesh(fieldNames,fieldBundle,mesh,init_value,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: fieldNames(:)
    TYPE(ESMF_mesh),INTENT(IN)            :: mesh
    REAL(ESMF_KIND_R8),INTENT(IN),OPTIONAL:: init_value

    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: fieldBundle
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: ii
    REAL(ESMF_KIND_R8)                    :: initial_value
    TYPE(ESMF_field)                      :: field
    REAL(ESMF_KIND_R8),POINTER            :: field_ptr(:) 

    rc = ESMF_FAILURE

    NULLIFY(field_ptr)

    IF (PRESENT(init_value)) THEN
       initial_value = init_value
    ELSE
       initial_value = 0.0
    END IF

    DO ii=1,SIZE(fieldNames)
       field = ESMF_FieldCreate(mesh, typekind=ESMF_TYPEKIND_R8, name=TRIM(fieldNames(ii)), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_FieldGet(field=field, localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       field_ptr = initial_value       
       CALL ESMF_FieldBundleAdd(fieldBundle, (/field/), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    NULLIFY(field_ptr)

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_populateFieldBundleOnMesh
  

  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_getStringListFromConfig(config,label,stringList,rc,returnCount)

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE,INTENT(INOUT) :: stringList(:)
    TYPE(ESMF_config),INTENT(INOUT)      :: config
    CHARACTER(len=ESMF_MAXSTR),INTENT(IN):: label
    INTEGER,INTENT(OUT)                  :: rc
    INTEGER,INTENT(OUT),OPTIONAL         :: returnCount
    
    CHARACTER(len=ESMF_MAXSTR)           :: dummyString
    INTEGER                              :: listCount,ii

    CHARACTER(len=ESMF_MAXSTR)           :: local_label

    rc = ESMF_FAILURE

    ! point config to start of list 
    CALL ESMF_ConfigFindLabel(config,TRIM(label),rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         RETURN
    
    ! how many items in list?
    listCount = 0
    DO WHILE (rc.EQ.0)
       CALL ESMF_ConfigGetAttribute(config, dummyString,rc=rc) 
       IF  (rc.EQ.ESMF_RC_NOT_FOUND) EXIT
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       listCount = listCount + 1
    END DO
    CALL ESMF_ConfigFindLabel(config, TRIM(label),rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! retrieve list items
    ALLOCATE(stringList(listCount))
    DO ii = 1,listCount
       CALL ESMF_ConfigGetAttribute(config, stringList(ii),rc=rc) 
       IF  (rc.EQ.ESMF_RC_NOT_FOUND) EXIT
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    IF (PRESENT(returnCount)) THEN
      returnCount = listCount
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_getStringListFromConfig


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_getIntegerListFromConfig(config,label,IntegerList,rc)

    INTEGER,ALLOCATABLE,INTENT(INOUT)    :: IntegerList(:)

    TYPE(ESMF_config),INTENT(INOUT)      :: config

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN):: label
    INTEGER,INTENT(OUT),OPTIONAL         :: rc

    CHARACTER(len=ESMF_MAXSTR)           :: dummyInteger
    INTEGER                              :: listCount,ii

    rc = ESMF_FAILURE

    ! point config to start of list 
    CALL ESMF_ConfigFindLabel(config,TRIM(label),rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         RETURN
    
    ! how many items in list?
    listCount = 0
    DO WHILE (rc.EQ.0)
       CALL ESMF_ConfigGetAttribute(config, dummyInteger,rc=rc) 
       IF  (rc.EQ.ESMF_RC_NOT_FOUND) EXIT
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       listCount = listCount + 1
    END DO

    ALLOCATE(IntegerList(listCount))

    ! go back to the start of the list
    CALL ESMF_ConfigFindLabel(config, TRIM(label),rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! retrieve list items
    DO ii = 1,listCount
       CALL ESMF_ConfigGetAttribute(config, IntegerList(ii),rc=rc) 
       IF  (rc.EQ.ESMF_RC_NOT_FOUND) EXIT
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_getIntegerListFromConfig

  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_getRealListFromConfig(config,label,RealList,rc)

    REAL(ESMF_KIND_R8),ALLOCATABLE,INTENT(INOUT) :: RealList(:)

    TYPE(ESMF_config),INTENT(INOUT)      :: config

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN):: label
    INTEGER,INTENT(OUT),OPTIONAL         :: rc

    CHARACTER(len=ESMF_MAXSTR)           :: dummyReal
    INTEGER                              :: listCount,ii

    rc = ESMF_FAILURE

    ! point config to start of list 
    CALL ESMF_ConfigFindLabel(config,TRIM(label),rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         RETURN
    
    ! how many items in list?
    listCount = 0
    DO WHILE (rc.EQ.0)
       CALL ESMF_ConfigGetAttribute(config, dummyReal,rc=rc) 
       IF  (rc.EQ.ESMF_RC_NOT_FOUND) EXIT
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       listCount = listCount + 1
    END DO

    ALLOCATE(RealList(listCount))

    ! go back to the start of the list
    CALL ESMF_ConfigFindLabel(config, TRIM(label),rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! retrieve list items
    DO ii = 1,listCount
       CALL ESMF_ConfigGetAttribute(config, RealList(ii),rc=rc) 
       IF  (rc.EQ.ESMF_RC_NOT_FOUND) EXIT
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_getRealListFromConfig



  !--------------------------------------------------------------------
  SUBROUTINE FISOC_FB2NC(filename,fieldBundle,FISOC_config)

    CHARACTER(len=ESMF_MAXSTR),INTENT(INOUT):: filename
    TYPE(ESMF_FieldBundle),INTENT(IN)       :: fieldBundle
    TYPE(ESMF_config),INTENT(INOUT)         :: FISOC_config
    
    TYPE(ESMF_FileStatus_Flag)              :: NC_status
    INTEGER                                 :: rc
    CHARACTER(len=ESMF_MAXSTR)              :: output_dir

    CALL ESMF_ConfigGetAttribute(FISOC_config, output_dir, label='output_dir:', rc=rc)
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
       output_dir = "./"
       msg = "WARNING: output directory not found, setting to current dir"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)
    ELSE
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
 
    filename = TRIM(output_dir)//'/'//TRIM(filename)

    IF (ESMF_IO_NETCDF_PRESENT) THEN
       NC_status=ESMF_FILESTATUS_REPLACE
       CALL  ESMF_FieldBundleWrite(fieldBundle, TRIM(filename),  overwrite=.FALSE., & 
            status=NC_status, iofmt=ESMF_IOFMT_NETCDF, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       msg = "FB2NC: written netcdf file"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    ELSE
       msg = "ERROR: trying to write NetCDF output but NetCDF "// &
            "not present in this ESMF build."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    END IF
    
  END SUBROUTINE FISOC_FB2NC



  !--------------------------------------------------------------------
  ! This is a wrapper for ESMF_FieldRegridStore.  This subroutine 
  ! exists because ESMF_FieldRegridStore does not preserve data in the 
  ! field used to create the route handle.  This wrapper copies the 
  ! data before calling ESMF_FieldRegridStore and writes it back in 
  ! afterwards.
  !
  SUBROUTINE FISOC_FieldRegridStore(vm, InField, OutField, regridmethod, &
       unmappedaction, extrapMethod, routeHandle, rc)

    TYPE(ESMF_VM),INTENT(IN)                          :: vm
    TYPE(ESMF_Field),INTENT(INOUT)                    :: InField, OutField
    TYPE(ESMF_RegridMethod_Flag),INTENT(IN),OPTIONAL  :: regridmethod
    TYPE(ESMF_UnmappedAction_Flag),INTENT(IN),OPTIONAL:: unmappedaction
    TYPE(ESMF_ExtrapMethod_Flag),INTENT(IN),OPTIONAL  :: extrapMethod

    TYPE(ESMF_RouteHandle),INTENT(OUT),OPTIONAL       :: routeHandle
    INTEGER,INTENT(OUT),OPTIONAL                      :: rc

    REAL(ESMF_KIND_R8),ALLOCATABLE   :: InFieldData1D_cp(:), OutFieldData1D_cp(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: InFieldData2D_cp(:,:), OutFieldData2D_cp(:,:)
    REAL(ESMF_KIND_R8),POINTER       :: InFieldData1D(:), OutFieldData1D(:)
    REAL(ESMF_KIND_R8),POINTER       :: InFieldData2D(:,:), OutFieldData2D(:,:)
    INTEGER                          :: InDims, OutDims

    rc = ESMF_FAILURE

    ! check dimensionality of fields
    CALL ESMF_FieldGet(field=InField, dimCount=InDims, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=OutField, dimCount=OutDims, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! take copy of fields
    SELECT CASE(InDims)
    CASE(1)
       CALL ESMF_FieldGet(field=InField, farrayPtr=InFieldData1D, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       ALLOCATE(InFieldData1D_cp(SIZE(InFieldData1D)))
       InFieldData1D_cp = InFieldData1D
    CASE(2)
       CALL ESMF_FieldGet(field=InField, farrayPtr=InFieldData2D, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       ALLOCATE(InFieldData2D_cp(SIZE(InFieldData2D,1),SIZE(InFieldData2D,2)))
       InFieldData2D_cp = InFieldData2D
    CASE DEFAULT
       msg = 'ERROR: field neither 1D nor 2D in regrid wrapper'
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    SELECT CASE(OutDims)
    CASE(1)
       CALL ESMF_FieldGet(field=OutField, farrayPtr=OutFieldData1D, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       ALLOCATE(OutFieldData1D_cp(SIZE(OutFieldData1D)))
       OutFieldData1D_cp = OutFieldData1D
    CASE(2)
       CALL ESMF_FieldGet(field=OutField, farrayPtr=OutFieldData2D, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       ALLOCATE(OutFieldData2D_cp(SIZE(OutFieldData2D,1),SIZE(OutFieldData2D,2)))
       OutFieldData2D_cp = OutFieldData2D       
    CASE DEFAULT
       msg = 'ERROR: field neither 1D nor 2D in regrid wrapper'
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! create the routehandle for regridding
    CALL ESMF_FieldRegridStore(InField, OutField, regridmethod=regridmethod, &
         unmappedaction=unmappedaction, routehandle=routeHandle,             &
         dstMaskValues=(/MASK_OPEN_OCEAN/),                                   &
         extrapMethod=extrapMethod, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    
    ! copy the field data back in, cleaning up as we go
    SELECT CASE(InDims)
    CASE(1)
       InFieldData1D=InFieldData1D_cp
       DEALLOCATE(InFieldData1D_cp)
       NULLIFY(InFieldData1D)
    CASE(2)
       InFieldData2D=InFieldData2D_cp
       DEALLOCATE(InFieldData2D_cp)
       NULLIFY(InFieldData2D)
    END SELECT

    SELECT CASE(OutDims)
    CASE(1)
       OutFieldData1D=OutFieldData1D_cp
       DEALLOCATE(OutFieldData1D_cp)
       NULLIFY(OutFieldData1D)
    CASE(2)
       OutFieldData2D=OutFieldData2D_cp
       DEALLOCATE(OutFieldData2D_cp)
       NULLIFY(OutFieldData2D)
    END SELECT
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_FieldRegridStore


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_shrinkInt(arrayIn,ignoreVal)
    
    INTEGER,INTENT(INOUT),ALLOCATABLE :: arrayIn(:)
    INTEGER,INTENT(IN)                :: ignoreVal

    INTEGER,ALLOCATABLE               :: holder(:) ! temporary holding array
    INTEGER :: ii,numIgnores,nn

    numIgnores = 0
    DO ii = 1,SIZE(arrayIn)
       IF ( arrayIn(ii) .EQ. ignoreVal ) numIgnores = numIgnores + 1
    END DO

    IF ( numIgnores .GT. SIZE(arrayIn) ) THEN
       msg= "FATAL: too many ignored vals in FISOC_shrinkInt"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! populate the temporary holding array
    ALLOCATE(holder(SIZE(arrayIn) - numIgnores))

    holder = ignoreVal
    nn = 0
    DO ii = 1,SIZE(arrayIn)
       IF ( arrayIn(ii) .NE. ignoreVal ) THEN
          nn = nn+1
          holder(nn) = arrayIn(ii)
       END IF
    END DO

    ! now we can resize the input array and use the temporary holding 
    ! array to repopulate the input array.
    DEALLOCATE(arrayIn)
    ALLOCATE(arrayIn(SIZE(holder)))
    arrayIn = holder

    ! tidy up
    DEALLOCATE(holder)
    
  END SUBROUTINE FISOC_shrinkInt




  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_shrinkReal(arrayIn,ignoreVal)
    
    REAL(ESMF_KIND_R8),INTENT(INOUT),ALLOCATABLE :: arrayIn(:)
    REAL(ESMF_KIND_R8),INTENT(IN)                :: ignoreVal

    REAL(ESMF_KIND_R8),ALLOCATABLE               :: holder(:) ! temporary holding array
    INTEGER :: ii,numIgnores,nn

    numIgnores = 0
    DO ii = 1,SIZE(arrayIn)
       IF ( arrayIn(ii) .EQ. ignoreVal ) numIgnores = numIgnores + 1
    END DO

    IF ( numIgnores .GE. SIZE(arrayIn) ) THEN
       msg= "FATAL: too many ignored vals in FISOC_shrinkReal"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! populate the temporary holding array
    ALLOCATE(holder(SIZE(arrayIn) - numIgnores))

    holder = ignoreVal
    nn = 0
    DO ii = 1,SIZE(arrayIn)
       IF ( arrayIn(ii) .NE. ignoreVal ) THEN
          nn = nn+1
          holder(nn) = arrayIn(ii)
       END IF
    END DO

    ! now we can resize the input array and use the temporary holding 
    ! array to repopulate the input array.
    DEALLOCATE(arrayIn)
    ALLOCATE(arrayIn(SIZE(holder)))
    arrayIn = holder

    ! tidy up
    DEALLOCATE(holder)
    
  END SUBROUTINE FISOC_shrinkReal



  ! Wrapper for ESMF_VMAllGather, but allows array sizes to vary over PETs.
  ! We assume Array_local is already allocated and populated, but 
  ! Array_all is not.  Array_all is allocated within this routine.
  SUBROUTINE FISOC_VMAllGather_Int(vm,Array_local,Array_all,rc)

    TYPE(ESMF_vm),INTENT(IN)        :: vm
    INTEGER,INTENT(IN)              :: Array_local(:) ! on the local PET
    INTEGER,ALLOCATABLE,INTENT(OUT) :: Array_all(:)   ! over all PETs
    INTEGER,INTENT(OUT),OPTIONAL    :: rc

    INTEGER,ALLOCATABLE  :: Array_local_holder(:)
    INTEGER :: ArrLocSize, ArrLocSizeMax, ArrAllSize, PETcount

    IF (ALLOCATED(Array_all)) DEALLOCATE(Array_all)

    CALL ESMF_VMGet(vm, petCount=PETcount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ArrLocSize = SIZE(Array_local)

    ! Get the max number of array elements per pet.
    CALL ESMF_VMAllFullReduce(vm, (/ArrLocSize/),ArrLocSizeMax, 1, ESMF_REDUCE_MAX, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Get the total number of array elements over all PETs.
    CALL ESMF_VMAllFullReduce(vm,(/ArrLocSize/),ArrAllSize,1,ESMF_REDUCE_SUM, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! The holding array is a copy of the local array but padded with missing values
    ALLOCATE(Array_local_holder(ArrLocSizeMax))
    Array_local_holder = FISOC_missing
    Array_local_holder(1:ArrLocSize) = Array_local

    ALLOCATE(Array_all(ArrLocSizeMax*PETcount))
    CALL ESMF_VMAllGather(vm,Array_local_holder,Array_all,ArrLocSizeMax,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_shrink(Array_all,FISOC_missing)

    IF (ArrAllSize.NE.SIZE(Array_all)) THEN
       msg = "Array size mismatch after AllGather and shrink. "
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

  END SUBROUTINE FISOC_VMAllGather_Int



  ! Wrapper for ESMF_VMAllGather, but allows array sizes to vary over PETs.
  ! We assume Array_local is already allocated and populated, but 
  ! Array_all is not.  Array_all is allocated within this routine.
  SUBROUTINE FISOC_VMAllGather_Real(vm,Array_local,Array_all,rc)

    TYPE(ESMF_vm),INTENT(IN)                   :: vm
    REAL(ESMF_KIND_R8),INTENT(IN)              :: Array_local(:) ! on the local PET
    REAL(ESMF_KIND_R8),ALLOCATABLE,INTENT(OUT) :: Array_all(:)   ! over all PETs
    INTEGER,INTENT(OUT),OPTIONAL               :: rc

    REAL(ESMF_KIND_R8),ALLOCATABLE  :: Array_local_holder(:)
    INTEGER :: ArrLocSize, ArrLocSizeMax, ArrAllSize, PETcount

    IF (ALLOCATED(Array_all)) DEALLOCATE(Array_all)

    CALL ESMF_VMGet(vm, petCount=PETcount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ArrLocSize = SIZE(Array_local)

    ! Get the max number of array elements per pet.
    CALL ESMF_VMAllFullReduce(vm, (/ArrLocSize/),ArrLocSizeMax, 1, ESMF_REDUCE_MAX, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Get the total number of array elements over all PETs.
    CALL ESMF_VMAllFullReduce(vm,(/ArrLocSize/),ArrAllSize,1,ESMF_REDUCE_SUM, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! The holding array is a copy of the local array but padded with missing values
    ALLOCATE(Array_local_holder(ArrLocSizeMax))
    Array_local_holder = FISOC_missing_R8
    Array_local_holder(1:ArrLocSize) = Array_local

    ALLOCATE(Array_all(ArrLocSizeMax*PETcount))
    CALL ESMF_VMAllGather(vm,Array_local_holder,Array_all,ArrLocSizeMax,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_shrink(Array_all,FISOC_missing_R8)

    IF (ArrAllSize.NE.SIZE(Array_all)) THEN
       msg = "Array size mismatch after AllGather and shrink. "
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

  END SUBROUTINE FISOC_VMAllGather_Real



  !------------------------------------------------------------------------------
  SUBROUTINE Unique1DArray_I4(Arr_a)
    ! Author: Kong, kinaxj@gmail.com
    IMPLICIT NONE
    INTEGER(ESMF_KIND_I4),DIMENSION(:),ALLOCATABLE ::Arr_a,Arr_b
    LOGICAL,DIMENSION(:),ALLOCATABLE            ::mask
    INTEGER,DIMENSION(:),ALLOCATABLE            ::index_vector,indexSos
    INTEGER                                     ::i,j,num
    
    num=SIZE(Arr_a);  ALLOCATE(mask(num)); mask = .TRUE.
    DO i=num,2,-1
       mask(i)=.NOT.(ANY(Arr_a(:i-1)==Arr_a(i)))
    END DO
    
    ! Make an index vector
    ALLOCATE(indexSos(SIZE(PACK([(i,i=1,num)],mask))))
    ALLOCATE(index_vector(SIZE(indexSos))); index_vector=PACK([(i,i=1,num)],mask)
    
    ! Now copy the unique elements of a into b
    ALLOCATE(Arr_b(SIZE(index_vector)))
    Arr_b=Arr_a(index_vector)
    CALL move_alloc(Arr_b,Arr_a)
    
  END SUBROUTINE  Unique1DArray_I4
  
  !------------------------------------------------------------------------------
  SUBROUTINE Unique1DArray_D(Arr_a)
    ! Author: Kong, kinaxj@gmail.com
    IMPLICIT NONE
    REAL(ESMF_KIND_R8),DIMENSION(:),ALLOCATABLE ::Arr_a,Arr_b
    LOGICAL,DIMENSION(:),ALLOCATABLE            ::mask
    INTEGER,DIMENSION(:),ALLOCATABLE            ::index_vector,indexSos
    INTEGER                                     ::i,j,num
    
    num=SIZE(Arr_a);  ALLOCATE(mask(num)); mask = .TRUE.
    DO i=num,2,-1
       mask(i)=.NOT.(ANY(Arr_a(:i-1)==Arr_a(i)))
    END DO
    
    ! Make an index vector
    ALLOCATE(indexSos(SIZE(PACK([(i,i=1,num)],mask))))
    ALLOCATE(index_vector(SIZE(indexSos))); index_vector=PACK([(i,i=1,num)],mask)
    
    ! Now copy the unique elements of a into b
    ALLOCATE(Arr_b(SIZE(index_vector)))
    Arr_b=Arr_a(index_vector)
    CALL move_alloc(Arr_b,Arr_a)
    
  END SUBROUTINE  Unique1DArray_D

  !--------------------------------------------------------------------
  ! Each gridded component should have either a mesh or a grid.
  ! If a gridded component has neither or both it is generally fatal.
  SUBROUTINE FISOC_OneGrid(fatal,grid,mesh)

    LOGICAL,INTENT(IN)                  :: fatal
    TYPE(ESMF_grid),OPTIONAL,INTENT(IN) :: grid
    TYPE(ESMF_mesh),OPTIONAL,INTENT(IN) :: mesh
    
    IF   (                                                        &
         ( (PRESENT(grid)) .AND. (PRESENT(mesh)) )                & 
         .OR.                                                     &
         ( (.NOT.(PRESENT(grid))) .AND. (.NOT.(PRESENT(mesh))) )  & 
         ) THEN

       msg = "Expecting one mesh OR one grid per gridded component"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)

       IF (fatal) CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END IF

  END SUBROUTINE FISOC_OneGrid


  ! --------------------------------------------------------------------------
  ! Get the rank (this is .le. dimension count) fo the first 
  ! field in a field bundle (assume all fields have same rank)
  SUBROUTINE FISOC_getFirstFieldRank(fieldBundle,rank,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT) :: fieldBundle
    INTEGER,INTENT(OUT)                  :: rank
    INTEGER,INTENT(OUT),OPTIONAL         :: rc

    TYPE(ESMF_Field)             :: field
    TYPE(ESMF_Field),ALLOCATABLE :: fieldList(:)
    INTEGER                      :: fieldCount

    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ALLOCATE(fieldList(fieldCount))

    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldGet(fieldList(1), rank=rank, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_getFirstFieldRank


  !--------------------------------------------------------------------
  ! Make a route handle for regridding operations.  Use the first 
  ! fields from the source and target field bundles to set up the 
  ! route handle.
  !
  SUBROUTINE FISOC_makeRHfromFB(sourceFB,targetFB, &
       Regrid_method,Extrap_method,verbose_coupling,regridRouteHandle,vm,rc)
    
    TYPE(ESMF_fieldBundle)                  :: sourceFB, targetFB
    TYPE(ESMF_RouteHandle),INTENT(OUT)      :: regridRouteHandle
    TYPE(ESMF_RegridMethod_Flag),INTENT(IN) :: Regrid_method
    TYPE(ESMF_ExtrapMethod_Flag),INTENT(IN) :: Extrap_method
    TYPE(ESMF_vm),INTENT(IN)                :: vm
    LOGICAL,INTENT(IN)                      :: verbose_coupling
    INTEGER,INTENT(OUT),OPTIONAL            :: rc

    INTEGER                            :: sourceFieldCount, targetFieldCount
    TYPE(ESMF_field),ALLOCATABLE       :: sourceFieldList(:), targetFieldList(:)

    rc = ESMF_FAILURE

    ! Extract source fields from bundle for regridding...
    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(sourceFB, fieldCount=sourceFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ! ... get list of fields from bundle.
    ALLOCATE(sourceFieldList(sourceFieldCount))
    CALL ESMF_FieldBundleGet(sourceFB, fieldList=sourceFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! Extract target fields from bundle...
    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(targetFB, fieldCount=targetFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ! ... get list of fields from bundle.
    ALLOCATE(targetFieldList(targetFieldCount))
    CALL ESMF_FieldBundleGet(targetFB, fieldList=targetFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF (verbose_coupling) THEN
       msg = "coupler extracted fields for creating routehandle"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
           
    IF (SIZE(sourceFieldList).LT.1) THEN
       msg = "source field list less than length 1"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    IF (SIZE(targetFieldList).LT.1) THEN
       msg = "target field list less than length 1"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! Create a route handle to add to the state object.  This will be used for regridding.
    CALL FISOC_FieldRegridStore(vm, sourceFieldList(1), targetFieldList(1), &
         regridmethod=Regrid_method, &
         extrapMethod=Extrap_method, &
         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         routehandle=regridRouteHandle, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_makeRHfromFB


  !--------------------------------------------------------------------
  SUBROUTINE FISOC_getGridFromFB(FB,grid,rc)

    TYPE(ESMF_FieldBundle),INTENT(INOUT)  :: FB
    TYPE(ESMF_grid),INTENT(OUT)           :: grid
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    TYPE(ESMF_Field),ALLOCATABLE  :: FieldList(:)
    INTEGER                       :: FieldCount     

    rc = ESMF_FAILURE

    ! How many fields?...
    CALL ESMF_FieldBundleGet(FB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(FieldList(FieldCount))
    CALL ESMF_FieldBundleGet(FB, fieldList=FieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Get grid from first field in list
    CALL ESMF_FieldGet(FieldList(1), grid=grid, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_getGridFromFB


  !--------------------------------------------------------------------
  SUBROUTINE FISOC_getMeshFromFB(FB,mesh,rc)

    TYPE(ESMF_FieldBundle),INTENT(INOUT)  :: FB
    TYPE(ESMF_mesh),INTENT(OUT)           :: mesh
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    TYPE(ESMF_Field),ALLOCATABLE  :: FieldList(:)
    INTEGER                       :: FieldCount     

    rc = ESMF_FAILURE

    ! How many fields?...
    CALL ESMF_FieldBundleGet(FB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(FieldList(FieldCount))
    CALL ESMF_FieldBundleGet(FB, fieldList=FieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Get grid from first field in list
    CALL ESMF_FieldGet(FieldList(1), mesh=mesh, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_getMeshFromFB



  !--------------------------------------------------------------------
  ! Loop over the field bundle using the same routehandle for each 
  ! regridding operation
  SUBROUTINE FISOC_regridFB(sourceFB,targetFB,regridRH,rc)
    TYPE(ESMF_FieldBundle),INTENT(INOUT) :: sourceFB,targetFB
    TYPE(ESMF_RouteHandle),INTENT(INOUT) :: regridRH
    INTEGER,INTENT(OUT),OPTIONAL         :: rc

    INTEGER                       :: fieldCount, targetFieldCount, ii
    TYPE(ESMF_field),ALLOCATABLE  :: sourceFieldList(:),targetFieldList(:)
    CHARACTER(len=ESMF_MAXSTR)    :: fieldName

    rc = ESMF_FAILURE

    ! How many fields?  And consistency check.
    CALL ESMF_FieldBundleGet(sourceFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldBundleGet(targetFB, fieldCount=targetFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (fieldCount.NE.targetFieldCount) THEN
       msg = "ERROR: mismatch in fieldbundle lengths"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! Get list of fields to be regridded
    ALLOCATE(sourceFieldList(fieldCount))
    CALL ESMF_FieldBundleGet(sourceFB, fieldList=sourceFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Get list of target fields for regridding
    ALLOCATE(targetFieldList(fieldCount))
    CALL ESMF_FieldBundleGet(targetFB, fieldList=targetFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,fieldCount

       CALL ESMF_FieldRegrid(sourceFieldList(ii),targetFieldList(ii), &
            routehandle=regridRH, zeroregion= ESMF_REGION_TOTAL, &
            checkflag=.TRUE.,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldGet(sourceFieldList(ii), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       msg = "Regridded field "//fieldName
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       
    END DO

    DEALLOCATE(targetFieldList)
    DEALLOCATE(sourceFieldList)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_regridFB


  !------------------------------------------------------------------------------
  ! localNodeIDs is an array of all node ids on the current partition.  This can
  ! include nodes at partition boundaries (or halo nodes) owned by neighbouring 
  ! partitions. 
  ! This subroutine returns an array of node ids owned by the current partition. 
  ! This is a subset of all nodes contained in localNodeIDs.
  SUBROUTINE FISOC_locallyOwnedNodes(localPet,localNodeIDS,nodeOwners,ownedNodeIDs)

    INTEGER,INTENT(IN)              :: localPet, localNodeIDs(:), nodeOwners(:)
    INTEGER,ALLOCATABLE,INTENT(OUT) :: ownedNodeIDs(:)

    INTEGER                         :: LON_count, ii

    ! TODO: fail fast if size of nodeOwners .ne. size of localnodeIDs
    
    ! How many of the nodes are locally owned?
    LON_count = 0
    DO ii = 1,SIZE(localNodeIDs)
       IF (nodeOwners(ii).EQ.localPet) THEN
          LON_count = LON_count + 1
       END IF
    END DO
    !TODO: fail fast if already allocated/associated
    ALLOCATE(ownedNodeIDs(LON_count))

    LON_count = 0
    DO ii = 1,SIZE(localNodeIDs)
       IF (nodeOwners(ii).EQ.localPet) THEN
          LON_count = LON_count + 1
          ownedNodeIDs(LON_count) = localNodeIDs(ii)
       END IF
    END DO

  END SUBROUTINE FISOC_locallyOwnedNodes


  !------------------------------------------------------------------------------
  ! The route handle is for array mappings between the ESMF fields defined on 
  ! meshes (with unique nodes) and component (e.g. Elmer) fields with some node 
  ! duplication across partition boundaries).
  !
  ! Variable naming here: "source" refers to the ESMF mesh and fields and "dest" 
  ! (short for destination) refers to the component (e.g. Elmer) mesh and fields. 
  ! Note that both are for holding component information, but the source ones are
  ! in ESMF run time structures and the dest are in native component structures.
  !
  ! Input args:
  ! sourceMesh - Arrays created on this mesh will not duplicate nodes 
  ! nodeIDs - global dest node IDs on the current partition, including nodes 
  ! duplicated across partition boundaries.
  !
  ! Output args:
  ! a routehandle and dest distgrid object to be used later to implement the one 
  ! to many mapping.
  !
  SUBROUTINE FISOC_CreateOneToManyRouteHandle(sourceMesh,nodeIDs,RH,distgridDest,vm)

    TYPE(ESMF_mesh), INTENT(IN) :: sourceMesh
    INTEGER,INTENT(IN)          :: nodeIDs(:) 
    TYPE(ESMF_vm),INTENT(IN)    :: vm

    TYPE(ESMF_RouteHandle),INTENT(OUT) :: RH 
    TYPE(ESMF_distgrid),INTENT(OUT)    :: distgridDest

    TYPE(ESMF_distgrid):: distgridSource
    TYPE(ESMF_array)   :: DummyArr_source  ! will not contain duplicate nodes
    TYPE(ESMF_array)   :: DummyArr_dest    ! may contain duplicate nodes
    INTEGER            :: rc
    INTEGER            :: localPET
    REAL(ESMF_KIND_R8),POINTER :: ptr_source(:),ptr_dest(:)


    ! Create dummy array on the distgrid for ESMF source mesh
    CALL ESMF_MeshGet(sourceMesh, nodalDistgrid=distgridSource, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    DummyArr_source = ESMF_ArrayCreate(distgridSource, ESMF_TYPEKIND_R8, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (localPet.EQ.0) THEN
       msg = "one2many routehandles: creating distrgrid"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    ! Create a distgrid containing sequence indices for the dest array, i.e. 
    ! containing the duplicate node IDs.  This can be used to create an array 
    ! including duplicates, like the dest fields.
    distgridDest  = ESMF_DistgridCreate(arbSeqIndexList=nodeIDs, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    DummyArr_dest = ESMF_ArrayCreate(distgridDest, ESMF_TYPEKIND_R8, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! initialilse arrays to zero, probably not needed
    CALL ESMF_ArrayGet(DummyArr_dest, farrayPtr=ptr_dest, rc=rc)
    CALL ESMF_ArrayGet(DummyArr_source,  farrayPtr=ptr_source, rc=rc)
    ptr_dest=0.0
    ptr_source=0.0
    IF (ASSOCIATED(ptr_dest)) NULLIFY(ptr_dest)
    IF (ASSOCIATED(ptr_source)) NULLIFY(ptr_source)
    
    ! Create the route handle for later use 
    IF (localPet.EQ.0) THEN
       msg = "one2many routehandles: store routehandle"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    CALL ESMF_ArrayRedistStore(DummyArr_source, DummyArr_dest, &
         RH, ignoreUnmatchedIndices=.TRUE., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Tidy up 
    IF (localPet.EQ.0) THEN
       msg = "one2many routehandles: tidy up"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    CALL ESMF_ArrayDestroy(DummyArr_dest, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_ArrayDestroy(DummyArr_source, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_CreateOneToManyRouteHandle

  !------------------------------------------------------------------------------
  ! Taking as input a source field,routehandle and distgrid for the destination
  ! field: extract the array from the source field, create a destination array 
  ! using the distgrid, and use these with the routehandle in an ESMF_ArrayRedist 
  ! operation.  Return a pointer to the the destination array.
  !
  ! Written for application of FISOC_CreateOneToManyRouteHandle.
  !
  SUBROUTINE FISOC_ArrayRedistFromField(RH,sourceField,distgridDest,destArrPtr)
    TYPE(ESMF_routeHandle),INTENT(INOUT) :: RH
    TYPE(ESMF_field),INTENT(IN)          :: sourceField
    TYPE(ESMF_distgrid),INTENT(IN)       :: distgridDest
    REAL(ESMF_KIND_R8),POINTER           :: destArrPtr(:)

    TYPE(ESMF_Array)                     :: sourceArray, destArray
    INTEGER                              :: rc
    
    CALL ESMF_FieldGet(sourceField, array=sourceArray, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! we need to create a dest array on the fly (and then convert from the 
    ! array to the dest variable after the redist)
    destArray = ESMF_ArrayCreate(distgridDest, ESMF_TYPEKIND_R8, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)       
    
    ! now redist the source array onto the dest array
    CALL ESMF_ArrayRedist(sourceArray, destArray, RH, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! ...and get a pointer to the data
    CALL ESMF_ArrayGet(destArray, farrayPtr=destArrPtr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_ArrayRedistFromField


  !------------------------------------------------------------------------------
  ! Routine copied from NSIDC, github repo here:
  ! https://github.com/nsidc/polarstereo-latlon-convert-fortran
  !
  !    DESCRIPTION:                                                            *
  !                                                                            *
  !    This subroutine converts from geodetic latitude and longitude to Polar  *
  !    Stereographic (X,Y) coordinates for the polar regions.  The equations   *
  !    are from Snyder, J. P., 1982,  Map Projections Used by the U.S.         *
  !    Geological Survey, Geological Survey Bulletin 1532, U.S. Government     *
  !    Printing Office.  See JPL Technical Memorandum 3349-85-101 for further  *
  !    details.                                                                *
  !                                                                            *
  !    ARGUMENTS:                                                              *
  !                                                                            *
  !    Variable    Type        I/O    Description                              *
  !                                                                            *
  !    ALAT       REAL*4        I     Geodetic Latitude (degrees, +90 to -90)  *
  !    ALON       REAL*4        I     Geodetic Longitude (degrees, 0 to 360)   *
  !    X          REAL*4        O     Polar Stereographic X Coordinate (km)    *
  !    Y          REAL*4        O     Polar Stereographic Y Coordinate (km)    *
  !                                                                            *
  !                  Written by C. S. Morris - April 29, 1985                  *
  !                  Revised by C. S. Morris - December 11, 1985               *
  !                                                                     
  !                  Revised by V. J. Troisi - January 1990                    *
  !                  SGN - provides hemisphere dependency (+/- 1)              *
  !		    Revised by Xiaoming Li - October 1996                      *
  !		    Corrected equation for RHO                                 *
  !*****************************************************************************
  SUBROUTINE FISOC_MAPLL (X,Y,ALAT,ALON)

    IMPLICIT NONE

    REAL(ESMF_KIND_R8) :: X,Y,ALAT,ALON
    REAL(ESMF_KIND_R8) :: SLAT,SGN,E,RE,RLAT,RLON ! "R" for radians
    REAL(ESMF_KIND_R8) :: E2,PI,MC,RHO,SL,T,TC

    !    REAL*4 X,Y,ALAT,ALONG,E,E2,CDR,PI,SLAT,MC
    !    !    Conversion constant from degrees to radians 
    !    CDR=57.29577951

    !*****************************************************************************
    !                                                                            *
    !    DEFINITION OF CONSTANTS:                                                *
    !                                                                            *
    E2 = .006693883 ! Eccentricity of the Hughes ellipsoid
    E =  sqrt(E2)
    PI=3.141592654
    SLAT=71.0       ! magnitude of latitude of equal area
    SGN=-1.0        ! minus 1 for southern hemisphere
    RE = 6378.273   ! Radius of the earth in kilometers.
    !*****************************************************************************

!    ! Convert longitude to positive degrees
!    if (alon.le.0.0) alon=alon+360.
!    if (alon.ge.360.0) alon=alon-360.

    ! Transform degrees to radians
    RLAT=ABS(ALAT)*PI/180.
    RLON=(ALON)*PI/180.

    ! Compute X and Y in grid coordinates.
    IF (ABS(RLAT).LT.PI/2.) THEN
       T=TAN(PI/4.-RLAT/2.)/((1.-E*SIN(RLAT))/(1.+E*SIN(RLAT)))**(E/2.)
       IF (ABS(90.-SLAT).LT.1.E-5) THEN
          RHO=2.*RE*T/((1.+E)**(1.+E)*(1.-E)**(1.-E))**(1/2.)
       ELSE
          SL=SLAT*PI/180.
          TC=TAN(PI/4.-SL/2.)/((1.-E*SIN(SL))/(1.+E*SIN(SL)))**(E/2.)
          MC=COS(SL)/SQRT(1.0-E2*(SIN(SL)**2))
          RHO=RE*MC*T/TC
       END IF
       Y=-RHO*SGN*COS(SGN*RLON)*1000.0
       X= RHO*SGN*SIN(SGN*RLON)*1000.0
       ! RMG: added the *1000 to go from km to m.
    ELSE
       X=0.0
       Y=0.0
    END IF

  END SUBROUTINE FISOC_MAPLL
  
END MODULE FISOC_utils_MOD
