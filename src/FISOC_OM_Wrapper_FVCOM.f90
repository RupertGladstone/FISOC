
!
! This is the ocean model spcific code for FISOC.  The main purpose is to transfer information 
! between ESMF structures and the OM's internal structures.
!

MODULE FISOC_OM_Wrapper
  
  USE ESMF
  
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  
  !FVCOM specific:
  USE mod_ocean_control
  USE ALL_VARS
  USE LIMS
  USE mod_par
  USE mod_wd
  USE mod_isf
 !  USE mod_main

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_OM_Wrapper_Init_Phase1,  FISOC_OM_Wrapper_Init_Phase2,  &
       FISOC_OM_Wrapper_Run, FISOC_OM_Wrapper_Finalize

  
  ! TODO: don't use mpic for anything other than passing to FVCOM.  All mpi operations 
  ! should be through the ESMF VM commands.
  INTEGER, SAVE                :: mpic
  INTEGER, ALLOCATABLE, SAVE   :: ownedNodeIDs(:) ! subset of local nodes owned by current partition

  ! Route handle for mapping from ESMF fields to FVCOM fields (one to many mapping)
  TYPE(ESMF_RouteHandle),SAVE  :: RH_ESMF2FVCOM

  ! nodal distgrid, an ESMF object holding information about the distribution of 
  ! FVCOM nodes across partitions, needed for the one to many mapping.
  TYPE(ESMF_distgrid),SAVE     :: distgridFVCOM

CONTAINS
  
  !--------------------------------------------------------------------------------------
  ! The first phase of initialisation is mainly to initialise the ocean model, and access 
  ! grid and variable initial information.
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase1(FISOC_config,vm,OM_ExpFB,OM_mesh,rc)
    
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)              :: vm ! ESMF virtual machine (parallel context)
    TYPE(ESMF_mesh),INTENT(OUT)           :: OM_mesh
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ExpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: localPet, petCount
    CHARACTER(len=ESMF_MAXSTR)            :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: OM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: OM_configFile, OM_stdoutFile
    LOGICAL                               :: verbose_coupling, first

    first = .TRUE.

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_stdoutFile, label='OM_stdoutFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    WRITE (OM_stdoutFile, "(a,I0)") TRIM(OM_stdoutFile), localPet
    OPEN(unit=OM_outputUnit, file=OM_stdoutFile, STATUS='REPLACE', ERR=101)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Initialising FVCOM"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_configFile, label='OM_configFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"*******************************************************************************"
       PRINT*,"**********       OM wrapper.  Init phase 1 method.         ********************"
       PRINT*,"********** Creating and registering ESMF FVCOM components. ********************"
       PRINT*,"*******************************************************************************"
       PRINT*,""
    END IF

    ! We will pass an MPI communicator to FVCOM
    CALL FISOC_VM_MPI_Comm_dup(vm,mpic,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (localPet.EQ.0) THEN
       WRITE (OM_outputUnit,*) 'FISOC is about to call FVCOM init method.'
    END IF
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FVCOM_initialize(casename_opt=OM_configFile, MPI_COMM_opt=mpic)
    IF (localPet.EQ.0) THEN
      WRITE (OM_outputUnit,*) 'FISOC has just called FVCOM init method.'
    END IF
    
    CALL FVCOM2ESMF_mesh(FISOC_config,OM_mesh,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    

    ! extract a list of required ocean variables from the configuration object
    label = 'FISOC_OM_ReqVars:' ! the FISOC names for the vars
    CALL FISOC_getListFromConfig(FISOC_config, label, OM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_populateFieldBundle(OM_ReqVarList,OM_ExpFB,OM_mesh, &
         init_value=FISOC_missingData,                             &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
 
    CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    RETURN
    
101 msg = "OM failed to open stdoutFile "//OM_stdoutFile
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase1
  
  
  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase2(FISOC_config,vm,OM_ImpFB,OM_ExpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ImpFB, OM_ExpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm
    
    LOGICAL   :: verbose_coupling, OM_initCavityFromISM, ISM2OM_init_vars
    INTEGER   :: localpet, urc

    rc = ESMF_FAILURE
    
    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********      OM wrapper.  Init phase 2 method.        *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we have access to the initialised ISM fields, just in case the OM needs "
       PRINT*,"to know about these in order to complete its initialisation."
       PRINT*,""
    END IF

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_initCavityFromISM, label='OM_initCavityFromISM:', rc=rc)
    IF  (rc.EQ.ESMF_RC_NOT_FOUND) THEN
       OM_initCavityFromISM = .FALSE.
    ELSE
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM2OM_init_vars, 'ISM2OM_init_vars',rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (ISM2OM_init_vars) THEN
       CALL sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    IF (OM_initCavityFromISM) THEN
      CALL CavityReset(OM_ImpFB,FISOC_config,rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) & 
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) & 
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase2
  
  
  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB,OM_ImpFB,rc_local)
    
    TYPE(ESMF_config),INTENT(INOUT)                :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL  :: OM_ExpFB, OM_ImpFB 
    TYPE(ESMF_VM),INTENT(IN)                       :: vm
    INTEGER,INTENT(OUT),OPTIONAL                   :: rc_local

    INTEGER                    :: localPet, rc
    LOGICAL                    :: verbose_coupling
    TYPE(ESMF_field)           :: ISM_dTdz_l0,ISM_z_l0, OM_bmb
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:,:), ISM_z_l0_ptr(:,:), OM_bmb_ptr(:,:)
    INTEGER                    :: OM_dt_sec
    REAL(ESMF_KIND_R8)         :: OM_dt_sec_float
    TYPE(ESMF_TimeInterval)    :: OM_dt
    TYPE(ESMF_time)            :: interval_startTime, interval_endTime
    CHARACTER(len=ESMF_MAXSTR) :: interval_startTime_char, interval_endTime_char

    rc_local = ESMF_FAILURE
    
    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (PRESENT(OM_ImpFB)) THEN       
       CALL sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! Get OM time interval (how long to run OM for) in seconds...
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_dt_sec, label='OM_dt_sec:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    OM_dt_sec_float = REAL(OM_dt_sec,ESMF_KIND_R8)

    ! ...convert to ESMF type...
    CALL ESMF_TimeIntervalSet(OM_dt, s=OM_dt_sec, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... use to create start and end times...
    interval_startTime = FISOC_time
    interval_endTime   = FISOC_time + OM_dt

    ! ... and convert these to character strings.
    CALL ESMF_TimeGet(interval_startTime, timeStringISOFrac=interval_startTime_char, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_TimeGet(interval_endTime, timeStringISOFrac=interval_endTime_char, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF (localPet.EQ.0) THEN
      WRITE (OM_outputUnit,*) 'FISOC is about to call FVCOM run method, period (sec): ',OM_dt_sec_float
      IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
        msg = "Calling FVCOM run method now."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
             line=__LINE__, file=__FILE__, rc=rc)
      END IF
    END IF

    CALL ESMF_VMBarrier(vm, rc=rc)
    CALL FVCOM_run(START_DATE_opt=interval_startTime_char,END_DATE_opt=interval_endTime_char)
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (localPet.EQ.0) THEN
      WRITE (OM_outputUnit,*) 'FISOC has just called FVCOM run method.'
    END IF
    
! TODO: is there an accessible exit flag or similar for FVCOM?
!    IF (exit_flag.NE.NoError) THEN
!      WRITE (msg, "(A,I0,A)") "ERROR: FVCOM has returned non-safe exit_flag=", &
!           exit_flag,", see FVCOM mod_scalars.f90 for exit flag meanings."
!      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
!           line=__LINE__, file=__FILE__, rc=rc)
!      RETURN
!    END IF
    
    IF (PRESENT(OM_ExpFB)) THEN
      CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************         OM wrapper.  Run method.           **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""

       IF ((PRESENT(OM_ExpFB)).AND.(.NOT.(PRESENT(OM_ImpFB)))) THEN
          PRINT*,"We have no new inputs for the OM from the ISM.  We need to call the OM "
          PRINT*,"and record its output in the OM export field bundle."
       END IF

       IF ((PRESENT(OM_ExpFB)).AND.(PRESENT(OM_ImpFB))) THEN
          PRINT*,"We have new inputs for the OM from the ISM in the OM import field bundle. "
          PRINT*,"We need to send these inputs to the OM, run one timestep, and record the OM "
          PRINT*,"outputs. "
       END IF
       
       IF ((.NOT.(PRESENT(OM_ExpFB))).AND.(PRESENT(OM_ImpFB))) THEN
          PRINT*,"We have new inputs for the OM from the ISM in the OM import field bundle. "
          PRINT*,"We need to send these inputs to the OM and run one timestep. We do not "
          PRINT*,"need to collect OM outputs.  Just run the OM one timestep."
       END IF

       IF ((.NOT.(PRESENT(OM_ExpFB))).AND.(.NOT.(PRESENT(OM_ImpFB)))) THEN
          PRINT*,"We have no new inputs for the OM from the ISM, and we do not need to "
          PRINT*,"collect OM outputs.  Just run the OM one timestep."
       END IF

    END IF

    rc_local = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_OM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Finalize(FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)           :: vm
    INTEGER,INTENT(OUT),OPTIONAL       :: rc

    INTEGER                            :: localPet
    LOGICAL                            :: verbose_coupling

    rc = ESMF_FAILURE


    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************       OM wrapper.  Finalise method.        **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"OM finalise method."
    END IF

    IF (localPet.EQ.0) THEN
       WRITE (OM_outputUnit,*) 'FISOC is about to call FVCOM finalize.'
    END IF
    CALL FVCOM_finalize()
    IF (localPet.EQ.0) THEN
       WRITE (OM_outputUnit,*) 'FISOC has just called FVCOM finalize.'
    END IF

!    CLOSE(unit=OM_outputUnit, ERR=102)
    
    rc = ESMF_SUCCESS
    
    RETURN
    
102 msg = "OM failed to close stdoutFile"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_OM_Wrapper_Finalize
  

  !--------------------------------------------------------------------------------------
  ! update the fields in the ocean export field bundle from the OM
  !--------------------------------------------------------------------------------------
  SUBROUTINE getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm

    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)
    INTEGER                               :: IstrR, IendR, JstrR, JendR ! tile start and end coords
    INTEGER                               :: ii, jj, nn
!    INTEGER                               :: LBi, UBi, LBj, UBj ! tile start and end coords including halo

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
    ! get a list of fields and their names from the OM export field bundle
    fieldCount = 0
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    fieldLoop: DO nn = 1,fieldCount
       
       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       ptr = FISOC_missingData
       
       SELECT CASE (TRIM(ADJUSTL(fieldName)))
         
       CASE ('OM_z_l0')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = -zisf(ownedNodeIds(ii))
         END DO
         
       CASE ('OM_bmb')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = melt_avg(ownedNodeIds(ii))
         END DO
         
       CASE DEFAULT
         msg = "ERROR: unknown variable: "//TRIM(ADJUSTL(fieldName))
         CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
              line=__LINE__, file=__FILE__, rc=rc)
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
         
       END SELECT
       
       IF (ASSOCIATED(ptr)) THEN
         NULLIFY(ptr)
       END IF
       
    END DO fieldLoop
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE GetFieldDataFromOM
  

  SUBROUTINE CavityReset(OM_ImpFB,FISOC_config,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ImpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)
    INTEGER                               :: IstrR, IendR, JstrR, JendR ! tile start and end coords
    INTEGER                               :: ii, jj, nn
!    INTEGER                               :: LBi, UBi, LBj, UBj ! tile start and end coords including halo

    rc = ESMF_FAILURE
        
    ! get a list of fields and their names from the OM export field bundle
     fieldCount = 0
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    fieldLoop: DO nn = 1,fieldCount
       
       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (TRIM(ADJUSTL(fieldName)).EQ.'ISM_z_l0') THEN

          CALL FISOC_ArrayRedistFromField(RH_ESMF2FVCOM,fieldList(nn),distgridFVCOM,ptr)

          msg = "FVCOM: use ISM_z_l0 to set ZISF (reset cavity)"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          DO ii  =  1, MT
             ZISF(ii)     =   -ptr(ii)
          END DO
          CALL     ISF_JUDGE
          CALL     WET_JUDGE
          D    =  H + EL-ZISF
          DT   =  H + ET-ZISF
          WHERE (D<=MIN_DEPTH)  D  =  MIN_DEPTH
          WHERE (DT<=MIN_DEPTH) DT =  MIN_DEPTH
          DTFA =  D
          CALL N2E2D(D,D1)
          CALL N2E2D(DT,DT1)
          CALL N2E2D(ZISF,ZISF1)       
       END IF
       
       IF (ASSOCIATED(ptr)) THEN
          NULLIFY(ptr)
       END IF
       
    END DO fieldLoop
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE CavityReset



  !--------------------------------------------------------------------------------------
  ! update the fields in the OM from the ocean import field bundle
  !--------------------------------------------------------------------------------------
  SUBROUTINE sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ImpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm

    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)
    INTEGER                               :: IstrR, IendR, JstrR, JendR ! tile start and end coords
    INTEGER                               :: ii, jj, nn

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
    ! get a list of fields and their names from the OM export field bundle
    fieldCount = 0
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!TODO: remove or move print statements to ESMF logs
!    PRINT*,'fieldCount is', fieldCount 

    fieldLoop: DO nn = 1,fieldCount
      
      CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      IF (FISOC_ISM2OM(fieldName,FISOC_config,rc=rc)) THEN

        CALL FISOC_ArrayRedistFromField(RH_ESMF2FVCOM,fieldList(nn),distgridFVCOM,ptr)
      
        SELECT CASE (TRIM(ADJUSTL(fieldName)))
        
        CASE ('ISM_dddt')
          DO ii = 1,SIZE(ptr)
            ISF_DDDT(ii) = ptr(ii)
         END DO
!         WRITE(msg,*) "ISF_DDDT max ", MAXVAL(ISF_DDDT) ," min ", MINVAL(ISF_DDDT)
!         CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
!              line=__LINE__, file=__FILE__, rc=rc)
!         WRITE(msg,*) "ISF ptr max ", MAXVAL(ptr) ," min ", MINVAL(ptr)
!         CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
!              line=__LINE__, file=__FILE__, rc=rc)
          
        CASE ('ISM_dsdt')
          DO ii = 1,SIZE(ptr)
            ISF_DSDT(ii) = ptr(ii)
          END DO
          
        CASE ('ISM_thick','ISM_z_l0','ISM_z_l0_previous','ISM_z_lts','ISM_z_lts_previous')
          msg = "WARNING: ignored variable: "//TRIM(ADJUSTL(fieldName))
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__, rc=rc)          
          
        CASE DEFAULT
          msg = "ERROR: unknown variable: "//TRIM(ADJUSTL(fieldName))
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          
        END SELECT
      
        IF (ASSOCIATED(ptr)) THEN
          NULLIFY(ptr)
        END IF

      END IF
        
    END DO fieldLoop
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE SendFieldDataToOM
  

  !--------------------------------------------------------------------------------
  ! Extract the FVCOM mesh information and use it to create an ESMF_mesh object
  SUBROUTINE FVCOM2ESMF_mesh(FISOC_config,ESMF_FVCOMmesh,vm,rc)
  
    TYPE(ESMF_config),INTENT(INOUT)  :: FISOC_config
    TYPE(ESMF_mesh),INTENT(INOUT)    :: ESMF_FVCOMMesh
    TYPE(ESMF_VM),INTENT(IN)         :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    INTEGER                          :: ii, nn, IERR
    CHARACTER(len=ESMF_MAXSTR)       :: subroutineName = "FVCOM2ESMF_mesh", OM_coords
    INTEGER,ALLOCATABLE              :: elemTypes(:), elemIds(:),elemConn(:), localNodeIDs(:)
    INTEGER,ALLOCATABLE              :: nodeIds(:),nodeOwners(:),nodeOwnersGL(:),nodeOwnersGL_recv(:)
    INTEGER,ALLOCATABLE              :: nodeMask(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:) 
    INTEGER                          :: localPet, petCount 
    INTEGER                          :: FVCOM_numNodes, FVCOM_numElems, OM_OPEN_OCEAN
    LOGICAL                          :: verbose_coupling
	
    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_coords, label='OM_coords:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_OPEN_OCEAN, label='OM_OPEN_OCEAN:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "FVCOM mesh: creating ESMF mesh structure"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
  
    !-------------------------------------------------------------------!
    ! Use fvcom grid variables directly, may need to use MODULE ALL_VARS
    ! and MODULE LIMS from FVCOM library, see mod_main.F for details
    ! and module par for retreating  EGID NGID
    FVCOM_numNodes        = MT
    FVCOM_numElems        = NT
!TODO: just use FVCOM vars directly, no need for intermediate vars, i.e. remove FVCOM_numNodes etc

    ALLOCATE(localNodeIds(M))

    ALLOCATE(nodeIds(FVCOM_numNodes))
    ALLOCATE(nodeCoords(FVCOM_numNodes*2))
    ALLOCATE(nodeOwners(FVCOM_numNodes))
    ALLOCATE(nodeMask(FVCOM_numNodes))
    
    ALLOCATE(elemIds(FVCOM_numElems))
    ALLOCATE(elemConn(FVCOM_numElems*3))
    ALLOCATE(elemTypes(FVCOM_numElems))

    ALLOCATE(nodeOwnersGL(MGL)) 
    ALLOCATE(nodeOwnersGL_recv(MGL)) 

    ! Setup mask to avoid regridding ice variables in the open ocean.
    nodeMask = MASK_ICE
    WHERE (ISISFN.EQ.OM_OPEN_OCEAN) nodeMask = MASK_OPEN_OCEAN

    WRITE (msg, "(A,I5)") "ISISFN size:   ", SIZE(ISISFN)
    CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)       
    WRITE (msg, "(A,I5,I5)") "ISISFN maxmin: ", MAXVAL(ISISFN), MINVAL(ISISFN)
    CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)       
    WRITE (msg, "(A,I5)") "nodeMask size: ", SIZE(nodeMask)
    CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)       
    WRITE (msg, "(A,I5)") "NGID_X size:   ", SIZE(NGID_X)
    CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)       
    WRITE (msg, "(A,I5)") "Num nodes:     ", FVCOM_NumNodes
    CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)       


    ! Construct a global array of node owners in which an arbitrary decision is 
    ! taken about which partition boundary nodes should belong to.
    nodeOwnersGL = -1
    DO nn=1,M
      nodeOwnersGL(NGID(nn))=localPET
      localNodeIds(nn) = nn
    END DO
    CALL MPI_Allreduce(nodeOwnersGL,nodeOwnersGL_recv,MGL,MPI_INT,MPI_MAX,mpic,IERR)

    ! sanity check (we initialised owners to -1, but PET count starts at 0, so if it 
    ! works then all nodes should have been assigned an owner .GE. 0)
    IF (MINVAL(nodeOwnersGL_recv).LT.0) THEN
      msg = "ERROR: Some nodes not assigned owners"
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
           line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    !populate nodeOwners from nodeOwnersGL
    DO ii = 1, FVCOM_numNodes
      nodeOwners(ii) = nodeOwnersGL_recv(NGID_X(ii))
    END DO
    elemTypes         =  ESMF_MESHELEMTYPE_TRI

    ! loop over to get nodeIds nodeCoords
    DO ii = 1, FVCOM_numNodes

       nn = (ii-1)*2
       nodeIds(ii)      = NGID_X(ii)
       
       SELECT CASE (OM_coords)

       CASE("spherical","Spherical","SPHERICAL")
          CALL FISOC_MAPLL(nodeCoords(nn+1),nodeCoords(nn+2),LAT(ii),LON(ii))

       CASE DEFAULT
          nodeCoords(nn+1) = XM(ii)
          nodeCoords(nn+2) = YM(ii)

       END SELECT

       IF ( (NGID_X(ii).EQ.3785) .OR. &
            (NGID_X(ii).EQ.5276) .OR. &
            (NGID_X(ii).EQ.6063) .OR. &
            (NGID_X(ii).EQ.7575) .OR. &
            (NGID_X(ii).EQ.7943) ) THEN
        
          WRITE (msg, "(A,I5)") "MASK_ICE:        ",MASK_ICE
          CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)
          WRITE (msg, "(A,I5)") "MASK_OPEN_OCEAN: ",MASK_OPEN_OCEAN
          CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)
          WRITE (msg, "(A,I5)") "OM_OPEN_OCEAN:   ",OM_OPEN_OCEAN
          CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)
          
          WRITE (msg, "(A,I5)") "NGID_X:   ", NGID_X(ii)
          CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)
          WRITE (msg, "(A,I5)") "ISISFN:   ", ISISFN(ii)
          CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)
          WRITE (msg, "(A,I5)") "nodeMask: ", nodeMask(ii)
          CALL ESMF_LogWrite(msg,logmsgFlag=ESMF_LOGMSG_INFO,line=__LINE__, file=__FILE__, rc=rc)

       END IF

    END DO

    ! loop over to get elemConn
    DO ii = 1, FVCOM_numElems
       nn = (ii-1)*3
       elemIds(ii)    = EGID_X(ii)
       elemConn(nn+1) = NV(ii,4)
       elemConn(nn+2) = NV(ii,3)
       elemConn(nn+3) = NV(ii,2)
    END DO 

    ! get a list od node ids for the subset of nodes that are locally owned
    ! (to be used in mapping FVCOM variables to ESMF fields)
    CALL FISOC_locallyOwnedNodes(localPet,localNodeIds,nodeOwners,ownedNodeIDs)


    !----------------------------------------------------------------!       
    ! Create Mesh structure in 1 step
    ESMF_FVCOMMesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
         nodeIds=nodeIds, nodeCoords=nodeCoords,                   &
         nodeOwners=nodeOwners,                                    &
         nodeMask=nodeMask, elementIds=elemIds,                    &
         elementTypes=elemTypes, elementConn=elemConn,             &
         coordSys=ESMF_COORDSYS_CART,                              &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
	
    CALL FISOC_CreateOneToManyRouteHandle(ESMF_FVCOMMesh,NGID_X(1:),RH_ESMF2FVCOM,distgridFVCOM,vm)

    DEALLOCATE(nodeIds)
    DEALLOCATE(nodeCoords)
    DEALLOCATE(nodeOwners)
    DEALLOCATE(nodeMask)
    
    DEALLOCATE(elemIds)
    DEALLOCATE(elemConn)
    DEALLOCATE(elemTypes)	

    !TODO: put this in PET logs
    print*,"FINISHED FVCOM MESH CREATION"

  END SUBROUTINE FVCOM2ESMF_mesh


END MODULE FISOC_OM_Wrapper
