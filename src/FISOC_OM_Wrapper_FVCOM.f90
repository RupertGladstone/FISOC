
!
! This is the ocean model spcific code for FISOC.  The main purpose is to transfer information 
! between ESMF structures and the OM's internal structures.
!

MODULE FISOC_OM_Wrapper
  
  USE ESMF
  
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  
  !FVCOM specific:
  USE Mod_driver
  USE NestingFVCOMMod, only : NestingFVCOM_register

***use something that gives access to the fvcom grid info and the variables we need 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_OM_Wrapper_Init_Phase1,  FISOC_OM_Wrapper_Init_Phase2,  &
       FISOC_OM_Wrapper_Run, FISOC_OM_Wrapper_Finalize

  ! Component, and State
  type(ESMF_GridComp) :: nestingFVCOMComp  ! the nesting FVCOM Component
  type(ESMF_State)    :: nestingFVCOMState ! the nesting FVCOM State
  

***need to check on timing stuff, how we modify this...  

  ! Clock, TimeInterval, and Times
  type(ESMF_Clock)        :: clock
  type(ESMF_TimeInterval) :: timeStep
  type(ESMF_Time)         :: startTime
  type(ESMF_Time)         :: stopTime
  
  ! Namelist and related variables
  integer :: fileunit
  
  ! Return codes for error checks
  integer :: rc, urc, petCount, pet_id, i
  
  character(len=80) ::CASENAME  
  
  CALL GET_CASENAME(CASENAME)
  
CONTAINS
  
  !--------------------------------------------------------------------------------------
  ! The first phase of initialisation is mainly to initialise the ocean model, and access 
  ! grid and variable initial information.
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase1(FISOC_config,vm,OM_ExpFB,OM_grid,rc)
    
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)              :: vm ! ESMF virtual machine (parallel context)
    TYPE(ESMF_grid),INTENT(OUT)           :: OM_grid
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ExpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: localPet, petCount
    INTEGER                               :: mpic
    CHARACTER(len=ESMF_MAXSTR)            :: label
    TYPE(ESMF_staggerLoc),ALLOCATABLE     :: OM_ReqVars_stagger(:)
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: OM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: OM_configFile, OM_stdoutFile
    LOGICAL                               :: verbose_coupling, first

    first = .TRUE.


    CALL ESMF_VMGet(vm, petCount=petCount, localPet=localPet, rc=rc)
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


    !
    ! Read in input file  TODO: remove this to subroutine in nesting code or FISOC_config file (duplication)
    !
    call ESMF_UtilIOUnitGet(unit=fileunit, rc=rc) ! get an available Fortran unit number
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    open(fileunit, status="old", file=TRIM(CASENAME)//"_input", &
      action="read", iostat=rc)
    if (rc .ne. 0) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN, &
        msg="Failed to open namelist file 'CASENAME_input'", &
        line=__LINE__, &
        file=__FILE__)
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif
    read(fileunit, domain,   end=20)
    read(fileunit, input,    end=20)
    read(fileunit, np_model, end=20)
 20 continue
    close(fileunit)

    !
    ! Check processor numbers according to input file.
    ! (np_total, np_domainX both read in from input file)
    ! TODO: move to subroutine (duplication with nesting code) 
    !
    if (petCount /= np_total) then
      if(pet_id == 0) write(*,*) petCount, np_total
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_BAD, &
        msg="The np_total is not equal to -np on the command.", &
        line=__LINE__, &
        file=__FILE__)
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    if (np_total /= np_domain1 + np_domain2) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_BAD, &
        msg="The np_total should equal to np_domain1 + np_domain2.", &
        line=__LINE__, &
        file=__FILE__)
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if  

    allocate(petlist(np_total))
    allocate(petlist1(np_domain1))
    allocate(petlist2(np_domain2))
    do i=1,np_total
      petlist(i)=i-1
      if (i<=np_domain1) then
        petlist1(i)=i-1
      else
        petlist2(i-np_domain1)=i-1
      end if
    end do

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"*******************************************************************************"
       PRINT*,"**********       OM wrapper.  Init phase 1 method.         ********************"
       PRINT*,"********** Creating and registering ESMF FVCOM components. ********************"
       PRINT*,"*******************************************************************************"
       PRINT*,""
       IF (pet_id == 0) then
         PRINT*,"Processor decomposition:"
         WRITE(*,'(A,100I5)') "Domain 1 uses : ", petlist1
         WRITE(*,'(A,100I5)') "Domain 2 uses : ", petlist2
         WRITE(*,'(A,100I5)') "Coupler  uses : ", petlist
       END IF
    END IF

    ! Create the top level Gridded Component and register services.
    nestingFVCOMComp = ESMF_GridCompCreate(name="Nesting FVCOM Application", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_GridCompSetServices(nestingFVCOMComp, NestingFVCOM_register, &
      userRc=urc, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    

TODO: sort out clocks.  Allow nesting code its own clocks, but check for clashes with FISOC clocks.
Overwrite end time according to FISOC_OM_DT

    IF (localPet.EQ.0) THEN
       WRITE (OM_outputUnit,*) 'FISOC is about to call FVCOM init method.'
    END IF
    msg = "Calling FVCOM initialisation"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ***    CALL ROMS_initialize(first,mpic,OM_configFile)
    msg = "Completed FVCOM initialisation"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)
    
    IF (localPet.EQ.0) THEN
      WRITE (OM_outputUnit,*) 'FISOC has just called FVCOM init method.'
    END IF
    
    ! extract a list of required ocean variables from the configuration object
    label = 'FISOC_OM_ReqVars:' ! the FISOC names for the vars
    CALL FISOC_getListFromConfig(FISOC_config, label, OM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_populateFieldBundle(OM_ReqVarList,OM_ExpFB,OM_grid, &
         init_value=FISOC_missingData,                             &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,ignoreAveragesOpt=.TRUE.,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    RETURN
    
101 msg = "OM failed to open stdoutFile"
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
    INTEGER   :: localpet
    
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
       CALL CavityReset(OM_ImpFB,FISOC_config,localPet,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) & 
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,ignoreAveragesOpt=.TRUE.,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) & 
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase2
  
  
  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB,OM_ImpFB,rc_local)
    
    use mod_grid , only : GRID

    TYPE(ESMF_config),INTENT(INOUT)                :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL  :: OM_ExpFB, OM_ImpFB 
    TYPE(ESMF_VM),INTENT(IN)                       :: vm
    INTEGER,INTENT(OUT),OPTIONAL                   :: rc_local

    INTEGER                    :: localPet, rc
    LOGICAL                    :: verbose_coupling
    TYPE(ESMF_field)           :: ISM_dTdz_l0,ISM_z_l0, OM_dBdt_l0
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:,:), ISM_z_l0_ptr(:,:), OM_dBdt_l0_ptr(:,:)
    INTEGER                    :: OM_dt_sec
    REAL(ESMF_KIND_R8)         :: OM_dt_sec_float


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

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_dt_sec, label='OM_dt_sec:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    OM_dt_sec_float = REAL(OM_dt_sec,ESMF_KIND_R8)

    IF (localPet.EQ.0) THEN
      WRITE (OM_outputUnit,*) 'FISOC is about to call FVCOM run method, period (sec): ',OM_dt_sec_float
      IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
        msg = "Calling FVCOM run method now."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
             line=__LINE__, file=__FILE__, rc=rc)
      END IF
    END IF
    CALL ESMF_VMBarrier(vm, rc=rc)
    CALL FVCOM_run(OM_dt_sec_float)
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (localPet.EQ.0) THEN
       WRITE (OM_outputUnit,*) 'FISOC has just called FVCOM run method.'
    END IF

    IF (exit_flag.NE.NoError) THEN
       WRITE (msg, "(A,I0,A)") "ERROR: FVCOM has returned non-safe exit_flag=", &
            exit_flag,", see FVCOM mod_scalars.f90 for exit flag meanings."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       RETURN
    END IF
    
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
    CALL FVCOM_finalize
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
  ! Use the cavity from the ISM first stage initialisation to set the OM cavity
  !--------------------------------------------------------------------------------------
  SUBROUTINE CavityReset(OM_ImpFB,FISOC_config,localPet,rc)

    USE mod_iceshelfvar, ONLY : ICESHELFVAR
    USE mod_param, ONLY : BOUNDS, Ngrids
    USE mod_grid , ONLY : GRID
    USE mod_stepping, ONLY : nnew, nstp
    USE set_depth_mod, ONLY : set_depth

    INTEGER,INTENT(IN)                    :: localPet
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ImpFB 
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config

    TYPE(ESMF_FIELD)                      :: ISM_z_l0
    REAL(ESMF_KIND_R8),POINTER            :: ISM_z_l0_ptr(:,:)
    INTEGER                               :: ii, jj
    INTEGER                               :: JstrR, JendR, IstrR, IendR
    INTEGER                               :: Jstr, Jend, Istr, Iend
    REAL(ESMF_KIND_R8)                    :: OM_WCmin
# ifdef FVCOM_DSDT
    TYPE(ESMF_FIELD)                      :: ISM_z_lts
    REAL(ESMF_KIND_R8),POINTER            :: ISM_z_lts_ptr(:,:)
# endif

    rc = ESMF_FAILURE

    IF (Ngrids.GT.1) THEN
       msg = "ERROR: FVCOM has nested grids, FISOC cannot yet handle this"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_WCmin, 'OM_WCmin',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldname='ISM_z_l0', field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(ISM_z_l0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

# ifdef FVCOM_DSDT
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldname='ISM_z_lts', field=ISM_z_lts, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(ISM_z_lts, farrayPtr=ISM_z_lts_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif

    IstrR=BOUNDS(Ngrids)%IstrR(localPet)
    IendR=BOUNDS(Ngrids)%IendR(localPet)
    JstrR=BOUNDS(Ngrids)%JstrR(localPet)
    JendR=BOUNDS(Ngrids)%JendR(localPet)

    Istr =BOUNDS(Ngrids)%Istr (localPet)
    Iend =BOUNDS(Ngrids)%Iend (localPet)
    Jstr =BOUNDS(Ngrids)%Jstr (localPet) 
    Jend =BOUNDS(Ngrids)%Jend (localPet)

!    CALL cp2bdry(ptr,JstrR,JendR,IstrR,IendR)
    DO jj = JstrR, JendR
       DO ii = IstrR, IendR
# ifdef FVCOM_DSDT
          GRID(1) % sice (ii, jj) = ISM_z_lts_ptr (ii, jj)
# endif
          IF ( ISM_z_l0_ptr(ii,jj) .GT. OM_WCmin-GRID(Ngrids)%h(ii,jj) ) THEN
             ICESHELFVAR(1) % iceshelf_draft(ii, jj, nstp) = ISM_z_l0_ptr (ii, jj)
             ICESHELFVAR(1) % iceshelf_draft(ii, jj, nnew) = ISM_z_l0_ptr (ii, jj)
             GRID(1) % zice (ii, jj) = ISM_z_l0_ptr (ii, jj)
          ELSE
             ICESHELFVAR(1) % iceshelf_draft(ii, jj, nstp) = OM_WCmin-GRID(Ngrids)%h(ii,jj)
             ICESHELFVAR(1) % iceshelf_draft(ii, jj, nnew) = OM_WCmin-GRID(Ngrids)%h(ii,jj)
             GRID(1) % zice (ii, jj) = OM_WCmin-GRID(Ngrids)%h(ii,jj)
          END IF
       END DO
    END DO

!    CALL set_depth(Ngrids,localPet)
!check draft is not below bedrock

# ifdef FVCOM_DSDT
    IF (ASSOCIATED(ISM_z_lts_ptr)) THEN
       NULLIFY(ISM_z_lts_ptr)
    END IF
# endif
    IF (ASSOCIATED(ISM_z_l0_ptr)) THEN
       NULLIFY(ISM_z_l0_ptr)
    END IF
    
    msg = "Ocean cavity reset.  Masks may be needed (NYI)"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
         line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE CavityReset


  !--------------------------------------------------------------------------------------
  ! update the fields in the ocean export field bundle from the OM
  !--------------------------------------------------------------------------------------
  SUBROUTINE getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,ignoreAveragesOpt,rc)

    USE mod_iceshelfvar, ONLY : ICESHELFVAR
    USE mod_param, ONLY       : BOUNDS, Ngrids
    USE mod_stepping, ONLY    : nnew
    USE mod_grid , ONLY       : GRID
    USE mod_average, ONLY     : AVERAGE

    IMPLICIT NONE

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm
    LOGICAL,INTENT(IN),OPTIONAL           :: ignoreAveragesOpt

    LOGICAL                               :: ignoreAverages
    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:,:)
    INTEGER                               :: IstrR, IendR, JstrR, JendR ! tile start and end coords
    INTEGER                               :: ii, jj, nn
!    INTEGER                               :: LBi, UBi, LBj, UBj ! tile start and end coords including halo

    rc = ESMF_FAILURE

    IF (PRESENT(ignoreAveragesOpt)) THEN
      ignoreAverages = ignoreAveragesOpt
    ELSE
      ignoreAverages = .FALSE.
    END IF

    IF (Ngrids .ne. 1) THEN
       msg = "ERROR: not expecting multiple FVCOM grids"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
  
    ! get the tile position from the OM (a tile is a rectangular domain decomposition element)
    IstrR=BOUNDS(Ngrids)%IstrR(localPet)
    IendR=BOUNDS(Ngrids)%IendR(localPet)
    JstrR=BOUNDS(Ngrids)%JstrR(localPet)
    JendR=BOUNDS(Ngrids)%JendR(localPet)
        
    ! get a list of fields and their names form the OM export field bundle
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
         
       CASE ('OM_dBdt_l0')
         DO jj = JstrR, JendR
           DO ii = IstrR, IendR
# if defined(FVCOM_AVERAGES)
             IF (ignoreAverages) THEN
               ptr(ii,jj) = ICESHELFVAR(1) % m(ii,jj)
             ELSE
               ptr(ii,jj) = AVERAGE(1) % avgismr(ii,jj)
             END IF
# else
             ptr(ii,jj) = ICESHELFVAR(1) % m(ii,jj)
# endif
           END DO
         END DO
!         CALL cp2bdry(ptr,JstrR,JendR,IstrR,IendR)
         
       CASE ('OM_temperature_l0')
         DO jj = JstrR, JendR
           DO ii = IstrR, IendR
# if defined(FVCOM_AVERAGES)
             IF (ignoreAverages) THEN
               ptr(ii,jj) = ICESHELFVAR(1) % Tb(ii,jj)
             ELSE
               ptr(ii,jj) = AVERAGE(1) % avgisTb(ii,jj)
             END IF
# else
             ptr(ii,jj) = ICESHELFVAR(1) % Tb(ii,jj)
# endif
           END DO
         END DO
         
       CASE ('OM_z_l0')
         DO jj = JstrR, JendR
           DO ii = IstrR, IendR
             ptr(ii,jj) = ICESHELFVAR(1) % iceshelf_draft(ii,jj,nnew(1))
           END DO
         END DO
         
       CASE ('OM_z_lts')
# if defined FVCOM_DSDT
         DO jj = JstrR, JendR
           DO ii = IstrR, IendR
             ptr(ii,jj) = GRID(1) % sice(ii,jj)
           END DO
         END DO
# else
         msg = "ERROR: cpp FVCOM_DSDT must be set in order to access OM_z_lts"
         CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
              line=__LINE__, file=__FILE__, rc=rc)
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
         
       CASE ('OM_bed')
         DO jj = JstrR, JendR
           DO ii = IstrR, IendR
             ptr(ii,jj) = -GRID(1) % h(ii,jj)
           END DO
         END DO
         
       CASE DEFAULT
         msg = "ERROR: unknown variable"
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
  

  !--------------------------------------------------------------------------------------
  ! Copy values to boundary cells from neihgbouring interior cells.
  SUBROUTINE cp2bdry(ptr,JstrR,JendR,IstrR,IendR)

    USE mod_param, ONLY : Lm, Mm, Ngrids

    REAL(ESMF_KIND_R8),POINTER,INTENT(IN) :: ptr(:,:)
    INTEGER,INTENT(IN)                    :: IstrR, IendR, JstrR, JendR

    INTEGER                               :: ii, jj, rc

    IF (Ngrids .ne. 1) THEN
       msg = "ERROR: not expecting multiple FVCOM grids"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    DO jj = JstrR, JendR
       IF (jj .EQ. 0) THEN
          ptr(:,jj) = ptr(:,jj+1)
       ELSE IF (jj .EQ. Mm(1)+1) THEN
          ptr(:,jj) = ptr(:,jj-1)
       END IF
    END DO
    
    DO ii = IstrR, IendR
       IF (ii .EQ. 0) THEN
          ptr(ii,:) = ptr(ii+1,:)
       ELSE IF (ii .EQ. Lm(1)+1) THEN
          ptr(ii,:) = ptr(ii-1,:)
       END IF
    END DO

  END SUBROUTINE cp2bdry


  !--------------------------------------------------------------------------------------
  SUBROUTINE sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc)

    USE mod_iceshelfvar, ONLY : ICESHELFVAR
    USE mod_param, ONLY       : BOUNDS, Ngrids
    USE mod_stepping, ONLY    : nnew, nstp
    USE mod_grid , ONLY       : GRID

    TYPE(ESMF_fieldBundle),INTENT(INOUT)     :: OM_ImpFB 
    TYPE(ESMF_config),INTENT(INOUT)          :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)                 :: vm
    INTEGER,INTENT(OUT),OPTIONAL             :: rc

    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:,:)
    INTEGER                               :: IstrR, IendR, JstrR, JendR ! tile start and end coords
    INTEGER                               :: ii, jj, nn

    rc = ESMF_FAILURE

    IF (Ngrids .ne. 1) THEN
       msg = "ERROR: not expecting multiple FVCOM grids"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! get the tile position from the OM (a tile is a rectangular domain decomposition element)
    IstrR=BOUNDS(Ngrids)%IstrR(localPet)
    IendR=BOUNDS(Ngrids)%IendR(localPet)
    JstrR=BOUNDS(Ngrids)%JstrR(localPet)
    JendR=BOUNDS(Ngrids)%JendR(localPet)
    
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
       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (FISOC_ISM2OM(fieldName,FISOC_config,rc=rc)) THEN

          SELECT CASE (TRIM(ADJUSTL(fieldName)))
             
          CASE ('ISM_dddt')
# ifdef FVCOM_DDDT
             DO jj = JstrR, JendR
                DO ii = IstrR, IendR
                   ICESHELFVAR(1) % iceshelf_dddt(ii,jj) = ptr(ii,jj)
                END DO
             END DO
# else
             msg = "Trying to pass DDDT to FVCOM but incompatible cpp"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
                  line=__LINE__, file=__FILE__, rc=rc)          
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
 
          CASE ('ISM_dsdt')
# ifdef FVCOM_DSDT
             DO jj = JstrR, JendR
                DO ii = IstrR, IendR
                   ICESHELFVAR(1) % iceshelf_dsdt(ii,jj) = ptr(ii,jj)
                END DO
             END DO
# else
             msg = "Trying to pass DSDT to FVCOM but incompatible cpp"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
                  line=__LINE__, file=__FILE__, rc=rc)          
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
             
          CASE ('ISM_z_l0','ISM_z_l0_linterp')
# ifdef FVCOM_DRAFT
             ! A note on the OM ICESHELFVAR(1) % iceshelf_draft:
             ! iceshelf_draft(:,:,nstp) is the previous draft and iceshelf_draft(:,:,nnew) is 
             ! the current draft.  We only update the new draft from the ISM.
             ! The OM var zice will be set internally by the OM based on the iceshelf_draft.
             CALL cp2bdry(ptr,JstrR,JendR,IstrR,IendR)
             DO jj = JstrR, JendR
                DO ii = IstrR, IendR
                   ICESHELFVAR(1) % iceshelf_draft(ii,jj,nnew) = ptr(ii,jj)
!                   GRID(1) % zice (ii, jj) = ptr (ii, jj)
                END DO
             END DO
# else
             msg = "Trying to pass draft to FVCOM but incompatible cpp"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
                  line=__LINE__, file=__FILE__, rc=rc)          
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
             
          CASE('ISM_temperature_l0', 'ISM_temperature_l1', 'ISM_z_l1', 'ISM_velocity_l0', 'ISM_z_l0_previous', 'ISM_dTdz_l0')
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
    
  END SUBROUTINE sendFieldDataToOM


  ! Extract the FVCOM mesh information and use it to create an ESMF_mesh object
  SUBROUTINE FVCOM2ESMF_mesh(FISOC_config,OM_mesh,vm,rc=rc)

    print*,"Tore and Qin to write this"

  END SUBROUTINE FVCOM2ESMF_mesh

END MODULE FISOC_OM_Wrapper
