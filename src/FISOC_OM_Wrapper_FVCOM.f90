
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

  ! TODO figure out which modules we need to get the FVCOM variables, assuming that we 
  ! can't get what we want through the nesting state.

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_OM_Wrapper_Init_Phase1,  FISOC_OM_Wrapper_Init_Phase2,  &
       FISOC_OM_Wrapper_Run, FISOC_OM_Wrapper_Finalize

  
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
    INTEGER                               :: mpic
    CHARACTER(len=ESMF_MAXSTR)            :: label
    TYPE(ESMF_staggerLoc),ALLOCATABLE     :: OM_ReqVars_stagger(:)
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
    CALL FVCOM_initialize
!    CALL FVCOM_run
! TODO: call FVCOM init
    IF (localPet.EQ.0) THEN
      WRITE (OM_outputUnit,*) 'FISOC has just called FVCOM init method.'
    END IF
    
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
 
! TODO: rewrite this subroutine for FVCOM   
!    CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,ignoreAveragesOpt=.TRUE.,rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
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

! TODO: write send vars subroutine
!    IF (ISM2OM_init_vars) THEN
!       CALL sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    END IF

    IF (OM_initCavityFromISM) THEN
      msg = "Cavity reset NYI for FVCOM"
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
           line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

! TODO: allow cavity reset
!       CALL CavityReset(OM_ImpFB,FISOC_config,localPet,rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) & 
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    END IF

! TODO: write this subroutine
!    CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,ignoreAveragesOpt=.TRUE.,rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) & 
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase2
  
  
  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB,OM_ImpFB,rc_local)
    
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

!TODO: write this
!    IF (PRESENT(OM_ImpFB)) THEN       
!       CALL sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    END IF


    ! FVCOM needs date strings that look like this:
    ! START_DATE      = '2000-01-01 00:00:00',
    ! END_DATE        = '2000-01-01 00:01:00',

    ! Get OM time interval (how long to run OM for) in seconds...
    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_dt_sec, label='OM_dt_sec:', rc=rc)
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
    CALL ESMF_TimeGet(startTime, timeStringISOFrac=interval_startTime_char, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_TimeGet(endTime, timeStringISOFrac=interval_endTime_char, rc=rc)
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
    
! TODO: wirte this subroutine
!    IF (PRESENT(OM_ExpFB)) THEN
!       CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    END IF

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
!TODO: FVCOM finalize
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
  

  ! Extract the FVCOM mesh information and use it to create an ESMF_mesh object
  SUBROUTINE FVCOM2ESMF_mesh(FISOC_config,ESMF_FVCOMmesh,vm,rc)
  
    TYPE(ESMF_config),INTENT(INOUT)  :: FISOC_config
    TYPE(ESMF_mesh),INTENT(INOUT)    :: ESMF_FVCOMMesh
    TYPE(ESMF_VM),INTENT(IN)         :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    INTEGER                          :: ii, nodeIndex
    CHARACTER(len=ESMF_MAXSTR)       :: subroutineName = "FVCOM2ESMF_mesh"
    INTEGER,ALLOCATABLE              :: elemTypes(:), elemIds(:),elemConn(:)
    INTEGER,ALLOCATABLE              :: nodeIds(:),nodeOwners(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:) 
    INTEGER                          :: localPet, petCount 
    INTEGER                          :: FVCOM_numNodes, FVCOM_numElems 

	
    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
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
    FVCOM_numNodes        = MGL                            
    FVCOM_numElems        = NGL	
	
    ALLOCATE(nodeIds(FVCOM_numNodes))
    ALLOCATE(nodeCoords(FVCOM_numNodes*2))
    ALLOCATE(nodeOwners(FVCOM_numNodes))
    
    ALLOCATE(elemIds(FVCOM_numElems))
    ALLOCATE(elemConn(FVCOM_numElems*3))
    ALLOCATE(elemTypes(FVCOM_numElems))
 
    nodeIds           = (/(ii, ii=1, FVCOM_numNodes, 1)/)
    nodeCoords(:,1)   =  XG                                 
    nodeCoords(:,2)   =  YG
    nodeOwners        =  partition
    
    elemIds           = (/(ii, ii=1, FVCOM_numElems, 1)/)
    elemConn          =  NV
    elemTypes         =  ESMF_MESHELEMTYPE_TRI
    
    !----------------------------------------------------------------!       
    ! Create Mesh structure in 1 step
    ESMF_FVCOMMesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
         nodeIds=nodeIds, nodeCoords=nodeCoords,            &
         nodeOwners=nodeOwners, elementIds=elemIds,      &
         elementTypes=elemTypes, elementConn=elemConn,        &
         coordSys=ESMF_COORDSYS_CART,                              &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
	
    DEALLOCATE(nodeIds)
    DEALLOCATE(nodeCoords)
    DEALLOCATE(nodeOwners)
    
    DEALLOCATE(elemIds)
    DEALLOCATE(elemConn)
    DEALLOCATE(elemTypes)	
    

  END SUBROUTINE FVCOM2ESMF_mesh


END MODULE FISOC_OM_Wrapper
