
MODULE FISOC_ISM_Wrapper

  USE ESMF
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  USE ElmerSolver_mod
  USE MainUtils
!  USE Messages, ONLY : MessageUnit
  USE Messages, ONLY : InfoOutUnit
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_ISM_Wrapper_Init_Phase1,  FISOC_ISM_Wrapper_Init_Phase2,  &
       FISOC_ISM_Wrapper_Run, FISOC_ISM_Wrapper_Finalize


  ! Note that Elmer's CurrentModel is shared through the Types module (via MainUtils)

  ! Elmer element types (not directly available from Elmer, though ideally they should be)
  INTEGER, PARAMETER :: ELMER_ELEMENT_NODAL            = 101
  INTEGER, PARAMETER :: ELMER_ELEMENT_LINE_LINEAR      = 202
  INTEGER, PARAMETER :: ELMER_ELEMENT_LINE_QUADRAT     = 203
  INTEGER, PARAMETER :: ELMER_ELEMENT_LINE_CUBIC       = 204
  INTEGER, PARAMETER :: ELMER_ELEMENT_TRIANGLE_LINEAR  = 303
  INTEGER, PARAMETER :: ELMER_ELEMENT_TRIANGLE_QUADRAT = 306
  INTEGER, PARAMETER :: ELMER_ELEMENT_TRIANGLE_CUBIC   = 310
  INTEGER, PARAMETER :: ELMER_ELEMENT_QUADRIL_BILINEAR = 404
  INTEGER, PARAMETER :: ELMER_ELEMENT_QUADRIL_QUADRAT  = 408
  INTEGER, PARAMETER :: ELMER_ELEMENT_QUADRIL_QUADRAT2 = 409
  INTEGER, PARAMETER :: ELMER_ELEMENT_QUADRIL_CUBIC    = 412
  INTEGER, PARAMETER :: ELMER_ELEMENT_TETRA_LINEAR     = 504
  INTEGER, PARAMETER :: ELMER_ELEMENT_TETRA_QUADRAT    = 510
  INTEGER, PARAMETER :: ELMER_ELEMENT_PYRAMID_LINEAR   = 605
  INTEGER, PARAMETER :: ELMER_ELEMENT_PYRAMID_QUADRAT  = 613
  INTEGER, PARAMETER :: ELMER_ELEMENT_WEDGE_LINEAR     = 706
  INTEGER, PARAMETER :: ELMER_ELEMENT_WEDGE_QUADRAT    = 715
  INTEGER, PARAMETER :: ELMER_ELEMENT_HEXAHED_TRILIN   = 808
  INTEGER, PARAMETER :: ELMER_ELEMENT_HEXAHED_QUADRAT  = 820
  INTEGER, PARAMETER :: ELMER_ELEMENT_HEXAHED_QUADRAT2 = 827

  ! These variable names should be given in the ISM_varNames list in the FISOC 
  ! config file.  If omitted, defaults are given here.  In either case the user 
  ! must ensure these names correspond to Elmer variable names.
  ! TODO: is this information in the manual?
  ! For Elmer SSA runs, FISOC config should contain something like this:
  !  FISOC_ISM_ReqVars:  ISM_z_l0 ISM_z_lts 
  !  ISM_varNames:       'Zb' 'Zs'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_mask           = 'groundedmask'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_bmb            = 'meltRate'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_temperature_l0 = 'oceanTemperature'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_temperature_l1 = 'oceanTemperature'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_velocity_l0    = 'Velocity'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_thick          = 'Depth'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_dddt           = 'dhdt'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_z_l0           = 'Coordinate 3'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_z_l1           = 'Coordinate 3'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_z_lts          = 'Coordinate 3'
  CHARACTER(len=ESMF_MAXSTR),SAVE :: EIname_H_l0           = 'depth'

  ! The following mesh related properties are calculated during mesh conversion 
  ! during initialisation, and are needed during variable transfer while 
  ! timestepping.

  ! EI_NodeIDs contains unique node ids used in transferring ocean field 
  ! data from esmf fields to Elmer variables
  INTEGER, ALLOCATABLE,SAVE :: EI_NodeIDs(:)

! TODO: some of these vars below might not need to be module scope, need revising
  INTEGER, ALLOCATABLE,SAVE :: EI_ElementIDs(:)
  INTEGER, ALLOCATABLE,SAVE :: nodeOwners(:)

  ! This is an array of Elmer node IDs for the current partition that are owned by 
  ! the local pet.  This is needed for mapping Elmer data on to the ESMF mesh.
  INTEGER, ALLOCATABLE,SAVE :: ownedNodeIDs(:)
  
  ! Global node numbering reference (Elmer uses local node numbering).  We 
  ! must add this number to the node id whenever converting from Elmer node 
  ! ids to ESMF node ids.  Same for element id.
  ! Edit: probably no longer needed now that full ID lists are kept in memory.
  INTEGER,SAVE  :: EI_firstNodeThisPET, EI_firstElemThisPET

  ! Route handles for switching between Elmer arrays (which include duplicated nodes 
  ! along partition boundaries) and corresponding ESMF arrays (which don't)
  TYPE(ESMF_RouteHandle),SAVE :: RH_ESMF2Elmer

  ! nodal distgrid, an ESMF object holding information about the distribution of 
  ! Elmer nodes across partitions, needed for the redist related operations.
  TYPE(ESMF_distgrid),SAVE :: distgridElmer

CONTAINS

  !--------------------------------------------------------------------------------------
  ! This initialisation wrapper aims to convert the Elmer mesh and required variables 
  ! to the ESMF formats.  It also performs simple sanity/consistency checks.
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1(FISOC_config,vm,ISM_ExpFB,ISM_mesh,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(INOUT)           :: vm
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_mesh),INTENT(OUT)           :: ISM_mesh
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR)            :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: ISM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: ISM_configFile_FISOC, ISM_stdoutFile
    CHARACTER(len=MAX_STRING_LEN)         :: ISM_configFile_Elmer
    LOGICAL                               :: verbose_coupling
    TYPE(Mesh_t)                          :: Elmer_Mesh
    REAL(ESMF_KIND_R8)                    :: Elmer_dt, FISOC_ISM_dt
    INTEGER                               :: localpet, ISM_BodyID
    LOGICAL                               :: UseFootprint

    rc = ESMF_FAILURE

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Initialising Elmer/Ice"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    label = 'FISOC_ISM_ReqVars:'
    CALL FISOC_getListFromConfig(FISOC_config, label, ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ElmerSetVarNames(FISOC_config)

    ! FISOC tells Elmer which .sif to use
    CALL ESMF_ConfigGetAttribute(FISOC_config, ISM_configFile_FISOC, label='ISM_configFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    ISM_configFile_Elmer = ""
    ISM_configFile_Elmer = ISM_configFile_FISOC 

    ! information to get the appropriate surface for the Elmer mesh.
    CALL ESMF_ConfigGetAttribute(FISOC_config, ISM_BodyID, label='ISM_BodyID:', rc=rc)
    UseFootprint = .FALSE.
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
       msg = "ISM_BodyID not found, expecting Elmer footprint."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)
       UseFootprint = .TRUE.
    ELSEIF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN       
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF


    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********      ISM wrapper.  Init phase 1 method.        *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
    END IF

    ! FISOC sets the unit for Elmer's standard messaging routines to use instead of stdout
    CALL ESMF_ConfigGetAttribute(FISOC_config, ISM_stdoutFile, label='ISM_stdoutFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    OPEN(unit=ISM_outputUnit, file=ISM_stdoutFile, STATUS='REPLACE', ERR=101)
!    MessageUnit = ISM_outputUnit
    InfoOutUnit = ISM_outputUnit
    
    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Initialising Elmer/Ice parallel environment"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    CALL Initialise_Elmer_ParEnv(vm,rc=rc)

    ! initialise Elmer, and pull out the footprint before extrusion if we don't 
    ! have a bodyID defined.
    IF (UseFootprint) THEN
          msg = "Elmer footprint mesh not currently supported; please define body id in FISOC config"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!       CALL ElmerSolver_init(meshFootprint=Elmer_Mesh, ParEnvInitialised=.TRUE., &
!            inputFileName=ISM_configFile_Elmer) 
!       CALL Elmer2ESMF_meshFootprint(Elmer_mesh,ISM_mesh,vm,rc=rc)
    ELSE
       IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
          msg = "Calling Elmer/ice independent initialisation"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
       END IF
       CALL ElmerSolver_init(ParEnvInitialised=.TRUE., &
            inputFileName=ISM_configFile_Elmer) 
       IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
          msg = "Converting Elmer/Ice mesh to ESMF"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
       END IF
       CALL Elmer2ESMF_mesh(FISOC_config,ISM_BodyID,ISM_mesh,vm,rc=rc)
    END IF

    ! now Elmer is initialised we can check for timestep consistency
    IF (localpet.EQ.0) THEN
       IF ( .NOT.TimeStepConsistent(FISOC_config) ) THEN
          WRITE (msg, "(A,I0,A,I0,A)") &
               "FATAL: Elmer/FISOC timestep inconsistency"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    END IF

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Initialising Elmer/Ice fields for coupling"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    CALL FISOC_populateFieldBundle(ISM_ReqVarList,ISM_ExpFB,ISM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! get hold of list of required variables from Elmer and convert them here from elmer to esmf type.
    CALL getFieldDataFromISM(ISM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Initialising Elmer/Ice: first phase complete"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    rc = ESMF_SUCCESS

    RETURN

101 msg = "ISM failed to open stdoutFile"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1
  

  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2(FISOC_config,vm,ISM_ImpFB,ISM_ExpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ImpFB, ISM_ExpFB
    TYPE(ESMF_VM),INTENT(IN)              :: vm
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    LOGICAL                      :: verbose_coupling
    INTEGER                      :: localpet

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL sendFieldDataToISM(ISM_ImpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********      ISM wrapper.  Init phase 2 method.        *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we have access to the re-initialised OM fields, just in case the ISM needs "
       PRINT*,"to know about these in order to complete its initialisation."
       PRINT*,""
    END IF

    CALL getFieldDataFromISM(ISM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2




  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Run(FISOC_config,vm,ISM_ExpFB,ISM_ImpFB,rc)

    TYPE(ESMF_fieldbundle),INTENT(INOUT),OPTIONAL :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)      :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)             :: vm
    INTEGER,INTENT(OUT),OPTIONAL         :: rc

    INTEGER                      :: localPet
    TYPE(ESMF_field)             :: OM_bmb, ISM_z_l0, ISM_z_l1
    REAL(ESMF_KIND_R8),POINTER   :: OM_bmb_ptr(:),ISM_z_l0_ptr(:),ISM_z_l1_ptr(:)
    LOGICAL                      :: verbose_coupling

    rc = ESMF_FAILURE

! TODO:
! make sure to run only one timestep
! do some time step consistency checks

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
       PRINT*,"************        ISM wrapper.  Run method.           **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
    END IF

    ! get hold of the elmer variables for receiving inputs, and convert them here from esmf to elmer type.
    CALL sendFieldDataToISM(ISM_ImpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    

    WRITE (ISM_outputUnit,*) 'FISOC is about to call Elmer/Ice run method.'
    CALL ESMF_VMBarrier(vm, rc=rc)
    CALL ElmerSolver_run()
    CALL ESMF_VMBarrier(vm, rc=rc)
    WRITE (ISM_outputUnit,*) 'FISOC has just called Elmer/Ice run method.'

    ! get hold of list of required variables from Elmer and convert them here from elmer to esmf type.
    CALL getFieldDataFromISM(ISM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Finalize(FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL       :: rc
    TYPE(ESMF_VM),INTENT(IN)           :: vm

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
       PRINT*,"************      ISM wrapper.  Finalise method.        **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"ISM finalise method."
    END IF

! TODO: fix Elmer finalise call, not working for some reason...
!    CALL ElmerSolver_finalize()
!    CALL ElmerSolver_finalize(PreserveParEnvOpt=.TRUE.)

    CLOSE(unit=ISM_outputUnit, ERR=102)

    rc = ESMF_SUCCESS

    RETURN

102 msg = "ISM failed to close stdoutFile"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_ISM_Wrapper_Finalize




  !--------------------------------------------------------------------------------------
  ! Get the Elmer timestep (assume it is in units of years) and convert it to seconds.
  ! This should match the FISOM ISM timestep.
  !
  ! An additional check is that only one timestep should be run per call to Elmer run.
  !
  LOGICAL FUNCTION TimeStepConsistent(FISOC_config)

    TYPE(ESMF_config),INTENT(INOUT)  :: FISOC_config

    REAL(ESMF_KIND_R8)    :: Elmer_dt
    INTEGER               :: ISM_dt_sec, Elmer_dt_sec, tol, Elmer_nt, rc

    TimeStepConsistent = .TRUE.

    Elmer_dt = ListGetConstReal( CurrentModel % Simulation, 'Timestep Sizes' )
    Elmer_dt_sec = INT(FISOC_secPerYear * Elmer_dt)

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_sec, 'ISM_dt_sec',rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    tol = 10 ! tolerance in seconds
    IF ( (Elmer_dt_sec.GT.(ISM_dt_sec+tol)) .OR. (Elmer_dt_sec.LT.(ISM_dt_sec-tol)) ) THEN
       TimeStepConsistent = .FALSE.
       WRITE (msg, "(A,I0,A,I0,A)") &
            "WARNING: Elmer/FISOC timestep inconsistency (", &
            Elmer_dt_sec, " and ", ISM_dt_sec, ")"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)          
    END IF
       
    Elmer_nt = ListGetInteger( CurrentModel % Simulation, 'Timestep Intervals' )
    IF ( Elmer_nt.NE.1) THEN
       msg = "ERROR: Elmer/Ice .sif must set Timestep Intervals to 1"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
  END FUNCTION TimeStepConsistent



  !------------------------------------------------------------------------------
  SUBROUTINE  ElmerSetVarNames(FISOC_config)
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: ISM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: ISM_varNames(:)
    CHARACTER(len=ESMF_MAXSTR)            :: label
    INTEGER                               :: ii, rc
    
    label = 'FISOC_ISM_ReqVars:'
    CALL FISOC_getListFromConfig(FISOC_config, label, ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    label = 'ISM_varNames:'
    CALL FISOC_getListFromConfig(FISOC_config, label, ISM_varNames,rc=rc)
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
       msg = "ISM_varNames not found, using hard coded Elmer defaults"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)
    ELSE IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ELSE
       DO ii=1,SIZE(ISM_ReqVarList)
          ! Some ugly hard coding here.  Would be better to write this 
          ! stuff to the config object, but ESMF doesn't allow much in 
          ! the way of runtime editing of the config object.  So there 
          ! isn't a good and obvious way of setting defaults that are 
          ! model-specific (i.e. can't be hard coded in FISOC_utils).
          SELECT CASE (ISM_ReqVarList(ii))
          CASE ('OM_bmb')
             EIname_bmb             = ISM_varNames(ii)
          CASE ('OM_temperature_l0')
             EIname_temperature_l0  = ISM_varNames(ii)
          CASE ('ISM_z_l0')   
             EIname_z_l0            = ISM_varNames(ii)
          CASE ('ISM_z_lts')
             EIname_z_lts           = ISM_varNames(ii)
          CASE ('ISM_mask')
             EIname_mask            = ISM_varNames(ii)
          CASE ('ISM_thick')
             EIname_thick           = ISM_varNames(ii)
          CASE ('ISM_dddt')
             EIname_dddt            = ISM_varNames(ii)
          CASE DEFAULT
             msg = "unknown varName "//ISM_ReqVarList(ii)
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
                  line=__LINE__, file=__FILE__, rc=rc)
             CALL ESMF_Finalize(endflag=ESMF_END_ABORT)                 
          END SELECT

       END DO
    END IF
    
    IF (ALLOCATED(ISM_ReqVarList)) DEALLOCATE(ISM_ReqVarList)
    IF (ALLOCATED(ISM_varNames)) DEALLOCATE(ISM_varNames)

  END SUBROUTINE ElmerSetVarNames


  !------------------------------------------------------------------------------
  SUBROUTINE Initialise_Elmer_ParEnv(vm,rc)

    TYPE(ESMF_VM),INTENT(INOUT)      :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    INTEGER                          :: mpic, mpic_dup ! mpi comm, duplicate from the ISM VM
    INTEGER                          :: localPet, petCount, peCount, ierr

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, mpiCommunicator=mpic, localPet=localPet, &
         peCount=peCount,petCount=petCount,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ! The returned MPI communicator spans the same MPI processes that the VM
    ! is defined on.

    CALL MPI_Comm_dup(mpic, mpic_dup, ierr)
    ! Duplicate the MPI communicator not to interfere with ESMF communications.
    ! The duplicate MPI communicator can be used in any MPI call in the user
    ! code. 

    IF (peCount.NE.petCount) THEN
       msg = "FATAL: expected peCount = petCount"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ParEnv % MyPE = localPet
    OutputPE = ParEnv % MyPe
    ParEnv % PEs  = petCount
    ELMER_COMM_WORLD = mpic_dup
    ParEnv % ActiveComm = ELMER_COMM_WORLD
!    ParEnv % ActiveComm = MPI_COMM_WORLD ! or mpic_dup
    Parenv % NumOfNeighbours = 0
    ParEnv % Initialized = .TRUE.

    rc = ESMF_SUCCESS

  END SUBROUTINE Initialise_Elmer_ParEnv
  


  !--------------------------------------------------------------------------------------
  SUBROUTINE sendFieldDataToISM(ISM_ImpFB,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)     :: ISM_ImpFB 
    TYPE(ESMF_config),INTENT(INOUT)          :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)                 :: vm
    INTEGER,INTENT(OUT),OPTIONAL             :: rc


    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    TYPE(ESMF_Array)                      :: ESMFarr, Elmerarr 
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:), ptr_Elmer(:)
    INTEGER                               :: ii, jj, nn
    INTEGER,POINTER                       :: EI_fieldPerm(:)
    REAL(KIND=dp),POINTER                 :: EI_fieldVals(:)
    TYPE(Variable_t),POINTER              :: EI_field

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! get a list of fields and their names from the ISM import field bundle
    fieldCount = 0
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    fieldLoop: DO nn = 1,fieldCount

       ! access the FISOC version of the current field, and extract the array 
       ! from the field (because we can't do a redist on a field)
       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_FieldGet(fieldList(nn), array=ESMFarr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ! we need to create an elmer array on the fly and then convert from the 
       ! array to the elmer variable after the redist...
       Elmerarr = ESMF_ArrayCreate(distgridElmer, ESMF_TYPEKIND_R8, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)       
       
       ! now redist the ESMF array onto the Elmer array (this is still not yet 
       ! in an Elmer variable)
       CALL ESMF_ArrayRedist(ESMFarr, Elmerarr, RH_ESMF2Elmer, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ! ...and get a pointer to the data
       CALL ESMF_ArrayGet(Elmerarr, farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL FISOC_ArrayRedistFromField(RH_ESMF2Elmer,fieldList(nn),distgridElmer,ptr)

       ! check if the current variable is in the list of variables to be passed 
       ! to Elmer
       IF (FISOC_OM2ISM(fieldName,FISOC_config,rc=rc)) THEN

          SELECT CASE (TRIM(ADJUSTL(fieldName)))
              
          CASE ('OM_bmb')

             ! access the Elmer/Ice variable for the current field
             msg = "Getting Elmer melt var, name: "//EIname_bmb
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                  line=__LINE__, file=__FILE__, rc=rc)
             EI_field => VariableGet( CurrentModel % Mesh % Variables, &
                  EIname_bmb, UnFoundFatal=.TRUE.)
             EI_fieldVals => EI_field % Values
             EI_fieldPerm => EI_field % Perm
             ! copy the data from the Elmer array (ESMF object) to the Elmer 
             ! variable (native Elmer object), converting from m/sec to m/a.
             DO ii = 1,SIZE(ptr)
                EI_fieldVals(EI_fieldPerm(EI_NodeIDs(ii))) = ptr(ii) * FISOC_secPerYear
             END DO
          
          CASE ('OM_temperature_l0')
! TODO:
!             EI_field => VariableGet( CurrentModel % Mesh % Variables, &
!                  EIname_temperature_l0, UnFoundFatal=.TRUE.)
            msg = "WARNING: temperature exchange not tested... "//TRIM(ADJUSTL(fieldName))
            CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
                 line=__LINE__, file=__FILE__, rc=rc)          
            
          CASE ('OM_turnips')
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
!TODO: destroy or nullify temporary arrays?
       END IF

    END DO fieldLoop

    rc = ESMF_SUCCESS
    
  END SUBROUTINE sendFieldDataToISM




  !--------------------------------------------------------------------------------------
  SUBROUTINE getFieldDataFromISM(ISM_ExpFB,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)     :: ISM_ExpFB 
    TYPE(ESMF_config),INTENT(INOUT)          :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)                 :: vm
    INTEGER,INTENT(OUT),OPTIONAL             :: rc


    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)
    INTEGER                               :: ii, jj, nn
    INTEGER,POINTER                       :: EI_fieldPerm(:)
    REAL(KIND=dp),POINTER                 :: EI_fieldVals(:)
    TYPE(Variable_t),POINTER              :: EI_field


    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! get a list of fields and their names from the ISM export field bundle
    fieldCount = 0
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    fieldLoop: DO nn = 1,fieldCount

       ! access the FISOC version of the current field 
       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ! access the Elmer/Ice version of the current field
       SELECT CASE (TRIM(ADJUSTL(fieldName)))

       CASE ('ISM_z_l0')
          EI_field => VariableGet( CurrentModel % Mesh % Variables, &
               EIname_z_l0, UnFoundFatal=.TRUE.)
          EI_fieldVals => EI_field % Values
          EI_fieldPerm => EI_field % Perm ! don't need perm for coords
          DO ii = 1,SIZE(ownedNodeIDs)
             IF (ASSOCIATED(EI_fieldPerm)) THEN
                ptr(ii) = EI_fieldVals(EI_fieldPerm(ownedNodeIds(ii)))
             ELSE
                ptr(ii) = EI_fieldVals(ownedNodeIds(ii))
             END IF
          END DO
          
       CASE ('ISM_thick')
          EI_field => VariableGet( CurrentModel % Mesh % Variables, &
               EIname_thick, UnFoundFatal=.TRUE.)
          EI_fieldVals => EI_field % Values
          EI_fieldPerm => EI_field % Perm ! don't need perm for coords
          DO ii = 1,SIZE(ownedNodeIDs)
             IF (ASSOCIATED(EI_fieldPerm)) THEN
                ptr(ii) = EI_fieldVals(EI_fieldPerm(ownedNodeIds(ii)))
             ELSE
                ptr(ii) = EI_fieldVals(ownedNodeIds(ii))
             END IF
          END DO
          
       CASE ('ISM_dddt')
          EI_field => VariableGet( CurrentModel % Mesh % Variables, &
               EIname_dddt, UnFoundFatal=.TRUE.)
          EI_fieldVals => EI_field % Values
          EI_fieldPerm => EI_field % Perm ! don't need perm for coords
          DO ii = 1,SIZE(ownedNodeIDs)
             IF (ASSOCIATED(EI_fieldPerm)) THEN
                ptr(ii) = EI_fieldVals(EI_fieldPerm(ownedNodeIds(ii)))
             ELSE
                ptr(ii) = EI_fieldVals(ownedNodeIds(ii))
             END IF
          END DO
          ptr = ptr / FISOC_secPerYear
          WRITE(msg,*) "ISM_DDDT max ", MAXVAL(ptr) ," min ", MINVAL(ptr)
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          
       CASE ('ISM_z_lts')
          EI_field => VariableGet( CurrentModel % Mesh % Variables, &
               EIname_z_lts, UnFoundFatal=.TRUE.)
          EI_fieldVals => EI_field % Values
          EI_fieldPerm => EI_field % Perm
          DO ii = 1,SIZE(ownedNodeIDs)
             IF (ASSOCIATED(EI_fieldPerm)) THEN
                ptr(ii) = EI_fieldVals(EI_fieldPerm(ownedNodeIds(ii)))
             ELSE
                ptr(ii) = EI_fieldVals(ownedNodeIds(ii))
             END IF
          END DO
          
       CASE ('ISM_mask')
          EI_field => VariableGet( CurrentModel % Mesh % Variables, &
               EIname_mask, UnFoundFatal=.TRUE.)
          EI_fieldVals => EI_field % Values
          EI_fieldPerm => EI_field % Perm
          DO ii = 1,SIZE(ownedNodeIDs)
            ptr(ii) = EI_fieldVals(EI_fieldPerm(ownedNodeIds(ii)))
            IF (ptr(ii).GT.0) ptr(ii) = MASK_GROUNDED_ICE
            IF (ptr(ii).EQ.0) ptr(ii) = MASK_GL
            IF (ptr(ii).LT.0) ptr(ii) = MASK_FLOATING_ICE
          END DO
           
       CASE ('ISM_temperature_l0','ISM_temperature_l1','ISM_velocity_l0','ISM_z_l1','ISM_z_l0_previous','ISM_z_lts_previous')
          msg = "WARNING: ignored variable: "//TRIM(ADJUSTL(fieldName))
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__, rc=rc)          
       
       CASE ('ISM_dTdz_l0','ISM_dsdt','ISM_z_l0_linterp')
          msg = "INFO: not extracting derived variable: "//TRIM(ADJUSTL(fieldName))
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
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
       
    END DO fieldLoop
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE getFieldDataFromISM
  

  !------------------------------------------------------------------------------
  !
  ! Convert an Elmer mesh to ESMF structures 
  ! 
  ! Only the surface defined by BodyID will be extracted.  It is expected that 
  ! elements on this surface will be triangles and quads.
  SUBROUTINE Elmer2ESMF_mesh(FISOC_config,BodyID,ESMF_ElmerMesh,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)  :: FISOC_config
    INTEGER,INTENT(IN)               :: BodyID
    TYPE(ESMF_mesh),INTENT(INOUT)    :: ESMF_ElmerMesh
    TYPE(ESMF_VM),INTENT(IN)         :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    CHARACTER(len=ESMF_MAXSTR)       :: label
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: ISM_ProjVector(:)
    TYPE(mesh_t)                     :: ElmerMesh
    INTEGER                          :: localPet, petCount
    INTEGER                          :: ii, nodeIndex
    INTEGER                          :: EI_numNodes 
    LOGICAL                          :: verbose_coupling
    INTEGER,ALLOCATABLE              :: ESMF_elemTypes(:)
    INTEGER,ALLOCATABLE              :: elemConn(:), ElementIDs_global(:)
    INTEGER,ALLOCATABLE              :: nodeIds(:), nodeIds_global(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:)
    INTEGER                          :: numQuadElems, numTriElems, numElems


    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Obtain basic parallel information (pet = persistent execution thread)
    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! We can access the main Elmer derived types from here.  We assume there is
    ! one mesh in this Elmer simulation.
    ElmerMesh = CurrentModel % Mesh
    CALL ElmerMeshSanityChecks(ElmerMesh)

    ! If a projection vector is not specified, the default is downwards.
    label = "ISM_ProjVector:"
    CALL FISOC_getListFromConfig(FISOC_config,label,ISM_ProjVector,rc)
    IF  (rc.EQ.ESMF_RC_NOT_FOUND) THEN
       msg = "ISM_ProjVector not found, assuming downwards."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)
       ALLOCATE(ISM_ProjVector(3))
       ISM_ProjVector = (/0.0_dp, 0.0_dp, -1.0_dp/)   
    ELSEIF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! We are only interested in the BodyID corresponding to the surface 
    ! on which information will be exchanged with the ocean (typically 
    ! the lower surface).
    ! Count the elements on the relevant surface.
    numQuadElems = countElements(ElmerMesh,BodyID,(/  &
         ELMER_ELEMENT_QUADRIL_BILINEAR, ELMER_ELEMENT_QUADRIL_QUADRAT, &
         ELMER_ELEMENT_QUADRIL_QUADRAT2,ELMER_ELEMENT_QUADRIL_CUBIC/))
    numTriElems = countElements(ElmerMesh,BodyID,(/  &
         ELMER_ELEMENT_TRIANGLE_LINEAR, ELMER_ELEMENT_TRIANGLE_QUADRAT, &
         ELMER_ELEMENT_TRIANGLE_CUBIC/))
    numElems    = numQuadElems+numTriElems
    
    ! Some explanation of Nodes variables:
    ! NodeIDs is a unique list of all the Elmer IDs on this partition/body
    ! EI_numNodes is the number of nodes on this partition/body

    ALLOCATE(elemConn(4*numQuadElems+3*numTriElems))
    ALLOCATE(ESMF_elemTypes(numElems))

    ALLOCATE(ElementIDs_global(numElems))

    ALLOCATE(EI_ElementIDs(numElems))
    ALLOCATE(EI_NodeIDs(4*numQuadElems+3*numTriElems))

    CALL findNodesAndElements(ElmerMesh,BodyID, (/ &
         ELMER_ELEMENT_QUADRIL_BILINEAR, ELMER_ELEMENT_QUADRIL_QUADRAT, &
         ELMER_ELEMENT_QUADRIL_QUADRAT2,ELMER_ELEMENT_QUADRIL_CUBIC, &
         ELMER_ELEMENT_TRIANGLE_LINEAR, ELMER_ELEMENT_TRIANGLE_QUADRAT, &
         ELMER_ELEMENT_TRIANGLE_CUBIC/))
    ! Note: the findNodes call will calculate a reallocated array for EI_NodeIDs, 
    ! smaller than the original array (due to uniqueness of node IDs).
    ! The call will also populate EI_ElementIDs.
    ! These are module-wide variables.

    EI_numNodes = SIZE(EI_NodeIDs)
    ALLOCATE(NodeIDs(EI_numNodes))
    ALLOCATE(NodeIDs_global(EI_numNodes))

    ! Calculate some properties dependent on other PETs, and 
    ! which will be needed for variable exchange.
    EI_firstNodeThisPET = firstItemThisPET(EI_numNodes,vm)
    EI_firstElemThisPET = firstItemThisPET(numElems,vm)

    ElementIDs_global = (/(ii, ii=EI_firstElemThisPET, EI_firstElemThisPET+numElems-1, 1)/)
    nodeIDs_global    = (/(ii, ii=EI_firstNodeThisPET, EI_firstNodeThisPET+EI_numNodes-1, 1)/)

    ALLOCATE(nodeOwners(EI_numNodes))
    nodeOwners=localPet

    ALLOCATE(nodeCoords(EI_numNodes*2))

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Elmer mesh: getting node coords"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    CALL getNodeCoords(nodeCoords,ElmerMesh,ISM_ProjVector,EI_numNodes)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Elmer mesh: making node ids unique"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    CALL uniquifyGlobalNodeIDs(nodeIDs_global,nodeCoords,vm)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Elmer mesh: identifing local nodes"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    CALL FISOC_locallyOwnedNodes(localPet,EI_nodeIDs,nodeOwners,ownedNodeIDs)

    ! TODO (probably not urgent, maybe not needed at all)
    ! Note that ElmerMesh%ParallelInfo%GlobalDOFs may already contain global unique node ids 
    ! (not just for the boundary of interest but for the whole mesh). It might make sense to 
    ! use these rather than construct the uniqueIDs here. 
    
    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Elmer mesh: building element connectivity"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    CALL buildElementConnectivity(elemConn, nodeIDs_global, BodyID, ESMF_elemTypes, ISM_ProjVector,ELmerMesh)

    ! Create Mesh structure in 1 step
    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Elmer mesh: creating ESMF mesh structure"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    ESMF_ElmerMesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
         nodeIds=nodeIds_global, nodeCoords=nodeCoords,            &
         nodeOwners=nodeOwners, elementIds=ElementIDs_global,      &
         elementTypes=ESMF_elemTypes, elementConn=elemConn,        &
         coordSys=ESMF_COORDSYS_CART,                              &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Elmer mesh: creating routehandles"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
!    CALL CreateArrayMappingRouteHandles(ESMF_ElmerMesh,nodeIDs_global,vm)
    CALL FISOC_CreateOneToManyRouteHandle(ESMF_ElmerMesh,nodeIDs_global,RH_ESMF2Elmer,distgridElmer,vm)

    IF (ALLOCATED(elemConn)) DEALLOCATE(elemConn)
    IF (ALLOCATED(nodeIds_global)) DEALLOCATE(nodeIds_global)
    IF (ALLOCATED(nodeCoords)) DEALLOCATE(nodeCoords)
    
  END SUBROUTINE Elmer2ESMF_mesh
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  ! These module-scope route handles are for array mappings between the ESMF 
  ! fields (with unique nodes) and Elmer fields (with some node duplication 
  ! across partition boundaries).
  !
  ! Input vars:
  ! ESMF_ElmerMesh - Arrays created on this mesh will not duplicate nodes 
  ! ElmerIDs - global node IDs including duplicates across partition boundaries
  !
  SUBROUTINE CreateArrayMappingRouteHandles(ESMF_ElmerMesh,ElmerIDs,vm)

    TYPE(ESMF_mesh), INTENT(IN) :: ESMF_ElmerMesh
    INTEGER,INTENT(IN)          :: ElmerIDs(:) 
    TYPE(ESMF_vm),INTENT(IN)    :: vm

    TYPE(ESMF_array)   :: DummyArr_ESMF  ! will not contain duplicate nodes
    TYPE(ESMF_array)   :: DummyArr_Elmer ! will contain duplicate nodes
    TYPE(ESMF_distgrid):: distgridESMF
    INTEGER            :: rc

    INTEGER,TARGET,ALLOCATABLE :: dg_Elmer(:), dg_ESMF(:)
    INTEGER :: localPET,PETcount,tmp
    REAL(ESMF_KIND_R8),POINTER :: ptr_ESMF(:),ptr_Elmer(:)

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create dummy array on the distgrid for ESMF Elmer mesh
    CALL ESMF_MeshGet(ESMF_ElmerMesh, nodalDistgrid=distgridESMF, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!distgridESMF  = ESMF_DistgridCreate(arbSeqIndexList=ElmerIDs, rc=rc)
!IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, file=__FILE__)) &
!     CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DummyArr_ESMF = ESMF_ArrayCreate(distgridESMF, ESMF_TYPEKIND_R8, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (localPet.EQ.0) THEN
       msg = "Elmer routehandles: creating distrgrid"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    ! Create a distgrid containing sequence indices for the Elmer array, i.e. 
    ! containing the duplicate node IDs.  This can be used to create an array 
    ! including duplicates, like the Elmer fields.
    distgridElmer  = ESMF_DistgridCreate(arbSeqIndexList=ElmerIDs, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    DummyArr_Elmer = ESMF_ArrayCreate(distgridElmer, ESMF_TYPEKIND_R8, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!get stuff for printing
!    CALL ESMF_DistGridGet(distgridElmer, 0,elementCount=ct_Elmer, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    ALLOCATE(dg_Elmer(ct_Elmer))
!    CALL ESMF_DistGridGet(distgridElmer, 0, seqIndexList=dg_Elmer, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    CALL ESMF_DistGridGet(distgridESMF, 0,elementCount=ct_ESMF, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    ALLOCATE(dg_ESMF(ct_ESMF))
!    CALL ESMF_DistGridGet(distgridESMF, 0, seqIndexList=dg_ESMF, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! initialilse arrays to zero, probably not needed
    CALL ESMF_ArrayGet(DummyArr_Elmer, farrayPtr=ptr_Elmer, rc=rc)
    CALL ESMF_ArrayGet(DummyArr_ESMF,  farrayPtr=ptr_ESMF, rc=rc)
    ptr_Elmer=0.0
    ptr_ESMF=0.0
    IF (ASSOCIATED(ptr_Elmer)) NULLIFY(ptr_Elmer)
    IF (ASSOCIATED(ptr_ESMF)) NULLIFY(ptr_ESMF)
    
    ! Create the route handles for later use 
    IF (localPet.EQ.0) THEN
       msg = "Elmer routehandles: store routehandle"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    CALL ESMF_ArrayRedistStore(DummyArr_ESMF, DummyArr_Elmer, &
         RH_ESMF2Elmer, ignoreUnmatchedIndices=.TRUE., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

! Note: many to one default behaviour is not suitable, so don't 
! use this method for Elmer to ESMF data exchange
!    CALL ESMF_ArrayRedistStore(DummyArr_Elmer, DummyArr_ESMF, &
!         RH_Elmer2ESMF, ignoreUnmatchedIndices=.TRUE., rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Tidy up 
    IF (localPet.EQ.0) THEN
       msg = "Elmer routehandles: tidy up"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    CALL ESMF_ArrayDestroy(DummyArr_Elmer, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_ArrayDestroy(DummyArr_ESMF, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    CALL ESMF_DistGridDestroy(distgridElmer, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    CALL ESMF_DistGridDestroy(distgridESMF, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE CreateArrayMappingRouteHandles



  !------------------------------------------------------------------------------
  SUBROUTINE ElmerMeshSanityChecks(ElmerMesh)
    
    TYPE(mesh_t),INTENT(IN) :: ElmerMesh
    
    ! some basic sanity checks
    IF (.NOT.ASSOCIATED(ElmerMesh % Elements)) THEN
       msg = "Elmer mesh elements not associated"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (.NOT.ASSOCIATED(ElmerMesh % Nodes)) THEN
       msg = "Elmer mesh nodes not associated"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
!    IF (SIZE(ElmerMesh % Elements).NE.ElmerMesh % NumberOfBulkElements + ElmerMesh % NumberOfBoundaryElements) THEN
!       msg = "Elmer mesh number of elements inconsistency"
!       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
!            line=__LINE__, file=__FILE__)
!       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    END IF
    IF (SIZE(ElmerMesh % Nodes % x) .NE. ElmerMesh % Nodes % NumberOfNodes) THEN
       msg = "Elmer mesh number of nodes inconsistency"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

  END SUBROUTINE ElmerMeshSanityChecks



  !------------------------------------------------------------------------------
  SUBROUTINE getNodeCoords(nodeCoords,ElmerMesh,ISM_ProjVector,EI_numNodes)
    
    REAL(ESMF_KIND_R8),DIMENSION(:),INTENT(OUT) :: nodeCoords
    TYPE(mesh_t), INTENT(IN)                    :: ElmerMesh
    REAL(ESMF_KIND_R8),DIMENSION(3),INTENT(IN)  :: ISM_ProjVector
    INTEGER,INTENT(IN)                          :: EI_numNodes 

    INTEGER  :: ii, nodeIndex

    ! loop over nodes to get coords
    DO ii = 1,EI_numNodes
       nodeIndex = (ii-1)*2+1
       nodeCoords(nodeIndex) = ElmerMesh % Nodes % x(EI_NodeIDs(ii))
       nodeIndex = (ii-1)*2+2
       nodeCoords(nodeIndex) = ElmerMesh % Nodes % y(EI_NodeIDs(ii))
    END DO
    IF ( ALL(ISM_ProjVector.NE. (/0.0_dp, 0.0_dp, -1.0_dp/)) ) THEN 
       msg = "Projection vector is not yet implemented - more code here please!"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

  END SUBROUTINE getNodeCoords

  
  !------------------------------------------------------------------------------
  ! Element connectivity: an ordered list of nodes for each element.
  ! Also convert Elmer element types to ESMF element types here.
  SUBROUTINE buildElementConnectivity(elemConn, nodeIDs_global, BodyID, ESMF_elemTypes, ISM_ProjVector, ELmerMesh)

    INTEGER,DIMENSION(:),INTENT(INOUT)  :: elemConn(:)
    INTEGER,INTENT(IN)                  :: nodeIDs_global(:)
    INTEGER,INTENT(IN)                  :: BodyID
    INTEGER,INTENT(OUT)                 :: ESMF_elemTypes(:)
    REAL(ESMF_KIND_R8),INTENT(IN)       :: ISM_ProjVector(3)
    TYPE(mesh_t), INTENT(IN)            :: ElmerMesh

    INTEGER                             :: ii, nn, cc, Direction, nNodes
    TYPE(element_t)                     :: Element
    LOGICAL                             :: ReOrderNodes = .TRUE.
    REAL(ESMF_KIND_R8),ALLOCATABLE      :: ssNO_coords(:,:)
    INTEGER,ALLOCATABLE,DIMENSION(:)    :: NO_IDs

    Direction = ANTI_CLOCKWISE

    cc = 0
    AllElements: DO ii = 1,SIZE(EI_ElementIDs)
       Element = ElmerMesh%Elements(EI_ElementIDs(ii))
       ESMF_elemTypes(ii) = get_ESMF_elementType(ElmerMesh%Elements(EI_ElementIDs(ii) )%TYPE%ElementCode)
       nNodes = SIZE(Element%NodeIndexes) ! number of nodes on this element
!       IF (ReOrderNodes) THEN
!          ALLOCATE(NO_coords(nNodes,3))
!          ALLOCATE(NO_IDs(nNodes))
!          NodesThisElement: DO nn = 1,nNodes
!***build node IDs and coords here
!          END DO NodesThisElement
!          CALL  NodeOrdering(Direction, ISM_ProjVector, NO_Ids, NO_coords,&
!               OrderedNodeIds, nNodes)
!          DEALLOCATE(NO_IDs);DEALLOCATE(NO_coords)
!       ELSE
          NodesThisElement: DO nn = 1,nNodes
             cc = cc+1
!             elemConn(cc) = nodeIDs_global(findIndex(Element%NodeIndexes(nn),EI_NodeIDs))
             elemConn(cc) = findIndex(Element%NodeIndexes(nn),EI_NodeIDs)
          END DO NodesThisElement
!       END IF
    END DO AllElements
!print*,"node ordering to go here if needed..."

  END SUBROUTINE buildElementConnectivity

  
  !------------------------------------------------------------------------------
  INTEGER FUNCTION findIndex(val,IDs)
    INTEGER,INTENT(IN)              :: val
    INTEGER,INTENT(IN),DIMENSION(:) :: IDs    
    INTEGER :: nn
    findIndex = -1
    DO nn=1,SIZE(IDs)
       IF (val .EQ. IDs(nn)) THEN
          findIndex = nn
          EXIT
       END IF
    END DO
    IF (findIndex.LT.0) THEN
       msg = "cannot find index..."
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    RETURN
  END FUNCTION findIndex



  !------------------------------------------------------------------------------
  ! Count the number of elements of the given type  on the required 
  ! Elmer body on the local pet.
  INTEGER FUNCTION countElements(ElmerMesh,BodyID,typeList)
    TYPE(mesh_t), INTENT(IN)      :: ElmerMesh
    INTEGER, INTENT(IN)           :: typeList(:)
    INTEGER, INTENT(IN)           :: BodyID
    INTEGER                       :: ii, tt
    countElements = 0
    AllElements: DO ii = 1,SIZE(ElmerMesh%Elements)
       IF (ASSOCIATED (ElmerMesh%Elements(ii)%TYPE)) then
          ElementTypes: DO tt = 1,SIZE(typeList)
             IF (typeList(tt) .EQ. ElmerMesh%Elements(ii)%TYPE%ElementCode) THEN
                IF (BodyID .EQ. ElmerMesh%Elements(ii)%BodyID) THEN 
                   countElements = countElements + 1
                END IF
             END IF
          END DO ElementTypes
       END IF
    END DO AllElements
    RETURN 
  END FUNCTION countElements



  !------------------------------------------------------------------------------
  ! Count the number of nodes on the required Elmer body on the local pet. 
  ! A unique array of (local) node IDs will be calculated (EI_NodeIDs).
  ! An array of (global) element IDs will be created (EI_ElementIDs).
  ! These are the Elmer IDs.
  SUBROUTINE findNodesAndElements(ElmerMesh,BodyID,elemTypeList) 

    TYPE(mesh_t), INTENT(IN)          :: ElmerMesh
    INTEGER, INTENT(IN)               :: elemTypeList(:),BodyID

    REAL(ESMF_KIND_R8), ALLOCATABLE   :: NodeIDs_real(:)
    INTEGER                           :: ii, tt, nn, cc, ee

    ALLOCATE(NodeIDs_real(SIZE(EI_NodeIDs)))

    ! find all nodes on this body via elements on the body (and populate 
    ! the element list while we're at it)
    cc = 0
    ee = 0
    AllElements: DO ii = 1,SIZE(ElmerMesh%Elements)
       IF (ASSOCIATED (ElmerMesh%Elements(ii)%TYPE)) then
          ElementTypes: DO tt = 1,SIZE(elemTypeList)
             IF (elemTypeList(tt) .EQ. ElmerMesh%Elements(ii)%TYPE%ElementCode) THEN
                IF (BodyID .EQ. ElmerMesh%Elements(ii)%BodyID) THEN 
                   ee = ee + 1
                   EI_ElementIDs(ee) = ElmerMesh % Elements(ii) % ElementIndex
                   DO nn = 1,SIZE(ElmerMesh%Elements(ii)%NodeIndexes)
                      cc = cc+1
                      NodeIDs_real(cc) = REAL(ElmerMesh%Elements(ii)%NodeIndexes(nn),ESMF_KIND_R8)
                   END DO
                END IF
             END IF
          END DO ElementTypes
       END IF
    END DO AllElements

    ! Some nodes will occur on multiple elements.  We only want a 
    ! unique list of nodes here.
    CALL  Unique1DArray(NodeIDs_real)
    IF (ALLOCATED(EI_NodeIDs))  DEALLOCATE(EI_NodeIDs)
    ALLOCATE(EI_NodeIDs(SIZE(NodeIDs_real)))
    EI_NodeIDs = NINT(NodeIDs_real)

  END SUBROUTINE findNodesAndElements


  
  !------------------------------------------------------------------------------
  INTEGER FUNCTION numElementsByType(Elmer_mesh,typeList,ESMF_elementTypeList,elementIDlist,elemConn)

    TYPE(mesh_t), INTENT(IN)      :: Elmer_mesh
    INTEGER, INTENT(IN)           :: typeList(:)
    INTEGER, INTENT(OUT),OPTIONAL :: ESMF_elementTypeList(:),elementIDlist(:),elemConn(:)!,maxNodes

    INTEGER                  :: numElems, ii, jj, kk, nn

    numElementsByType = 0

    numElems = SIZE(Elmer_mesh % Elements)

    kk = 1

    DO ii = 1,numElems
       DO jj = 1,SIZE(typeList)
          IF (typeList(jj) .EQ. Elmer_mesh % Elements(ii) % TYPE % ElementCode) THEN
             numElementsByType = numElementsByType + 1          
!             IF (PRESENT(maxNodes)) THEN
!                maxNodes = MAX(getNumNodes(Elmer_mesh % Elements(ii) % TYPE % ElementCode),maxNodes)
!             END IF

             ! convert element type codes from Elmer to ESMF
             IF (PRESENT(ESMF_elementTypeList)) THEN
                ESMF_elementTypeList(ii) = get_ESMF_elementType(Elmer_mesh % Elements(ii) % TYPE % ElementCode)
             END IF
             IF (PRESENT(elementIDlist)) THEN
                elementIDlist(ii) = Elmer_mesh % Elements(ii) % ElementIndex
                IF (elementIDlist(ii).LT.-1) THEN
                   msg = "WARNING: Elmer element index less than zero"
                   CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                        line=__LINE__, file=__FILE__)
                   CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
                END IF
             END IF
             IF (PRESENT(elemConn)) THEN
                DO nn = 1,SIZE(Elmer_mesh % Elements(ii) % NodeIndexes)
                   elemConn(kk) = Elmer_mesh % Elements(ii) % NodeIndexes(nn)
                   kk=kk+1
                END DO
             END IF
          END IF
       END DO
    END DO
    
    RETURN

  END FUNCTION numElementsByType


  ! 
  SUBROUTINE findMatch(matching_nodes_indices,nodeCoords_all,nodeCoords,startNode,toler)
    INTEGER,ALLOCATABLE,INTENT(OUT)  :: matching_nodes_indices(:)
    REAL(ESMF_KIND_R8),INTENT(IN)    :: nodeCoords_all(:) ! node coordinates over all PETs
    REAL(ESMF_KIND_R8),INTENT(IN)    :: nodeCoords(2)     ! coordinates for this node
    INTEGER,INTENT(IN)               :: startNode
    REAL(ESMF_KIND_R8),INTENT(IN)    :: toler             ! tolerance
    
    INTEGER :: nNodes2test, nNodes, matches, ii

    IF (ALLOCATED(matching_nodes_indices)) DEALLOCATE (matching_nodes_indices) 
    
    nNodes = SIZE(nodeCoords_all)/2
    nNodes2test = nNodes + 1 - startNode
    ALLOCATE(matching_nodes_indices(nNodes2test))
    matching_nodes_indices = FISOC_missing

    matches = 0
    DO ii = startNode,nNodes
       IF ((                                                        &
            (nodeCoords(1) .GE. nodeCoords_all((ii-1)*2+1)-toler)   &
            .AND.                                                   &
            (nodeCoords(1) .LE. nodeCoords_all((ii-1)*2+1)+toler) ) &
            .AND. (                                                 &
            (nodeCoords(2) .GE. nodeCoords_all((ii-1)*2+2)-toler)   &
            .AND.                                                   &
            (nodeCoords(2) .LE. nodeCoords_all((ii-1)*2+2)+toler) ) &
            ) THEN
          matches=matches+1
          matching_nodes_indices(matches) = ii
       END IF
    END DO

    CALL FISOC_shrink(matching_nodes_indices,FISOC_missing)

  END SUBROUTINE findMatch



  !------------------------------------------------------------------------------
  ! This subroutine deserves a silly name on account of being a pain in the *** 
  ! to write.
  ! The main purpose is to make sure that where a node exists on two different 
  ! partitions but at the same physical location it should have the same global 
  ! ID on both partitions.
  ! It also makes sure that each node has one owner.
  !
  SUBROUTINE uniquifyGlobalNodeIDs(nodeIDs_global,nodeCoords,vm)
 
    INTEGER,INTENT(INOUT),ALLOCATABLE :: nodeIDs_global(:) ! global IDs on local PET
    REAL(ESMF_KIND_R8),INTENT(INOUT):: nodeCoords(:)     ! node coordinates on local PET
    TYPE(ESMF_vm),INTENT(IN)        :: vm

    INTEGER  :: ii, jj, nn, pp, index_all, lastUniqIndex
    INTEGER              :: nNodes       ! number of nodes this PET
    INTEGER              :: nNodes_all   ! total number over all PETs
    INTEGER,ALLOCATABLE  :: nNodes_arr(:)! array of nNodes on all PETs
    INTEGER              :: MAXnNodes    ! highest nNodes of all PETs
    INTEGER              :: nNodes_prevPETs ! total number over previous PETs (while looping over PETs)
    INTEGER              :: nNodes_PETpp ! nNodes for PET number pp
    INTEGER              :: PETcount, localPET, rc
    INTEGER,ALLOCATABLE  :: matching_nodes_indices(:)
    INTEGER,ALLOCATABLE  :: IDs_global_all(:) ! global IDs over all PETs
    INTEGER,ALLOCATABLE  :: EI_nodeIDs_all(:) ! Elmer local node IDs over all PETs
    INTEGER,ALLOCATABLE  :: EI_nodeIDs_all_temp(:) ! temporary holder for EI_nodeIDs_all 
    INTEGER,ALLOCATABLE  :: nodeOwners_all(:) ! node owners over all PET
    INTEGER,ALLOCATABLE  :: sendOffsets(:)    ! offsets for scattering _all array to PETs
    REAL(ESMF_KIND_R8),ALLOCATABLE  :: nodeCoords_all(:) ! node coordinates over all PETs

    CALL ESMF_VMGet(vm, localPet=localPET, petCount=PETcount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Get the number of nodes on this PET
    nNodes = SIZE(nodeIDs_global)

    ! Get the max number of nodes per pet.
    CALL ESMF_VMAllFullReduce(vm, (/nNodes/), MAXnNodes, 1, ESMF_REDUCE_MAX, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Get the total number of nodes over all PETs.
    CALL ESMF_VMAllFullReduce(vm,(/nNodes/),nNodes_all,1,ESMF_REDUCE_SUM, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Get an array of the number of nodes on each PET.
    ALLOCATE(nNodes_arr(PETcount))
    CALL ESMF_VMGather(vm, (/nNodes/), nNodes_arr, 1, 0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ALLOCATE(IDs_global_all(nNodes_all))
    IDs_global_all = FISOC_missing
    

    ! Gather the info we need from individual PETs onto arrays spanning 
    ! all PETs.

    CALL FISOC_VMAllGather(vm,EI_NodeIDs,EI_nodeIDs_all,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(EI_nodeIDs_all_temp(SIZE(EI_nodeIDs_all)))
    EI_nodeIDs_all_temp = EI_nodeIDs_all

    CALL FISOC_VMAllGather(vm,nodeOwners,nodeOwners_all,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_VMAllGather(vm,nodeCoords,nodeCoords_all,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! The following code all takes place in serial.  It could probably be made more efficient.
    ! We loop over PETs whilst on PET 0
    IF (localPET.EQ.0) THEN
       nNodes_prevPETs = 0
       PETloop: DO pp = 0,PETcount-1
          nNodes_PETpp = nNodes_arr(pp+1) ! the number of ndoes on PET pp
          ! loop over the nodes on each PET, finding matches on other PETs, and updating 
          ! the IDs accordingly (and the EI IDs too)
          NodesLoop: DO nn = 1,nNodes_PETpp
             index_all = nNodes_prevPETs+nn
             IF (index_all.EQ.1) THEN
                IDs_global_all(index_all) = 1
                lastUniqIndex = 1
             ELSE
                ! Only populate missing values.  If it aint missing it has been populated 
                ! as a duplicate node already.
                IF (IDs_global_all(index_all).EQ.FISOC_missing) THEN
                   IDs_global_all(index_all) = IDs_global_all(lastUniqIndex)+1
                   lastUniqIndex = index_all
                END IF
                CALL findMatch(matching_nodes_indices,nodeCoords_all, &
                     nodeCoords_all((index_all-1)*2+1:(index_all-1)*2+2), index_all+1, 0.1_dp )
                matchingNodes: DO ii = 1,SIZE(matching_nodes_indices)
                   IF (IDs_global_all(matching_nodes_indices(ii)).EQ.FISOC_missing) THEN
                      IDs_global_all(matching_nodes_indices(ii)) = IDs_global_all(index_all)
                      nodeOwners_all(matching_nodes_indices(ii)) = pp
!                      EI_nodeIDs_all(matching_nodes_indices(ii)) = EI_nodeIDs_all_temp(index_all)
                   END IF
                END DO matchingNodes
                IF (ALLOCATED(matching_nodes_indices)) DEALLOCATE(matching_nodes_indices)
             END IF
          END DO NodesLoop
          nNodes_prevPETs = nNodes_prevPETs + nNodes_PETpp
       END DO PETloop
    END IF
    
    ! scatter the three important arrays from PET 0 back to all PETs.

    ALLOCATE(sendOffsets(PETcount))
    sendOffsets(1) = 0
    IF (PETcount.GT.1) THEN
       DO ii = 2,PETcount
          sendOffsets(ii) = sendOffsets(ii-1)+nNodes_arr(ii-1)
       END DO
    END IF

    CALL  ESMF_VMScatterV(vm, IDs_global_all, nNodes_arr, &
         sendOffsets, nodeIDs_global, nNodes, 0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL  ESMF_VMScatterV(vm, nodeOwners_all, nNodes_arr, &
         sendOffsets, nodeOwners, nNodes, 0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL  ESMF_VMScatterV(vm, EI_nodeIDs_all, nNodes_arr, &
         sendOffsets, EI_nodeIDs, nNodes, 0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (ALLOCATED(sendOffsets))         DEALLOCATE(sendOffsets)
    IF (ALLOCATED(nNodes_arr))          DEALLOCATE(nNodes_arr)
    IF (ALLOCATED(IDs_global_all))      DEALLOCATE(IDs_global_all)
    IF (ALLOCATED(EI_nodeIDs_all))      DEALLOCATE(EI_nodeIDs_all)
    IF (ALLOCATED(EI_nodeIDs_all_temp)) DEALLOCATE(EI_nodeIDs_all_temp)
    IF (ALLOCATED(nodeOwners_all))      DEALLOCATE(nodeOwners_all)
    IF (ALLOCATED(nodeCoords_all))      DEALLOCATE(nodeCoords_all)

  END SUBROUTINE uniquifyGlobalNodeIDs


  !------------------------------------------------------------------------------
  !
  ! Elmer/Ice partitions have local node ids (i.e. node ids start at 1 for each
  ! partition).  But ESMF needs global node ids. So here we calculate a starting 
  ! id for each partition.
  !
  ! Same for elements.
  !
  ! In this function an "item" refers to a node or element
  ! 
  INTEGER FUNCTION firstItemThisPET(numItems,vm)

    TYPE(ESMF_VM),INTENT(IN)   :: vm
    INTEGER, INTENT(IN)        :: numItems

    INTEGER                    :: localPET, PETcount, ii, rc
    INTEGER, ALLOCATABLE       :: firstItemID(:)
    INTEGER, ALLOCATABLE       :: numItemsAllPETS(:), numItemsArr(:)


    CALL ESMF_VMGet(vm, localPet=localPET, petCount=PETcount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) RETURN

    ALLOCATE(firstItemID(PETcount))
    ALLOCATE(numItemsAllPETS(PETcount))
    ALLOCATE(numItemsArr(1))
    numItemsArr = numItems

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! gather the number of items for each PET to PET zero
    CALL ESMF_VMGather(vm, sendData=numItemsArr, recvData=numItemsAllPETS, count=1, rootPet=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (localPET.EQ.0) THEN
       firstItemID(1) = 1 ! first item id on first PET 
       
       ! the starting item for each PET should be the total number of items 
       ! on all previous items plus 1.
       DO ii = 2, PETcount
          firstItemID(ii) = firstItemID(ii-1) + numItemsAllPETS(ii-1)
       END DO

    END IF
       
    CALL ESMF_VMBroadcast(vm, bcstData=firstItemID, count=PETcount, rootPet=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    firstItemThisPET = firstItemID(localPET+1)

    IF (ALLOCATED(firstItemID))     DEALLOCATE(firstItemID)
    IF (ALLOCATED(numItemsAllPETS)) DEALLOCATE(numItemsAllPETS)
    IF (ALLOCATED(numItemsArr))     DEALLOCATE(numItemsArr)

  END FUNCTION firstItemThisPET


  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_ESMF_elementType(Elmer_ElementCode)

    INTEGER, INTENT(IN)          :: Elmer_ElementCode

    SELECT CASE(Elmer_ElementCode)

    CASE(ELMER_ELEMENT_TRIANGLE_LINEAR,ELMER_ELEMENT_TRIANGLE_QUADRAT,ELMER_ELEMENT_TRIANGLE_CUBIC)

       get_ESMF_elementType = ESMF_MESHELEMTYPE_TRI

    CASE(ELMER_ELEMENT_QUADRIL_BILINEAR,ELMER_ELEMENT_QUADRIL_QUADRAT, &
         ELMER_ELEMENT_QUADRIL_QUADRAT2,ELMER_ELEMENT_QUADRIL_CUBIC)

       get_ESMF_elementType = ESMF_MESHELEMTYPE_QUAD

    CASE DEFAULT

       msg = "ESMF element type not found for this Elmer element type "//CHAR(Elmer_ElementCode)
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END SELECT

    RETURN

  END FUNCTION get_ESMF_elementType


  !------------------------------------------------------------------------------
  !
  ! Convert an Elmer mesh to ESMF structures 
  !
  ! Note: this subroutine expects to recieve a 2D Elmer mesh containing triangles 
  ! or and quads.
  ! 
  ! *** This subroutine should be superceded by the more generic ***
  ! *** Elmer2ESMF_mesh as of Nov 2016                           ***
  !
  SUBROUTINE Elmer2ESMF_meshFootprint(Elmer_mesh,ESMF_ElmerMesh,vm,rc)

    TYPE(ESMF_mesh),INTENT(INOUT)    :: ESMF_ElmerMesh
    TYPE(Mesh_t),INTENT(IN)          :: Elmer_Mesh
    TYPE(ESMF_VM),INTENT(IN)         :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    INTEGER                          :: ii, nodeIndex
    CHARACTER(len=ESMF_MAXSTR)       :: subroutineName = "Elmer2ESMF_meshFootprint"
    INTEGER,ALLOCATABLE              :: ESMF_elementTypeList(:),elementIDlist(:),elementIDlist_global(:)
    INTEGER,ALLOCATABLE              :: elemConn(:), nodeIds(:), nodeIds_global(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:) 
    INTEGER                          :: numQuadElems, numTriElems, numTotElems
    INTEGER                          :: localPet, petCount, numElems
    INTEGER                          :: EI_numNodes 


    rc = ESMF_FAILURE

    msg = "Elmer to ESMF mesh format conversion"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__)

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) RETURN

    ! some basic sanity checks
    IF (Elmer_mesh % MeshDim.NE.2) THEN
       msg = "Elmer mesh dimension not equal to 2"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (.NOT.ASSOCIATED(Elmer_mesh % Elements)) THEN
       msg = "Elmer mesh elements not associated"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (.NOT.ASSOCIATED(Elmer_mesh % Nodes)) THEN
       msg = "Elmer mesh nodes not associated"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (SIZE(Elmer_mesh % Elements).NE.Elmer_mesh % NumberOfBulkElements + Elmer_mesh % NumberOfBoundaryElements) THEN
       msg = "Elmer mesh number of elements inconsistency"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (SIZE(Elmer_mesh % Nodes % x) .NE. Elmer_mesh % Nodes % NumberOfNodes) THEN
       msg = "Elmer mesh number of nodes inconsistency"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF


    ! count the numbers of triangles and quadrilaterals in the Elmer mesh, as these are the element types 
    ! we intend to use. Note: this code could be made more efficient if we find it takes up much time. 
    numTriElems = numElementsByType(Elmer_mesh,(/ELMER_ELEMENT_TRIANGLE_LINEAR,ELMER_ELEMENT_TRIANGLE_QUADRAT,&
         ELMER_ELEMENT_TRIANGLE_CUBIC/))

    numQuadElems = numElementsByType(Elmer_mesh,(/ELMER_ELEMENT_QUADRIL_BILINEAR,ELMER_ELEMENT_QUADRIL_QUADRAT,&
         ELMER_ELEMENT_QUADRIL_QUADRAT2,ELMER_ELEMENT_QUADRIL_CUBIC/))


    ! Use number of elements and nodes on this PET, and some collective functions, to assign  
    ! global element and node identifiers.
    ! For now we assume that all nodes are used.  May not always be true.
    EI_numNodes         = Elmer_mesh % Nodes % NumberOfNodes 
    numElems            = numQuadElems+numTriElems
    EI_firstNodeThisPET = firstItemThisPET(EI_numNodes,vm)
    EI_firstElemThisPET = firstItemThisPET(numElems,vm)

    ALLOCATE(nodeIds(EI_numNodes))
    ALLOCATE(nodeIds_global(EI_numNodes))
    ALLOCATE(nodeCoords(EI_numNodes*Elmer_mesh % MeshDim))
    ALLOCATE(nodeOwners(EI_numNodes))

    ALLOCATE(ESMF_elementTypeList(numElems))
    ALLOCATE(elementIDlist(numElems))
    ALLOCATE(elementIDlist_global(numElems))

    ALLOCATE(elemConn(4*numQuadElems+3*numTriElems))

    nodeOwners=localPet

    numTotElems = numElementsByType(Elmer_mesh,(/ELMER_ELEMENT_QUADRIL_BILINEAR,ELMER_ELEMENT_QUADRIL_QUADRAT,&
         ELMER_ELEMENT_QUADRIL_QUADRAT2,ELMER_ELEMENT_QUADRIL_CUBIC,ELMER_ELEMENT_TRIANGLE_LINEAR,&
         ELMER_ELEMENT_TRIANGLE_QUADRAT,ELMER_ELEMENT_TRIANGLE_CUBIC/), &
         ESMF_elementTypeList=ESMF_elementTypeList,elementIDlist=elementIDlist, &
         elemConn = elemConn)

    IF (numTotElems .NE. numQuadElems+numTriElems) THEN
       msg = "Elmer mesh total number of viable elements inconsistency"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    nodeIds        = (/(ii, ii=1, EI_numNodes, 1)/)
    nodeIds_global = (/(ii, ii=EI_firstNodeThisPET, EI_firstNodeThisPET+EI_numNodes-1, 1)/)


    ! loop over nodes to get coords
    DO ii = 1,EI_numNodes
       nodeIndex = (ii-1)*Elmer_mesh%MeshDim+1
       nodeCoords(nodeIndex) = Elmer_mesh % Nodes % x(ii)
       nodeIndex = (ii-1)*Elmer_mesh%MeshDim+2
       nodeCoords(nodeIndex) = Elmer_mesh % Nodes % y(ii)
    END DO


    ! Create Mesh structure in 1 step
    ESMF_ElmerMesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
         nodeIds=nodeIds_global, nodeCoords=nodeCoords, &
         nodeOwners=nodeOwners, elementIds=elementIDlist,&
         elementTypes=ESMF_elementTypeList, elementConn=elemConn, &
         coordSys=ESMF_COORDSYS_CART,                         &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF (ALLOCATED(ESMF_elementTypeList)) DEALLOCATE(ESMF_elementTypeList)
    IF (ALLOCATED(elementIDlist)) DEALLOCATE(elementIDlist)
    IF (ALLOCATED(elementIDlist_global)) DEALLOCATE(elementIDlist_global)
    IF (ALLOCATED(elemConn)) DEALLOCATE(elemConn)
    IF (ALLOCATED(nodeIds)) DEALLOCATE(nodeIds)
    IF (ALLOCATED(nodeIds_global)) DEALLOCATE(nodeIds_global)
    IF (ALLOCATED(nodeCoords)) DEALLOCATE(nodeCoords)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE Elmer2ESMF_meshFootprint
  

  
END MODULE FISOC_ISM_Wrapper

