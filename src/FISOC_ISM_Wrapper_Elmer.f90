
MODULE FISOC_ISM_Wrapper

  USE ESMF
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  USE ElmerSolver_mod
  USE MainUtils
  USE Messages, ONLY : MessageUnit

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_ISM_Wrapper_Init_Phase1,  FISOC_ISM_Wrapper_Init_Phase2,  &
       FISOC_ISM_Wrapper_Run, FISOC_ISM_Wrapper_Finalize

  ! Note that CurrentModel is shared through the Types module (via MainUtils)

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

  REAL(ESMF_KIND_R8), PARAMETER :: Elmer_secpyr = 365.25_dp*24.0_dp*60.0_dp*60.0_dp

  ! These variable names are hard coded to match the names given in the Elmer/Ice .sif
  ! The user must ensure the names correspond.
  ! TODO: is this information in the manual?
  CHARACTER(len=ESMF_MAXSTR), PARAMETER :: EIname_dBdt_l0        = 'meltRate'
  CHARACTER(len=ESMF_MAXSTR), PARAMETER :: EIname_temperature_l0 = 'oceanTemperature'
  CHARACTER(len=ESMF_MAXSTR), PARAMETER :: EIname_temperature_l1 = 'oceanTemperature'
  CHARACTER(len=ESMF_MAXSTR), PARAMETER :: EIname_velocity_l0    = 'Velocity'
  CHARACTER(len=ESMF_MAXSTR), PARAMETER :: EIname_z_l0           = 'Coordinate 3'
  CHARACTER(len=ESMF_MAXSTR), PARAMETER :: EIname_z_l1           = 'Coordinate 3'

  INTEGER  :: EI_numNodesAtBed ! how many nodes on the curret PET at the lower surface of the mesh

  ! global node numbering reference (Elmer uses local node numbering).  We must add this number 
  ! to the node id whenever converting from Elmer node ids to ESMF node ids.  Same for element id.
  INTEGER  :: EI_firstNodeThisPET, EI_firstElemThisPET


CONTAINS

  !--------------------------------------------------------------------------------------
  ! This initialisation wrapper aims to convert the Elmer mesh and required variables 
  ! to the ESMF formats.  It also performs simple sanity/consistency checks.
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1(ISM_ReqVarList,ISM_ExpFB,ISM_mesh,&
       FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: ISM_ReqVarList(:)
    TYPE(ESMF_VM),INTENT(INOUT)           :: vm

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_mesh),INTENT(OUT)           :: ISM_mesh
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR)            :: ISM_configFile_FISOC, ISM_stdoutFile
    CHARACTER(len=MAX_STRING_LEN)         :: ISM_configFile_Elmer

    TYPE(Mesh_t)                          :: Elmer_Mesh
    REAL(ESMF_KIND_R8)                    :: Elmer_dt, FISOC_ISM_dt
    INTEGER                               :: localpet

! TODO:
! -double check that the sif really does specify to extrude the mesh, and that the non-extruded mesh really is 2d.
! -do some consistency checks between elmer and fisoc config files: time stepping mainly
! -get elmer variables list, recieve esmf required variables list
! -check the elmer contains those variables (also get Elmer names list from config for variable name checking
! -convert required variables to esmf format (a new subroutine fr this, to be called in init and run)
! -store required vars with mesh in export state (higher level wrapper can manage the states, here we just care about field bundles)

! note: variables to be input to Elmer (basal melt rate) should be defined (perhaps as exported vars) in the 
! sif.  These also to be checked for their presence against a list of required vars from ESMF

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! FISOC tells Elmer which .sif to use
    CALL ESMF_ConfigGetAttribute(FISOC_config, ISM_configFile_FISOC, label='ISM_configFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_configFile_Elmer = ""
    ISM_configFile_Elmer = ISM_configFile_FISOC 

    ! FISOC sets the unit for Elmer's standard messaging routines to use instead of stdout
    CALL ESMF_ConfigGetAttribute(FISOC_config, ISM_stdoutFile, label='ISM_stdoutFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    OPEN(unit=ISM_outputUnit, file=ISM_stdoutFile, STATUS='REPLACE', ERR=101)
    MessageUnit = ISM_outputUnit

    CALL Initialise_Elmer_ParEnv(vm,rc=rc)

    CALL ElmerSolver_init(Elmer_Mesh,.TRUE.,ISM_configFile_Elmer) 
    ! It is intended that ElmerSolver_init should return the mesh prior to Elmer's internal extrusion 
    
    IF (localpet.EQ.0) THEN
       IF ( .NOT.TimeStepConsistent(FISOC_config) ) THEN
          WRITE (msg, "(A,I0,A,I0,A)") &
               "FATAL: Elmer/FISOC timestep inconsistency"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    END IF

    CALL Elmer2ESMF_mesh(Elmer_mesh,ISM_mesh,vm,rc=rc)

    CALL FISOC_populateFieldBundle(ISM_ReqVarList,ISM_ExpFB,ISM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! get hold of list of required variables from Elmer and convert them here from elmer to esmf type.
    CALL getFieldDataFromISM(ISM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

    RETURN

101 msg = "ISM failed to open stdoutFile"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1
  

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
    Elmer_dt_sec = INT(Elmer_secpyr * Elmer_dt)

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

  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2(ISM_ImpFB,ISM_ExpFB,FISOC_config,vm,rc)

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
       PRINT*,"**********      OM wrapper.  Init phase 2 method.        *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we have access to the initialised ISM fields, just in case the OM needs "
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

    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_config)      :: FISOC_config
    TYPE(ESMF_VM)          :: vm

    INTEGER,INTENT(OUT),OPTIONAL :: rc

    INTEGER                      :: localPet
    TYPE(ESMF_field)             :: OM_dBdt_l0, ISM_z_l0, ISM_z_l1
    REAL(ESMF_KIND_R8),POINTER   :: OM_dBdt_l0_ptr(:),ISM_z_l0_ptr(:),ISM_z_l1_ptr(:)
    LOGICAL                      :: verbose_coupling

    rc = ESMF_FAILURE

! TODO:
! make the interface to this have correspondingarguments to the OM side (e.g. pass vm not just localpet)
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

    IF (verbose_coupling) THEN
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
  SUBROUTINE FISOC_ISM_Wrapper_Finalize(FISOC_config,localPet,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL       :: rc
    INTEGER,INTENT(IN)                 :: localPet

    rc = ESMF_FAILURE

    CALL ElmerSolver_finalize()

    CLOSE(unit=ISM_outputUnit, ERR=102)

    rc = ESMF_SUCCESS

    RETURN

102 msg = "ISM failed to close stdoutFile"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_ISM_Wrapper_Finalize



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

    DEALLOCATE(firstItemID)
    DEALLOCATE(numItemsAllPETS)
    DEALLOCATE(numItemsArr)

  END FUNCTION firstItemThisPET

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
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END SELECT

    RETURN

  END FUNCTION get_ESMF_elementType


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
    ParEnv % ActiveComm = MPI_COMM_WORLD ! or mpic_dup
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
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)
    INTEGER                               :: ii, jj, nn, numNodes
    INTEGER,ALLOCATABLE                   :: nodeIds(:)
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
    
    ! set up list of node ids.
    ! TODO: resolve code duplication with elmer2esmf mesh routine for setting up nodeIds array (and partner routine to this one).
!    numNodes =  CurrentModel % mesh % Nodes % NumberOfNodes 
    numNodes = EI_numNodesAtBed
    ALLOCATE(nodeIds(numNodes))
    nodeIds = (/(ii, ii=1, numNodes, 1)/)

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
       IF (FISOC_OM2ISM(fieldName,FISOC_config,rc=rc)) THEN

          SELECT CASE (TRIM(ADJUSTL(fieldName)))
              
          CASE ('OM_dBdt_l0')

             EI_field => VariableGet( CurrentModel % Mesh % Variables, &
                  EIname_dBdt_l0, UnFoundFatal=.TRUE.)
             EI_fieldVals => EI_field % Values
             EI_fieldPerm => EI_field % Perm
             DO ii = 1,numNodes
!                EI_fieldVals(EI_fieldPerm(nodeIds(ii))) = ptr(EI_firstNodeThisPET+ii-1) * FISOC_secPerYear
                EI_fieldVals(EI_fieldPerm(nodeIds(ii))) = ptr(ii) * FISOC_secPerYear
             END DO

          CASE ('OM_temperature_l0')
             EI_field => VariableGet( CurrentModel % Mesh % Variables, &
                  EIname_temperature_l0, UnFoundFatal=.TRUE.)
             EI_fieldVals => EI_field % Values
             EI_fieldPerm => EI_field % Perm
             DO ii = 1,numNodes
!                EI_fieldVals(EI_fieldPerm(nodeIds(ii))) = ptr(EI_firstNodeThisPET+ii-1)
                EI_fieldVals(EI_fieldPerm(nodeIds(ii))) = ptr(ii)
             END DO

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

       END IF

    END DO fieldLoop

    DEALLOCATE(nodeIds)

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
    INTEGER                               :: ii, jj, nn, numNodes
    INTEGER,ALLOCATABLE                   :: nodeIds(:)
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
    
    ! set up list of node ids.
    ! TODO: resolve code duplication with elmer2esmf mesh routine for setting up nodeIds array (and partner routine to this one).
!    numNodes =  CurrentModel % mesh % Nodes % NumberOfNodes 
    numNodes = EI_numNodesAtBed
    ALLOCATE(nodeIds(numNodes))
    nodeIds = (/(ii, ii=1, numNodes, 1)/)

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
          DO ii = 1,numNodes
!              ptr(EI_firstNodeThisPET+ii-1) = EI_fieldVals(nodeIds(ii))
              ptr(ii) = EI_fieldVals(nodeIds(ii))
          END DO
          
       CASE ('ISM_temperature_l0','ISM_temperature_l1','ISM_velocity_l0','ISM_z_l1','ISM_z_l0_previous')
          msg = "WARNING: ignored variable: "//TRIM(ADJUSTL(fieldName))
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__, rc=rc)          
          
       CASE ('ISM_dTdz_l0','ISM_dddt')
          msg = "INFO: not exporting derived variable: "//TRIM(ADJUSTL(fieldName))
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
    
    DEALLOCATE(nodeIds)
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE getFieldDataFromISM
  
  !------------------------------------------------------------------------------
  !
  ! Convert an Elmer mesh to ESMF structures 
  !
  ! Note: this subroutine expects to recieve a 2D Elmer mesh containing triangles 
  ! or and quads.
  !
  SUBROUTINE Elmer2ESMF_mesh(Elmer_mesh,ESMF_ElmerMesh,vm,rc)

    TYPE(ESMF_mesh),INTENT(INOUT)    :: ESMF_ElmerMesh
    TYPE(Mesh_t),INTENT(IN)          :: Elmer_Mesh
    TYPE(ESMF_VM),INTENT(IN)         :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    INTEGER                          :: ii, nodeIndex
    CHARACTER(len=ESMF_MAXSTR)       :: subroutineName = "Elmer2ESMF_mesh"
    INTEGER,ALLOCATABLE              :: ESMF_elementTypeList(:),elementIDlist(:),elementIDlist_global(:)
    INTEGER,ALLOCATABLE              :: nodeOwners(:)
    INTEGER,ALLOCATABLE              :: elemConn(:), nodeIds(:), nodeIds_global(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:) 
    INTEGER                          :: numNodes, numQuadElems, numTriElems, numTotElems
    INTEGER                          :: localPet, petCount, numElems

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
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (.NOT.ASSOCIATED(Elmer_mesh % Elements)) THEN
       msg = "Elmer mesh elements not associated"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (.NOT.ASSOCIATED(Elmer_mesh % Nodes)) THEN
       msg = "Elmer mesh nodes not associated"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (SIZE(Elmer_mesh % Elements).NE.Elmer_mesh % NumberOfBulkElements + Elmer_mesh % NumberOfBoundaryElements) THEN
       msg = "Elmer mesh number of elements inconsistency"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (SIZE(Elmer_mesh % Nodes % x) .NE. Elmer_mesh % Nodes % NumberOfNodes) THEN
       msg = "Elmer mesh number of nodes inconsistency"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF


    ! count the numbers of triangles and quadrilaterals in the Elmer mesh, as these are the element types 
    ! we intend to use. Note: this code could be made more efficient if we find it takes up much time. 
    numTriElems = numElementsByType(Elmer_mesh,(/ELMER_ELEMENT_TRIANGLE_LINEAR,ELMER_ELEMENT_TRIANGLE_QUADRAT,&
         ELMER_ELEMENT_TRIANGLE_CUBIC/))

    numQuadElems = numElementsByType(Elmer_mesh,(/ELMER_ELEMENT_QUADRIL_BILINEAR,ELMER_ELEMENT_QUADRIL_QUADRAT,&
         ELMER_ELEMENT_QUADRIL_QUADRAT2,ELMER_ELEMENT_QUADRIL_CUBIC/))


    ! Use number of elements and ndoes on this PET, and some collective functions, to assign  
    ! global element and node identifiers.
    ! For now we assume that all nodes are used.  May not always be true.
    numNodes            = Elmer_mesh % Nodes % NumberOfNodes 
    numElems            = numQuadElems+numTriElems
    EI_numNodesAtBed    = numNodes
    EI_firstNodeThisPET = firstItemThisPET(numNodes,vm)
    EI_firstElemThisPET = firstItemThisPET(numElems,vm)

    ALLOCATE(nodeIds(numNodes))
    ALLOCATE(nodeIds_global(numNodes))
    ALLOCATE(nodeCoords(numNodes*Elmer_mesh % MeshDim))
    ALLOCATE(nodeOwners(numNodes))

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
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    nodeIds        = (/(ii, ii=1, numNodes, 1)/)
    nodeIds_global = (/(ii, ii=EI_firstNodeThisPET, EI_firstNodeThisPET+numNodes-1, 1)/)


    ! loop over nodes to get coords
    DO ii = 1,numNodes
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
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    DEALLOCATE(ESMF_elementTypeList)
    DEALLOCATE(elementIDlist)
    DEALLOCATE(elementIDlist_global)
    DEALLOCATE(elemConn)
    DEALLOCATE(nodeIds)
    DEALLOCATE(nodeIds_global)
    DEALLOCATE(nodeCoords)
    DEALLOCATE(nodeOwners)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE Elmer2ESMF_mesh


END MODULE FISOC_ISM_Wrapper

