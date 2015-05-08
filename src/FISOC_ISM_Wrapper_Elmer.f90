
MODULE FISOC_ISM_Wrapper

  USE ESMF
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  USE ElmerSolver_mod
  USE MainUtils

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_ISM_Wrapper_Init_Phase1,  FISOC_ISM_Wrapper_Init_Phase2,  &
       FISOC_ISM_Wrapper_Run, FISOC_ISM_Wrapper_Finalize

  ! Note that CurrentModel is shared through the Types module (via MainUtils)

CONTAINS

  !--------------------------------------------------------------------------------------
  ! This initialisation wrapper aims to convert the Elmer mesh and required variables 
  ! to the ESMF formats.  It also performs simple sanity/consistency checks.
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1(ISM_ReqVarList,ISM_ExpFB,ISM_mesh,&
       FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(IN)          :: FISOC_config
    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: ISM_ReqVarList(:)
    TYPE(ESMF_VM),INTENT(INOUT)           :: vm

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_mesh),INTENT(OUT)           :: ISM_mesh
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    TYPE(Mesh_t)                          :: Elmer_Mesh

! TODO:
! -get access to the sif
! -double check that the sif really does specify to extrude the mesh, and that the non-extruded mesh really is 2d.
! -do some consistency checks between elmer and fisoc config files: time stepping mainly
! -get elmer variables list, recieve esmf required variables list
! -check the elmer contains those variables (also get Elmer names list from config for variable name checking
! -convert required variables to esmf format (a new subroutine fr this, to be called in init and run)
! -store required vars with mesh in export state (higher level wrapper can manage the states, here we just care about field bundles)

! note: variables tobe input to Elmer (basal melt rate) should be defined (perhaps as exported vars) in the 
! sif.  These also to be checked for their presence against a list of required vars from ESMF

    rc = ESMF_FAILURE

    CALL Initialise_Elmer_ParEnv(vm,rc=rc)

    CALL ElmerSolver_init(Elmer_Mesh,.TRUE.) 
    ! It is intended that ElmerSolver_init should return the mesh prior to extrusion 

    CALL Elmer2ESMF_mesh(Elmer_mesh,ISM_mesh,rc=rc)

    CALL FISOC_populateFieldBundle(ISM_ReqVarList,ISM_ExpFB,ISM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1
  


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2(ISM_ImpFB,FISOC_config,localPet,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ImpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    INTEGER,INTENT(IN)                    :: localPet

    LOGICAL                               :: verbose_coupling

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Run(FISOC_config,localPet,ISM_ImpFB,ISM_ExpFB,rc)

    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_config)      :: FISOC_config
    INTEGER,INTENT(IN)     :: localPet

    INTEGER,INTENT(OUT),OPTIONAL :: rc

    TYPE(ESMF_field)             :: OM_dBdt_l0, ISM_z_l0, ISM_z_l1
    REAL(ESMF_KIND_R8),POINTER   :: OM_dBdt_l0_ptr(:),ISM_z_l0_ptr(:),ISM_z_l1_ptr(:)
    LOGICAL                      :: verbose_coupling

    rc = ESMF_FAILURE

! make sure to run only one timestep

! get hold of the elmer variables for receiving inputs, and convert them here from esmf to elmer type.

    CALL ElmerSolver_run()
! get hold of list of required variables from Elmer and convert them here from elmer to esmf type.

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Finalize(FISOC_config,localPet,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL       :: rc
    INTEGER,INTENT(IN)                 :: localPet

    rc = ESMF_FAILURE

    CALL ElmerSolver_finalize()

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Finalize


  !------------------------------------------------------------------------------
  !
  ! Convert an Elmer mesh to ESMF structures 
  !
  ! Note: this subroutine expects to recieve a 2D Elmer mesh containing triangles and 
  ! quads.
  !
  SUBROUTINE Elmer2ESMF_mesh(Elmer_mesh,ESMF_ElmerMesh,rc)
    TYPE(ESMF_mesh),INTENT(INOUT)    :: ESMF_ElmerMesh
    TYPE(Mesh_t),INTENT(IN)          :: Elmer_Mesh
    INTEGER,INTENT(OUT),OPTIONAL     :: rc

    INTEGER                          :: ii, nodeIndex
    CHARACTER(len=ESMF_MAXSTR)       :: subroutineName = "Elmer2ESMF_mesh"
    INTEGER,ALLOCATABLE              :: ESMF_elementTypeList(:),elementIDlist(:)

    ! ESMF mesh vars
    INTEGER,ALLOCATABLE              :: nodeOwners(:),elemIds(:), elemTypes(:)
    INTEGER,ALLOCATABLE              :: elemConn(:), nodeIds(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:) 
    INTEGER                          :: numNodes, numQuadElems, numTriElems, numTotElems

    rc = ESMF_FAILURE

    msg = "Elmer to ESMF mesh format conversion"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__)

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

    ALLOCATE(ESMF_elementTypeList(numQuadElems+numTriElems))
    ALLOCATE(elementIDlist(numQuadElems+numTriElems))
    ALLOCATE(elemConn(4*numQuadElems+3*numTriElems))

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

    !for now we assume that all nodes are used.  may not always be true.
    numNodes =  Elmer_mesh % Nodes % NumberOfNodes 
    ALLOCATE(nodeIds(numNodes))
    ALLOCATE(nodeCoords(numNodes*Elmer_mesh % MeshDim))
    ALLOCATE(nodeOwners(numNodes))

    nodeOwners=0 ! everything on PET 0

    nodeIds = (/(ii, ii=1, numNodes, 1)/)

    ! loop over nodes to get coords
    DO ii = 1,numNodes
       nodeIndex = (ii-1)*Elmer_mesh%MeshDim+1
       nodeCoords(nodeIndex) = Elmer_mesh % Nodes % x(ii)
       nodeIndex = (ii-1)*Elmer_mesh%MeshDim+2
       nodeCoords(nodeIndex) = Elmer_mesh % Nodes % y(ii)
    END DO

    ! Create Mesh structure in 1 step
    ESMF_ElmerMesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
         nodeIds=nodeIds, nodeCoords=nodeCoords, &
         nodeOwners=nodeOwners, elementIds=elementIDlist,&
         elementTypes=ESMF_elementTypeList, elementConn=elemConn, &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    
    DEALLOCATE(ESMF_elementTypeList)
    DEALLOCATE(elementIDlist)
    DEALLOCATE(elemConn)
    DEALLOCATE(nodeIds)
    DEALLOCATE(nodeCoords)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE Elmer2ESMF_mesh


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
  
END MODULE FISOC_ISM_Wrapper

