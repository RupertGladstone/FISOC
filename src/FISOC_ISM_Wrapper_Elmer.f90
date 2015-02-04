
MODULE FISOC_ISM_Wrapper

  USE FISOC_Elmer_types
  USE ESMF
  USE ElmerSolver_mod
  USE MainUtils

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_ISM_Wrapper_Init,  FISOC_ISM_Wrapper_Run, FISOC_ISM_Wrapper_Finalize

  CHARACTER(len=ESMF_MAXSTR) :: msg

  INTEGER                    :: rc

  ! Note that CurrentModel is shared through the Types module

CONTAINS

  ! This initialisation wrapper aims to convert the Elmer mesh and required variables 
  ! to the ESMF formats.  It also performs simple sanity/consistency checks.
  SUBROUTINE FISOC_ISM_Wrapper_Init(FISOC_ElmerFieldBundle,FISOC_ElmerMesh,FISOC_config)

    TYPE(ESMF_fieldBundle),INTENT(INOUT) :: FISOC_ElmerFieldBundle
    TYPE(ESMF_mesh),INTENT(INOUT)        :: FISOC_ElmerMesh
    TYPE(ESMF_config),INTENT(IN)         :: FISOC_config

    TYPE(Mesh_t)                         :: Elmer_Mesh


!double check that the sif really does specify to extrude the mesh, and that the non-extruded mesh really is 2d.

!do some consistency checks between elmer and fisoc config files: time stepping mainly

! get elmer input params, especially time stepping information

! get elmer variables list, recieve esmf required variables list (from esmf_config?), 
! check the elmer contains those variables
! (use a look up table? or lookup table embedded in esmf config object?).

!convert required variables to esmf format

!store required vars with mesh in export state 

!sanity checks: timestep size, number of time intervals to run for

    CALL ElmerSolver_init(Elmer_Mesh) ! Intended to return the mesh prior to extrusion (need to add check for this mesh)

    CALL Elmer2ESMF_mesh(Elmer_mesh,FISOC_ElmerMesh)

!note: variables tobe input (basal melt rate) should be defined (perhaps as exported vars) in the 
! sif.  These also to be checked for their presence against a list of required vars from ESMF

  END SUBROUTINE FISOC_ISM_Wrapper_Init
  
  SUBROUTINE FISOC_ISM_Wrapper_Run()

! get hold of the elmer variables for receiving inputs, and convert them here from esmf to elmer type.
    CALL ElmerSolver_run()
! get hold of list of required variables from Elmer and convert them here from elmer to esmf type.

  END SUBROUTINE FISOC_ISM_Wrapper_Run

  SUBROUTINE FISOC_ISM_Wrapper_Finalize()

    CALL ElmerSolver_finalize()

  END SUBROUTINE FISOC_ISM_Wrapper_Finalize

  !------------------------------------------------------------------------------
  !
  ! Convert an Elmer mesh to ESMF structures 
  !
  ! Note: this subroutine expects to recieve a 2D Elmer mesh containing triangles and 
  ! quads.
  !
  SUBROUTINE Elmer2ESMF_mesh(Elmer_mesh,ESMF_ElmerMesh)
    TYPE(ESMF_mesh),INTENT(INOUT)    :: ESMF_ElmerMesh
    TYPE(Mesh_t),INTENT(IN)          :: Elmer_Mesh

    INTEGER                          :: ii, nodeIndex
    CHARACTER(len=ESMF_MAXSTR)       :: subroutineName = "Elmer2ESMF_mesh"
    INTEGER,ALLOCATABLE              :: ESMF_elementTypeList(:),elementIDlist(:)

    ! ESMF mesh vars
    INTEGER,ALLOCATABLE              :: nodeOwners(:),elemIds(:), elemTypes(:)
    INTEGER,ALLOCATABLE              :: elemConn(:), nodeIds(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:) 
    INTEGER                          :: numNodes, numQuadElems, numTriElems, numTotElems


    msg = "Starting Elmer to ESMF mesh format conversion"
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

END MODULE FISOC_ISM_Wrapper

