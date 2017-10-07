
PROGRAM TestNodeOrderingCode
  
  IMPLICIT NONE

  INTEGER, PARAMETER    :: CLOCKWISE=1, ANTI_CLOCKWISE=2
  INTEGER               :: Direction   ! The direction of desired ordering
  INTEGER               :: NumNodes    ! How many nodes this element contains
  INTEGER,ALLOCATABLE   :: NodeIds(:),OrderedNodeIds(:)
  REAL,ALLOCATABLE      :: NodeCoords(:,:)
  REAL                  :: ProjectionVector(3)

  NumNodes = 3
  ALLOCATE(NodeCoords(3,NumNodes))
  ALLOCATE(NodeIds(NumNodes))
  ALLOCATE(OrderedNodeIds(NumNodes))
  Direction=ANTI_CLOCKWISE
  NodeIds = (/10, 12, 13/)

  NodeCoords = RESHAPE((/ &
       1, 1, 1, &
       2, 1, 1, &
       0, 0, 3 /), SHAPE(NodeCoords))
  ProjectionVector = (/1, 0, 0/)

  CALL  NodeOrdering(Direction, ProjectionVector, NodeIds, NodeCoords, OrderedNodeIds)

  Print*,"Final ordering is ",OrderedNodeIds
  
CONTAINS

  
  SUBROUTINE NodeOrdering(Direction, ProjectionVector, NodeIds, NodeCoords, OrderedNodeIds)

    INTEGER,INTENT(IN)      :: Direction   ! The direction of desired ordering
    REAL,INTENT(IN)         :: NodeCoords(:,:)
    REAL,INTENT(IN)         :: ProjectionVector(3)
    INTEGER,INTENT(IN)      :: NodeIds(:)
    INTEGER,INTENT(OUT)     :: OrderedNodeIds(:)
    

    Print*,"First node has id ",NodeIds(1)," and cooords ",NodeCoords(1:3,1)
    Print*,"Second node has id ",NodeIds(2)," and cooords ",NodeCoords(1:3,2)
    Print*,"etc..."

    OrderedNodeIds = (/10, 13, 12/)

  END SUBROUTINE NodeOrdering
  
END PROGRAM TestNodeOrderingCode

