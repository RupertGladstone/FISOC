
PROGRAM TestNodeOrderingCode
  
  IMPLICIT NONE

  INTEGER, PARAMETER    :: CLOCKWISE=1, ANTI_CLOCKWISE=2
  INTEGER               :: Direction   ! The direction of desired ordering
  INTEGER               :: NumNodes    ! How many nodes this element contains
  INTEGER,ALLOCATABLE   :: NodeIds(:),OrderedNodeIds(:)
  REAL,ALLOCATABLE      :: NodeCoords(:,:)
  REAL                  :: ProjectionVector(3)

  NumNodes = 3
  !ALLOCATE(NodeCoords(3,NumNodes))
  ALLOCATE(NodeCoords(NumNodes,3))
  ALLOCATE(NodeIds(NumNodes))
  ALLOCATE(OrderedNodeIds(NumNodes))
  Direction=ANTI_CLOCKWISE
  Direction=CLOCKWISE
  !NodeIds = (/10, 12, 13/)
  NodeIds = (/9, 12, 13/)
  ProjectionVector = (/1, 0, 0/)
  !Direction=1

!  NodeCoords = RESHAPE((/ &
!       1, 1, 1, &
!       2, 1, 1, &
!       0, 0, 3 & !/), SHAPE(NodeCoords))
!     ,-1, 0,-1/), SHAPE(NodeCoords))

!  NodeCoords = RESHAPE( (/ &
!   1.0,  1.0,  1.0,  1.0, &
!  +1.0, +0.8, -1.4, -1.2, &
!  +1.0, +1.0, -0.9, -1.0 /), SHAPE(NodeCoords))

  NodeCoords = RESHAPE( (/ &
   1.0, -1.0, -1.0,  &
  +1.0, +0.8, -1.2,  &
  +1.0, +1.0, -0.9  /), SHAPE(NodeCoords))

! print*,NodeCoords(4,1:3)

!  NodeCoords = RESHAPE( (/ &
!   1.0,  1.0,  1.0,  1.0, &
!  +1.0, +0.8, -1.4, -1.2, &
!  +1.0, +1.0, -0.9, -1.0 /), SHAPE(NodeCoords))

  CALL  NodeOrdering(Direction, ProjectionVector, NodeIds, NodeCoords,&
                     OrderedNodeIds, NumNodes)

  Print*,"Final ordering is ",OrderedNodeIds
  !call flip(OrderedNodeIds,NumNodes)
  !Print*,"Final ordering is ",OrderedNodeIds  

CONTAINS

  !----------------------------------------------------------------------------------------------
  ! Node ordering routine.  Determines node ordering for triangular and quadrilateral elements.
  !
  ! The ProjectionVector defines the normal to a plane on which to project the node coords.
  ! This plane is assumed to pass through the origin.
  ! It is assumed that the ProjectionVector is a unit vector. TODO: add a check for this
  !  
  SUBROUTINE NodeOrdering(Direction, ProjectionVector, NodeIds, NodeCoords, OrderedNodeIds, nn)

    INTEGER,INTENT(IN)      :: nn !!! NumNodes, length of NodeIds
    INTEGER,INTENT(IN)      :: Direction   ! The direction of desired ordering
    REAL,INTENT(IN)         :: NodeCoords(nn,3)
    REAL,INTENT(IN)         :: ProjectionVector(3)
    INTEGER,INTENT(IN)      :: NodeIds(nn)
    INTEGER,INTENT(OUT)     :: OrderedNodeIds(nn)
    
    !Print*,"First node has id ",NodeIds(1)," and cooords ",NodeCoords(1:3,1)
    !Print*,"Second node has id ",NodeIds(2)," and cooords ",NodeCoords(1:3,2)
    !Print*,"etc..."
    !OrderedNodeIds = (/10, 13, 12/)

    integer          :: sz,ii,jj
    real,allocatable :: theta(:),phi(:),lambda(:),projectedNodes(:,:)
    real,dimension(3):: crs,mcrs,center
    real             :: dd
    integer          :: sect1len,sect2len,sect3len,sect1(99),sect2(99),sect3(99)
    real             :: ang1(99),ang2(99),ang3(99)

    !write(*,"('NodeID :',3I4)") NodeIds
    !write(*,*) NodeCoords

    ALLOCATE(projectedNodes(nn,3))
    ALLOCATE(theta(nn))
    ALLOCATE(phi(nn))    
    ALLOCATE(lambda(nn))

    ! The dot product of the NodeCoord vector and the ProjectionVector gives the magnitude of 
    ! the component of the NodeCoord vector normal to the plane. 
    ! This dot_product*ProjectionVector is node i projected on the ProjectionVector.
    ! projectedNodes(i,:) is node i projected on the plane whose normal vector is ProjectionVector
    ! and which passes through the origin.
    DO ii=1,nn
       projectedNodes(ii,:) = NodeCoords(ii,:) - dot( NodeCoords(ii,:), ProjectionVector ) * ProjectionVector 
    END DO

    !!! Then calculate the center point of these nodes.
    do ii=1,3
      center(ii)=average( projectedNodes(:,ii), nn )
    enddo
    !write(*,*) "center:"
    !write(*,*) center

    !!! projectedNode(ii) minus center
    DO ii=1,nn
       projectedNodes(ii,:) = projectedNodes(ii,:) - center
    END DO    

    call cross( ProjectionVector, projectedNodes(1,:), crs )
    mcrs=-1*crs; !!! crs is first_node x ProjectionVector('x' means cross product)
    !!! mcrs is minus crs

    do ii=1,nn
       call cos_angle( projectedNodes(ii,:), projectedNodes(1,:), theta(ii)  )
       call cos_angle( projectedNodes(ii,:), crs                , phi(ii)    )
       call cos_angle( projectedNodes(ii,:), mcrs               , lambda(ii) )
    enddo

    write(*,"( 'theta:' ,4F10.3 )") theta
    write(*,"( 'phi:'   ,4F10.3 )") phi
    write(*,"( 'lambda:',4F10.3 )") lambda

    !!! To get the order, split the plane into 3 sector,
    sect1len=0 !!! length of sector1. 'sect1' contains indexes of nodes in this sector.
    sect2len=0
    sect3len=0

    do ii=2,nn
      if( theta(ii).ge.0 .and. phi(ii).ge.0 )then

        sect1len=sect1len+1
        sect1( sect1len )=ii
        ang1 ( sect1len )=theta(ii)

      else if( theta(ii).ge.0 .and. lambda(ii).ge.0 )then

        sect3len=sect3len+1
        sect3( sect3len )=ii
        ang3 ( sect3len )=lambda(ii)

      else

        sect2len=sect2len+1
        sect2( sect2len )=ii
        ang2 ( sect2len )=phi(ii)

      endif
    enddo
    
    !!! sort nodes in each sector, based on their angles
    call sort( sect1, ang1, sect1len ) 
    call sort( sect2, ang2, sect2len )
    call sort( sect3, ang3, sect3len )

    !!! Get the array of ordered node IDs.
    OrderedNodeIds(1)=NodeIds(1)
    jj=2
    if(sect1len.ge.1)then
      do ii=1,sect1len
        OrderedNodeIds(jj)=NodeIds( sect1(ii) )
        !write(*,"( 'Id(',I2,')=',I4 )") sect1(i), NodeIds( sect1(ii) )
        write(*,"( 'Sector1 : ',I2 )") sect1(ii)
        jj=jj+1
      enddo
    endif

    if(sect2len.ge.1)then
      do ii=1,sect2len
        OrderedNodeIds(jj)=NodeIds( sect2(ii) )
        !write(*,"( 'Id(',I2,')=',I4 )") sect2(ii), NodeIds( sect2(ii) )
        write(*,"( 'Sector2 : ',I2 )") sect2(ii)
        jj=jj+1
      enddo
    endif

    if(sect3len.ge.1)then
      do ii=1,sect3len
        OrderedNodeIds(jj)=NodeIds( sect3(ii) )
        !write(*,"( 'Id(',I2,')=',I4 )") sect3(ii), NodeIds( sect3(ii) )
        write(*,"( 'Sector3 : ',I2 )") sect3(ii)
        jj=jj+1
      enddo
    endif


    !!! Here we get the order of clockwise. If we want anti-clockwise, flip it!
!    if(Direction.ne.1)then !!! Direction == 1, is clockwise
!      call flip(OrderedNodeIds(2:n),n-1)
!    endif

    SELECT CASE (Direction)
    CASE (CLOCKWISE)
    CASE (ANTI_CLOCKWISE)
      CALL flip(OrderedNodeIds(2:nn),nn-1)
    CASE DEFAULT
       PRINT*,"Direction not recognised, add fatal call here"
    END SELECT
    
  END SUBROUTINE NodeOrdering

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function func(i) result(j)
    integer, intent(in) :: i ! input
    integer             :: j ! output
    j = i**2 + i**3
  end function func


  function average(arr,n) result(r)
    integer,intent(in)::n 
    real,intent(in),dimension(n)::arr 
    real::r
    integer::i
    r=0
    do i=1,n
      r=r+arr(i)
    enddo
    r=r/n
  end function average
  

  REAL FUNCTION dot(u,v)
    real,intent(in)::u(3)
    real,intent(in)::v(3)    

    dot=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
  END FUNCTION dot


  SUBROUTINE cross(u,v,r)
    real,intent(in)::u(3)
    real,intent(in)::v(3)    
    real,intent(out)::r(3)
    r=(/ u(2)*v(3)-u(3)*v(2),&
         u(3)*v(1)-u(1)*v(3),&
         u(1)*v(2)-u(2)*v(1) /)
  END SUBROUTINE cross


  SUBROUTINE cos_angle(u,v,r) !!! calculate the cos angle between two vector.
    real,intent(in)::u(3)
    real,intent(in)::v(3)    
    real,intent(out)::r
    real::len1,len2,d

    len1=sqrt( u(1)**2 + u(2)**2 + u(3)**2 )
    len2=sqrt( v(1)**2 + v(2)**2 + v(3)**2 )
    d = dot(u,v)
    r=d/(len1*len2)
  END SUBROUTINE cos_angle


  SUBROUTINE sort(x,y,n)
    integer,intent(in)::n
    integer,intent(inout),dimension(n)::x
    real,intent(inout),dimension(n)::y    
    !real,intent(out),dimension(n)::r
    integer::i,j,tmp1
    real::tmp

    if( n.ge.2 )then

      do i=1, n-1
        do j=1+i, n
          if( y(j).gt.y(i) )then
            tmp=y(i)
            y(i)=y(j)
            y(j)=tmp
                
            tmp1=x(i)
            x(i)=x(j)
            x(j)=tmp1
          endif
        end do  
      end do

    end if
    !r=x

  END SUBROUTINE sort


  subroutine flip(x,n)
    integer,intent(in)::n
    integer,intent(inout)::x(n)
    integer tmp,i

    if(mod(n,2)==0)then

      do i=1,n/2
        tmp=x(i)
        x(i)=x(n-i+1)
        x(n-i+1)=tmp
      enddo

    else if(mod(n,2)==1)then

      do i=1,(n-1)/2
        tmp=x(i)
        x(i)=x(n-i+1)
        x(n-i+1)=tmp
      enddo

    endif

  end subroutine flip 


END PROGRAM TestNodeOrderingCode














