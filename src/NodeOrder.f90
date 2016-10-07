
PROGRAM TestNodeOrderingCode
  
  IMPLICIT NONE

  INTEGER, PARAMETER    :: CLOCKWISE=1, ANTI_CLOCKWISE=2
  INTEGER               :: Direction   ! The direction of desired ordering
  INTEGER               :: NumNodes    ! How many nodes this element contains
  INTEGER,ALLOCATABLE   :: NodeIds(:),OrderedNodeIds(:)
  REAL,ALLOCATABLE      :: NodeCoords(:,:)
  REAL                  :: ProjectionVector(3)

  NumNodes = 4
  !ALLOCATE(NodeCoords(3,NumNodes))
  ALLOCATE(NodeCoords(NumNodes,3))
  ALLOCATE(NodeIds(NumNodes))
  ALLOCATE(OrderedNodeIds(NumNodes))
  Direction=ANTI_CLOCKWISE
  !NodeIds = (/10, 12, 13/)
  NodeIds = (/9, 12, 13, 20/)
  ProjectionVector = (/1, 0, 0/)
  !Direction=1

!  NodeCoords = RESHAPE((/ &
!       1, 1, 1, &
!       2, 1, 1, &
!       0, 0, 3 & !/), SHAPE(NodeCoords))
!     ,-1, 0,-1/), SHAPE(NodeCoords))

  NodeCoords = RESHAPE( (/ &
   1.0,  1.0,  1.0,  1.0, &
  +1.0, +0.8, -1.4, -1.2, &
  +1.0, +1.0, -0.9, -1.0 /), SHAPE(NodeCoords))

  !write(*,*) NodeCoords

  CALL  NodeOrdering(Direction, ProjectionVector, NodeIds, NodeCoords,&
                     OrderedNodeIds, NumNodes)

  Print*,"Final ordering is ",OrderedNodeIds
  !call flip(OrderedNodeIds,NumNodes)
  !Print*,"Final ordering is ",OrderedNodeIds  

CONTAINS

  
  SUBROUTINE NodeOrdering(Direction, ProjectionVector, NodeIds, NodeCoords, OrderedNodeIds,n)

    integer,intent(in)      :: n !!! NumNodes, length of NodeIds
    INTEGER,INTENT(IN)      :: Direction   ! The direction of desired ordering
    REAL,INTENT(IN)         :: NodeCoords(n,3)
    REAL,INTENT(IN)         :: ProjectionVector(3)
    INTEGER,INTENT(IN)      :: NodeIds(n)
    INTEGER,INTENT(OUT)     :: OrderedNodeIds(n)
    
    !Print*,"First node has id ",NodeIds(1)," and cooords ",NodeCoords(1:3,1)
    !Print*,"Second node has id ",NodeIds(2)," and cooords ",NodeCoords(1:3,2)
    !Print*,"etc..."
    !OrderedNodeIds = (/10, 13, 12/)

    integer::sz,i,j
    real,allocatable::theta(:),phi(:),lambda(:),projectedNodes(:,:)
    real,dimension(3)::crs,mcrs
    real d
    integer sect1len,sect2len,sect3len,sect1(99),sect2(99),sect3(99)
    real ang1(99),ang2(99),ang3(99)

    !write(*,"('NodeID :',3I4)") NodeIds
    !write(*,*) NodeCoords
    !sz=size(NodeIds)

    allocate(projectedNodes(n,3))
    allocate(theta(n))
    allocate(phi(n))    
    allocate(lambda(n))

    do i=1,n
      call dot(NodeCoords(i,:),ProjectionVector,d) !!! d is result of dot product
      projectedNodes(i,:)=NodeCoords(i,:)-d*ProjectionVector 
      !!! d*ProjectionVector is projected node_i on ProjectionVector
      !!! projectedNodes(i,:) is projected node_i on the plane whose normal vector is ProjectionVector
    enddo

    call cross( ProjectionVector, projectedNodes(1,:), crs )
    mcrs=-1*crs; !!! crs is first_node x ProjectionVector('x' means cross product)
    !!! mcrs is minus crs

    do i=1,n
      call angle( projectedNodes(i,:), projectedNodes(1,:), theta(i) )
      call angle( projectedNodes(i,:), crs,                 phi(i) )
      call angle( projectedNodes(i,:), mcrs,                lambda(i) )
    enddo

    write(*,"( 'theta:' ,4F10.3 )") theta
    write(*,"( 'phi:'   ,4F10.3 )") phi
    write(*,"( 'lambda:',4F10.3 )") lambda

    !!! To get the order, split the plane into 3 sector,
    sect1len=1 !!! length of sector1. 'sect1' contains indexes of nodes in this sector.
    sect2len=1
    sect3len=1

    do i=2,n
      if( theta(i).ge.0 .and. phi(i).ge.0 )then

        sect1( sect1len )=i
        ang1 ( sect1len )=theta(i)
        sect1len=sect1len+1

      else if( theta(i).ge.0 .and. lambda(i).ge.0 )then

        sect3( sect3len )=i
        ang3 ( sect3len )=lambda(i)
        sect3len=sect3len+1

      else

        sect2( sect2len )=i
        ang2 ( sect2len )=phi(i)
        sect2len=sect2len+1

      endif
    enddo
    
    !!! sort nodes in each sector, based on their angles
    call sort( sect1, ang1, sect1len-1 ) 
    call sort( sect2, ang2, sect2len-1 )
    call sort( sect3, ang3, sect3len-1 )

    !!! Get the array of ordered node IDs.
    OrderedNodeIds(1)=NodeIds(1)
    j=2
    if(sect1len.ge.2)then
      do i=1,sect1len-1
        OrderedNodeIds(j)=NodeIds( sect1(i) )
        !write(*,"( 'Id(',I2,')=',I4 )") sect1(i), NodeIds( sect1(i) )
        write(*,"( 'Sector1 : ',I2 )") sect1(i)
        j=j+1
      enddo
    endif

    if(sect2len.ge.2)then
      do i=1,sect2len-1
        OrderedNodeIds(j)=NodeIds( sect2(i) )
        !write(*,"( 'Id(',I2,')=',I4 )") sect2(i), NodeIds( sect2(i) )
        write(*,"( 'Sector2 : ',I2 )") sect2(i)
        j=j+1
      enddo
    endif

    if(sect3len.ge.2)then
      do i=1,sect3len-1
        OrderedNodeIds(j)=NodeIds( sect3(i) )
        !write(*,"( 'Id(',I2,')=',I4 )") sect3(i), NodeIds( sect3(i) )
        write(*,"( 'Sector3 : ',I2 )") sect3(i)
        j=j+1
      enddo
    endif


    !!! Here we get the order of clockwise. If we want anti-clockwise, flip it!
    if(Direction.ne.1)then !!! Direction == 1, is clockwise
      call flip(OrderedNodeIds(2:n),n-1)
    endif


  END SUBROUTINE NodeOrdering

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE dot(u,v,r)
    real,intent(in)::u(3)
    real,intent(in)::v(3)    
    real,intent(out)::r
    !write(*,*) u
    r=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
  END SUBROUTINE dot
  

  SUBROUTINE cross(u,v,r)
    real,intent(in)::u(3)
    real,intent(in)::v(3)    
    real,intent(out)::r(3)
    r=(/ u(2)*v(3)-u(3)*v(2),&
         u(3)*v(1)-u(1)*v(3),&
         u(1)*v(2)-u(2)*v(1) /)
  END SUBROUTINE cross


  SUBROUTINE angle(u,v,r)
    real,intent(in)::u(3)
    real,intent(in)::v(3)    
    real,intent(out)::r
    real::len1,len2,d

    len1=sqrt( u(1)**2 + u(2)**2 + u(3)**2 )
    len2=sqrt( v(1)**2 + v(2)**2 + v(3)**2 )
    call dot(u,v,d)
    r=d/(len1*len2)
    !r=d
  END SUBROUTINE angle


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














