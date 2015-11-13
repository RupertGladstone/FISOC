
module FISh_MOD

  implicit none
  
  !      flowline ice shelf model
  
  double precision :: secpyr, A, C, L, tend, dt, rhoi
  double precision :: n, m, rhos, grav
  integer :: maxx, fluxBC
  
  !      flow parameter ice and basal sliding
  parameter (secpyr = 365 * 24 * 3600)
  parameter (A = 1.d-25 * secpyr, C = 1.d7)
  
  !      flux condition ocean: fluxBC=1 (water pressure); fluxBC=0 (u=0)
  parameter (fluxBC = 0)
  
  !      domain length
  parameter (L = 200e3)
  
  parameter (maxx = 5, tend = 10.0, dt = .1)
  parameter (rhoi = 900., rhos = 1000., grav = 9.8)
  parameter (n = 3., m = 1./3.)

  !      h = ice thickness
  !      hb = lower boundary of ice sheet/ice shelf
  !      b = bedrock elevation (bathymetry)
  !      x = horizontal distance
  !      s = surface elevation: s = h + hb
  
  double precision,save :: x(maxx), b(maxx), h(maxx), h2(maxx)
  double precision,save :: hb(maxx), s(maxx), hstag(maxx)
  double precision,save :: xstag(maxx), slope(maxx), taud(maxx)
  double precision,save :: gridx, ub(maxx), u(maxx), up(maxx)
  double precision,save :: con(maxx), dtdx2, f(maxx), g(maxx)
  double precision,save :: cen(maxx), txx(maxx), haf(maxx)
  double precision,save :: bstag(maxx), dn(maxx), h0(maxx), mb(maxx)
  real,save :: totH, time
  integer,save :: grl(maxx)
  integer,save :: tc, tl, grlj
  
contains
  
  subroutine FISh_initialize

    integer :: i
    
    !---------   Initialization  -----------------

    mb = 0.
    
    gridx = L / (maxx-1)
    do i=1,maxx
       x(i) = (i-1)*gridx
       b(i) = -900.
       h(i) = -rhos * b(i) / rhoi + 0.1
       if (x(i) > L/2.) then
          h(i) = h(1) - (x(i) - L / 2.) / 1.e2
       end if
       h2(i) = h(i)
       h0(i) = h2(i)
    enddo
    
    tl = nint(tend / dt) + 1

    tc = 0

  end subroutine FISh_initialize


  !----------	  Loop in time ----------------------
  subroutine FISh_run()
    
    integer :: k, i

    tc = tc + 1
print*,"fishRun", tc

    time = (tc - 1) * dt
    !          floating condition for ice sheet geometry
    do i=1,maxx
       haf(i) = b(i) + h(i) * rhoi / rhos
       if (haf(i)<0) then
          hb(i) = (-rhoi) * h(i) / rhos
       else
          hb(i) = b(i)
       end if
       s(i) = hb(i) + h(i)
    enddo
    
    !          variables and floating condition on u-grid
    do i=1,maxx-1
       xstag(i) = (x(i+1) + x(i)) / 2.
       hstag(i) = (h(i+1) + h(i)) / 2.
       bstag(i) = (b(i+1) + b(i)) / 2.
       haf(i) = bstag(i) + hstag(i) * rhoi / rhos
       if (haf(i)>0) then
          grl(i) = 0
       else
          grl(i) = 2
       end if
    enddo
    grl(maxx) = 2
    
    !          determining the grounding line
    grlj = maxx;
    do i=1,maxx-2
       if ((grl(i)==0).and.(grl(i+1)==2)) then
          grl(i) = 1
          grlj = i
       end if
    enddo
    
    !          longitudinal and driving stresses
    do i=1,maxx-1
       slope(i) = (s(i+1) - s(i)) / gridx
       taud(i) = (-rhoi) * grav * hstag(i) * slope(i)
       txx(i) = .25 * rhoi * grav * h(i) * (1. - rhoi / rhos)
    enddo
    txx(maxx) = 0.25 * rhoi * grav * h(maxx) * (1. - rhoi / rhos)
    
    !          initialization of velocity field based on SIA
    if (tc==1) then
       do i=1,maxx-1
          ub(i) = C**(-1./m) * abs(taud(i)**(1./m-1.)) * taud(i) * secpyr
          u(i) = ub(i) + 2./(n+2.) * A * hstag(i) * abs(taud(i))**(n-1.) * taud(i)
       enddo
    end if
    
       !          iteration on implicit scheme (u=f(h2))
    do k=1,3
       if (k==1) then
          if (tc==1) then
!             call ShelfU(u,h,grl,taud,txx,A,C,gridx,secpyr,maxx,m,n,fluxBC)
              call ShelfU(h)
          endif
       else
          call ShelfU(h2)
!          call ShelfU(u,h2,grl,taud,txx,A,C,gridx,secpyr,maxx,m,n,fluxBC)
       end if
       
       !   		   arrange staggered u-grid
       dtdx2 = dt / (2. * gridx);
       do i=2,maxx-1
          up(i) = (-u(i)) * dtdx2
          dn(i) = u(i-1) * dtdx2
          cen(i) = 1. - up(i) - dn(i)
          con(i) = h(i) + mb(i) * dt
       enddo

       
       !              tridiagonal matrix solution
       f(1) = 0
       g(1) = h(2)
       do i=2,maxx-1
          f(i) = up(i) / (cen(i) - dn(i) * f(i-1))
          g(i) = (con(i) + dn(i) * g(i-1)) / (cen(i) - dn(i) * f(i-1))
       enddo
       
       !              boundary conditions
       if (grl(maxx)==2) then
          h2(maxx) = h(maxx-1)
       else
          h2(maxx) = 0
       end if

       do i=maxx-1,2,-1
          h2(i) = g(i) + f(i) * h2(i+1)
       enddo
       h2(1) = h2(3)
    enddo
    
    !          copy new ice thickness on old ones
    totH = 0
    do i=1,maxx
       h(i) = h2(i)
       totH = totH + h(i)
    enddo
    

  end subroutine FISh_run
  !-------------- end of loop in time -------------
  
  
  subroutine FISh_finalize()
  end subroutine FISh_finalize

  
!  subroutine ShelfU(u,h,grl,taud,txx,A,c,gridx,secpyr,maxx,m,n,fluxBC)
   subroutine ShelfU(h_loc)
    
!    integer :: maxx, fluxBC
!    integer :: grl(maxx), k, i
!    double precision :: u(maxx), h(maxx), taud(maxx)
!    double precision :: secpyr, m, n, A, c, gridx, txx(maxx)

    double precision :: h_loc(maxx)
    
    integer :: kk, ii
!    double precision :: up(maxx), dn(maxx), cen(maxx)
!    double precision :: g(maxx)
!    double precision :: f(maxx), con(maxx)
    double precision ::  eeff(maxx), mu(maxx), beta(maxx)
    
    do kk=1,5
       !          ice shelf viscosity
       do ii=2,maxx
          eeff(ii) = ((u(ii) - u(ii-1)) / gridx)**2.
          if (eeff(ii)<1.d-30) then
             eeff(ii) = 1.d-30
          endif
       enddo
       eeff(1) = eeff(3)
       eeff(maxx) = eeff(maxx-1)
       do ii=1,maxx
          mu(ii) = 0.5 * h_loc(ii) * A**(-1./n) * eeff(ii)**((1.-n)/(2.*n))
       enddo
       
       !          calculation of basal friction beta
       do ii=1,maxx
          if (grl(ii)<2) then
             beta(ii) = c * abs(u(ii))**(m-1.)/(secpyr**m)
          else
             beta(ii) = 0.
          end if
       enddo
       beta(1)=beta(2)
       
       do ii=2,maxx-1
          dn(ii) = (-4.) * mu(ii) / (gridx**2.)
          cen(ii) = (-4.) * (mu(ii) + mu(ii+1)) / (gridx**2.) - beta(ii)
          up(ii)= (-4.) * mu(ii+1) / (gridx**2.)
          con(ii) = (-taud(ii))
       enddo
       
       !          tridiagonal matrix solution
       if (fluxBC==1) then
          f(maxx) = 1.
          g(maxx) = (A * txx(maxx)**n) * gridx
       else
          f(maxx) = 0
          g(maxx) = 0
       end if
       do ii=maxx-1,2,-1
          f(ii) = dn(ii) / (cen(ii) - up(ii) * f(ii+1))
          g(ii) = (con(ii) + up(ii) * g(ii+1)) / (cen(ii) - up(ii) * f(ii+1))
       enddo
       u(1) = (-g(2)) / (1. + f(2))
       do ii=2,maxx
          u(ii) = g(ii) + f(ii) * u(ii-1)
       enddo
    enddo
    
    return
    
  end subroutine ShelfU
  
end module FISh_MOD


	   
