       program IShelf
	   
c      flowline ice shelf model
	   
       double precision :: secpyr, A, C, L, tend, dt, rhoi
       double precision :: n, m, mb, rhos, grav
       integer :: maxx, fluxBC
	   
c      flow parameter ice and basal sliding
       parameter (secpyr = 365 * 25 * 3600)
       parameter (A = 1.d-25 * secpyr, C = 1.d7)
	   
c      flux condition ocean: fluxBC=1 (water pressure); fluxBC=0 (u=0)
       parameter (fluxBC = 0)
	   
c      domain length
       parameter (L = 200e3)
	   
       parameter (maxx = 500, tend = 500.0, dt = .1)
       parameter (rhoi = 900., rhos = 1000., grav = 9.8)
       parameter (n = 3., mb = 0., m = 1./3.)

c      h = ice thickness
c      hb = lower boundary of ice sheet/ice shelf
c      b = bedrock elevation (bathymetry)
c      x = horizontal distance
c      s = surface elevation: s = h + hb

       double precision :: x(maxx), b(maxx), h(maxx), h2(maxx)
       double precision :: hb(maxx), s(maxx), hstag(maxx)
       double precision :: xstag(maxx), slope(maxx), taud(maxx)
       double precision :: gridx, ub(maxx), u(maxx), up(maxx)
       double precision :: con(maxx), dtdx2, f(maxx), g(maxx)
       double precision :: cen(maxx), txx(maxx), haf(maxx)
       double precision :: bstag(maxx), dn(maxx), h0(maxx)
       real :: totH, time
       integer :: grl(maxx)
       integer :: tl, tc, grlj, i
	   
c---------   Initialization  -----------------
       
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
	   
c----------	  Loop in time ----------------------

       do tc=1,tl
           time = (tc - 1) * dt
c          floating condition for ice sheet geometry
           do i=1,maxx
               haf(i) = b(i) + h(i) * rhoi / rhos
               if (haf(i)<0) then
                   hb(i) = (-rhoi) * h(i) / rhos
               else
                   hb(i) = b(i)
               end if
               s(i) = hb(i) + h(i)
           enddo
		   
c          variables and floating condition on u-grid
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
		   
c          determining the grounding line
           grlj = maxx;
           do i=1,maxx-2
               if ((grl(i)==0).and.(grl(i+1)==2)) then
                   grl(i) = 1
                   grlj = i
               end if
           enddo
		   
c          longitudinal and driving stresses
           do i=1,maxx-1
               slope(i) = (s(i+1) - s(i)) / gridx
               taud(i) = (-rhoi) * grav * hstag(i) * slope(i)
               txx(i) = .25 * rhoi * grav * h(i) * (1. - rhoi
     .			      / rhos)
           enddo
           txx(maxx) = 0.25 * rhoi * grav * h(maxx) * (1. -
     .        rhoi / rhos)
	 
c          initialization of velocity field based on SIA
           if (tc==1) then
               do i=1,maxx-1
                   ub(i) = C**(-1./m) * abs(taud(i)**(1./m-1.)) *
     .                taud(i) * secpyr
                   u(i) = ub(i) + 2./(n+2.) * A * hstag(i) *
     .                abs(taud(i))**(n-1.) * taud(i)
               enddo
           end if
		   
c          iteration on implicit scheme (u=f(h2))
           do k=1,3
               if (k==1) then
                   if (tc==1) then
                       call ShelfU(u,h,grl,taud,txx,A,C,gridx,
     .                    secpyr,maxx,m,n,fluxBC)
                   endif
               else
                   call ShelfU(u,h2,grl,taud,txx,A,C,gridx,secpyr,
     .                maxx,m,n,fluxBC)
               end if
		   
c   		   arrange staggered u-grid
               dtdx2 = dt / (2. * gridx);
               do i=2,maxx-1
                   up(i) = (-u(i)) * dtdx2
                   dn(i) = u(i-1) * dtdx2
                   cen(i) = 1. - up(i) - dn(i)
                   con(i) = h(i) + mb * dt
               enddo
           
c              tridiagonal matrix solution
               f(1) = 0
               g(1) = h(2)
               do i=2,maxx-1
                   f(i) = up(i) / (cen(i) - dn(i) * f(i-1))
                   g(i) = (con(i) + dn(i) * g(i-1)) /
     .                (cen(i) - dn(i) * f(i-1))
               enddo
		   
c              boundary conditions
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
		   
c          copy new ice thickness on old ones
           totH = 0
           do i=1,maxx
               h(i) = h2(i)
               totH = totH + h(i)
               if (tc==tl) print *, i, h(i), u(i)
           enddo
           print *, totH
       enddo
	   
c-------------- end of loop in time -------------

       end program IShelf
	   
	   
	   
	   
c-----------  SUBROUTINES ----------------------
	   
       subroutine ShelfU(u,h,grl,taud,txx,A,c,gridx,secpyr,
     .    maxx,m,n,fluxBC)
	
          integer :: maxx, fluxBC
          integer :: grl(maxx)
          double precision :: u(maxx), h(maxx), taud(maxx)
          double precision :: secpyr, m, n, A, c, gridx, txx(maxx)
	   
       double precision :: up(maxx), dn(maxx), cen(maxx)
       double precision :: g(maxx), eeff(maxx), mu(maxx)
       double precision :: f(maxx), con(maxx), beta(maxx)
	   
       do k=1,5
c          ice shelf viscosity
           do i=2,maxx
               eeff(i) = ((u(i) - u(i-1)) / gridx)**2.
			   if (eeff(i)<1.d-30) then
                   eeff(i) = 1.d-30
               endif
           enddo
           eeff(1) = eeff(3)
           eeff(maxx) = eeff(maxx-1)
           do i=1,maxx
               mu(i) = 0.5 * h(i) * A**(-1./n) * eeff(i)**
     .           ((1.-n)/(2.*n))
           enddo
		   
c          calculation of basal friction beta
           do i=1,maxx
               if (grl(i)<2) then
                   beta(i) = c * abs(u(i))**(m-1.)/(secpyr**m)
               else
                   beta(i) = 0.
               end if
           enddo
           beta(1)=beta(2)
		   
           do i=2,maxx-1
               dn(i) = (-4.) * mu(i) / (gridx**2.)
               cen(i) = (-4.) * (mu(i) + mu(i+1)) / (gridx**2.)
     .            - beta(i)
               up(i)= (-4.) * mu(i+1) / (gridx**2.)
               con(i) = (-taud(i))
           enddo
	   
c          tridiagonal matrix solution
           if (fluxBC==1) then
               f(maxx) = 1.
               g(maxx) = (A * txx(maxx)**n) * gridx
           else
               f(maxx) = 0
               g(maxx) = 0
           end if
           do i=maxx-1,2,-1
               f(i) = dn(i) / (cen(i) - up(i) * f(i+1))
               g(i) = (con(i) + up(i) * g(i+1)) / (cen(i) -
     .            up(i) * f(i+1))
           enddo
           u(1) = (-g(2)) / (1. + f(2))
           do i=2,maxx
               u(i) = g(i) + f(i) * u(i-1)
           enddo
       enddo
       return
       end
	   
	   
	   
	   