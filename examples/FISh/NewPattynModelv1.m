% SSA ice-shelf model
% Staggered grid for u
%
% Frank PATTYN, Laboratoire de Glaciologie, ULB, 2013
% Version 1.0

function [grlj,h] = NewPattynModelv1()

clear all;
close all;

tic;

A=1e-25;
C=1e7;
L=200e3;
maxx=500;
time_end=500;
dt=0.1;

fluxBC=0;

rho_ice=900.; % Ice density
rho_sea=1000.; % Sea water density

% Vectors
gridx=L/(maxx-1); % grid spacing (m)
x=linspace(0,L,maxx)'; % horizontal distance
b=-900-x*0;
h=-rho_sea*b/rho_ice+0.1;
h(x>=L/2)=h(1)-(x(x>=L/2)-L/2)/1.e2;
h2=h;

up=zeros(maxx,1);
dn=zeros(maxx,1);
cen=zeros(maxx,1);
con=zeros(maxx,1);
f=zeros(maxx,1);
g=zeros(maxx,1);
u=zeros(maxx,1);
hstag=zeros(maxx,1);
xstag=zeros(maxx,1);
bstag=zeros(maxx,1);

% Constants
secperyear=365*24*3600;
A=A*secperyear; % Ice flow parameter
grav=9.8; % Gravitation constant
n=3; % Glen index
mb=0.; % Surface mass balance
sea_level=0; % height of eustatic sea level
m=1./3.; % basal sliding exponent

% time stepping

time_lapse = round(time_end/dt)+1;

figure;

for time_count=1:time_lapse
    
    % floating condition for ice sheet geometry determination
    haf=b-sea_level+h*rho_ice/rho_sea; % height above floating
    hb=b;
    hb(haf<0)=sea_level-rho_ice*h(haf<0)/rho_sea;
	s=hb+h;
    
    if time_count==1
        h0=h;
        hb0=hb;
    end
    
    % geometry variables on u-grid
    for j=1:maxx-1
        hstag(j,1)=(h(j+1)+h(j))/2.; % ice thickness on u-grid
        xstag(j,1)=(x(j+1)+x(j))/2.; % distance on u-grid
        bstag(j,1)=(b(j+1)+b(j))/2.;
    end
    hstag(maxx)=NaN;
    xstag(maxx)=NaN;
    bstag(maxx)=NaN;
    
    % floating condition on u-grid
    haf=bstag-sea_level+hstag*rho_ice/rho_sea; % height above floating
    grl=ones(maxx,1)*2; % initialize grl-vector as floating
    grl(haf>0)=0; % grounded
    
     % Subgrid determination of grounding line position on u-grid
    grlj=maxx;
    for j=1:maxx-2
        if grl(j)<.5 && grl(j+1)>1.5
            grl(j)=1; % last grounded grid point grlj
            grlj=j;
        end
    end
    
    % Longitudinal and driving stresses
    slope=diff(s)/gridx; % slope on u-grid
    slope(maxx)=NaN;
    taud=-rho_ice*grav.*hstag.*slope; % driving stress on u-grid
    txx=.25*rho_ice*grav*h.*(1.-rho_ice/rho_sea); % on h-grid
    
   if time_count==1 % initialization of velocity field with SIA
       ub=C^(-1./m)*abs(taud.^(1./m-1.)).*taud*secperyear; % basal velocity (m/year) on u-grid
       u=ub+2./(n+2.)*A*hstag.*abs(taud).^(n-1.).*taud; % ice velocity in ice sheet
   end
   
   for kk=1:3 % Iteration on implicit scheme (u=f(h2))
       if kk==1
           if time_count==1
               u=NewShelfyStream(u,h,grl,taud,txx,A,C,gridx,secperyear,maxx, ...
                   m,n,fluxBC);
           end
       else
           u=NewShelfyStream(u,h2,grl,taud,txx,A,C,gridx,secperyear,maxx, ...
               m,n,fluxBC);
       end
    
       % Arrange staggered grid 
       dtdx2=dt/(2.*gridx);
       for j=2:maxx-1
           up(j)=-u(j)*dtdx2;
           dn(j)=u(j-1)*dtdx2;
           cen(j)=1.-up(j)-dn(j);
           con(j)=h(j)+mb*dt;
       end

       % Tridiagonal matrix solution    
       f(1)=0;
       g(1)=h(2); % boundary condition ice divide (symmetric)
       for j=2:maxx-1
           f(j)=up(j)/(cen(j)-dn(j)*f(j-1));
           g(j)=(con(j)+dn(j)*g(j-1))/(cen(j)-dn(j)*f(j-1));
       end
    
       % Boundary conditions
       if grl(maxx)>1.5
           h2(maxx)=h(maxx-1); % ice shelf
       else
           h2(maxx)=0; % ice sheet
       end
       for j=maxx-1:-1:2
           h2(j)=g(j)+f(j)*h2(j+1);
       end
       h2(1)=h2(3);
   end
   h=h2;
    
    % plotting
    
    if rem(time_count,2)==1
        subplot(1,2,1)
        plot(x,h0+hb0,'-r','linewidth',2); hold on;
        plot(x,hb0,'-r','linewidth',2);
        plot(x,h+hb);
        plot(x,b,'-r','linewidth',3);
        plot(x,hb); hold off;
        grid on;
        pause(0.001);
        subplot(1,2,2);
        plot(x,u); % hold on;
        grid on;
        disp([time_count sum(h)]);
    end
   
    
end %% end of time stepping

toc;

save toto;

end


%%%%% Subroutines


function u = NewShelfyStream(u,h,grl,taud,txx,A,c,gridx,secperyear, ...
    maxx,m,n,fluxBC)

up=zeros(maxx,1);
dn=zeros(maxx,1);
cen=zeros(maxx,1);
con=zeros(maxx,1);
f=zeros(maxx,1);
g=zeros(maxx,1);
eeff=zeros(maxx,1);

for i=1:5
    eeff(2:maxx,1)=max((diff(u)/gridx).^2,1e-30); % eeff on h-grid
    eeff(1)=eeff(3);
    eeff(maxx)=eeff(maxx-1);
    mu=0.5*h*A^(-1./n).*eeff.^((1.-n)/(2.*n)); % mu on h-grid
    
    % Calculation of beta
    beta=zeros(maxx,1);
    beta(grl<2)=c*abs(u(grl<2)).^(m-1.)/(secperyear^m);
   
    for j=2:maxx-1
        dn(j)=-4.*mu(j)/(gridx^2.);
        cen(j)=-4.*(mu(j)+mu(j+1))/(gridx^2.)-beta(j);
        up(j)=-4.*mu(j+1)/(gridx^2.);
        con(j)=-taud(j);
    end
    beta(1)=beta(2);
    
    % Tridiagonal matrix solution   
    if fluxBC==1
        f(maxx)=1.;
        g(maxx)=(A*txx(maxx)^n)*gridx; % boundary condition on u-grid
    else
        f(maxx)=0;
        g(maxx)=0;
    end
    for j=maxx-1:-1:2
        f(j)=dn(j)/(cen(j)-up(j)*f(j+1));
        g(j)=(con(j)+up(j)*g(j+1))/(cen(j)-up(j)*f(j+1));
    end
    u(1)=-g(2)/(1.+f(2));
    for j=2:maxx
        u(j)=g(j)+f(j)*u(j-1);
    end
end

end


%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%


