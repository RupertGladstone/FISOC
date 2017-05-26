function [x,y]=polar_stereo_deluxe(alat,alon,xcentre,ycentre,gridsize,ref_lat)

%common map, mflag, slat, slon, x0, y0, re, e, e2, pi, cdr

%slat = -71.0
slat=ref_lat;
slon = 0.0;
%x0 = 3400
x0=0.;
y0=0.;
%y0 = 3400
re = 6378.137;
e2 = 6.694379852e-3;
e = sqrt(e2);
pi = 3.141592654;
cdr = pi/180.0;
mflag = -1;

%;
%function mapll, alat, alon

%common map, mflag, slat, slon, x0, y0, re, e, e2, pi, cdr
%common mapll, llflag, tc, mc, rho

%if (n_elements(mflag) eq 0) then map_init
%if (n_elements(llflag) eq 0) then begin
   if (abs(slat) ==  90)
      rho = 2.*re./((1+e).^(1+e).*(1-e).^(1-e)).^(e/2);
   else
      sl = abs(slat).*cdr;
      tc = tan(pi/4-sl/2)./((1-e.*sin(sl))./(1+e.*sin(sl))).^(e./2);
      mc = cos(sl)./sqrt(1-e2.*(sin(sl).^2));
      rho = re.*mc./tc;
   end
%   llflag = -1
%endif

lat = abs(alat).*cdr;
t = tan(pi/4-lat/2)./((1-e.*sin(lat))./(1+e.*sin(lat))).^(e/2);
lon = -(alon-slon).*cdr;
x = (x0 - (rho.*t.*sin(lon)))./gridsize+xcentre;
y = (y0 + (rho.*t.*cos(lon)))./gridsize+ycentre;
