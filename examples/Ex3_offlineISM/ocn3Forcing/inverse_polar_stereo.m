%  INVERSE_POLAR_STEREO
% [alat,alon] = inverse_polar_stereo(x,y,xcentre,ycentre,gridsize,ref_lat)
% this m-function calculates the latitude and longitude of a point
% from its x,y coordinates in the polar stereographic coordinate system
% with standard reference latitude -71.0. It should be the inverse to
% polar_stereo_deluxe.m.  Like that function it is based on the 
% IDL function map_loc.pro
% 20070715 Roland Warner


function [alat,alon] = inverse_polar_stereo(x,y,xcentre,ycentre,gridsize,ref_lat)


%slat = -71.0
slat=ref_lat;
slon = 0.0;
re = 6378.137;
e2 = 6.694379852e-3;
pi = 3.141592654;
e = sqrt(e2);

cdr = pi./180.0;
mflag = -1;


if (abs(slat) ==  90) then 
      % rho = 2.*re./((1+e).^(1+e).*(1-e).^(1-e)).^(e/2);
      % PUZZLE
      % NOTE according to the mapxy function in map_loc.pro 
      rho = 2.*re./sqrt((1+e).^(1+e).*(1-e).^(1-e));
      % for this case 
      % which looks like it could be a typo for the definition in the mapll
      % since the definitions below are identical.
      
   else 
      sl = abs(slat).*cdr;
      tc = tan(pi/4-sl/2)./((1-e.*sin(sl))./(1+e.*sin(sl))).^(e./2);
      mc = cos(sl)./sqrt(1-e2.*(sin(sl).^2));
      rho = re.*mc./tc;
end
   

   a1 =  5./24.*e2.^2;
   a2 =  1./12.*e2.^3;
   a3 =  7./48.*e2.^2;
   a4 = 29./240.*e2.^3;
   a5 =  7./120.*e2.^3;


t = sqrt((x-xcentre).^2+(y-ycentre).^2)./rho;
chi = (pi./2) - 2.*atan(t);
alat = chi + ((e2./2)+a1+a2).*sin(2.*chi) + (a3+a4).*sin(4.*chi) + a5.*sin(6.*chi);
alat = -alat./cdr;
alon = -atan2(-(x-xcentre),(y-ycentre))./cdr + slon;
