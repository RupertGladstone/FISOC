
clear all
!rm isomip_plus_ocean1.nc
InpGrd = 'Ocean1_input_geom_v1.01.nc'
GrdName = 'isomip_plus_ocean1.nc'


addpath(genpath('/home/elmeruser/Source/roms_matlab'))
%addpath(genpath('/home/elmeruser/Source/roms_wilkin'))
%addpath(genpath('roms_matlab'))

icedraft = ncread(InpGrd,'lowerSurface');
h = -ncread(InpGrd,'bedrockTopography');
xload = ncread(InpGrd,'x');
yload = ncread(InpGrd,'y');
[Yload,Xload]=meshgrid(yload,xload);

%% define ISOMIP_PLUS arbitrary cartesian grid
dHorz = 2e3; %m
x_min =320e3; 
x_max =800e3;
y_min =0;
y_max =80e3;
x = [x_min+dHorz/2:dHorz:x_max-dHorz/2]; %grid cell centers
y = [y_min+dHorz/2:dHorz:y_max-dHorz/2]; %grid cell centers
% alternatively, load x and y from input file.
%x = ncread(InpGrd,'x');
%y = ncread(InpGrd,'y');
[s.y s.x ] = meshgrid(y,x); % yep... they reversed x and y. :|
y = fliplr(y); % and they made positive y go to the left. Go figure.


% smoothing ?
wct = h + icedraft;
wcts = smooth_bath(h+icedraft,(icedraft<0).*double(wct>0).*~isnan(wct),4,0.1); %smooth all wct where ocean exists
icedraft_orig = icedraft; icedraft = -smooth_bath(-icedraft,double(wct>0),8,0.1); %smooth ice front 
h_orig=h; h = h_orig+(wcts-wct); %shift bathy up/down so that it's smoother


if 1
% Because this script is shit, insert this hack to increase the size of h,zice,x,y
% since we chop: (1:end-1,2:end), add a NaN row at end, and a nan col at start.
h(end+1,:)=NaN; h(:,2:end+1)=h; h(:,1)=NaN;
icedraft(end+1,:)=NaN; icedraft(:,2:end+1)=icedraft; icedraft(:,1)=NaN;
%now add shitty x,y increments, in the same way as above, to s.x and s.y
Xload(end+1,:)=NaN; Xload(:,2:end+1)=Xload; Xload(:,1)=NaN; Xload = inpaint_nans(Xload,1);
Yload(end+1,:)=NaN; Yload(:,2:end+1)=Yload; Yload(:,1)=NaN; Yload = inpaint_nans(Yload,1);
s.x(end+1,:)=NaN; s.x(:,2:end+1)=s.x; s.x(:,1)=NaN;
s.y(end+1,:)=NaN; s.y(:,2:end+1)=s.y; s.y(:,1)=NaN;
increment_hack=1;
end


% interpolate original grid to ocean grid
H = griddata(Xload,Yload,h,s.x,s.y);
Z = griddata(Xload,Yload,icedraft,s.x,s.y);
wct = H+Z;

% calculate corresponding geographic grid...
referenceLat = -75;
[refLatInPSG_x0,refLatInPSG_y0]=polar_stereo_deluxe(referenceLat,0,0,0,1e-3,-71); %in km

sex = reshape( s.x(:,2).' ,1,numel(s.x(:,2)));
sex(find(isnan(sex)))= [];
PSG_y0 = s.x(:,2)' + refLatInPSG_y0 - median(sex); % ref to median of grid.

%PSG_y0 = s.x(:,2)' + refLatInPSG_y0 - nanmedian(s.x(:,2)'); % ref to median of grid.
[lat0,lon0] = inverse_polar_stereo(zeros(size(PSG_y0))/1000,PSG_y0/1000,0,0,1,-71);
geogrid_lat = repmat(lat0,[size(h,2) 1]);
%dlon = dHorz * 360 ./(2*pi*earthRadius('m')*cos(deg2rad(lat0)));
dlon = dHorz * 360 ./(2*pi*6371000*cos(deg2rad(lat0)));
geogrid_lon(1,:)   = lon0;
for iii=2:size(h,2)
geogrid_lon(iii,:) = geogrid_lon(iii-1,:)-dlon;
end
[PSG_x,PSG_y]=polar_stereo_deluxe(geogrid_lat,geogrid_lon,0,0,1e-3,-71); %in km



% calculate masks
min_wct = 20; %set minimum water column thickness
max_wct = 7000; %set max WCT?
mask = wct;
mask(wct > min_wct) = 1;
mask(wct <= min_wct) = 0;

% manually remove those ugly lobes.
if increment_hack
mask(98:118,2:6)=0;
mask(98:118,37:41)=0;
%mask(101,8)=0;
else
mask(98:118,1:5)=0;
mask(98:118,36:40)=0;
%mask(101,7)=0;
end

imask = Z;

% add calving
imask(abs(imask)<10)=0;

%finalise mask
imask(abs(imask)>0) = 1;
imask = imask.*mask;



% use s.m* for internal stuff, and mask/imask for saving to netcdf
s.mw = mask;   
s.mi = imask; 

rmask = logical(s.mw(1:end-1,2:end));
imask = logical(s.mi(1:end-1,2:end));  % change this to be consistent..


s.dx = ones(size(s.x)).*dHorz;
s.dy = ones(size(s.y)).*dHorz;


[m,n] = size(s.x);

%geogrid_lon = s.lon;
%geogrid_lat = s.lat;
geometry{1} = s.dx;
geometry{2} = s.dy;


if ~isequal(size(rmask), size(s.x)-1)
    if ~isempty(rmask)
        disp(' ## Wrong size mask.')
    end
    rmask = zeros(m-1, n-1);
end

if ~isequal(size(imask), size(s.x)-1)
    if ~isempty(imask)
        disp(' ## Wrong size mask.')
    end
    imask = zeros(m-1, n-1);
end


%% Checks of masking and depths:
bathymetry = H; %rename to these so old names are archived.
ice_draft = Z;
%ang = s.angle;   % ROMS needs radians. %need to work out angle for grid.
ang = zeros(size(s.x));
min_depth = min_wct;
max_depth = max_wct;
ice_draft = ice_draft.*s.mi;  %ensure ice is masked!
ice_draft(ice_draft > 0) = 0; %ensure ice draft is not positive
bathymetry(find(isnan(bathymetry))) = min_depth; %change any nans to min_depth
bathymetry(bathymetry<min_depth) = min_depth;    %enforce minimum depth to min_depth
bathymetry((bathymetry + ice_draft) > max_depth) = max_depth; %enforce max depth
s.mwNaN = s.mw; % Ensure WCT=0 over land
s.mwNaN(s.mw == 0) = NaN;
wct = (bathymetry + ice_draft).*s.mwNaN;
wct(isnan(wct)) = 0; 
bathymetry(wct<min_depth) = -ice_draft(wct<min_depth)+min_depth; %any bathy where wct is too thin is set to ice draft+min_depth
%ii = find(wct < min_depth);		%find any water thinner than min_wct
%bathymetry(ii) = -ice_draft(ii) + min_depth; %set bathy where wct<min_wct = ice + min_wct
%wct = bathymetry + ice_draft; 		%recalc wct.
%bathymetry(bathymetry<min_depth) = min_depth; %check again to enforce minimum depth of bathy...

if 1 %manual nearest interp
clear xtmp8 ytmp8 xtmp ytmp
for indi=1:size(s.x,1)
  xtmp8(2*indi-1,:) = s.x(indi,:);
  xtmp8(2*indi,:) = s.x(indi,:);
  ytmp8(2*indi-1,:) = s.y(indi,:);
  ytmp8(2*indi,:) = s.y(indi,:);
end
for indi=1:size(s.x,2)
    xtmp(:,2*indi-1) = xtmp8(:,indi);
    xtmp(:,2*indi) = xtmp8(:,indi);
    ytmp(:,2*indi-1) = ytmp8(:,indi);
    ytmp(:,2*indi) = ytmp8(:,indi);
end
%xtmp(:,1)=[];
%ytmp(:,1)=[];
%xtmp(end,:)=[];  %Bug fix: 'end' instead of '1'
%ytmp(end,:)=[];  %Bug fix: 'end' instead of '1'

xtmp = xtmp(1:end-1,2:end); %or just do this...
ytmp = ytmp(1:end-1,2:end);

grid_x=xtmp;
grid_y=ytmp;

s.lat = geogrid_lat';
s.lon = geogrid_lon';
clear lattmp8 lattmp lontmp lontmp8
for indi=1:size(s.x,1)
lattmp8(2*indi-1,:) = s.lat(indi,:);
lattmp8(2*indi,:) = s.lat(indi,:);
lontmp8(2*indi-1,:) = s.lon(indi,:);
lontmp8(2*indi,:) = s.lon(indi,:);
end
for indi=1:size(s.x,2)
lattmp(:,2*indi-1) = lattmp8(:,indi);
lattmp(:,2*indi) = lattmp8(:,indi);
lontmp(:,2*indi-1) = lontmp8(:,indi);
lontmp(:,2*indi) = lontmp8(:,indi);
end
%lattmp(:,1)=[];
%lattmp(end,:)=[];  %Bug fix: 'end' instead of '1'
%lontmp(:,1)=[];
%lontmp(end,:)=[];  %Bug fix: 'end' instead of '1'
lontmp = lontmp(1:end-1,2:end); %or just do this...
lattmp = lattmp(1:end-1,2:end); %or just do this...

geogrid_lat = lattmp;
geogrid_lon = lontmp;
elseif 0
theInterpFcn = 'interp2';
theInterpMethod = 'nearest';
grid_x = feval(theInterpFcn, s.x, 1, theInterpMethod);
grid_y = feval(theInterpFcn, s.y, 1, theInterpMethod);
geogrid_lon = feval(theInterpFcn, geogrid_lon, 1, theInterpMethod);
geogrid_lat = feval(theInterpFcn, geogrid_lat, 1, theInterpMethod);
end
[n, m] = size(grid_x);


FLIPPING = 0;
if FLIPPING
    grid_x = flipud(grid_x);
    grid_y = flipud(grid_y);
    geogrid_lon = flipud(geogrid_lon);
    geogrid_lat = flipud(geogrid_lat);
    geometry{1} = flipud(geometry{1});
    geometry{2} = flipud(geometry{2});
    mask = flipud(mask);
    imask = flipud(imask);
    bathymetry = flipud(bathymetry);
    ice_draft = flipud(ice_draft);
    ang = flipud(ang);
end

TRANSPOSE = 0;
if TRANSPOSE
    grid_x = grid_x';
    grid_y = grid_y';
    geogrid_lon = geogrid_lon';
    geogrid_lat = geogrid_lat';
    geometry{1} = geometry{1}';
    geometry{2} = geometry{2}';
    mask = mask';
    imask = imask';
    bathymetry = bathymetry';
    ice_draft = ice_draft';
    ang = ang';
end


xl = max(grid_x(:)) - min(grid_x(:));
el = max(grid_y(:)) - min(grid_y(:));
if ~rem(m, 2), m = m-1; end   % m, n must be odd.
if ~rem(n, 2), n = n-1; end

i_rho = 2:2:m-1; j_rho = 2:2:n-1;
i_psi = 3:2:m-2; j_psi = 3:2:n-2;
i_u   = 3:2:m-2; j_u   = 2:2:n-1;
i_v   = 2:2:m-1; j_v   = 3:2:n-2;

% The xi direction (left-right):
LP = (m-1)/2;   % The rho dimension.
L = LP-1;       % The psi dimension.

% The eta direction (up-down):
MP = (n-1)/2;   % The rho dimension.
M = MP-1;       % The psi dimension.
f = 2.*7.2921e-5.*sin(geogrid_lat(j_rho, i_rho).*pi./180);
f(isnan(f))=0;  % doesn't like f to be NaN

% Calculate other masking arrays.

umask = zeros(size(rmask));
vmask = zeros(size(rmask));
pmask = zeros(size(rmask));

for i = 2:LP
    for j = 1:MP
        umask(j, i-1) = rmask(j, i) * rmask(j, i-1);
    end
end

for i = 1:LP
    for j = 2:MP
        vmask(j-1, i) = rmask(j, i) * rmask(j-1, i);
    end
end

for i = 2:LP
    for j = 2:MP
        pmask(j-1, i-1) = rmask(j, i) * rmask(j, i-1) * rmask(j-1, i) * rmask(j-1, i-1);
    end
end



% Average angle -- We should do this via (x, y) components.
temp = ang;
ang = zeros(n, m);
ang(1:2:end-1,2:2:end) = temp(1:end-1,2:end);
gx = geometry{1};   % Spherical distances in meters.
gy = geometry{2};

sx = 0.5*(gx(1:end-1, :) + gx(2:end, :));
sy = 0.5*(gy(:, 1:end-1) + gy(:, 2:end));

pm = 1 ./ sx(:,2:end);
pn = 1 ./ sy(2:end,:);

% pm and pn cannot be Inf, even if on land, so if values
% are Inf, set to an arbitrary non-zero value
pm(isinf(pm))=0.999e-3;
pn(isinf(pn))=0.999e-3;
pm(isnan(pm))=0.999e-3;
pn(isnan(pn))=0.999e-3;
dmde = zeros(size(pm));
dndx = zeros(size(pn));

dmde(2:end-1, :) = 0.5*(1./pm(3:end, :) - 1./pm(1:end-2, :));
dndx(:, 2:end-1) = 0.5*(1./pn(:, 3:end) - 1./pn(:, 1:end-2));
dmde(isinf(dmde))=0;
dndx(isinf(dndx))=0;
dmde(isnan(dmde))=0;
dndx(isnan(dndx))=0;



%% Final write
addpath(genpath('roms_matlab'))
c_grid(LP,MP,GrdName);

nc_write(GrdName,'xl', xl);
nc_write(GrdName,'el', el);
nc_write(GrdName,'f',f');
nc_write(GrdName,'x_rho', grid_x(j_rho, i_rho)');
nc_write(GrdName,'y_rho', grid_y(j_rho, i_rho)');
nc_write(GrdName,'x_psi', grid_x(j_psi, i_psi)');
nc_write(GrdName,'y_psi', grid_y(j_psi, i_psi)');
nc_write(GrdName,'x_u', grid_x(j_u, i_u)');
nc_write(GrdName,'y_u', grid_y(j_u, i_u)');
nc_write(GrdName,'x_v', grid_x(j_v, i_v)');
nc_write(GrdName,'y_v', grid_y(j_v, i_v)');
nc_write(GrdName,'lon_rho', geogrid_lon(j_rho, i_rho)');
nc_write(GrdName,'lat_rho', geogrid_lat(j_rho, i_rho)');
nc_write(GrdName,'lon_psi', geogrid_lon(j_psi, i_psi)');
nc_write(GrdName,'lat_psi', geogrid_lat(j_psi, i_psi)');
nc_write(GrdName,'lon_u', geogrid_lon(j_u, i_u)');
nc_write(GrdName,'lat_u', geogrid_lat(j_u, i_u)');
nc_write(GrdName,'lon_v', geogrid_lon(j_v, i_v)');
nc_write(GrdName,'lat_v', geogrid_lat(j_v, i_v)');
nc_write(GrdName,'h', bathymetry(1:end-1,2:end)');
nc_write(GrdName,'zice', ice_draft(1:end-1,2:end)');
nc_write(GrdName,'mask_rho', double(rmask'));
nc_write(GrdName,'mask_zice', double(imask'));
nc_write(GrdName,'mask_psi', double(pmask(1:end-1, 1:end-1)'));
nc_write(GrdName,'mask_u', double(umask(1:end, 1:end-1)'));
nc_write(GrdName,'mask_v', double(vmask(1:end-1, 1:end)'));
nc_write(GrdName,'angle', ang(j_rho, i_rho)');
nc_write(GrdName,'pm', pm');
nc_write(GrdName,'pn', pn');
nc_write(GrdName,'dmde', dmde');
nc_write(GrdName,'dndx', dndx');

