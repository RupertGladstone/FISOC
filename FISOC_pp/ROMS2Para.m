
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to prepare ROMS output netcdf file for reading into
% Paraview. Involves creation of two (x,y) 1-dimensional coordinate
% variables. Some assumptions about the ROMS output file are
% hard-coded here.  Also, converts ROMS hybrid vertical coord to
% cartesian (assumes that some ROMS matlab functionality is
% available, specifically set_depth and stretching).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [returnCode] = ROMS2Para(inFileName)

returnCode = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Using high level netcdf interface for variable reading'])
x_rho = ncread(inFileName,'x_rho');
y_rho = ncread(inFileName,'y_rho');
time  = squeeze(ncread(inFileName,'ocean_time'));
zeta  = ncread(inFileName,'zeta');
draft = ncread(inFileName,'draft');
w_vel = ncread(inFileName,'w');
% assume coord data is replicated across the grid, i.e. we're using
% a regular cartesian grid here, so that we can simply use 1D coord
% vars. 
xi_rho_data  = squeeze(x_rho(:,1));
eta_rho_data = squeeze(y_rho(1,:));
% read meta data to use in the coordinate mapping from ROMS hybrid
% vertical coord to z coord.
Vtransform  = ncread(inFileName,'Vtransform');
Vstretching = ncread(inFileName,'Vstretching');
theta_s     = ncread(inFileName,'theta_s');
theta_b     = ncread(inFileName,'theta_b');
hc          = ncread(inFileName,'hc');
h           = double(ncread(inFileName,'h'));
igrid_rho = 1; igrid_w   = 5; % see set_depth for grid
                              % stagger info
report = 0; % verbose or not
nt = length(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['using low level netcdf interface for writing to file'])
ncid = netcdf.open(inFileName,"WRITE");

% first get dimension ids
xi_rho_dimid  = netcdf.inqDimID(ncid,'xi_rho');
eta_rho_dimid = netcdf.inqDimID(ncid,'eta_rho');
s_rho_dimid   = netcdf.inqDimID(ncid,'s_rho');
s_w_dimid     = netcdf.inqDimID(ncid,'s_w');
time_dimid    = netcdf.inqDimID(ncid,'ocean_time');
[~,N]  = netcdf.inqDim(ncid,s_rho_dimid); % number of levels in the
                                          % vertical (rho points)
[~,Nw] = netcdf.inqDim(ncid,s_w_dimid);   % number of levels in the
                                          % vertical (w points)
[~,nx] = netcdf.inqDim(ncid,xi_rho_dimid);
[~,ny] = netcdf.inqDim(ncid,eta_rho_dimid);
 
% enter define mode to define the new 1D coordinate variables
netcdf.reDef(ncid);
xi_rho_varid  = netcdf.defVar(ncid,'xi_rho','double',xi_rho_dimid);
eta_rho_varid = netcdf.defVar(ncid,'eta_rho','double',eta_rho_dimid);
z_rho_varid   = netcdf.defVar(ncid,'z_rho','double',[xi_rho_dimid,eta_rho_dimid,s_rho_dimid,time_dimid]);
w_vel_varid   = netcdf.defVar(ncid,'w_rho','double',[xi_rho_dimid,eta_rho_dimid,s_rho_dimid,time_dimid]);

z_rho_data = zeros(nx,ny,N,nt);  z_w_data = zeros(nx,ny,Nw,nt);
w_rho_data = zeros(nx,ny,N,nt); % we want w vel on rho points

disp(['process z coord '])
for tt = 1:nt
    disp([tt])
    dpz = zeta(:,:,tt) + draft(:,:,tt);
    z_rho_data(:,:,:,tt)=set_depth(Vtransform, Vstretching,theta_s, ...
                                                      theta_b, hc, N,igrid_rho, ...
                                                      h, dpz, report);
    z_w_data(:,:,:,tt)=set_depth(Vtransform, Vstretching,theta_s, ...
                                                    theta_b, hc, N,igrid_w, ...
                                                    h, dpz, report);
end
%disp(['adjust for ice draft'])
%for ll = 1:N
%    z_rho_data(:,:,ll,1:nt) = z_rho_data(:,:,ll,1:nt) + reshape(draft(:,:,1:nt),nx,ny,1,nt);
%end
%for ll = 1:Nw
%    z_w_data(:,:,ll,1:nt)   = z_w_data(:,:,ll,1:nt)   + reshape(draft(:,:,1:nt),nx,ny,1,nt);
%end

disp(['interpolate w velocity on to rho points'])
for tt = 1:nt
    disp(['processing time step ',num2str(tt)])
    for xx = 1:nx
        for yy = 1:ny
            w_rho_data(xx,yy,:,tt) = interp1(squeeze(z_w_data(xx,yy,:,tt)), ...
                                             squeeze(w_vel(xx,yy,:,tt)), ...
                                             squeeze(z_rho_data(xx,yy,:,tt)) );
        end
    end
end

disp(['leave define mode to populate coord variables'])
netcdf.endDef(ncid);
netcdf.putVar(ncid,xi_rho_varid,xi_rho_data);
netcdf.putVar(ncid,eta_rho_varid,eta_rho_data);
netcdf.putVar(ncid,z_rho_varid,z_rho_data);
netcdf.putVar(ncid,w_vel_varid,w_rho_data);

% Paraview can more easily match the time stepping of Elmer if we
% set the time coord to give just integer increments.
time_varid = netcdf.inqVarID(ncid,'ocean_time');
timeSteps = linspace(0,nt-1,nt);

netcdf.putVar(ncid,time_varid,timeSteps);

disp(['all done, closing netcdf file'])
netcdf.close(ncid)

returnCode = 0;

return