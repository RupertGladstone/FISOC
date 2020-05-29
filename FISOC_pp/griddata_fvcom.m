function [xg,yg,data]=griddata_fvcom(ncfile,var,resg)
% read out variables from fvcom output and interpolate them
% onto 1km structured grid
%
% input:
%   ncfile - fvcom output netcdf results file
%   resg    - resolution of structured grid in meters
%   var     - fvcom output variable
%     - 'u','v'   - 3D velocity field in m/s
%     - 'meltrate' - melting rate in m/s
%     - 'zisf '    - ice draft in meters
%     -  'h'       - bathymetry
%     - 'temp'     - temperature in Celcus degrees
%     - 'salinity'     - salinity in PSU


% output
%   - 'xg','yg' - structured grid in 1 km resolution
%   - data      - fvcom variables interpolated onto xg,yg
%----------------------------------%


x = ncread(ncfile,'x');
y = ncread(ncfile,'y');
xc = ncread(ncfile,'xc');
yc = ncread(ncfile,'yc');

if strcmp(var,'u') | strcmp(var,'v')
    xg = [min(xc):resg:max(xc)];
    yg = [min(yc):resg:max(yc)];
    [xg yg] = meshgrid(xg,yg);
    tmp = ncread(ncfile,var);
    nt  = size(tmp,3);  % time dimension
    for j = 1:nt
      for k = 1:size(tmp,2) % vertical layers
        data(:,:,k,j) = griddata(xc,yc,tmp(:,k,j),xg,yg);
      end
    end

elseif strcmp (var,'zisf') | strcmp(var,'meltrate')| strcmp(var,'h')
    xg = [min(x):resg:max(x)];
    yg = [min(y):resg:max(y)];
    [xg yg] = meshgrid(xg,yg);
    tmp = ncread(ncfile,var);
    nt  = size(tmp,2);  % time dimension
    for j = 1:nt
        data(:,:,j) = griddata(x,y,tmp(:,j),xg,yg);
    end
    
else  strcmp(var,'temp') | strcmp(var,'salinity')
    xg = [min(x):resg:max(x)];
    yg = [min(y):resg:max(y)];
    [xg yg] = meshgrid(xg,yg);
    tmp = ncread(ncfile,var);
    nt  = size(tmp,3);  % time dimension
    for j = 1:nt
      for k = 1:size(tmp,2) % vertical layers
        data(:,:,k,j) = griddata(xc,yc,tmp(:,k,j),xg,yg);
      end
    end
    
end
end