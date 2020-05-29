

ncfile = '/media/sf_VBshare/FISOC_Ex5_bil2/ocean_his_select.nc';

xr = ncread(ncfile,'x_rho');
yr = ncread(ncfile,'y_rho');

ur = ncread(ncfile,'u_eastward');
vr = ncread(ncfile,'v_northward');
urb = ncread(ncfile,'ubar_eastward');
vrb = ncread(ncfile,'vbar_northward');

meltr = ncread(ncfile,'m');

size(xr);

