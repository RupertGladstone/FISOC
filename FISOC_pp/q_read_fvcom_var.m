
%ncfile = '/global/work/qin/FVCOM_results/ISFWD.Ex5/output/Ex5_0001.nc';
ncfile = '/media/sf_VBshare/FISOC_Ex5_EF/Ex5_fvcom_0003.nc';


%var  = 'zisf';
var  = 'u';        [xg ,yg,uu]   = griddata_fvcom(ncfile,var,res);
var  = 'v';        [xg ,yg,vv]   = griddata_fvcom(ncfile,var,res);
var  = 'meltrate'; [xgr,ygr,melt]= griddata_fvcom(ncfile,var,res);


%%
% load M
% fzi = ncread(ncfile,'zisf');
% fu  = ncread(ncfile,'u');
% nt  = 10;
% k   =1 ;
% 
% figure
% 
% subplot(1,2,1)
% 
% % plot_field(Mobj,fzi(:,nt));
% scatter(Mobj.xc,Mobj.yc,20,squeeze(fu(:,k,nt)),'filled');
% colorbar
% set(gca,'linewidth',2,'fontsize',20);
% xlabel('Distance (m)')
% ylabel('Distance (m)')
% axis([min(Mobj.x) max(Mobj.x) min(Mobj.y) max(Mobj.y)])
%  
% subplot(1,2,2)
% %pcolor(xg,yg,squeeze(data(:,:,nt))); shading interp
% pcolor(xg,yg,squeeze(data(:,:,k,nt))); shading interp
% 
% 
% colorbar
% set(gca,'linewidth',2,'fontsize',20);
% xlabel('Distance (m)')
% ylabel('Distance (m)')
% 
% colorbar