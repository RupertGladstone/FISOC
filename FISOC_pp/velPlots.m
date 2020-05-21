
%clear all
%close all

nn = 2; % downsampling stride

% box size for plots
x1 = -2.0; x2 = 32.0;
y1 = 10.0; y2 = 102.0;

squish = 0.1; barWidth = 0.03; % manipulating relative position
                               % colorbar and subplots

res  = 1000; % regridding resolution for fvcom

cutoff = 0.17; % remove fvcom melt rates above this value (grounded
              % region seems to have high values for some reason...)

fs = 9 ; % font size


read_from_roms;
xr = xr'; yr = yr';
meltr = meltr(:,:,900); meltr = meltr'
meltr_yr = meltr * 360.0 * 24.0 * 3600.0;
xr_km = xr/1000.0; yr_km = yr/1000.0;

% subsample ROMS for quiver plot
ur2d = ur(:,:,11,900); vr2d = vr(:,:,11,900);
%ur2d = urb(:,:,900); vr2d = vrb(:,:,900);
ur2d_s = ur2d(1:nn:end,1:nn:end); vr2d_s = vr2d(1:nn:end,1:nn:end);
xr_km_s = xr_km(1:nn:end,1:nn:end); yr_km_s = yr_km(1:nn:end,1:nn:end);


%q_read_fvcom_var;
melt_yr = melt * 360.0 * 24.0 * 3600.0;
xg_km  = xg/1000.0  ; yg_km  = yg/1000.0 ;
xgr_km = xgr/1000.0 ; ygr_km = ygr/1000.0 ;
melt_yr(melt_yr>cutoff) = 0.0;

% subsample fvcom for quiver plot
uu2d = uu(:,:,1,106); vv2d = vv(:,:,1,106);
uu2d_s = uu2d(1:nn:end,1:nn:end); vv2d_s = vv2d(1:nn:end,1:nn:end);
xg_km_s = xg_km(1:nn:end,1:nn:end); yg_km_s = yg_km(1:nn:end,1:nn:end);



figure(1); clf;

subplot(1,2,1); hold on
ap1 = get(gca,'position');
ap1(3) = ap1(3) - squish/2.0;

subplot(1,2,2); hold on
ap2 = get(gca,'position');
ap2(1) = ap2(1) - squish/2.0;
ap2(3) = ap2(3) - squish/2.0;
set(gca,'position',ap2);

%set(gca,'position',[ 0.5703-LeftShift 0.1100 0.3347+LeftShift 0.8150])
set(gca,'FontSize',fs,'FontName','Helvetica');

imagesc(xgr_km(1,:),ygr_km(12:end,1),melt_yr(12:end,:,106),[-0.2 0.2])  
colormap(bluewhitered(256))
colorbar('Position', [ap2(1)+ap2(3)+squish/2.0 ap2(2) barWidth ap2(4)])
q = quiver(xg_km_s,yg_km_s,uu2d_s,vv2d_s,1.5,'k','linewidth',1);
xlim([x1 x2])
ylim([y1 y2])
xlabel(['Distance, km'])
ylabel(['Distance, km'])


%set(gcf,'PaperPositionMode','manual')
%set(gcf,'PaperUnits','centimeters')
set(gca,'box','on')


subplot(1,2,1); hold on
set(gca,'position',ap1)
%get(gca,'position')
%set(gca,'position',[ 0.1300 0.1100 0.3347-LeftShift 0.8150])
set(gca,'FontSize',fs,'FontName','Helvetica');

imagesc(xr_km(1,:),yr_km(13:end,1),meltr_yr(13:end,:),[-0.2 0.2]);
%ph = pcolor(xr_km(1,:),yr_km(:,1),meltr_yr);
%set(ph,'MarkerEdgeColor','none')
%set(ph,'MeshStyle','none')      
%colormap(bluewhitered(256)), colorbar
q = quiver(xr_km_s,yr_km_s,ur2d_s',vr2d_s',1.5,'k','linewidth',1);
xlim([x1 x2])
ylim([y1 y2])
xlabel(['Distance, km'])
ylabel(['Distance, km'])

set(gca,'box','on')

%uistack(ph,'bottom')
%ax = gca();
%uistack(ax,'top')         

%subplot(1,3,3); hold on; colorbar



%set(gcf,'PaperPositionMode','manual')
%set(gcf,'PaperUnits','centimeters')
%set(gcf,'PaperPosition',[0.25 1.25 17 15])
saveas(gcf,'VelComparison.png')
