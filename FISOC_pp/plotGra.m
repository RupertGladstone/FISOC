
EF5_rowArea = 53.6; ER5_rowArea = 30.0; off = -0.0


% Experiment 5 (VE2 in the paper) 

% bilinear interp, Elmer/ROMS Ex5
iab = load('/media/sf_VBshare/FISOC_Ex5_bil2/graOverTime_ER5.asc');
oab = load('/media/sf_VBshare/FISOC_Ex5_bil2/ocean_his_select_gra.asc');

% Elmer/FVCOM Ex5
iaf = load('/media/sf_VBshare/FISOC_Ex5_EF/graOverTime_EF5.asc');
%oaf_all = load('/media/sf_VBshare/FISOC_Ex5_EF/totaldryarea_Ex5.txt');
oaf_all = load('/media/sf_VBshare/FISOC_Ex5_EF/dryarea_fvcom.txt');
time_oaf = oaf_all(:,1); oaf = oaf_all(:,2);

% convert from metres squared to km squared
msq2kmsq = 1000000.0;
iab = iab/msq2kmsq; oab = oab/msq2kmsq;
iaf = iaf/msq2kmsq; oaf = oaf/msq2kmsq; 

nt = 1440;

time = linspace(0,40,nt);
l1_ER5 = linspace(0,0,nt);
l2_ER5 = linspace(ER5_rowArea,ER5_rowArea,nt);
l1_EF5 = linspace(off,off,nt);
l2_EF5 = linspace(off+1*EF5_rowArea,off+1*EF5_rowArea,nt);
l3_EF5 = linspace(off+2*EF5_rowArea,off+2*EF5_rowArea,nt);
l4_EF5 = linspace(off+3*EF5_rowArea,off+3*EF5_rowArea,nt);
l5_EF5 = linspace(off+4*EF5_rowArea,off+4*EF5_rowArea,nt);

oaf_ti = interp1(time_oaf,oaf,time);
%    Vq = INTERP1(X,V,Xq) interpolates to find Vq, the values of the
%    underlying function V=F(X) at the query points Xq. 
 
%    X must be a vector. The length of X is equal to N.
%    If V is a vector, V must have length N, and Vq is the same size as Xq.

%dab = oab(2:nt+1)-iab(1:nt);
dab = oab(1:nt)-iab(1:nt);
daf = oaf(1:nt)-iaf(1:nt); 


figure(2) ; clf ;

subplot(2,1,1)
%plot(time_oaf,oaf,'r-'); hold on
plot(time,oaf_ti,'k-'); hold on
plot(time,iaf(1:nt),'k--'); 
xlim([0 40])
ylim([0.0e8 12.0e2])
xlabel(['Time, a'])
ylabel(['Area, km^2'])
%plot(time,oab(1:nt),'k-'); hold on
%plot(time,iab(1:nt),'k--'); 
legend('Ocean dry cell area', 'Ice sheet grounded area')

set(gca,'FontSize',7,'FontName','Helvetica');
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperUnits','centimeters')

subplot(2,1,2)
plot(time,daf(1:nt),'k-');hold on
%plot(time,l1_EF5,'g-')
plot(time,l2_EF5,'g-')
plot(time,l3_EF5,'g-')
plot(time,l4_EF5,'g-')
plot(time,l5_EF5,'g-')
xlim([0 40])
ylim([0.2e2 1.3e2])
xlabel(['Time, a'])
ylabel(['Area, km^2'])
legend( 'Ocean minus ice grounded area')

set(gca,'FontSize',7,'FontName','Helvetica');
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0.25 1.25 8 8])

saveas(gcf,'GroundedArea_EF5.png')


figure(1) ; clf ;

subplot(2,1,1)
plot(time,oab(1:nt),'k-'); hold on
plot(time,iab(1:nt),'k--'); 
xlim([0 40])
ylim([0.0e2 12.0e2])
xlabel(['Time, a'])
ylabel(['Area, km^2'])
legend('Ocean dry cell area', 'Ice sheet grounded area')

set(gca,'FontSize',7,'FontName','Helvetica');
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperUnits','centimeters')

subplot(2,1,2)
plot(time,dab(1:nt),'k-');hold on
plot(time,l1_ER5,'g-')
plot(time,l2_ER5,'g-')
xlim([0 40])
ylim([-0.1e2 0.5e2])
xlabel(['Time, a'])
ylabel(['Area, km^2'])
legend( 'Ocean minus ice grounded area')

set(gca,'FontSize',7,'FontName','Helvetica');
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0.25 1.25 8 8])

saveas(gcf,'GroundedArea_ER5.png')

