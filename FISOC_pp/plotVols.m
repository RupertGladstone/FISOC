
start = 3;

% Experiment 4 (VE1 in the paper) 

% bilinear interp
ovb = load('/media/sf_VBshare/FISOC_Ex4_interp/ocean_his_select_vol.asc'); % ocean
ivb = load('/media/sf_VBshare/FISOC_Ex4_interp/metrics_Ex4.dat'); % ice

% nearest neighbour regridding
ivn = load('/media/sf_VBshare/FISOC_Ex4_cp/metrics_Ex4.dat'); % ice
ovn = load('/media/sf_VBshare/FISOC_Ex4_cp/ocean_his_select_vol.asc'); % ocean

ivo = load('/media/sf_VBshare/FISOC_Ex4_EO/metrics_Ex4.dat'); % ice only

% Elmer/FVCOM coupling (others are Elmer/ROMS)
ivf = load('/media/sf_VBshare/FISOC_Ex4_EF/metrics_Ex4.dat');
ovf_all = load('/media/sf_VBshare/FISOC_Ex4_EF/totalmass_Ex4.txt');
time_ovf = ovf_all(:,1); omf = ovf_all(:,2);

%scp gladston@puhti.csc.fi:/scratch/project_2000339/gladston/EO_FX4/*dat ./elmerVol_Ex4_IO.asc
%scp gladston@puhti.csc.fi:/scratch/project_2000339/gladston/ER_FX4_interp/*dat ./elmerVolEx4Bil.asc
%ivo = load('elmerVolEx4_IO.asc');    % ice only
%ivb = load('elmerVolEx4Bil.asc');    % bilinear interp, ice
%ovb = load('ocean_his_bil_vol.asc'); % bilinear interp, ocean

rho_o = 1027.0;
rho_i = 910.0;

% convert volume to mass
omb = ovb*rho_o; imb = ivb*rho_i; % coupled bilinear interp ER5
omn = ovn*rho_o; imn = ivn*rho_i; % coupled nearest neighbour
imo = ivo*rho_i;                  % ice only
imf = ivf*rho_i; 

time = linspace(0,100,3600);
time_if = linspace(0,48.986,1763);
time_io = linspace(0,100,1201);

omf_ti = interp1(time_ovf,omf,time_if');

tot = ovb+ivb;
totmb = omb+imb; % total mass, bilinear interp
totmn = omn+imn; % total mess, nearest neighbour
totmf = omf_ti+imf; % total mass, FVCOM coupling


%ovb(end)-ovb(1);
%ivb(end)-ivb(1);
%tot(end)-tot(1);
%omb(end)-omb(1)
%imb(end)-imb(1);
%totm(end)-totm(1);

% nodemalised volume trends
%ovbn  = ovb -   ovb(3);
%ivbn  = ivb -   ivb(3);
%ivon  = ivo -   ivo(3);
%totn  = tot -   tot(3);

% normalised mass trends
omb_n  = omb - omb(start); imb_n  = imb - imb(start); totmb_n = totmb - totmb(start);
omn_n  = omn - omn(start); imn_n  = imn - imn(start); totmn_n = totmn - totmn(start);
omf_n  = omf - omf(start); imf_n  = imf - imf(start); totmf_n = totmf - totmf(start);
imo_n  = imo - imo(start);

%figure(1) ; clf ;
%subplot(3,1,1); hold on ;
%plot(time, tot)
%xlim([0 100]); ylim([1.54e12 1.55e12]);
%subplot(3,1,2); hold on ;
%plot(time, ivb)
%xlim([0 100]); ylim([8.44e11 8.54e11]);
%subplot(3,1,3); hold on ;
%plot(time, ovb)
%xlim([0 100]); ylim([6.9e11 7.0e11]);

%figure(2) ; clf ;
%plot(time,ovbn); hold on
%plot(time,ivbn)
%%plot(time_io,ivon)
%plot(time,totn)

figure(3) ; clf ;
plot(time,omb_n,'k--'); hold on
plot(time,imb_n,'k-.')
plot(time,totmb_n,'k-')
plot(time_ovf,omf_n,'r--')
plot(time_if,imf_n,'r-.')
plot(time_if,totmf_n,'r-')
xlim([0 40])
ylim([-1.3e12 2.3e12])

xlabel(['Time, a'])
ylabel(['Mass, kg'])

legend('VE1\_ER Ocean',  ['VE1\_ER Ice shelf'], 'VE1\_ER Total',  ...  
       'VE1\_EF Ocean', ['VE1\_EF Ice shelf'], ['VE1\_EF Total'], ...
       'Location','northeast')   
 
set(gca,'FontSize',7,'FontName','Helvetica');
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0.25 1.25 8 8])
saveas(gcf,'MassCons.png')

