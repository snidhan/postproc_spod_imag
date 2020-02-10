
%% Written by Sheel Nidhan
%  Plotting the SPOD modes for a specific azimuthal wavenumber for the 3dplane SPOD
clear;
clc;
%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7000;
mode  = 1;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
nr = 354;
nplanes = 22;
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
x_sampled =  [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
Nrows = numvar*nr*nplanes*Nblk;
Nrows_permode = numvar*nr*nplanes;
dir2 = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/3dplane_spod/eigenmodes/';
disp(dir2);

%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

%% Reading the time file

time = importdata ('/home/sheel/Work/projects/spod_re5e4/post/frinf/time_stamps/time_stamp_1892600_2613200_uniform.txt');
time_spod = time(1:N,1);
dt   = time_spod(2:end,1) - time_spod(1:end-1,1);

%%  Partitioning the time data in different blocks

Nblk = floor((N-Novlp)/(Nfreq-Novlp));

qstart = zeros(Nblk,1);
qend  = zeros(Nblk,1);

for i = 1:Nblk
    qstart(i,1) = (i-1)*(Nfreq-Novlp) + 1;
    qend(i,1) = qstart(i) + Nfreq - 1;
end

%% Importing Karu's datafiles

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_TKE.dat';
LK_TKE   = importdata(filename);
LK_TKE_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_TKE(:,1)));
    LK_TKE_loc_planes(i,2) = LK_TKE(idx,4);
    LK_TKE_loc_planes(i,1) = LK_TKE(idx,1);

end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_WMEAN.dat';
LK_mean   = importdata(filename);
LK_mean_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_mean(:,1)));
    LK_mean_loc_planes(i,2) = LK_mean(idx,4);
    LK_mean_loc_planes(i,1) = LK_mean(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_centerline.dat';
TKE_Centerline   = importdata(filename);
TKE_centerline_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-TKE_Centerline(:,1)));
    TKE_centerline_loc_planes(i,2) = TKE_Centerline(idx,2);
    TKE_centerline_loc_planes(i,1) = TKE_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/UX_rms_centerline.dat';
ux_Centerline   = importdata(filename);
ux_centerline_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ux_Centerline(:,1)));
    ux_centerline_loc_planes(i,2) = ux_Centerline(idx,2);
    ux_centerline_loc_planes(i,1) = ux_Centerline(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Defect_centerline.dat';
ud_Centerline   = importdata(filename);
ud_centerline_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ud_Centerline(:,1)));
    ud_centerline_loc_planes(i,2) = ud_Centerline(idx,2);
    ud_centerline_loc_planes(i,1) = ud_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/MKE_areaI_FINF.dat';
mke_area   = importdata(filename);
mke_area_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-mke_area(:,1)));
    mke_area_loc_planes(i,2) = mke_area(idx,2);
    mke_area_loc_planes(i,1) = mke_area(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_areaI_FINF.dat';
tke_area   = importdata(filename);
tke_area_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-tke_area(:,1)));
    tke_area_loc_planes(i,2) = tke_area(idx,2);
    tke_area_loc_planes(i,1) = tke_area(idx,1);
end

%% Separating time data into blocks

time_blk = zeros(Nblk,Nfreq);

for i = 1:Nblk
    time_blk(i,:) = time_spod(qstart(i,1):qend(i,1),1)';
end

%% Fixing the frequency of SPOD spectrum

f = (0:Nfreq-1)/dt(1)/Nfreq;

if mod(Nfreq,2) == 0
    f(Nfreq/2 + 1:end) = f(Nfreq/2 + 1:end)-1/dt(1);
else
    f((Nfreq+1)/2 + 1:end) = f((Nfreq+1)/2 + 1:end) - 1/dt(1);
end

%% Reading SPOD eigenvalues files

eigmodes = zeros(Nrows_permode,1);
m = sprintf('%03d',mode);
freq =  sprintf('%04d',6);  
filename = strcat(dir2, 'eigenmode_freq_',freq,'_',m,'.mod');
A = importdata(filename);
eigmodes(:,1) = A(:,1) + sqrt(-1).*A(:,2);


%% Rearranging the eigmodes for visualization

eigenmodes_arranged(:,:) = reshape(eigmodes(:,1), [nr*numvar, nplanes]);

for i = 1:nplanes
    eigenmodes_separated(:,:,i) = reshape(eigenmodes_arranged(:,i), [nr, numvar]);
end

u_eigenmode = squeeze(eigenmodes_separated(:,1,:));
v_eigenmode = squeeze(eigenmodes_separated(:,2,:));
w_eigenmode = squeeze(eigenmodes_separated(:,3,:));


%% Normalizing eigenmodes


% u_eigenmode = u_eigenmode/(norm(u_eigenmode, 'fro'));
% v_eigenmode = v_eigenmode/(norm(v_eigenmode, 'fro'));
% w_eigenmode = w_eigenmode/(norm(w_eigenmode, 'fro'));

% for i = 1:nplanes
%     u_eigenmode(:,i) = u_eigenmode(:,i)/norm(u_eigenmode(:,i),2);
%     v_eigenmode(:,i) = v_eigenmode(:,i)/norm(v_eigenmode(:,i),2);
%     w_eigenmode(:,i) = w_eigenmode(:,i)/norm(w_eigenmode(:,i),2);
% end


%% Plotting the eigenmodes 

% figure;
% hold on;
% h1 = plot(rc./half_width_tke(1,1), abs(w_eigenmode(:,4)).^2, 'b-', 'Linewidth', 2);
% h2 = plot(rc./half_width_tke(2,1), abs(w_eigenmode(:,8)).^2, 'r-', 'Linewidth', 2);
% h3 = plot(rc./half_width_tke(3,1), abs(w_eigenmode(:,12)).^2, 'k-', 'Linewidth', 2);
% h4 = plot(rc./half_width_tke(4,1), abs(w_eigenmode(:,16)).^2, 'k--', 'Linewidth', 2);
% h5 = plot(rc./half_width_tke(5,1), abs(w_eigenmode(:,20)).^2, 'r--', 'Linewidth', 2);
% h6 = plot(rc./half_width_tke(6,1), abs(w_eigenmode(:,22)).^2, 'b--', 'Linewidth', 2);
% 
% hXLabel = xlabel('$r/L_{Kmean}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{x}|^{2}$','interpreter','latex','fontsize',15);
% % %hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);
% 
% hLegend = legend([h1, h2, h3, h4, h5, h6], 'x/D = 20', 'x/D = 40','x/D = 60', 'x/D = 80', 'x/D = 100', 'x/D = 120');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'phi_w_m2_st0_spodmode1_3dplane_tkebased.png','-dpng','-r600');  
% print(gcf,'phi_w_m2_st0_spodmode1_3dplane_tkebased.eps','-depsc','-r600');


% colormap('hot');
% colorbar;
% caxis([min(min(squeeze(real_ur_mode_contour(:,:,ct)))) max(max(squeeze(real_ur_mode_contour(:,:,ct))))])
% axis equal;
%% Three dimensional slices of real part of ur m=1 St=0.135

ur_mode = squeeze(u_eigenmode);
ux_mode = squeeze(w_eigenmode);

theta = linspace(0,2*pi,256);

% Fix the sign of mode based on the x/D = 10 mode
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ux_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         disp(i);
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end

figure; hold on;

for i = 2:2:20
    plot(rc, real(ur_mode(:,i)), '--', 'Linewidth', 2);
end

for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta)); %#ok<*SAGROW>
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:2:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ur_mode_contour(:,:,ct)),100);
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');
colorbar;
caxis([-1 1]);
%caxis([min(min(min(squeeze(real_ur_mode_contour(:,:,:))))) max(max(max(squeeze(real_ur_mode_contour(:,:,:)))))]);
axis equal;

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1]);


% dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ur_real_m1st0135_3d_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ur_real_m1st0135_3d_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of real part of ux m=1 St=0.135

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode);
ux_mode = squeeze(w_eigenmode);

% Fix the sign of mode based on the x/D = 10 mode

real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ur_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, real(ux_mode(:,i)), '--', 'Linewidth', 2);
end
close;

for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

end


set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ux_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');
colorbar;
caxis([-0.01 0.01]);
% caxis([min(min(min(squeeze(real_ur_mode_contour(:,:,:))))) max(max(max(squeeze(real_ur_mode_contour(:,:,:)))))]);
axis equal;

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1]);

% dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ux_real_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ux_real_m1st0135_x_D_slices.eps'),'-depsc2','-r600');
