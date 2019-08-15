%% Written by Sheel Nidhan
%  Plotting the SPOD eigenvalues for a given x/D
clear; clc; close all;
Nfreq             = 512;
Novlp             = 256;
N                 = 7000;
stride            = 100;
nstart            = 1892600;
nend              = nstart + (N-1)*stride;
mode_sampled      = linspace(0,10,11)';
x_sampled         = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];

% half_width        = [1.6378; 2.1689; 2.6075; 2.74025; 2.7666; 2.8444; 2.8634; 2.8839; 3.0012];
% half_width_tke    = [1.9557; 2.5796; 2.9654; 3.1111; 3.24953; 3.2983; 3.3487; 3.4004; 3.4227];
% tke_centerline    = [0.0036; 0.001308; 0.0007469; 0.0005866; 0.000493; 0.0004610; 0.0004283; 0.00040009; 0.000379];
% integrated_tke    = [0.04867; 0.03149; 0.02356; 0.01914; 0.01605; 0.01373];
% defect_centerline = [0.0497; 0.02575; 0.017019; 0.015231; 0.01476; 0.013879; 0.013694; 0.013363; 0.01249];
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nblk_sampled = 26; 
Nf_sampled = 100;
dir_modes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/matlab_files/spectrum/';
disp(dir_modes);

%% Load the 5p grid files
x_5p_grid = importdata('/home/sheel/Work/projects/spod_re5e4/grid/frinf/5p_grid_points.txt');

for i = 1:size(x_sampled,1)
   [val, idx] = min(abs(x_5p_grid-x_sampled(i,1)));
   x_sample_closest(i,1) =  x_5p_grid(idx,1);
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

%% Loading the matlabfiles for SPOD modes

eigenspectra_allm = zeros(Nf_sampled, Nblk_sampled, size(mode_sampled,1), size(x_sampled,1));
for x_d = 1:size(x_sampled,1)
for i_mode = 1:size(mode_sampled,1)
        filename = strcat(dir_modes, 'spectrum_x_D_', int2str(x_sampled(x_d,1)), '.mat');
        disp(filename);
        load(filename);
        eigenspectra_allm(:,:,i_mode,x_d) = eigvalue(1:Nf_sampled,1:Nblk_sampled,i_mode);
end
end


%% Importing Karu's datafiles

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_TKE.dat';
LK_TKE   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_TKE(:,1)));
    LK_TKE_loc_planes(i,2) = LK_TKE(idx,4);
    LK_TKE_loc_planes(i,1) = LK_TKE(idx,1);

end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_WMEAN.dat';
LK_mean   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_mean(:,1)));
    LK_mean_loc_planes(i,2) = LK_mean(idx,4);
    LK_mean_loc_planes(i,1) = LK_mean(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_centerline.dat';
TKE_Centerline   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-TKE_Centerline(:,1)));
    TKE_centerline_loc_planes(i,2) = TKE_Centerline(idx,2);
    TKE_centerline_loc_planes(i,1) = TKE_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/UX_rms_centerline.dat';
ux_Centerline   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ux_Centerline(:,1)));
    ux_centerline_loc_planes(i,2) = ux_Centerline(idx,2);
    ux_centerline_loc_planes(i,1) = ux_Centerline(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Defect_centerline.dat';
ud_Centerline   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ud_Centerline(:,1)));
    ud_centerline_loc_planes(i,2) = ud_Centerline(idx,2);
    ud_centerline_loc_planes(i,1) = ud_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Defect_centerline.dat';
ud_Centerline   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ud_Centerline(:,1)));
    ud_centerline_loc_planes(i,2) = ud_Centerline(idx,2);
    ud_centerline_loc_planes(i,1) = ud_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/MKE_areaI_FINF.dat';
mke_area   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-mke_area(:,1)));
    mke_area_loc_planes(i,2) = mke_area(idx,2);
    mke_area_loc_planes(i,1) = mke_area(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_areaI_FINF.dat';
tke_area   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-tke_area(:,1)));
    tke_area_loc_planes(i,2) = tke_area(idx,2);
    tke_area_loc_planes(i,1) = tke_area(idx,1);
end



%% Best fit line for the SPOD mode m=1, St = 0.136

log_eigmode_m1_st0136 = log(squeeze(eigenspectra_allm(6,1,2,:)));
log_x_sampled = log(x_sample_closest);

% Scaling till x/D = 60
log_eigmode_m1_st0136_60 = log_eigmode_m1_st0136(1:11,1);
log_x_sampled_60 = log_x_sampled(1:11,1);
[coeffs_60, S_60] = polyfit(log_x_sampled_60, log_eigmode_m1_st0136_60, 1);
y_60_fitted = polyval(coeffs_60, log_x_sampled_60);
CI_60 = polyparci(coeffs_60, S_60, 0.99);
% Scaling from x/D = 60 to end

log_eigmode_m1_st0136_100 = log_eigmode_m1_st0136(14:20,1);
log_x_sampled_100 = log_x_sampled(14:20,1);
[coeffs_100, S_100] = polyfit(log_x_sampled_100, log_eigmode_m1_st0136_100, 1);
y_100_fitted = polyval(coeffs_100, log_x_sampled_100);
CI_100 = polyparci(coeffs_100, S_100, 0.99);


area_integrated_mke = (ud_centerline_loc_planes(:,2).*LK_mean_loc_planes(:,2)).^2;
log_mke = log(area_integrated_mke(1:11,1));
log_x_sampled_mke = log(x_sample_closest(1:11,1));
[coeffs_mke, S_mke] = polyfit(log_x_sampled_mke, log_mke, 1);
y_mke_fitted = polyval(coeffs_mke, log_x_sampled_mke);
CI_mke = polyparci(coeffs_mke, S_mke, 0.99);
% Scaling from x/D = 60 to end

figure;
loglog((x_sample_closest(1:11,1)), exp(y_mke_fitted) , 'r-', 'LineWidth', 2);
hold on;
loglog(x_sample_closest(1:end-2,1), area_integrated_mke(1:end-2,1), 'ko', 'MarkerSize',7);
h=text(10, 5*10^-3,'slope=-1.13','interpreter','latex','FontSize', 15);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$(U_{d}L_{mean})^{2}$','interpreter','latex','fontsize',15);

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'udlmean_square_streamwsie_variation.png','-dpng','-r600');  
% print(gcf,'udlmean_square_streamwise_variation.eps','-depsc','-r600');

figure;
h1 = loglog(x_sample_closest(1:end-2), squeeze(eigenspectra_allm(6,1,2,1:end-2)), 'ko', 'MarkerSize',7);
hold on;
h2 = loglog(x_sample_closest(1:end-2,1), area_integrated_mke(1:end-2,1), 'ro', 'MarkerSize',7);
hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);

hLegend = legend([ h1, h2], 'm=1, St=0.136, SPOD1', '$(U_{d}L_{mean})^{2}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

%set(gcf, 'PaperPositionMode', 'auto');
%print(gcf,'udlmean_square_m1_st0.136_streamwsie_variation.png','-dpng','-r600');  
% print(gcf,'udlmean_square_m1_st0.136_streamwise_variation.eps','-depsc','-r600');

%% Best fit line for the SPOD mode m=2, St=0

log_eigmode_m2_st0 = log(squeeze(eigenspectra_allm(1,1,3,:)));
log_x_sampled = log(x_sample_closest);

% % Scaling till x/D = 60
% log_eigmode_m2_st0_60 = log_eigmode_m2_st0(1:12,1);
% log_x_sampled_60 = log_x_sampled(1:12,1);
% [coeffs_m2_60, S_m2_60] = polyfit(log_x_sampled_60, log_eigmode_m2_st0_60, 1);
% y_m2_60_fitted = polyval(coeffs_m2_60, log_x_sampled_60);
% CI_m2_60 = polyparci(coeffs_m2_60, S_m2_60, 0.99);
% mse_m2_60 = immse(y_m2_60_fitted, log_eigmode_m2_st0_60);
% % Scaling from x/D = 60 to end
% 
% log_eigmode_m2_st0_100 = log_eigmode_m2_st0(13:20,1);
% log_x_sampled_100 = log_x_sampled(13:20,1);
% [coeffs_m2_100, S_m2_100] = polyfit(log_x_sampled_100, log_eigmode_m2_st0_100, 1);
% y_m2_100_fitted = polyval(coeffs_m2_100, log_x_sampled_100);
% CI_m2_100 = polyparci(coeffs_m2_100, S_m2_100, 0.99);
% mse_m2_100 = immse(y_m2_100_fitted, log_eigmode_m2_st0_100);

log_eigmode_m2_st0 = log_eigmode_m2_st0(8:end-2,1);
[coeffs_m2, S_m2] = polyfit(log_x_sampled(8:end-2), log_eigmode_m2_st0, 1);
y_m2_fitted = polyval(coeffs_m2, log_x_sampled(8:end-2));
CI_m2 = polyparci(coeffs_m2, S_m2, 0.99);
mse_m2 = immse(y_m2_fitted, log_eigmode_m2_st0);


area_integrated_tke = ((TKE_centerline_loc_planes(:,2).^(0.5)).*LK_TKE_loc_planes(:,2)).^2;
log_tke = log(area_integrated_tke(8:end-2,1));
log_x_sampled_tke = log(x_sample_closest(8:end-2,1));
[coeffs_tke, S_tke] = polyfit(log_x_sampled_tke, log_tke, 1);
y_tke_fitted = polyval(coeffs_tke, log_x_sampled_tke);
CI_tke = polyparci(coeffs_tke, S_tke, 0.99);
% Scaling from x/D = 60 to end

figure;
loglog((x_sample_closest(8:end-2,1)), exp(y_tke_fitted) , 'r-', 'LineWidth', 2);
hold on;
loglog(x_sample_closest(1:end-2,1), area_integrated_tke(1:end-2,1), 'ko', 'MarkerSize',7);
h=text(10, 5*10^-3,'slope=-1.13','interpreter','latex','FontSize', 15);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$(K^{1/2}L_{k})^{2}$','interpreter','latex','fontsize',15);

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'udlmean_square_streamwsie_variation.png','-dpng','-r600');  
% print(gcf,'udlmean_square_streamwise_variation.eps','-depsc','-r600');

figure;
h1 = loglog(x_sample_closest(1:end-2), squeeze(eigenspectra_allm(1,1,3,1:end-2)), 'ko', 'MarkerSize',7);
hold on;
h2 = loglog(x_sample_closest(1:end-2,1), area_integrated_tke(1:end-2,1), 'ro', 'MarkerSize',7);
hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);


hLegend = legend([ h1, h2], 'm=2, St=0, SPOD1', '$(K^{1/2}L_{k})^{2}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

%set(gcf, 'PaperPositionMode', 'auto');
%print(gcf,'k05cent_square_m2_st0_streamwsie_variation.png','-dpng','-r600');  
%print(gcf,'k05cent_square_m2_st0_streamwise_variation.eps','-depsc','-r600');
%% Plotting the eigenspectra of interest
figure;
h1 = loglog(x_sample_closest(1:end-2), squeeze(eigenspectra_allm(6,1,2,1:end-2)), 'ko', 'MarkerSize',7);
hold on;

h2 = loglog(x_sample_closest(1:11,1), exp(y_60_fitted), 'r-', 'LineWidth', 2);
h=text(10, 5*10^-3,'slope=-1.13','interpreter','latex','FontSize', 15);

h3 = loglog(x_sample_closest(14:20,1), exp(y_100_fitted), 'b-', 'LineWidth', 2);
h=text(25, 3*10^-4,'slope=-1.57','interpreter','latex','FontSize', 15);

h4 = loglog(x_sample_closest(1:end-2), squeeze(eigenspectra_allm(1,1,3,1:end-2)), 'ks', 'MarkerSize',7);

h6 = loglog(x_sample_closest(2:end), exp(y_m2_fitted), 'g-.', 'LineWidth', 2);
h=text(5, 10^-3,'slope=-0.59','interpreter','latex', 'FontSize', 15);


xlim([3 200]);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda^{(1)}(m=1, St=0.136; m=2, St=0)$','interpreter','latex','fontsize',15);
%hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

hLegend = legend([ h1, h4], 'm=1, St = 0.136', 'm=2, St=0');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 10;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'eigenspectra_m12_mode1_streamwsie_variation.png','-dpng','-r600');  
% print(gcf,'eigenspectra_m12_mode1_streamwise_variation.eps','-depsc','-r600');