%% Written by Sheel Nidhan
%  Plotting the SPOD eigenvalues for a given x/D for similarity analysis

clear; clc; close all;
Nfreq             = 512;
Novlp             = 256;
N                 = 7000;
stride            = 100;
nstart            = 1892600;
nend              = nstart + (N-1)*stride;
mode_sampled      = linspace(0,10,11)';
x_sampled         =  [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
loc_planes        = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
half_width        = [1.6378; 2.1689; 2.6075; 2.74025; 2.7666; 2.8444; 2.8634; 2.8839; 3.0012];
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nblk_sampled = 27; 
Nf_sampled = 100;
dir_modes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/matlab_files/spectrum/';
disp(dir_modes);

%% Reading the time file

time = importdata ('/home/sheel/Work/projects/spod_re5e4/post/frinf/time_stamps/time_stamp_1892600_2613200_uniform.txt');
time_spod = time(1:N,1);
dt   = time_spod(2:end,1) - time_spod(1:end-1,1);

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



%% Plotting eigenspectra of interest
% 
% dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';
% save(strcat(dirout, 'eigvalues_similarity_diff_loc.mat'), 'f', 'eigenspectra_allm', 'TKE_centerline_loc_planes', 'LK_TKE_loc_planes', ...
%                     'ud_centerline_loc_planes', 'LK_mean_loc_planes', 'mke_area_loc_planes', 'tke_area_loc_planes');

% m=2, normalized by (K^0.5/2*L_k)^2
figure;
hold on;
h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,8)/((TKE_centerline_loc_planes(8,2)^0.5)*LK_TKE_loc_planes(8,2))^2,     'bo','MarkerSize',7); %#ok<*NASGU>
h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,9)/((TKE_centerline_loc_planes(9,2)^0.5)*LK_TKE_loc_planes(9,2))^2,     'r*','MarkerSize',7);
h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,10)/((TKE_centerline_loc_planes(10,2)^0.5)*LK_TKE_loc_planes(10,2))^2,  'kd','MarkerSize',7);
h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,11)/((TKE_centerline_loc_planes(11,2)^0.5)*LK_TKE_loc_planes(11,2))^2,  'c>','MarkerSize',7);
h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,12)/((TKE_centerline_loc_planes(12,2)^0.5)*LK_TKE_loc_planes(12,2))^2,  'g<','MarkerSize',7);
h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,13)/((TKE_centerline_loc_planes(13,2)^0.5)*LK_TKE_loc_planes(13,2))^2,  'b+','MarkerSize',7);
h7 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,14)/((TKE_centerline_loc_planes(14,2)^0.5)*LK_TKE_loc_planes(14,2))^2,  'r^','MarkerSize',7);
h8 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,15)/((TKE_centerline_loc_planes(15,2)^0.5)*LK_TKE_loc_planes(15,2))^2,  'kv','MarkerSize',7);
h9 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,16)/((TKE_centerline_loc_planes(16,2)^0.5)*LK_TKE_loc_planes(16,2))^2,  'gp','MarkerSize',7);
h10 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,17)/((TKE_centerline_loc_planes(17,2)^0.5)*LK_TKE_loc_planes(18,2))^2, 'ch','MarkerSize',7);
h11 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,18)/((TKE_centerline_loc_planes(18,2)^0.5)*LK_TKE_loc_planes(18,2))^2, 'bx','MarkerSize',7);
h12 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,19)/((TKE_centerline_loc_planes(19,2)^0.5)*LK_TKE_loc_planes(19,2))^2, 'b<','MarkerSize',7);
h13 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,20)/((TKE_centerline_loc_planes(20,2)^0.5)*LK_TKE_loc_planes(20,2))^2, 'mo','MarkerSize',7);
h14 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,21)/((TKE_centerline_loc_planes(21,2)^0.5)*LK_TKE_loc_planes(21,2))^2, 'md','MarkerSize',7);
h15 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,5,22)/((TKE_centerline_loc_planes(22,2)^0.5)*LK_TKE_loc_planes(22,2))^2, 'm^','MarkerSize',7);
xlim([0 0.4]);
ylim([0 0.05]);

% m=0, normalized by (K^0.5/2*L_k)^2
% figure;
% hold on;
% h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,8)/((TKE_centerline_loc_planes(8,2)^0.5)*LK_TKE_loc_planes(8,2))^2,     'bo','MarkerSize',7);
% h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,9)/((TKE_centerline_loc_planes(9,2)^0.5)*LK_TKE_loc_planes(9,2))^2,     'r*','MarkerSize',7);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,10)/((TKE_centerline_loc_planes(10,2)^0.5)*LK_TKE_loc_planes(10,2))^2,  'kd','MarkerSize',7);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,11)/((TKE_centerline_loc_planes(11,2)^0.5)*LK_TKE_loc_planes(11,2))^2,  'c>','MarkerSize',7);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,12)/((TKE_centerline_loc_planes(12,2)^0.5)*LK_TKE_loc_planes(12,2))^2,  'g<','MarkerSize',7);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,13)/((TKE_centerline_loc_planes(13,2)^0.5)*LK_TKE_loc_planes(13,2))^2,  'b+','MarkerSize',7);
% h7 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,14)/((TKE_centerline_loc_planes(14,2)^0.5)*LK_TKE_loc_planes(14,2))^2,  'r^','MarkerSize',7);
% h8 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,15)/((TKE_centerline_loc_planes(15,2)^0.5)*LK_TKE_loc_planes(15,2))^2,  'kv','MarkerSize',7);
% h9 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,16)/((TKE_centerline_loc_planes(16,2)^0.5)*LK_TKE_loc_planes(16,2))^2,  'gp','MarkerSize',7);
% h10 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,17)/((TKE_centerline_loc_planes(17,2)^0.5)*LK_TKE_loc_planes(18,2))^2, 'ch','MarkerSize',7);
% h11 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,18)/((TKE_centerline_loc_planes(18,2)^0.5)*LK_TKE_loc_planes(18,2))^2, 'bx','MarkerSize',7);
% h12 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,19)/((TKE_centerline_loc_planes(19,2)^0.5)*LK_TKE_loc_planes(19,2))^2, 'b<','MarkerSize',7);
% h13 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,20)/((TKE_centerline_loc_planes(20,2)^0.5)*LK_TKE_loc_planes(20,2))^2, 'mo','MarkerSize',7);
% h14 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,21)/((TKE_centerline_loc_planes(21,2)^0.5)*LK_TKE_loc_planes(21,2))^2, 'md','MarkerSize',7);
% h15 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,22)/((TKE_centerline_loc_planes(22,2)^0.5)*LK_TKE_loc_planes(22,2))^2, 'm^','MarkerSize',7);
% xlim([0 0.4]);
% ylim([0 0.3]);

% m=1 normalized by (U_dL_mean)^2
% figure;
% hold on
% h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,12)/((ud_centerline_loc_planes(1,2)^1)*LK_mean_loc_planes(12,2))^2,     'bo','MarkerSize',7);
% h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,13)/((ud_centerline_loc_planes(2,2)^1)*LK_mean_loc_planes(13,2))^2,     'r*','MarkerSize',7);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,14)/((ud_centerline_loc_planes(3,2)^1)*LK_mean_loc_planes(14,2))^2,  'kd','MarkerSize',7);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,15)/((ud_centerline_loc_planes(4,2)^1)*LK_mean_loc_planes(15,2))^2,  'c>','MarkerSize',7);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,16)/((ud_centerline_loc_planes(5,2)^1)*LK_mean_loc_planes(16,2))^2,  'g<','MarkerSize',7);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,17)/((ud_centerline_loc_planes(6,2)^1)*LK_mean_loc_planes(17,2))^2,  'b+','MarkerSize',7);
% h7 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,18)/((ud_centerline_loc_planes(7,2)^1)*LK_mean_loc_planes(18,2))^2,  'r^','MarkerSize',7);
% h8 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,19)/((ud_centerline_loc_planes(8,2)^1)*LK_mean_loc_planes(19,2))^2,  'kv','MarkerSize',7);
% h9 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,20)/((ud_centerline_loc_planes(9,2)^1)*LK_mean_loc_planes(20,2))^2,  'gp','MarkerSize',7);
% h10 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,21)/((ud_centerline_loc_planes(10,2)^1)*LK_mean_loc_planes(21,2))^2, 'ch','MarkerSize',7);
% h11 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,11)/((ud_centerline_loc_planes(11,2)^1)*LK_mean_loc_planes(11,2))^2, 'bx','MarkerSize',7);
% h12 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,12)/((ud_centerline_loc_planes(12,2)^1)*LK_mean_loc_planes(12,2))^2, 'b<','MarkerSize',7);
% h13 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,13)/((ud_centerline_loc_planes(13,2)^1)*LK_mean_loc_planes(13,2))^2, 'mo','MarkerSize',7);
% h14 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,14)/((ud_centerline_loc_planes(14,2)^1)*LK_mean_loc_planes(14,2))^2, 'md','MarkerSize',7);
% h15 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,15)/((ud_centerline_loc_planes(15,2)^1)*LK_mean_loc_planes(15,2))^2, 'm^','MarkerSize',7);
% xlim([0 0.4]);
% ylim([0 0.3]);


% m=1 normalized by (K^0.5/2*L_k)^2
% figure;
% hold on;
% h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,1)/((TKE_centerline_loc_planes(1,2)^0.5)*LK_TKE_loc_planes(1,2))^2,     'bo','MarkerSize',7);
% h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,2)/((TKE_centerline_loc_planes(2,2)^0.5)*LK_TKE_loc_planes(2,2))^2,     'r*','MarkerSize',7);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,3)/((TKE_centerline_loc_planes(3,2)^0.5)*LK_TKE_loc_planes(3,2))^2,  'kd','MarkerSize',7);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,4)/((TKE_centerline_loc_planes(4,2)^0.5)*LK_TKE_loc_planes(4,2))^2,  'c>','MarkerSize',7);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,5)/((TKE_centerline_loc_planes(5,2)^0.5)*LK_TKE_loc_planes(5,2))^2,  'g<','MarkerSize',7);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,6)/((TKE_centerline_loc_planes(6,2)^0.5)*LK_TKE_loc_planes(6,2))^2,  'b+','MarkerSize',7);
% h7 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,7)/((TKE_centerline_loc_planes(7,2)^0.5)*LK_TKE_loc_planes(7,2))^2,  'r^','MarkerSize',7);
% h8 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,8)/((TKE_centerline_loc_planes(8,2)^0.5)*LK_TKE_loc_planes(8,2))^2,  'kv','MarkerSize',7);
% h9 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,9)/((TKE_centerline_loc_planes(9,2)^0.5)*LK_TKE_loc_planes(9,2))^2,  'gp','MarkerSize',7);
% h10 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,10)/((TKE_centerline_loc_planes(10,2)^0.5)*LK_TKE_loc_planes(10,2))^2, 'ch','MarkerSize',7);
% h11 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,11)/((TKE_centerline_loc_planes(11,2)^0.5)*LK_TKE_loc_planes(11,2))^2, 'bx','MarkerSize',7);
% h12 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,12)/((TKE_centerline_loc_planes(12,2)^0.5)*LK_TKE_loc_planes(12,2))^2, 'b<','MarkerSize',7);
% h13 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,13)/((TKE_centerline_loc_planes(13,2)^0.5)*LK_TKE_loc_planes(13,2))^2, 'mo','MarkerSize',7);
% h14 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,14)/((TKE_centerline_loc_planes(14,2)^0.5)*LK_TKE_loc_planes(14,2))^2, 'md','MarkerSize',7);
% h15 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,15)/((TKE_centerline_loc_planes(15,2)^0.5)*LK_TKE_loc_planes(15,2))^2, 'm^','MarkerSize',7);
% xlim([0 0.4]);
% ylim([0 0.3]);


% m=1 unnormalized
% figure;
% hold on;
% h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,1),     'bo','MarkerSize',7);
% h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,2),     'r*','MarkerSize',7);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,3),  'kd','MarkerSize',7);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,4),  'c>','MarkerSize',7);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,5),  'g<','MarkerSize',7);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,6),  'b+','MarkerSize',7);
% h7 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,7),  'r^','MarkerSize',7);
% h8 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,8),  'kv','MarkerSize',7);
% h9 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,9),  'gp','MarkerSize',7);
% h10 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,10), 'ch','MarkerSize',7);
% h11 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,11), 'bx','MarkerSize',7);
% h12 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,12), 'b<','MarkerSize',7);
% h13 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,13), 'mo','MarkerSize',7);
% h14 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,14), 'md','MarkerSize',7);
% h15 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,15), 'm^','MarkerSize',7);
% xlim([0 0.4]);
% ylim([0 0.005]);

% figure;
% hold on
% h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,4)/((tke_centerline(4,1)^0.5)*half_width_tke(4,1)^1.5)^2, 'bo','MarkerSize',7);
% h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,5)/((tke_centerline(5,1)^0.5)*half_width_tke(5,1)^1.5)^2, 'rs','MarkerSize',7);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,6)/((tke_centerline(6,1)^0.5)*half_width_tke(6,1)^1.5)^2, 'kd','MarkerSize',7);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,7)/((tke_centerline(7,1)^0.5)*half_width_tke(7,1)^1.5)^2, 'k>','MarkerSize',7);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,8)/((tke_centerline(8,1)^0.5)*half_width_tke(8,1)^1.5)^2, 'r<','MarkerSize',7);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,1,9)/((tke_centerline(9,1)^0.5)*half_width_tke(9,1)^1.5)^2, 'bp','MarkerSize',7);
% xlim([0 0.4]);


% % h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,1)/(defect_centerline(1,1)*half_width(1,1)), 'b-','LineWidth',2);
% % h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,2)/(defect_centerline(2,1)*half_width(2,1)),'r-','LineWidth',2);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,3)/(defect_centerline(3,1)*half_width(3,1)),'k-','LineWidth',2);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,4)/(defect_centerline(4,1)*half_width(4,1)),'k--','LineWidth',2);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,5)/(defect_centerline(5,1)*half_width(5,1)),'r--','LineWidth',2);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,6)/(defect_centerline(6,1)*half_width(6,1)),'b--','LineWidth',2);


% % h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,1)/(integrated_tke(1,1)), 'b-','LineWidth',2);
% % h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,2)/(integrated_tke(2,1)),'r-','LineWidth',2);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,3)/(integrated_tke(3,1)),'k-','LineWidth',2);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,4)/(integrated_tke(4,1)),'k--','LineWidth',2);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,5)/(integrated_tke(5,1)),'r--','LineWidth',2);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,6)/(integrated_tke(6,1)),'b--','LineWidth',2);
 
% 
% figure
% hold on
% h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,4), 'bo','MarkerSize',7);
% h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,5), 'rs','MarkerSize',7);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,6), 'kd','MarkerSize',7);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,7), 'k>','MarkerSize',7);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,8), 'r<','MarkerSize',7);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,2,9), 'bp','MarkerSize',7);
% xlim([0 0.4]);


% h1 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,1)/sum(sum(eigenspectra_allm(:,:,2,1))), 'b-','LineWidth',2);
% h2 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,2)/sum(sum(eigenspectra_allm(:,:,2,2))),'r-','LineWidth',2);
% h3 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,3)/sum(sum(eigenspectra_allm(:,:,2,3))),'k-','LineWidth',2);
% h4 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,4)/sum(sum(eigenspectra_allm(:,:,2,4))),'k--','LineWidth',2);
% h5 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,5)/sum(sum(eigenspectra_allm(:,:,2,5))),'r--','LineWidth',2);
% h6 = plot(f(1:Nf_sampled)', eigenspectra_allm(:,1,3,6)/sum(sum(eigenspectra_allm(:,:,2,6))),'b--','LineWidth',2);


% hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$\lambda^{(1)}(m=1,St=0.136)/(U_{d}L_{d})^{2}$','interpreter','latex','fontsize',15);
% % %hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);
% 
% hLegend = legend([ h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15], 'x/D = 40', 'x/D = 45', 'x/D = 50', 'x/D = 55', 'x/D = 60', 'x/D = 65', ...
%             'x/D = 70', 'x/D = 75', 'x/D = 80', 'x/D = 85', 'x/D = 90', 'x/D = 95', 'x/D = 100', 'x/D = 110', 'x/D = 120');
        
% hLegend = legend([ h1, h2, h3, h4, h5, h6, h7, h8, h9, h10], 'x/D = 5', 'x/D = 10', 'x/D = 15', 'x/D = 20', 'x/D = 25', 'x/D = 30', ...
%                                                            'x/D = 35', 'x/D = 40', 'x/D = 45', 'x/D = 50');
                                                       
                                                               
% hLegend = legend([ h1, h2, h3, h4, h5, h6, h7, h8, h9, h10], 'x/D = 5', 'x/D = 10', 'x/D = 15', 'x/D = 20', 'x/D = 25', 'x/D = 30', ...
%                                                            'x/D = 35', 'x/D = 40', 'x/D = 45', 'x/D = 50');
%         
        
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 10;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'eigenspectra_m1_mode1_unnormalized_udefectld_x_D_5_50.png','-dpng','-r600');  
% print(gcf,'eigenspectra_m1_mode1_normalized_udefectld_x_D_5_50.eps','-depsc','-r600');