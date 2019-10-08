%% Written by Sheel Nidhan
%  Plotting the SPOD modes for a given x/D

Nfreq = 512;
Novlp = 256;
N     = 7000;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
mode_sampled = [0; 1; 2; 3; 4];
x_sampled =  [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];

nr = 354;
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nrows = numvar*nr*Nblk;
Nrows_permode = numvar*nr;
Nblk_sampled = 5; 
Nf_sampled = 50;

dir_modes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/matlab_files/modes/';
disp(dir_modes);

%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
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

u_eigenmode_allm = zeros(nr, Nblk_sampled, Nf_sampled, size(mode_sampled,1), size(x_sampled,1));
w_eigenmode_allm = zeros(nr, Nblk_sampled, Nf_sampled, size(mode_sampled,1), size(x_sampled,1));
v_eigenmode_allm = zeros(nr, Nblk_sampled, Nf_sampled, size(mode_sampled,1), size(x_sampled,1));

for x_d = 1:size(x_sampled,1)
for i_mode = 1:size(mode_sampled,1)
        filename = strcat(dir_modes, 'eigenmodes_x_D_', int2str(x_sampled(x_d,1)), '_m', int2str(mode_sampled(i_mode,1)), '.mat');
        disp(filename);
        load(filename);
        u_eigenmode_allm(:,:,:,i_mode, x_d) = u_eigenmode(:,:,:);
        v_eigenmode_allm(:,:,:,i_mode, x_d) = v_eigenmode(:,:,:);
        w_eigenmode_allm(:,:,:,i_mode, x_d) = w_eigenmode(:,:,:);
end
end
%% Normalizing eigenmodes

for x = 1:size(x_sampled,1)
for m = 1:size(mode_sampled,1)
    for fn = 1:Nf_sampled
        for Nb = 1:Nblk_sampled
            spod_mode = u_eigenmode_allm(:,Nb,fn,m,x);
%             spod_mode_mag = spod_mode.*conj(spod_mode);
%             spod_mode_mag_radial = spod_mode_mag.*rc;
            spod_mode_mag = norm(spod_mode, inf);
            alpha = spod_mode_mag;
%             alpha = trapz(rc, spod_mode_mag_radial);
            u_eigenmode_allm(:,Nb,fn,m,x) = u_eigenmode_allm(:,Nb,fn,m,x)/(alpha);
        end
    end
end
end

for x = 1:size(x_sampled,1)
for m = 1:size(mode_sampled,1)
    for fn = 1:Nf_sampled
        for Nb = 1:Nblk_sampled
            spod_mode = v_eigenmode_allm(:,Nb,fn,m,x);
%             spod_mode_mag = spod_mode.*conj(spod_mode);
%             spod_mode_mag_radial = spod_mode_mag.*rc;
            spod_mode_mag = norm(spod_mode, inf);
            alpha = spod_mode_mag;
%             alpha = trapz(rc, spod_mode_mag_radial);
            v_eigenmode_allm(:,Nb,fn,m,x) = v_eigenmode_allm(:,Nb,fn,m,x)/(alpha);
        end
    end
end
end

for x = 1:size(x_sampled,1)
for m = 1:size(mode_sampled,1)
    for fn = 1:Nf_sampled
        for Nb = 1:Nblk_sampled
            spod_mode = w_eigenmode_allm(:,Nb,fn,m,x);
%            spod_mode_mag = spod_mode.*conj(spod_mode);
%            spod_mode_mag_radial = spod_mode_mag.*rc;
            spod_mode_mag = norm(spod_mode, inf);
            alpha = spod_mode_mag;
%             alpha = trapz(rc, spod_mode_mag_radial);
            w_eigenmode_allm(:,Nb,fn,m,x) = w_eigenmode_allm(:,Nb,fn,m,x)/(alpha);
        end
    end
end
end

%% Find the peaks of interest in St = 0.136, m = 1 mode


% for  i = 1:size(loc_planes,1)
%     [pks, indexes] = findpeaks(abs(w_eigenmode_allm(:,1,6,2,i)).^2);
%     [pks_max, indx_max] = max(pks);
%     m1_st0136_local_max_loc(i,1) = rc(indexes(indx_max));
%     
% end
% 
% log_x_sampled = log(LK_TKE_loc_planes(:,1));
% log_m1_st0136_local_max_loc = log(m1_st0136_local_max_loc(1:end,1));
% [coeffs_m2, S_m2] = polyfit(log_x_sampled, log_m1_st0136_local_max_loc, 1);
% y_m2_fitted = polyval(coeffs_m2, log_x_sampled);
% CI_m2 = polyparci(coeffs_m2, S_m2, 0.99);
% mse_m2 = immse(y_m2_fitted, log_m1_st0136_local_max_loc);

% loading the xlsx file of manual peak finding

% A = importdata('/home/sheel/Dropbox/spod_re5e4/spod_analysis_v2.0/m1_st0136_mode_maximas.xlsx');
% 
% log_x_sampled = log(LK_TKE_loc_planes(6:end-1,1));
% log_m1_st0136_local_max_loc = log(A(2,6:end-1))';
% [coeffs_m2, S_m2] = polyfit(log_x_sampled, log_m1_st0136_local_max_loc, 1);
% y_m2_fitted = polyval(coeffs_m2, log_x_sampled);
% CI_m2 = polyparci(coeffs_m2, S_m2, 0.99);
% mse_m2 = immse(y_m2_fitted, log_m1_st0136_local_max_loc);
% h1 = loglog(LK_TKE_loc_planes(6:end-1,1), A(2,6:end-1)', 'ko', 'MarkerSize',7);
% hold on;
% 
% h2 = loglog(LK_TKE_loc_planes(6:end-1,1), exp(y_m2_fitted), 'r-', 'LineWidth', 2);

%% Plotting w eigenmodes of interest

% dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';
% save(strcat(dirout, 'eigenmodes_similarity_diff_loc.mat'), 'f', 'rc', 'w_eigenmode_allm', 'u_eigenmode_allm', 'v_eigenmode_allm', ...
%                     'TKE_centerline_loc_planes', 'LK_TKE_loc_planes', 'ud_centerline_loc_planes', 'LK_mean_loc_planes', ...
%                     'mke_area_loc_planes', 'tke_area_loc_planes');

% close all;

figure;
hold on
h1 = plot(rc, abs(w_eigenmode_allm(:,1,6,2,4)).^2,'bo','LineWidth',2); 
h2 = plot(rc, abs(w_eigenmode_allm(:,1,1,3,8)).^2,'r-','LineWidth',2);
h3 = plot(rc, abs(w_eigenmode_allm(:,1,1,3,16)).^2,'k-','LineWidth',2);
%h4 = plot(rc, abs(w_eigenmode_allm(:,1,1,4,20)).^2,'c-','LineWidth',2);
% h5 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,6)).^2,'m-','LineWidth',2);
% h6 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,7)).^2,'g-','LineWidth',2);
% h7 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,8)).^2,'y-','LineWidth',2);
% h8 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,9)).^2,'y--','LineWidth',2);
% h9 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,10)).^2,'g--','LineWidth',2);
% h10 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,11)).^2,'m--','LineWidth',2);
% h11 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,12)).^2,'c--','LineWidth',2);
% h12 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,13)).^2,'k--','LineWidth',2);
% h13 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,14)).^2,'r--','LineWidth',2);
% h14 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,15)).^2,'b--','LineWidth',2); %#ok<*NASGU>
% h15 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,16)).^2,'g-','LineWidth',2); %#ok<*NASGU>
% h16 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,17)).^2,'m-','LineWidth',2); %#ok<*NASGU>
% h17 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,18)).^2,'c-','LineWidth',2); %#ok<*NASGU>
% h18 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,19)).^2,'k-','LineWidth',2); %#ok<*NASGU>
% h19 = plot(rc, abs(w_eigenmode_allm(:,1,6,4,20)).^2,'r-','LineWidth',2);
% h20 = plot(rc, abs(w_eigenmode_allm(:,1,6,4 ,21)).^2,'b-','LineWidth',2);

ylim([0 1.2]);
% xlim([0 4]);

figure;
hold on
h1 = plot(rc, abs(u_eigenmode_allm(:,1,6,2,4)).^2,'bo','LineWidth',2); 
h2 = plot(rc, abs(u_eigenmode_allm(:,1,1,3,8)).^2,'ro','LineWidth',2);
h3 = plot(rc, abs(u_eigenmode_allm(:,1,1,3,16)).^2,'ko','LineWidth',2);
% h4 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,22)).^2,'c-','LineWidth',2);
% h5 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,6)).^2,'m-','LineWidth',2);
% h6 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,7)).^2,'g-','LineWidth',2);
% h7 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,8)).^2,'y-','LineWidth',2);
% h8 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,9)).^2,'y--','LineWidth',2);
% h9 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,10)).^2,'g--','LineWidth',2);
% h10 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,11)).^2,'m--','LineWidth',2);
% h11 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,12)).^2,'c--','LineWidth',2);
% h12 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,13)).^2,'k--','LineWidth',2);
% h13 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,14)).^2,'r--','LineWidth',2);
% h14 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,15)).^2,'b--','LineWidth',2); %#ok<*NASGU>
% h15 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,16)).^2,'g-','LineWidth',2); %#ok<*NASGU>
% h16 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,17)).^2,'m-','LineWidth',2); %#ok<*NASGU>
% h17 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,18)).^2,'c-','LineWidth',2); %#ok<*NASGU>
% h18 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,19)).^2,'k-','LineWidth',2); %#ok<*NASGU>
% h19 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,20)).^2,'r-','LineWidth',2);
% h20 = plot(rc, abs(u_eigenmode_allm(:,1,6,4,21)).^2,'b-','LineWidth',2);

ylim([0 1.2]);


% figure;
% hold on
% h1 = plot(rc/LK_TKE_loc_planes(2,2), abs(w_eigenmode_allm(:,1,6,2,2)).^2,'b-','LineWidth',2); 
% h2 = plot(rc/LK_TKE_loc_planes(3,2), abs(w_eigenmode_allm(:,1,6,2,3)).^2,'r-','LineWidth',2);
% h3 = plot(rc/LK_TKE_loc_planes(4,2), abs(w_eigenmode_allm(:,1,6,2,4)).^2,'k-','LineWidth',2);
% h4 = plot(rc/LK_TKE_loc_planes(5,2), abs(w_eigenmode_allm(:,1,6,2,5)).^2,'c-','LineWidth',2);
% h5 = plot(rc/LK_TKE_loc_planes(6,2), abs(w_eigenmode_allm(:,1,6,2,6)).^2,'m-','LineWidth',2);
% h6 = plot(rc/LK_TKE_loc_planes(7,2), abs(w_eigenmode_allm(:,1,6,2,7)).^2,'g-','LineWidth',2);
% h7 = plot(rc/LK_TKE_loc_planes(8,2), abs(w_eigenmode_allm(:,1,6,2,8)).^2,'y-','LineWidth',2);
% h8 = plot(rc/LK_TKE_loc_planes(9,2), abs(w_eigenmode_allm(:,1,6,2,9)).^2,'y--','LineWidth',2);
% h9 = plot(rc/LK_TKE_loc_planes(10,2), abs(w_eigenmode_allm(:,1,6,2,10)).^2,'g--','LineWidth',2);
% h10 = plot(rc/LK_TKE_loc_planes(11,2), abs(w_eigenmode_allm(:,1,6,2,11)).^2,'m--','LineWidth',2);
% h11 = plot(rc/LK_TKE_loc_planes(12,2), abs(w_eigenmode_allm(:,1,6,2,12)).^2,'c--','LineWidth',2);
% h12 = plot(rc/LK_TKE_loc_planes(13,2), abs(w_eigenmode_allm(:,1,6,2,13)).^2,'k--','LineWidth',2);
% h13 = plot(rc/LK_TKE_loc_planes(14,2), abs(w_eigenmode_allm(:,1,6,2,14)).^2,'r--','LineWidth',2);
% h14 = plot(rc/LK_TKE_loc_planes(15,2), abs(w_eigenmode_allm(:,1,6,2,15)).^2,'b--','LineWidth',2); %#ok<*NASGU>
% h15 = plot(rc/LK_TKE_loc_planes(16,2), abs(w_eigenmode_allm(:,1,6,2,16)).^2,'g-','LineWidth',2); %#ok<*NASGU>
% h16 = plot(rc/LK_TKE_loc_planes(17,2), abs(w_eigenmode_allm(:,1,6,2,17)).^2,'m-','LineWidth',2); %#ok<*NASGU>
% h17 = plot(rc/LK_TKE_loc_planes(18,2), abs(w_eigenmode_allm(:,1,6,2,18)).^2,'c-','LineWidth',2); %#ok<*NASGU>
% h18 = plot(rc/LK_TKE_loc_planes(19,2), abs(w_eigenmode_allm(:,1,6,2,19)).^2,'k-','LineWidth',2); %#ok<*NASGU>
% h19 = plot(rc/LK_TKE_loc_planes(20,2), abs(w_eigenmode_allm(:,1,6,2,20)).^2,'r-','LineWidth',2);
% h20 = plot(rc/LK_TKE_loc_planes(21,2), abs(w_eigenmode_allm(:,1,6,2,21)).^2,'b-','LineWidth',2);
% 
% ylim([0 1.2]);
% xlim([0 5]);

figure;
hold on
h1 = plot(rc/LK_TKE_loc_planes(8,2), abs(u_eigenmode_allm(:,1,6,4,8)),'b-','LineWidth',2); 
h2 = plot(rc/LK_TKE_loc_planes(9,2), abs(u_eigenmode_allm(:,1,6,4,9)),'r-','LineWidth',2);
h3 = plot(rc/LK_TKE_loc_planes(10,2), abs(u_eigenmode_allm(:,1,6,4,10)),'k-','LineWidth',2);
h4 = plot(rc/LK_TKE_loc_planes(11,2), abs(u_eigenmode_allm(:,1,6,4,11)),'c-','LineWidth',2);
h5 = plot(rc/LK_TKE_loc_planes(12,2), abs(u_eigenmode_allm(:,1,6,4,12)),'m-','LineWidth',2);
h6 = plot(rc/LK_TKE_loc_planes(13,2), abs(u_eigenmode_allm(:,1,6,4,13)),'g-','LineWidth',2);
h7 = plot(rc/LK_TKE_loc_planes(14,2), abs(u_eigenmode_allm(:,1,6,4,14)),'y-','LineWidth',2);
h8 = plot(rc/LK_TKE_loc_planes(15,2), abs(u_eigenmode_allm(:,1,6,4,15)),'y--','LineWidth',2);
h9 = plot(rc/LK_TKE_loc_planes(16,2), abs(u_eigenmode_allm(:,1,6,4,16)),'g--','LineWidth',2);
h10 = plot(rc/LK_TKE_loc_planes(17,2), abs(u_eigenmode_allm(:,1,6,4,17)),'m--','LineWidth',2);
h11 = plot(rc/LK_TKE_loc_planes(18,2), abs(u_eigenmode_allm(:,1,6,4,18)),'c--','LineWidth',2);
h12 = plot(rc/LK_TKE_loc_planes(19,2), abs(u_eigenmode_allm(:,1,6,4,19)),'k--','LineWidth',2);
h13 = plot(rc/LK_TKE_loc_planes(20,2), abs(u_eigenmode_allm(:,1,6,4,20)),'r--','LineWidth',2);
h14 = plot(rc/LK_TKE_loc_planes(21,2), abs(u_eigenmode_allm(:,1,6,4,21)),'b--','LineWidth',2);
ylim([0 1.2]);
xlim([0 5]);

%% Plotting the contourmap plot for ux and ur modes for different modes
theta = linspace(0, 2*pi, 100);
mode  = real(w_eigenmode_allm(:,1,1,3,8))*sin(2*theta);
[R THETA] = meshgrid(rc, theta);

figure;
[C h] = polarcont(rc,theta,mode);
axis equal;
colorbar;
set(h,'EdgeColor','none')


%% Putting the axes
hXLabel = xlabel('$r$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|^{2}$/$|\Phi_{x}|_{max}$(m = 3, St = 0, Mode1)','interpreter','latex','fontsize',15);
%hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3], 'x/D = 20', 'x/D = 40', 'x/D = 80');


% hLegend = legend([h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14], 'x/D = 40', 'x/D = 45','x/D = 50',  ...
%    'x/D = 55', 'x/D = 60', 'x/D = 65', 'x/D = 70', 'x/D = 75', 'x/D = 80', 'x/D = 85', 'x/D = 90', 'x/D = 95', 'x/D = 100', 'x/D = 110');

%hLegend = legend([h7, h8, h9, h10, h11, h12, h13, h14], 'x/D = 70', 'x/D = 75', 'x/D = 80', 'x/D = 85', 'x/D = 90', 'x/D = 95', 'x/D = 100', 'x/D = 110');


% hLegend = legend([h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14], 'x/D = 10','x/D = 15',  ...
%                                    'x/D = 20', 'x/D = 25', 'x/D = 30', 'x/D = 35', 'x/D = 40', 'x/D = 45', 'x/D = 50', 'x/D = 55', 'x/D = 60', 'x/D = 65', 'x/D = 70', 'x/D = 75');

% hLegend = legend([h9, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20], 'x/D = 50', 'x/D = 55', 'x/D = 60', 'x/D = 65', 'x/D = 70', 'x/D = 75', 'x/D = 80', 'x/D = 85', ...
%                                     'x/D  = 90', 'x/D = 95', 'x/D = 100', 'x/D = 110');

hLegend.Interpreter = 'Latex';
hLegend.FontSize = 10;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'phi_ux_m3_st0_spodmode1_x_D_20_40_80.png','-dpng','-r600');  
% print(gcf,'phi_ux_m3_st0136_spodmode1_lk_50_110.eps','-depsc','-r600');

%% Plotting a contour map 

% m = 2;
% theta = linspace(0,2*pi,100);
% f_theta =   exp(sqrt(-1)*m.*theta);

% Plotting St = 0 and m = 2 plot at different locations
% for i = 1:size(LK_mean_loc_planes,1)
%     u_eigenmode_2d_m2(:,:,i) = squeeze(u_eigenmode_allm(:,1,1,3,i))*f_theta;
%     w_eigenmode_2d_m2(:,:,i) = squeeze(w_eigenmode_allm(:,1,1,3,i))*f_theta;
% end
% 
% 
% for j = 1:size(rc,1)
%     for k = 1:size(theta)
%         y(j,k) = r(j)*cos(theta(k));
%         z(j,k) = r(j)*sin(theta(k));
%     end
% end

% Y = repmat(y, 1, 1, 22);
% Z = repmat(z, 1, 1, 22);
% 
% for i  = 1:size(LK_mean_loc_planes,1)
%     X(:,:,i) = LK_mean_loc_planes(i,1)*ones(354,100);
% end
% 
% xslice = LK_mean_loc_planes(:,1)';   
% yslice = [];
% zslice = [];
% contourslice(X,Y,Z,real(u_eigenmode_2d_m2),xslice,yslice,zslice)
% view(3)
% grid on
%axis equal;