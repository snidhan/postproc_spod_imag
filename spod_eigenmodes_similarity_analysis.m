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

%half_width = [1.6378; 2.1689; 2.6075; 2.7666; 3.0012; 3.1090];
%half_width_tke = [1.9557; 2.5796; 2.9654; 3.24593; 3.4227; 3.6530];

nr = 354;
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nrows = numvar*nr*Nblk;
Nrows_permode = numvar*nr;
Nblk_sampled = 5; 
Nf_sampled = 40;

dir_modes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/matlab_files/modes/';
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
            spod_mode_mag = norm(spod_mode);
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
             spod_mode_mag = norm(spod_mode);
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
             spod_mode_mag = norm(spod_mode);
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

A = importdata('/home/sheel/Dropbox/spod_re5e4/spod_analysis_v2.0/m1_st0136_mode_maximas.xlsx');

log_x_sampled = log(LK_TKE_loc_planes(6:end-1,1));
log_m1_st0136_local_max_loc = log(A(2,6:end-1))';
[coeffs_m2, S_m2] = polyfit(log_x_sampled, log_m1_st0136_local_max_loc, 1);
y_m2_fitted = polyval(coeffs_m2, log_x_sampled);
CI_m2 = polyparci(coeffs_m2, S_m2, 0.99);
mse_m2 = immse(y_m2_fitted, log_m1_st0136_local_max_loc);
h1 = loglog(LK_TKE_loc_planes(6:end-1,1), A(2,6:end-1)', 'ko', 'MarkerSize',7);
hold on;

h2 = loglog(LK_TKE_loc_planes(6:end-1,1), exp(y_m2_fitted), 'r-', 'LineWidth', 2);
%% Plotting w eigenmodes of interest

% close all;

figure;
hold on;
h1 = plot(rc/LK_TKE_loc_planes(14,2), abs(v_eigenmode_allm(:,1,6,2,8)).^2,'b-','LineWidth',2);
h2 = plot(rc/LK_TKE_loc_planes(15,2), abs(v_eigenmode_allm(:,1,6,2,9)).^2,'r-','LineWidth',2);
h3 = plot(rc/LK_TKE_loc_planes(16,2), abs(v_eigenmode_allm(:,1,6,2,10)).^2,'k-','LineWidth',2);
h4 = plot(rc/LK_TKE_loc_planes(17,2), abs(v_eigenmode_allm(:,1,6,2,11)).^2,'c-','LineWidth',2);
h5 = plot(rc/LK_TKE_loc_planes(18,2), abs(v_eigenmode_allm(:,1,6,2,12)).^2,'m-','LineWidth',2);
h6 = plot(rc/LK_TKE_loc_planes(19,2), abs(v_eigenmode_allm(:,1,6,2,13)).^2,'m--','LineWidth',2);
h7 = plot(rc/LK_TKE_loc_planes(20,2), abs(v_eigenmode_allm(:,1,6,2,14)).^2,'c--','LineWidth',2);
h8 = plot(rc/LK_TKE_loc_planes(21,2), abs(v_eigenmode_allm(:,1,6,2,15)).^2,'k--','LineWidth',2);
h9 = plot(rc/LK_TKE_loc_planes(22,2), abs(v_eigenmode_allm(:,1,6,2,16)).^2,'r--','LineWidth',2);
% 
figure;
hold on;
h1 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,8)).^2,'b-','LineWidth',2);
h2 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,9)).^2,'r-','LineWidth',2);
h3 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,10)).^2,'k-','LineWidth',2);
h4 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,11)).^2,'c-','LineWidth',2);
h5 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,12)).^2,'m-','LineWidth',2);
h6 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,13)).^2,'m--','LineWidth',2);
h7 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,14)).^2,'c--','LineWidth',2);
h8 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,15)).^2,'k--','LineWidth',2);
h9 = plot(rc, abs(v_eigenmode_allm(:,1,6,2,16)).^2,'r--','LineWidth',2);


ylim([0 0.02]);
xlim([0 10]);

%% Putting the axes
hXLabel = xlabel('$r/L_{Kmean}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|^{2}$(m = 1, St = 0.136, Mode1)','interpreter','latex','fontsize',15);
% %hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6, h7, h8, h9], 'x/D = 70', 'x/D = 75','x/D = 80',  ...
                                  'x/D = 85', 'x/D = 90', 'x/D = 95', 'x/D = 100', 'x/D = 110', 'x/D = 120');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

%%
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'phi_u_m1_st0136_spodmode1_tkebased_70_120.png','-dpng','-r600');  
print(gcf,'phi_u_m1_st0136_spodmode1_tkebased_70_120.eps','-depsc','-r600');