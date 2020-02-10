%% Written by Sheel Nidhan
clear; clc; close all;
%% SPOD Parameters
Nfreq  = 512;
Novlp  = 256;
N      = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
mode_sampled = [0; 1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%x_sampled =  [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
x_sampled = [40];
loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];

nr = 354;
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nrows = numvar*nr*Nblk;
Nrows_permode = numvar*nr;
Nblk_sampled = 20; 
Nf_sampled = 200;

dir_modes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/matlab_files/modes/';
dir_spec  = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/matlab_files/spectrum/';

disp(dir_modes);
disp(dir_spec);
%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

% Weights of radial direction
nothetar = length(rc);
weight_thetar = zeros(nr,1);
weight_thetar(1) = pi*( rc(1) + (rc(2)-rc(1))/2)^2 - pi*(rc(1))^2;

for i=2:nothetar-1
    weight_thetar(i) = pi*( rc(i) + (rc(i+1)-rc(i))/2 )^2 - pi*( rc(i) - (rc(i)-rc(i-1))/2 )^2;
end

weight_thetar(nothetar) = pi*rc(end)^2 - pi*( rc(end) - (rc(end)-rc(end-1))/2 )^2;
%% Importing Karu's datafiles

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_TKE.dat';
LK_TKE   = importdata(filename);
LK_TKE_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_TKE(:,1))); %#ok<*NCOMMA>
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

%% Loading the eigenvalues for SPOD modes

eigenspectra_allm = zeros(Nf_sampled, Nblk_sampled, size(mode_sampled,1), size(x_sampled,1));
for x_d = 1:size(x_sampled,1)
    for i_mode = 1:size(mode_sampled,1)
        filename = strcat(dir_spec, 'spectrum_x_D_', int2str(x_sampled(x_d,1)), '.mat');
        disp(filename);
        load(filename);
        eigenspectra_allm(:,:,i_mode,x_d) = eigvalue(1:Nf_sampled,1:Nblk_sampled,i_mode);
    end
end


%% Renormalizing the eigenmodes by multiplying with their respective sqrt(lambda)

for x_d = 1:size(x_sampled,1)
    for i_mode = 1:size(mode_sampled,1)
        for f_mode = 1:Nf_sampled
            for Nb = 1:Nblk_sampled
                u_eigenmode_allm(:,Nb,f_mode,i_mode,x_d) = sqrt(eigenspectra_allm(f_mode,Nb,i_mode,x_d))*u_eigenmode_allm(:,Nb,f_mode,i_mode,x_d);
                v_eigenmode_allm(:,Nb,f_mode,i_mode,x_d) = sqrt(eigenspectra_allm(f_mode,Nb,i_mode,x_d))*v_eigenmode_allm(:,Nb,f_mode,i_mode,x_d);
                w_eigenmode_allm(:,Nb,f_mode,i_mode,x_d) = sqrt(eigenspectra_allm(f_mode,Nb,i_mode,x_d))*w_eigenmode_allm(:,Nb,f_mode,i_mode,x_d);
            end
        end
    end
end

%% Normalizing eigenmodes

for x = 1:size(x_sampled,1)
for m = 1:size(mode_sampled,1)
    for fn = 1:Nf_sampled
        for Nb = 1:Nblk_sampled
             spod_mode_u = u_eigenmode_allm(:,Nb,fn,m,x);
             spod_mode_v = v_eigenmode_allm(:,Nb,fn,m,x);
             spod_mode_w = w_eigenmode_allm(:,Nb,fn,m,x);
%             full_spod_mode = [spod_mode_u; spod_mode_v; spod_mode_w];
%             spod_mode_u = spod_mode_u/norm(full_spod_mode);
%             spod_mode_v = spod_mode_v/norm(full_spod_mode);
%             spod_mode_w = spod_mode_w/norm(full_spod_mode);
             spod_mode_mag = spod_mode_u.*conj(spod_mode_u) + spod_mode_v.*conj(spod_mode_v) + spod_mode_w.*conj(spod_mode_w);
             alpha = dot(spod_mode_mag,weight_thetar);
             u_eigenmode_allm(:,Nb,fn,m,x) = spod_mode_u/sqrt(alpha);
             v_eigenmode_allm(:,Nb,fn,m,x) = spod_mode_v/sqrt(alpha);
             w_eigenmode_allm(:,Nb,fn,m,x) = spod_mode_w/sqrt(alpha);
        end
    end
end
end

%% Calculating Reynolds stresses from all m and f sampled

nf_sampled = 20; nblk_sampled = 20; nm_sampled = 3;
for k = 1:size(x_sampled,1)
    for i = 1:nf_sampled
        for j = 1:nm_sampled
            for l = 1:nblk_sampled
                reystress_uw_mode(:,i,j,l,k) = squeeze(eigenspectra_allm(i,l,j,k)*real((u_eigenmode_allm(:,l,i,j,k).*conj(w_eigenmode_allm(:,l,i,j,k)))));
                reystress_uu_mode(:,i,j,l,k) = squeeze(eigenspectra_allm(i,l,j,k)*real((u_eigenmode_allm(:,l,i,j,k).*conj(u_eigenmode_allm(:,l,i,j,k)))));
                reystress_ww_mode(:,i,j,l,k) = squeeze(eigenspectra_allm(i,l,j,k)*real((w_eigenmode_allm(:,l,i,j,k).*conj(w_eigenmode_allm(:,l,i,j,k)))));
                reystress_vv_mode(:,i,j,l,k) = squeeze(eigenspectra_allm(i,l,j,k)*real((v_eigenmode_allm(:,l,i,j,k).*conj(v_eigenmode_allm(:,l,i,j,k)))));
                tke_mode(:,i,j,l,k)          = reystress_uu_mode(:,i,j,l,k) + reystress_ww_mode(:,i,j,l,k) + reystress_vv_mode(:,i,j,l,k);
%               reystress_uw_mode(:,i,j,l,k) = real((u_eigenmode_allm(:,l,i,j,k).*conj(w_eigenmode_allm(:,l,i,j,k))));
%               reystress_uu_mode(:,i,j,l,k) = real((u_eigenmode_allm(:,l,i,j,k).*conj(u_eigenmode_allm(:,l,i,j,k))));
%               reystress_ww_mode(:,i,j,l,k) = real((w_eigenmode_allm(:,l,i,j,k).*conj(w_eigenmode_allm(:,l,i,j,k))));
%               reystress_vv_mode(:,i,j,l,k) = real((v_eigenmode_allm(:,l,i,j,k).*conj(v_eigenmode_allm(:,l,i,j,k))));
%               tke_mode(:,i,j,l,k)          = reystress_uu_mode(:,i,j,l,k) + reystress_ww_mode(:,i,j,l,k) + reystress_vv_mode(:,i,j,l,k);
            end
        end
    end
end


for k = 1:size(x_sampled,1)  %% only taking the leading order POD mode at each (m,f)
    reystress_uw_combined(:,:,k) = squeeze(sum(squeeze(sum(reystress_uw_mode(:,:,:,1,k),2)),3));
    reystress_ww_combined(:,:,k) = squeeze(sum(squeeze(sum(reystress_ww_mode(:,:,:,1,k),2)),3));
    reystress_uu_combined(:,:,k) = squeeze(sum(squeeze(sum(reystress_uu_mode(:,:,:,1,k),2)),3));
    reystress_vv_combined(:,:,k) = squeeze(sum(squeeze(sum(reystress_vv_mode(:,:,:,1,k),2)),3));
    tke_combined(:,:,k)          = squeeze(sum(squeeze(sum(tke_mode(:,:,:,1,k),2)),3));
    tke_combined_allnblk(:,:,k)  = squeeze(sum(squeeze(sum(squeeze(sum(tke_mode(:,:,:,:,k),4)),2)),3));
    reystress_uw_combined_allnblk(:,:,k) = squeeze(sum(squeeze(sum(squeeze(sum(reystress_uw_mode(:,:,:,:,k),4)),2)),3));
end

total_tke_modes = squeeze(sum(tke_combined_allnblk,2));
reystress_uw_modes = squeeze(sum(reystress_uw_combined_allnblk,2));

%% Importing the actual reynolds stresses 
nr = 354;
ntheta = 256;
dir_in_planes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/reystresses/run_2.0/';

loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
reystress_uw_2d = zeros(nr,ntheta,size(loc_planes,1));
reystress_uu_2d = zeros(nr,ntheta,size(loc_planes,1));
reystress_ww_2d = zeros(nr,ntheta,size(loc_planes,1));
reystress_vv_2d = zeros(nr,ntheta,size(loc_planes,1));

reystress_uw_1d = zeros(nr,size(loc_planes,1));
reystress_uu_1d = zeros(nr,size(loc_planes,1));
reystress_ww_1d = zeros(nr,size(loc_planes,1));
reystress_vv_1d = zeros(nr,size(loc_planes,1));
%% Reading the Reynolds stresses from data

for x_loc_planes = 1:size(loc_planes,1)

    filename = strcat(dir_in_planes, 'reystress_x_D_', int2str(loc_planes(x_loc_planes,1)), '.mat');
    disp(filename);
    load(filename);
    
    reystress_uw_2d(:,:,x_loc_planes) =  reystress_uw_av;
    reystress_ww_2d(:,:,x_loc_planes) =  reystress_ww_av;
    reystress_uu_2d(:,:,x_loc_planes) =  reystress_uu_av;
    reystress_vv_2d(:,:,x_loc_planes) =  reystress_vv_av;
    
    reystress_uw_1d(:,x_loc_planes) =  reystress_uw_av_th;
    reystress_ww_1d(:,x_loc_planes) =  reystress_ww_av_th;
    reystress_uu_1d(:,x_loc_planes) =  reystress_uu_av_th;
    reystress_vv_1d(:,x_loc_planes) =  reystress_vv_av_th; %#ok<*SAGROW>

end

total_averaged_tke = reystress_uu_1d + reystress_vv_1d + reystress_ww_1d;
total_averaged_tke_2d = reystress_uu_2d + reystress_vv_2d + reystress_ww_2d;
total_averaged_tke_1d = squeeze(mean(total_averaged_tke_2d,2));


for x_loc_planes = 1:size(loc_planes,1)
    total_integrated_tke(x_loc_planes,1) = (2*pi/256)*trapz(trapz(rc,rc.*total_averaged_tke_2d(:,:,x_loc_planes),1));
end

%% Calculating the area integrated TKE at each location

% for x_loc_planes = 1:size(loc_planes,1)
%     total_integrated_tke(x_loc_planes,1) = (2*pi/256)*trapz(trapz(rc,rc.*total_averaged_tke_2d(:,:,x_loc_planes),1));
% end
% 
% area_integrated_tke = total_integrated_tke;
% x_sample_closest = LK_mean_loc_planes(:,1);
% 
% % loglog(TKE_centerline_loc_planes(:,1), total_integrated_tke, 'ro');
% % hold on;
% % loglog(x_sample_closest(1:end), squeeze(eigenspectra_allm(1,1,3,1:end)), 'ko', 'MarkerSize',7);
% 
% log_tke = log(area_integrated_tke(8:end-2,1));
% log_x_sampled_tke = log(x_sample_closest(8:end-2,1));
% [coeffs_tke, S_tke] = polyfit(log_x_sampled_tke, log_tke, 1);
% y_tke_fitted = polyval(coeffs_tke, log_x_sampled_tke);
% CI_tke = polyparci(coeffs_tke, S_tke, 0.99);

%% Self-similarity of the Reynolds stresses for m modes

% nr = 354; ntheta = 256;
% figure; 
% hold on;
% C = {'k','b','r','g','m','c','k--','b--','r--','g--','m--','c--','k-.','b-.','r-.','g-.','m-.','c-.'}; % Cell array of colros.
% count = 1;
% 
% Legend = cell(7,1);
% Legend{1} = 'x/D = 70';
% Legend{2} = 'x/D = 75';
% Legend{3} = 'x/D = 80';
% Legend{4} = 'x/D = 85';
% Legend{5} = 'x/D = 90';
% Legend{6} = 'x/D = 95';
% Legend{7} = 'x/D = 100';
% Legend{8} = 'x/D = 40';
% Legend{9} = 'x/D = 45';
% Legend{10} = 'x/D = 50';

% for i = 14:20
%     disp(i);
%     plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_combined(:,2,i)/max(abs(reystress_uw_combined(:,2,i))), C{count}, 'Linewidth',2);
%     count = count + 1;
% end

% ylim([0 1.2])
% xlim([0 3]);
% 
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/max($<u_{x}u_{r}>_{r}$) ($m=1$)','interpreter','latex','fontsize',15);
% %hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 10;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

%set(gcf, 'PaperPositionMode', 'auto');
%print(gcf,strcat('unnormalized_uw_x_D_70_100_m1_lk', '.png'),'-dpng','-r600');  
%print(gcf,strcat('unnormalized_uw_x_D_70_100_m1_lk', '.eps'),'-depsc','-r600');

% close;

%% Scaling for Reynolds stresses in streamwise direction

% figure; 
% hold on;
% C = {'k','b','r','g','m','c','k--','b--','r--','g--','m--','c--','k-.','b-.','r-.','g-.','m-.','c-.'}; % Cell array of colros.
% 
% Legend = cell(7,1);
% Legend{1} = 'x/D = 40';
% Legend{2} = 'x/D = 45';
% Legend{3} = 'x/D = 50';
% Legend{4} = 'x/D = 55';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 65';
% Legend{7} = 'x/D = 70';
% Legend{8} = 'x/D = 75';
% Legend{9} = 'x/D = 80';
% Legend{10} = 'x/D = 85';

% Equation 7.10 of Dairay et al. 2015 to find d\delta/dx
% wake_width_mean = LK_mean_loc_planes(8:17,2);
% loc = LK_mean_loc_planes(8:17,1);
% 
% Scaling till x/D = 50
% log_wake_width_50 = log(wake_width_mean);
% log_loc =  log(loc);
% [coeffs, S] = polyfit(log_loc, log_wake_width_50, 1);
% y_fitted = polyval(coeffs, log_wake_width_50);
% CI = polyparci(coeffs, S, 0.99);

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)

% ddelta_dx = coeffs(1)*wake_width_mean./loc;
% scaling_factor = ud_centerline_loc_planes(8:17,2).*ddelta_dx;
% count = 1;

% for i = 8:17
%     disp(i);
%     plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i)/(ux_centerline_loc_planes(i,2)^2), C{count}, 'Linewidth',2);
%     plot(rc/LK_mean_loc_planes(i,2), -4*reystress_uw_combined(:,3,i)/(scaling_factor(count)), C{count}, 'Linewidth',2);
%     count = count + 1;
% end

% ylim([0 0.3])
% xlim([0 3]);


% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{\infty}U_{d}\frac{dL_{d}}{dx}$','interpreter','latex','fontsize',15);
%hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 10;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

%set(gcf, 'PaperPositionMode', 'auto');
%print(gcf,strcat('udddelta_uw_x_D_40_85_ld_m2', '.png'),'-dpng','-r600');  
%print(gcf,strcat('udddelta_uw_x_D_40_85_ld_m2', '.eps'),'-depsc','-r600');
close;

%% Reconstruction of RS from the modal structure

% dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';
% save(strcat(dirout, 'reystress_in_ranseq_construct_similarity_diff_loc.mat'), ...
%             'f', 'rc', 'reystress_uw_combined', 'reystress_ww_combined', ...
%             'reystress_ww_1d', 'reystress_uw_1d','LK_TKE_loc_planes', ...
%             'LK_mean_loc_planes', 'total_averaged_tke_1d', 'tke_combined');
close all;
i = 8;
figure;
hold on;
% h1 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'm-', 'Linewidth',2);
% h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
% h3 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);
%h5 = plot(rc/LK_TKE_loc_planes(i,2), total_tke_modes(:,1), 'k--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), 4*sum(tke_combined_allnblk(:,2:end),2)+ 2*squeeze(tke_combined_allnblk(:,1)),'k--', 'Linewidth',2);
xlim([0.03 3]);
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);
% hLegend = legend([h1, h2, h3], '$m=1$ contribution', '$m=2$ contribution', 'From simulation');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
%set(gcf, 'PaperPositionMode', 'auto');
%print(gcf,strcat('uw_x_D_80_reconstructed', '.png'),'-dpng','-r600');  
% print(gcf,strcat('uw_x_D_20_reconstructed', '.eps'),'-depsc','-r600');

%% Reconstruction of RS from the modal structure

% dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';
% save(strcat(dirout, 'reystress_in_ranseq_construct_similarity_diff_loc.mat'), ...
%             'f', 'rc', 'reystress_uw_combined', 'reystress_ww_combined', ...
%             'reystress_ww_1d', 'reystress_uw_1d','LK_TKE_loc_planes', 'LK_mean_loc_planes');
close all;
i = 8;
figure;
hold on;
% h1 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined_allnblk(:,2,i), 'm-', 'Linewidth',2);
% h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined_allnblk(:,3,i), 'r-', 'Linewidth',2);
% h3 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined_allnblk(:,1,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);
% h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_modes (:,1), 'k--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -4*sum(reystress_uw_combined_allnblk(:,2:end),2)-2*squeeze(reystress_uw_combined_allnblk(:,1)),'k--', 'Linewidth',2);

% h6 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined_allnblk(:,3,i)-4*reystress_uw_combined_allnblk(:,2,i)-2*reystress_uw_combined_allnblk(:,1,i),'k--', 'Linewidth',2);
% xlim([0 3]);
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);
% hLegend = legend([h1, h2, h3], '$m=1$ contribution', '$m=3$ contribution', 'From simulation');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
%set(gcf, 'PaperPositionMode', 'auto');
%print(gcf,strcat('uw_x_D_80_reconstructed', '.png'),'-dpng','-r600');  
% print(gcf,strcat('uw_x_D_20_reconstructed', '.eps'),'-depsc','-r600');



%% Reynolds stress from the m=2 mode 

% Scaling from the maximum value for reynolds stress

% max_ruw = max(abs(squeeze(reystress_uw_combined(:,2,:))),[],1)';
% loglog(LK_TKE_loc_planes(:,1), max_ruw, 'ko');
% log_max_ruw = log(max_ruw);
% log_x_sampled_max_uw = log(LK_TKE_loc_planes(:,1));
% [coeffs_max_uw, S_max_uw] = polyfit(log_x_sampled_max_uw, log_max_ruw, 1);
% y_max_uw_fitted = polyval(coeffs_max_uw, log_x_sampled_max_uw);
% CI_max_uw = polyparci(coeffs_max_uw, S_max_uw, 0.99);
% 

% figure;
% hold on
% h1  = plot(rc/LK_TKE_loc_planes(8,2),  -reystress_uw_mode(:,1,3,1,8)/TKE_centerline_loc_planes(8,2),'b-','LineWidth',2); 
% h2  = plot(rc/LK_TKE_loc_planes(9,2),  -reystress_uw_mode(:,1,3,1,9)/TKE_centerline_loc_planes(9,2),'r-','LineWidth',2);
% h3  = plot(rc/LK_TKE_loc_planes(10,2), -reystress_uw_mode(:,1,3,1,10)/TKE_centerline_loc_planes(10,2),'k-','LineWidth',2);
% h4  = plot(rc/LK_TKE_loc_planes(11,2), -reystress_uw_mode(:,1,3,1,11)/TKE_centerline_loc_planes(11,2),'c-','LineWidth',2);
% h5  = plot(rc/LK_TKE_loc_planes(12,2), -reystress_uw_mode(:,1,3,1,12)/TKE_centerline_loc_planes(12,2),'m-','LineWidth',2);
% h6  = plot(rc/LK_TKE_loc_planes(13,2), -reystress_uw_mode(:,1,3,1,13)/TKE_centerline_loc_planes(13,2),'g-','LineWidth',2);
% h7  = plot(rc/LK_TKE_loc_planes(14,2), -reystress_uw_mode(:,1,3,1,14)/TKE_centerline_loc_planes(14,2),'y-','LineWidth',2);
% h8  = plot(rc/LK_TKE_loc_planes(15,2), -reystress_uw_mode(:,1,3,1,15)/TKE_centerline_loc_planes(15,2),'y--','LineWidth',2);
% h9  = plot(rc/LK_TKE_loc_planes(16,2), -reystress_uw_mode(:,1,3,1,16)/TKE_centerline_loc_planes(16,2),'g--','LineWidth',2);
% h10 = plot(rc/LK_TKE_loc_planes(17,2), -reystress_uw_mode(:,1,3,1,17)/TKE_centerline_loc_planes(17,2),'m--','LineWidth',2);
% h11 = plot(rc/LK_TKE_loc_planes(18,2), -reystress_uw_mode(:,1,3,1,18)/TKE_centerline_loc_planes(18,2),'c--','LineWidth',2);
% h12 = plot(rc/LK_TKE_loc_planes(19,2), -reystress_uw_mode(:,1,3,1,19)/TKE_centerline_loc_planes(19,2),'k--','LineWidth',2);
% h13 = plot(rc/LK_TKE_loc_planes(20,2), -reystress_uw_mode(:,1,3,1,20)/TKE_centerline_loc_planes(20,2),'r--','LineWidth',2);
% h14 = plot(rc/LK_TKE_loc_planes(21,2), -reystress_uw_mode(:,1,3,1,21)/TKE_centerline_loc_planes(21,2),'b--','LineWidth',2);



% figure;
% hold on
% h1  = plot(rc,  reystress_uu_mode(:,1,3,1,8),'b-','LineWidth',2); 
% h2  = plot(rc,  reystress_uu_mode(:,1,3,1,9),'r-','LineWidth',2);
% h3  = plot(rc, reystress_uu_mode(:,1,3,1,10),'k-','LineWidth',2);
% h4  = plot(rc, reystress_uu_mode(:,1,3,1,11),'c-','LineWidth',2);
% h5  = plot(rc, reystress_uu_mode(:,1,3,1,12),'m-','LineWidth',2);
% h6  = plot(rc, reystress_uu_mode(:,1,3,1,13),'g-','LineWidth',2);
% h7  = plot(rc, reystress_uu_mode(:,1,3,1,14),'y-','LineWidth',2);
% h8  = plot(rc, reystress_uu_mode(:,1,3,1,15),'y--','LineWidth',2);
% h9  = plot(rc, reystress_uu_mode(:,1,3,1,16),'g--','LineWidth',2);
% h10 = plot(rc, reystress_uu_mode(:,1,3,1,17),'m--','LineWidth',2);
% h11 = plot(rc, reystress_uu_mode(:,1,3,1,18),'c--','LineWidth',2);
% h12 = plot(rc, reystress_uu_mode(:,1,3,1,19),'k--','LineWidth',2);
% h13 = plot(rc, reystress_uu_mode(:,1,3,1,20),'r--','LineWidth',2);
% h14 = plot(rc, reystress_uu_mode(:,1,3,1,21),'b--','LineWidth',2);
% 
% 
% %ylim([0 0.5]);
% xlim([0 10]);
% 
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$<u_{x}u_{r}>$($m=2$)','interpreter','latex','fontsize',15);
% % %hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);
% 
% hLegend = legend([h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14], 'x/D = 40', 'x/D = 45','x/D = 50',  ...
%    'x/D = 55', 'x/D = 60', 'x/D = 65', 'x/D = 70', 'x/D = 75', 'x/D = 80', 'x/D = 85', 'x/D = 90', 'x/D = 95', 'x/D = 100', 'x/D = 110');
% 
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 10;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
