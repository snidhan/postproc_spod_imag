%% Written by Sheel Nidhan
clear; clc; close all;
%% SPOD Parameters
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
dir_spec  = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/matlab_files/spectrum/';

disp(dir_modes);
disp(dir_spec);
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

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Defect_centerline.dat';
ud_Centerline   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ud_Centerline(:,1)));
    ud_centerline_loc_planes(i,2) = ud_Centerline(idx,2);
    ud_centerline_loc_planes(i,1) = ud_Centerline(idx,1);
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

%% Normalizing eigenmodes

for x = 1:size(x_sampled,1)
for m = 1:size(mode_sampled,1)
    for fn = 1:Nf_sampled
        for Nb = 1:Nblk_sampled
             spod_mode_u = u_eigenmode_allm(:,Nb,fn,m,x);
             spod_mode_v = v_eigenmode_allm(:,Nb,fn,m,x);
             spod_mode_w = w_eigenmode_allm(:,Nb,fn,m,x);
             spod_mode_mag = spod_mode_u.*conj(spod_mode_u) + spod_mode_v.*conj(spod_mode_v) ... 
                             + spod_mode_w.*conj(spod_mode_w);
             spod_mode_mag_radial = spod_mode_mag.*rc;
%            spod_mode_mag = norm(spod_mode);
%            alpha = spod_mode_mag;
             alpha = trapz(rc, spod_mode_mag_radial);
             u_eigenmode_allm(:,Nb,fn,m,x) = u_eigenmode_allm(:,Nb,fn,m,x)/sqrt(alpha);
             v_eigenmode_allm(:,Nb,fn,m,x) = v_eigenmode_allm(:,Nb,fn,m,x)/sqrt(alpha);
             w_eigenmode_allm(:,Nb,fn,m,x) = w_eigenmode_allm(:,Nb,fn,m,x)/sqrt(alpha);
        end
    end
end
end

% for x = 1:size(x_sampled,1)
% for m = 1:size(mode_sampled,1)
%     for fn = 1:Nf_sampled
%         for Nb = 1:Nblk_sampled
%              spod_mode = v_eigenmode_allm(:,Nb,fn,m,x);
%              spod_mode_mag = spod_mode.*conj(spod_mode);
%              spod_mode_mag_radial = spod_mode_mag.*rc;
% %            spod_mode_mag = norm(spod_mode);
% %            alpha = spod_mode_mag;
%              alpha = trapz(rc, spod_mode_mag_radial);
%              v_eigenmode_allm(:,Nb,fn,m,x) = v_eigenmode_allm(:,Nb,fn,m,x)/sqrt(alpha);
%         end
%     end
% end
% end
% 
% for x = 1:size(x_sampled,1)
% for m = 1:size(mode_sampled,1)
%     for fn = 1:Nf_sampled
%         for Nb = 1:Nblk_sampled
%              spod_mode = w_eigenmode_allm(:,Nb,fn,m,x);
%              spod_mode_mag = spod_mode.*conj(spod_mode);
%              spod_mode_mag_radial = spod_mode_mag.*rc;
% %            spod_mode_mag = norm(spod_mode);
% %            alpha = spod_mode_mag;
%              alpha = trapz(rc, spod_mode_mag_radial);
%              w_eigenmode_allm(:,Nb,fn,m,x) = w_eigenmode_allm(:,Nb,fn,m,x)/sqrt(alpha);
%         end
%     end
% end
% end

%% Calculating Reynolds stresses of specific modes

for i = 1:size(x_sampled,1)
    reystress_uw_mode2_st0(:,i)    = squeeze(eigenspectra_allm(1,1,3,i)*real((u_eigenmode_allm(:,1,1,3,i).*conj(w_eigenmode_allm(:,1,1,3,i)))));
    reystress_uw_mode1_st0136(:,i) = squeeze(eigenspectra_allm(6,1,2,i)*real((u_eigenmode_allm(:,1,6,2,i).*conj(w_eigenmode_allm(:,1,6,2,i)))));
    reystress_uw_mode0_st0(:,i)    = squeeze(eigenspectra_allm(1,1,1,i)*real((u_eigenmode_allm(:,1,1,1,i).*conj(w_eigenmode_allm(:,1,1,1,i)))));
    reystress_uu_mode2_st0(:,i)    = squeeze(eigenspectra_allm(1,1,3,i)*real((u_eigenmode_allm(:,1,1,3,i).*conj(u_eigenmode_allm(:,1,1,3,i)))));
    reystress_uu_mode1_st0136(:,i) = squeeze(eigenspectra_allm(6,1,2,i)*real((u_eigenmode_allm(:,1,6,2,i).*conj(u_eigenmode_allm(:,1,6,2,i)))));
    reystress_uu_mode0_st0(:,i)    = squeeze(eigenspectra_allm(1,1,1,i)*real((u_eigenmode_allm(:,1,1,1,i).*conj(u_eigenmode_allm(:,1,1,1,i)))));
    reystress_ww_mode2_st0(:,i)    = squeeze(eigenspectra_allm(1,1,3,i)*real((w_eigenmode_allm(:,1,1,3,i).*conj(w_eigenmode_allm(:,1,1,3,i)))));
    reystress_ww_mode1_st0136(:,i) = squeeze(eigenspectra_allm(6,1,2,i)*real((w_eigenmode_allm(:,1,6,2,i).*conj(w_eigenmode_allm(:,1,6,2,i)))));
    reystress_ww_mode0_st0(:,i)    = squeeze(eigenspectra_allm(1,1,1,i)*real((w_eigenmode_allm(:,1,1,1,i).*conj(w_eigenmode_allm(:,1,1,1,i)))));
    reystress_vv_mode2_st0(:,i)    = squeeze(eigenspectra_allm(1,1,3,i)*real((v_eigenmode_allm(:,1,1,3,i).*conj(v_eigenmode_allm(:,1,1,3,i)))));
    reystress_vv_mode1_st0136(:,i) = squeeze(eigenspectra_allm(6,1,2,i)*real((v_eigenmode_allm(:,1,6,2,i).*conj(v_eigenmode_allm(:,1,6,2,i)))));
    reystress_vv_mode0_st0(:,i)    = squeeze(eigenspectra_allm(1,1,1,i)*real((v_eigenmode_allm(:,1,1,1,i).*conj(v_eigenmode_allm(:,1,1,1,i)))));

end


%% Calculating Reynolds stresses from all m and f sampled

clear reystress_uw;

for k = 1:size(x_sampled,1)
    for i = 1:6
        for j = 1:5
            reystress_ww(:,i,j,k) = squeeze(eigenspectra_allm(i,1,j,k)*real((u_eigenmode_allm(:,1,i,j,k).*conj(w_eigenmode_allm(:,1,i,j,k)))));
        end
    end
end

for k = 1:size(x_sampled,1)
    reystress_ww_combined(:,k) = sum(squeeze(sum((reystress_ww(:,:,:,k)),2)),2);
end

hold on;

plot(rc, reystress_uw_combined(:,10),'k-')
plot(rc, reystress_uw_1d(1:end-2,10),'r-')
%% Importing the actual reynolds stresses 
nr = 356;
ntheta = 258;
dir_in_planes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/reystresses/';

loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
reystress_uw_2d = zeros(nr,ntheta,size(loc_planes,1));
reystress_uu_2d = zeros(nr,ntheta,size(loc_planes,1));
reystress_ww_2d = zeros(nr,ntheta,size(loc_planes,1));

reystress_uw_1d = zeros(nr,size(loc_planes,1));
reystress_uu_1d = zeros(nr,size(loc_planes,1));
reystress_ww_1d = zeros(nr,size(loc_planes,1));

%% Reading the Reynolds stresses
for x_loc_planes = 1:size(loc_planes,1)

    filename = strcat(dir_in_planes, 'reystress_x_D_', int2str(loc_planes(x_loc_planes,1)), '.mat');
    disp(filename);
    load(filename);
    
    reystress_uw_2d(:,:,x_loc_planes) =  reystress_uw_av;
    reystress_ww_2d(:,:,x_loc_planes) =  reystress_ww_av;
    reystress_uu_2d(:,:,x_loc_planes) =  reystress_uu_av;
    
    reystress_uw_1d(:,x_loc_planes) =  reystress_uw_av_th;
    reystress_ww_1d(:,x_loc_planes) =  reystress_ww_av_th;
    reystress_uu_1d(:,x_loc_planes) =  reystress_uu_av_th;
end



%% Plotting the Reynolds stresses versus actual data

for i = 1:size(x_sampled,1)
    figure;
    hold on;    
    h1 = plot(rc, (reystress_uw_mode2_st0(:,i)), 'r-', 'Linewidth',2);
    h2 = plot(rc, (reystress_uw_mode1_st0136(:,i)), 'k-', 'Linewidth',2);
    h3 = plot(rc, (reystress_uw_mode0_st0(:,i)), 'c-', 'Linewidth',2);
    h4 = plot(rc, (reystress_uw_mode2_st0(:,i) + reystress_ww_mode1_st0136(:,i) + reystress_ww_mode0_st0(:,i)), 'b-', 'Linewidth',2);
    h5 = plot(rc, (reystress_uw_1d(1:nr-2,i)), 'm-', 'Linewidth',2);
    %plot(rc, reystress_ww_mode2_st0(:,4)+reystress_ww_mode1_st0136(:,4)+ reystress_ww_mode0_st0(:,4), 'k-', 'Linewidth',2);
    hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',15);
    hYLabel = ylabel('$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);
    hTitle = title(strcat('Radial variation of ', ' $<u_{x}u_{x}>$', 'at x/D = ', int2str(x_sampled(i,1))),'interpreter','latex','fontsize',15);
    hLegend = legend([ h1, h2, h3, h4, h5], '(1) $m=2$, $St=0$', '(2) $m=1$, $St=0.136$', '(3) $m=0$, $St=0$', '(1) + (2) + (3)', 'Measured');         
    hLegend.Interpreter = 'Latex';
    hLegend.FontSize = 15;
    hLegend.FontWeight = 'bold';
    hLegend.Position = [0 0 1 1];
    %ylim([0 5*10^-4])
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf,strcat('unnormalized_ww_x_D', int2str(x_sampled(i,1)), '.png'),'-dpng','-r600');  
    %print(gcf,strcat('unnormalized_ww_x_D', int2str(x_sampled(i,1)), '.png'),'-depsc','-r600');
    close all;
end