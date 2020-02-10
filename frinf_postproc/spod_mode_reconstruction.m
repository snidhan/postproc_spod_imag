%% Written by Sheel Nidhan
%  Used to project the eigenmodes obtained on the actual flow snapshots and see the evolution of coherent structures

clear; clc; close all;
%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;
x = 10;

dir_modes = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/x_D_', int2str(x), '/eigenmodes/');
dir_spectrum = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/matlab_files/spectrum/';
dir_vel   = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');


%% Reading the time file
        
time = importdata ('/home/sheel/Work/projects/spod_re5e4/post/frinf/time_stamps/time_stamp_1892600_2613200_uniform.txt');
time_spod = time(1:N,1);
dt   = time_spod(2:end,1) - time_spod(1:end-1,1);
%% Read the data files of velocity and subtract the mean

nr = 356; ntheta = 258;

% Reading the data files of radial velocity
var1 = 'up';
var2 = 'vp';
var3 = 'wp';

dir = dir_vel;
u = zeros(nr,ntheta,N);
for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var1, '/', var1, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            u(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

disp('centering the velocities');
u_centered = 0.5*(u(2:nr-1,2:ntheta-1,:) + u(1:nr-2,2:ntheta-1,:));
u_mean = squeeze(mean(u_centered,3));
u_fluc = u_centered - u_mean;
u_fluc_sampled = u_fluc(:,:,1:Nfreq); %#ok<*NASGU>
clear u_fluc; clear u_centered; clear u;
% Reading the files of azimuthal velocity

v = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var2, '/', var2, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);              
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            v(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

disp('centering the velocities');
v_centered = 0.5*(v(2:nr-1,2:ntheta-1,:) + v(2:nr-1,1:ntheta-2,:));
v_mean = squeeze(mean(v_centered,3));
v_fluc = v_centered - v_mean;
v_fluc_sampled = v_fluc(:,:,1:Nfreq);
clear v_fluc; clear v_centered; clear v;

% Reading the data files of streamwise velocity

w = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var3, '/', var3, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            w(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

disp('centering the velocities');
w_centered = w(2:nr-1,2:ntheta-1,:);
w_mean = squeeze(mean(w_centered,3));
w_fluc = w_centered - w_mean;
w_fluc_sampled = w_fluc(:,:,1:Nfreq);
clear w_fluc; clear w_centered; clear w;


%% Reading SPOD modes files

Nblk = floor((N-Novlp)/(Nfreq-Novlp));
nr  = 354;
numvar = 3;
Nrows = numvar*nr*Nblk;
Nrows_permode = numvar*nr;
Nf_sampled = 512;
Nblk_sampled = 27;
mode = [0; 1; 2];
eigmodes = zeros(Nrows, Nf_sampled, size(mode,1));

% for i = 1:Nf_sampled
%     for j = 1:size(mode,1)
%         freq = sprintf('%04d',i);   
%         m = sprintf('%03d',mode(j,1));        
%         filename = strcat(dir_modes, 'eigenmode_freq_',freq,'_',m,'.mod');
%         disp(filename);
%         A = importdata(filename);
%         eigmodes(:,i,j) = A(:,1) + sqrt(-1).*A(:,2);
%     end
% end
% 
%save ('eigmodes_m012_nf512_x_D_40.mat', 'eigmodes', 'mode', 'Nf_sampled');
load('eigmodes_m012_nf512_x_D_40.mat');
%% Rearranging the eigmodes for visualization

for i = 1:Nf_sampled
    disp(i);
    for j = 1:size(mode,1)
        eigenmodes_arranged(:,:,i,j) = reshape(eigmodes(:,i,j), [Nrows_permode, Nblk]); %#ok<*AGROW>
    end
end
eigenmodes_proper_ranked = flip(eigenmodes_arranged, 2);
clear eigenmodes_arranged;
disp('Ranked the eigenmodes');

for i = 1:Nf_sampled
    disp(i);
    for j = 1:Nblk_sampled
        for k = 1:size(mode,1)
            eigenmodes_separated(:,:,j,i,k) = reshape(eigenmodes_proper_ranked(:,j,i,k), [nr, numvar]);
        end
    end
end

disp('Separated the eigenmodes')
clear eigenmodes_proper_ranked;
clear eigmodes;

u_eigenmode = squeeze(eigenmodes_separated(:,1,:,:,:));
v_eigenmode = squeeze(eigenmodes_separated(:,2,:,:,:));
w_eigenmode = squeeze(eigenmodes_separated(:,3,:,:,:));


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

% Weights in azimuthal direction
weight_theta = (2*pi/ntheta)*ones(ntheta,1);   % Check once again SHEEL NIDHAN

if numvar == 4
    weight_rtheta = weight_thetar*weight_theta';
    weight_rtheta_column = weight_rtheta(:);
elseif numvar == 3
    weight_rtheta = weight_thetar;
    weight_rtheta_column = weight_rtheta(:);
end

%% Loading the matlabfiles for SPOD modes

mode_sampled = [0 ;1; 2];
eigenspectra_allm = zeros(Nf_sampled, Nblk_sampled, size(mode_sampled,1));
for i_mode = 1:size(mode_sampled,1)
        filename = strcat(dir_spectrum, 'spectrum_x_D_', int2str(x), '.mat');
        disp(filename);
        load(filename);
        eigenspectra_allm(:,:,i_mode) = eigvalue(1:Nf_sampled,1:Nblk_sampled,i_mode);
end

%% Normalizing eigenmodes as per Moin and Moser eqn. 3.3

mode  = [0; 1; 2];
for m = 1:size(mode,1)
    for fn = 1:Nf_sampled
        for Nb = 1:Nblk_sampled
             spod_mode_u = u_eigenmode(:,Nb,fn,m);
             spod_mode_v = v_eigenmode(:,Nb,fn,m);
             spod_mode_w = w_eigenmode(:,Nb,fn,m);
             spod_mode_mag = spod_mode_u.*conj(spod_mode_u) + spod_mode_v.*conj(spod_mode_v) + spod_mode_w.*conj(spod_mode_w);
             alpha = dot(weight_thetar,spod_mode_mag);
             u_eigenmode(:,Nb,fn,m) = spod_mode_u/sqrt(alpha);
             v_eigenmode(:,Nb,fn,m) = spod_mode_v/sqrt(alpha);
             w_eigenmode(:,Nb,fn,m) = spod_mode_w/sqrt(alpha);
        end
    end
end

%% Fourier decomposition of the velocity fields in azimuthal and time 

nr = 354;
ntheta = 256;
for i = 1:Nf_sampled
    for j = 1:nr
            u_fluc_fft(j,:,i) = (1/ntheta)*fft(u_fluc_sampled(j,:,i)); %#ok<*SAGROW> %% Azimuthal FFT at each radial location
            v_fluc_fft(j,:,i) = (1/ntheta)*fft(v_fluc_sampled(j,:,i)); %% Azimuthal FFT at each radial location
            w_fluc_fft(j,:,i) = (1/ntheta)*fft(w_fluc_sampled(j,:,i)); %% Azimuthal FFT at each radial location
    end
end

window = zeros(Nfreq,1);
for i = 1:Nfreq
  window(i,1) = 0.54 - 0.46*cos((2*pi*(i-1))/(Nfreq-1));
end

winweight = 1.587845072201937; 

for j = 1:ntheta
    disp(j);
    for i = 1:nr
        u_fluc_fft2(i,j,:) = (winweight/Nf_sampled)*fft(window.*squeeze(u_fluc_fft(i,j,:)));  %% Temporal FFT at each grid point
        v_fluc_fft2(i,j,:) = (winweight/Nf_sampled)*fft(window.*squeeze(v_fluc_fft(i,j,:)));  %% Temporal FFT at each grid point
        w_fluc_fft2(i,j,:) = (winweight/Nf_sampled)*fft(window.*squeeze(w_fluc_fft(i,j,:)));  %% Temporal FFT at each grid point
    end
end


%% Calculating the coefficients 

c = zeros(Nblk_sampled, size(mode,1), Nf_sampled);

for i = 1:Nblk_sampled
    for j = 1:size(mode,1)
        for k = 1:Nf_sampled
            product_u = u_fluc_sampled(:,j,k).*conj(u_eigenmode(:,i,k,j));
            product_v = v_fluc_sampled(:,j,k).*conj(v_eigenmode(:,i,k,j));
            product_w = w_fluc_sampled(:,j,k).*conj(w_eigenmode(:,i,k,j));
            total_product = product_u + product_v + product_w;
            c(i,j,k) = dot(weight_thetar, total_product);
        end
    end
end

%% Reconstruction of the velocity back  from leading SPOD modes 

u_vel_reconstructed = zeros(nr, Nblk_sampled, size(mode,1), Nf_sampled);
v_vel_reconstructed = zeros(nr, Nblk_sampled, size(mode,1), Nf_sampled);
w_vel_reconstructed = zeros(nr, Nblk_sampled, size(mode,1), Nf_sampled);

for i = 1:Nblk_sampled
    for j = 1:size(mode,1)
        for k = 1:Nf_sampled
            u_vel_reconstructed(:,i,j,k) = c(i,j,k)*u_eigenmode(:,i,k,j);
            v_vel_reconstructed(:,i,j,k) = c(i,j,k)*v_eigenmode(:,i,k,j);
            w_vel_reconstructed(:,i,j,k) = c(i,j,k)*w_eigenmode(:,i,k,j);
        end
    end
end

%% Inverse fourier FFT in time to bring back the data to time domain

for i = 1:nr
    for j = 1:Nblk_sampled
        for k = 1:size(mode,1)
            u_vel_reconstructed(i,j,k,:) = Nf_sampled*ifft(u_vel_reconstructed(i,j,k,:));
            v_vel_reconstructed(i,j,k,:) = Nf_sampled*ifft(v_vel_reconstructed(i,j,k,:));
            w_vel_reconstructed(i,j,k,:) = Nf_sampled*ifft(w_vel_reconstructed(i,j,k,:));
        end
    end
end

%% Plotting certain quantities

time_sampled = time_spod(1:Nf_sampled,1);
[NR T] = meshgrid(time_sampled, rc);

figure;
h1 = contourf(NR, T, abs(squeeze(v_vel_reconstructed(:,1,2,:))), 'Linestyle', 'none');
colorbar;
figure;
h2 = contourf(NR, T, abs(squeeze(v_fluc_fft(:,2,:))), 'Linestyle', 'none');
colorbar;