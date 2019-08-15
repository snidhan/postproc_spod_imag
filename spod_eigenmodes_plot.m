%% Written by Sheel Nidhan
%  Plotting the SPOD modes for a specific azimuthal wavenumber
function [] = spod_eigenmodes_plot(x, m, Nf_sampled, Nblk_sampled)
%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7000;
mode  = m;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
nr = 354;
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nrows = numvar*nr*Nblk;
Nrows_permode = numvar*nr;

dir2 = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/x_D_', int2str(x), '/eigenmodes/');
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

eigmodes = zeros(Nrows, Nf_sampled);
m = sprintf('%03d',mode);

for i = 1:Nf_sampled
    
    freq = sprintf('%04d',i);   
    filename = strcat(dir2, 'eigenmode_freq_',freq,'_',m,'.mod');
    disp(filename);
    A = importdata(filename);
    eigmodes(:,i) = A(:,1) + sqrt(-1).*A(:,2);
   
end

%% Rearranging the eigmodes for visualization

for i = 1:Nf_sampled
    eigenmodes_arranged(:,:,i) = reshape(eigmodes(:,1), [Nrows_permode, Nblk]);
end

eigenmodes_proper_ranked = flip(eigenmodes_arranged, 2);
clear eigenmodes_arranged;

for i = 1:Nf_sampled
    for j = 1:Nblk_sampled
        eigenmodes_separated(:,:,j,i) = reshape(eigenmodes_proper_ranked(:,j,i), [nr, numvar]);
    end
end

clear eigenmodes_proper_ranked;
clear eigmodes;

u_eigenmode = squeeze(eigenmodes_separated(:,1,:,:));
v_eigenmode = squeeze(eigenmodes_separated(:,2,:,:));
w_eigenmode = squeeze(eigenmodes_separated(:,3,:,:));

matfile = strcat('eigenmodes_x_D_', int2str(x), '_m', int2str(mode), '.mat');
disp(matfile);
save(matfile);
%% Plotting the eigenmodes 

% figure;
% spod_mode = 1;
% wake_width = 2.5;
% h1 = plot(rc./wake_width, real(u_normalize(:,spod_mode,8)), 'r-', 'Linewidth', 2);
% hold on;
% h2 = plot(rc./wake_width, real(v_normalize(:,spod_mode,8)), 'b-', 'Linewidth', 2);
% h3 = plot(rc./wake_width, real(w_normalize(:,spod_mode,8)), 'k-', 'Linewidth', 2);
% 
% hXLabel = xlabel('$r/\L_{Hk}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$\Phi$','interpreter','latex','fontsize',15);
% hTitle = title('SPOD mode 4 for $m=0$ at $x/D = 40$','interpreter','latex','fontsize',15);
% 
% hLegend = legend([h1,h2,h3],'$u_{r}$', '$u_{\theta}$', '$u_{x}$');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

%% Saving the image file

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'spod_mode4_m0_x_D_40.png','-dpng','-r300');  