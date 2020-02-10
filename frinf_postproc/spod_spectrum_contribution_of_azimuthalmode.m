%% Written by Sheel Nidhan
%  Equation 4.1 of Johansson and George 2006b
function [FREQ, MODE, eigvalue, mode, eigenvalue_fraction_contri] = spod_spectrum_contribution_of_azimuthalmode(x)

close all;

%% SPOD Parameters
Nfreq = 512;
Novlp = 256;
N     = 7200;
mode  = linspace(0,11,12)';
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
dir2 = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/x_D_', int2str(x), '/eigenspectrum/');
mat_file = strcat('spectrum_x_D_', int2str(x), '.mat');
disp(dir2);

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

eigvalue = zeros(Nfreq,Nblk,size(mode,1));

for k = 1:size(mode,1)
    for i = 1:Nfreq
        m = sprintf('%03d',mode(k,1));
        freq = sprintf('%04d',i);   
        filename = strcat(dir2, 'eigenvalues_freq_',freq,'_',m,'.txt');
        disp(filename);
        A = importdata(filename);
        eigvalue(i,:,k) = A(:,1);
    end
end

for k = 1:size(mode,1)
    for i = 1:Nfreq
        eigvalue(i,:,k) = sort(eigvalue(i,:,k), 'descend');
    end
end

%% Integrating over all frequencies for specific azimuthal modes

freq_t = f(1:Nfreq/2)';
integrated_eigenvalues = zeros(Nblk,1);

for j = 1:size(mode,1)
    for i = 1:Nblk              % Contribution of modes integrated over all frequencies
        integrated_eigenvalues(i,j) = sum(eigvalue(:,i,j)); 
    end
end

%% Importing the actual reynolds stresses 
nr = 354;
ntheta = 256;
dir_in_planes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/reystresses/run_2.0/';

reystress_uu_2d = zeros(nr,ntheta,1);
reystress_ww_2d = zeros(nr,ntheta,1);
reystress_vv_2d = zeros(nr,ntheta,1);

%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

%% Reading the Reynolds stresses from data

filename = strcat(dir_in_planes, 'reystress_x_D_', int2str(x), '.mat');
disp(filename);
load(filename);
    
reystress_ww_2d(:,:,1) =  reystress_ww_av;
reystress_uu_2d(:,:,1) =  reystress_uu_av;
reystress_vv_2d(:,:,1) =  reystress_vv_av;

total_averaged_tke_2d = reystress_uu_2d + reystress_vv_2d + reystress_ww_2d;
total_integrated_tke(1,1) = (2*pi/256)*trapz(trapz(rc,rc.*total_averaged_tke_2d(:,:,1),1));

%% Each SPOD mode summed over all azimuthal modes

totE = sum(sum(integrated_eigenvalues));

for j = 1:size(mode,1)
    for i = 1:Nblk
        eigenvalue_fraction_contri(i,j) = integrated_eigenvalues(i,j)/(total_integrated_tke)*100; %#ok<*SAGROW,*AGROW>
    end
end

eigenvalue_fraction_contri(:,2:end) = 2*eigenvalue_fraction_contri(:,2:end);
%% Cummulative energy 

for j = 1:size(mode,1)
    for i = 1:Nblk
        eigenvalue_cum(i,j) = sum(eigenvalue_fraction_contri(1:i,j));
    end
end

%% Save mat file

%save(mat_file);

%% Plotting the eigenspectrum integrated over frequency as function of m
% 
% figure;
% h1 = bar(mode, eigenvalue_fraction_contri(1:3,:)','grouped');
% grid on;
% ylim([0 35]);
% 
% labels = {'POD Mode 1','POD Mode 2','POD Mode 3'};
% hLegend = legend(labels,'Location','NorthEast');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 10;
% hLegend.FontWeight = 'bold';
% 
% hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$\xi^{(n)}$','interpreter','latex','fontsize',15);
% hTitle  = title(strcat('Eigenspectrum integrated over frequency at $x/D=$', int2str(x)),'interpreter','latex','fontsize',15);
% 
% set(gcf, 'PaperPositionMode', 'auto');  
% print(gcf,strcat('integrated_eigenspectrm_over_f_x_D_',sprintf('%03d',x),'_spod.png'),'-dpng','-r600');
% print(gcf,strcat('integrated_eigenspectrm_over_f_x_D_',sprintf('%03d',x),'_spod.eps'),'-depsc','-r600');

%% Plotting the eigenspectrum integrated over frequency as function of m

% figure;
freq_sampled = f(1:50);
[FREQ MODE] = meshgrid(mode(1:6), freq_sampled);
% h1 = contourf(FREQ, MODE, (squeeze(eigvalue(1:50,1,1:6))), 'Linestyle', 'none'); %#ok<*NASGU>
% xticks([0 1 2 3 4 5]);
% colormap jet;
% colorbar

% hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$St$','interpreter','latex','fontsize',15);

% set(gcf, 'PaperPositionMode', 'auto');  
% print(gcf,strcat('contourf_allmf_x_D_',sprintf('%03d',x),'_spod.png'),'-dpng','-r600');
% print(gcf,strcat('contourf_allmf_f_x_D_',sprintf('%03d',x),'_spod.eps'),'-depsc','-r600');