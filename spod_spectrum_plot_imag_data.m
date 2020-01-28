function [eigvalue, f] = spod_spectrum_plot_imag_data(x, m)
%% Written by Sheel Nidhan
%  Plotting the SPOD spectrum for a specific azimuthal wavenumber
close all;
%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
mode  = m;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
dir2 = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/x_D_', int2str(x), '/eigenspectrum/');
%dir2 = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/';
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

eigvalue = zeros(Nfreq,Nblk);
m = sprintf('%03d',mode);

for i = 1:Nfreq
    freq = sprintf('%04d',i);   
    filename = strcat(dir2, 'eigenvalues_freq_',freq,'_',m,'.txt');
    disp(filename);
    A = importdata(filename);
    eigvalue(i,:) = A(:,1);
   
end

for i  = 1:Nfreq
    eigvalue(i,:) = sort(eigvalue(i,:),'descend');
end
%% Plotting SPOD eigenvalues

Nplot = 20;   % No. of modes to plot

C = repmat(linspace(1,0.1,Nplot).',1,3);

figure;

for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1(i) =  loglog(f',eigvalue(1:end,i), 'LineWidth',2,'Color',grey);
    hold on;
    grid on;
end

ylim([10^-14 1*10^-1]);
xlim([0 2]);

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
hYLabel = ylabel('SPOD mode energy','interpreter','latex','fontsize',15);
%% Saving images

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('spod_spectra_',int2str(mode),'x_D_',int2str(x),'.png'),'-dpng','-r600');  
% print(gcf,strcat('spod_spectra_',int2str(mode),'x_D_',int2str(x),'.eps'),'-depsc','-r600');  

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('spod_spectra_',int2str(mode),'_3dplane_','.png'),'-dpng','-r600');  
% print(gcf,strcat('spod_spectra_',int2str(mode),'_3dplane_','.eps'),'-depsc','-r600');  

