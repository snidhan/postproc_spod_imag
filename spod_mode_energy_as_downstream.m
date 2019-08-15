%% Written by Sheel Nidhan
%  Fig. 3 of Johannson and George 2006b
%  Evaluating the change in the energy of specific POD modes as a function of x/D 

clear; close all;
%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7000;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;

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

%% Generating the locations at which the eigenvalues are plotted

x_D = linspace(5,100,20)';
x_D(21,1) = 110; x_D(22,1) = 120;

%% Loading the mat files

dir = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/matlab_files/spectrum/';
mf = 11; Nf = 50;                             %% Number of frequencies and modes for saving low-dimensional data
eigvalue_spod1 = zeros(Nf,mf,size(x_D,1));
eigvalue_spod2 = zeros(Nf,mf,size(x_D,1));
eigvalue_spod3 = zeros(Nf,mf,size(x_D,1));
eigvalue_spod4 = zeros(Nf,mf,size(x_D,1));
eigvalue_spod5 = zeros(Nf,mf,size(x_D,1));
for z = 1:size(x_D,1)
    disp(z);
    filename = strcat('spectrum_x_D_', num2str(x_D(z,1)), '.mat');
    fullfile = strcat(dir, filename);
    disp(fullfile);
    load(fullfile);
    eigvalue_spod1(:,:,z) = squeeze(eigvalue(1:50,1,1:mf));
    eigvalue_spod2(:,:,z) = squeeze(eigvalue(1:50,2,1:mf));
    eigvalue_spod3(:,:,z) = squeeze(eigvalue(1:50,3,1:mf));
    eigvalue_spod4(:,:,z) = squeeze(eigvalue(1:50,4,1:mf));
    eigvalue_spod5(:,:,z) = squeeze(eigvalue(1:50,5,1:mf));
end


%% Plotting the variation of eigenvalues as a function of distance

hold on;
h1 = plot(x_D, squeeze(eigvalue_spod1(6,2,:)), 'ko','MarkerSize', 7);
h2 = plot(x_D, squeeze(eigvalue_spod1(1,3,:)), 'bs','MarkerSize', 7);
h3 = plot(x_D, squeeze(eigvalue_spod1(1,1,:)), 'r*','MarkerSize', 7);
h4 = plot(x_D, squeeze(eigvalue_spod1(8,1,:)), 'cd','MarkerSize', 7);
h5 = plot(x_D, squeeze(eigvalue_spod1(3,1,:)), 'g^','MarkerSize', 7);
%% Saving the plot

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda^{1}$','interpreter','latex','fontsize',15);
%hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

hLegend = legend([h1,h2,h3,h4,h5],'$m=1$, $St = 0.136$', '$m=2$, $St = 0$', ...
    '$m=0$, $St = 0$', '$m=0$, $St = 0.19$', '$m=0$, $St = 0.05$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

%ylim([0 8*10^-3]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'spod_mode4_x_D.png','-dpng','-r300');  
print(gcf,'spod_mode4_x_D.eps','-depsc','-r600');