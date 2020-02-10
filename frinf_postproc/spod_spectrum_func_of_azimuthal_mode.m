%%  Written by Sheel Nidhan
% % The comparison between SPOD modes for different azimuthal numbers

clc; clear all;
%% SPOD Parameters

Nfreq =  512;
Novlp =  256;
N     =  7000;
mode  = linspace(0,10,11)';
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;

dir2 = '/home/sheel/Work/projects/spod_re5e4/post/frinf/spod_data/3dplane_spod/eigenspectra/';
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

%% Saving the workspace to a mat file

%save('eigenspectra_x_D_65_m0_m5.mat');

%% Plotting SPOD eigenvalues for azimuthal wavenumbers

m = linspace(0,4,5)';

figure;

h1 = semilogy(f(1,1:Nfreq/2), eigvalue(1:Nfreq/2,1,1),'linewidth',2,'color','k');
hold on;
h2 = semilogy(f(1,1:Nfreq/2), eigvalue(1:Nfreq/2,1,2), 'linewidth', 2, 'color', 'r');
h3 = semilogy(f(1,1:Nfreq/2), eigvalue(1:Nfreq/2,1,3), 'linewidth', 2, 'color', 'b');
h4 = semilogy(f(1,1:Nfreq/2), eigvalue(1:Nfreq/2,1,4), 'linewidth', 2, 'color', 'g');
h5 = semilogy(f(1,1:Nfreq/2), eigvalue(1:Nfreq/2,1,5), 'linewidth', 2, 'color', 'm');
h6 = semilogy(f(1,1:Nfreq/2), eigvalue(1:Nfreq/2,1,6), 'linewidth', 2, 'color', 'c');

grid on;

xlim([0 0.5]);
ylim([10^-14 10^-2]);

% hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('SPOD mode energy','interpreter','latex','fontsize',15);
% hTitle  = title('$1^{st}$ SPOD mode for different $m$','interpreter','latex','fontsize',15);

hLegend = legend([h1,h2,h3,h4,h5,h6],'$m=0$', '$m=1$','$m=2$','$m=3$','$m=4$', '$m=5$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 20;
hLegend.FontWeight = 'bold';
hLegend.Location = 'southeast';

%% Plotting the SPOd modes

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,'first_spod_mode_m0_5_3dplanes_spod.png','-dpng','-r600');