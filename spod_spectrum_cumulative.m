%% Written by Sheel Nidhan
%  Plotting the cumulative SPOD spectrum integrated over frequency for a specific azimuthal wavenumber

clc; clear;

%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7000;
mode  = 5;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;

dir2 = '/home/sheel/Work/projects/spod_re5e4/post/frinf/spod_data/x_D_60/eigenspectra/';
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

%% Percentage contribution of modes

percent_contri = zeros(Nfreq/2,Nblk);

for i = 1:Nfreq/2
    s = sum(eigvalue(i,:));
    percent_contri(i,:) = 100*(eigvalue(i,:)/s); 
end

percent_cumulative = zeros(Nfreq/2,Nblk);

for i = 1:Nfreq/2
    for j = 1:Nblk
        percent_cumulative(i,j) = sum(percent_contri(i,1:j)); 
    end
end

%% Summing up energy of modes from all frequencies

freq_t = f(1:Nfreq/2)';
integrated_eigenvalues = zeros(Nblk,1);

for i = 1:Nblk              % Contribution of modes integrated over all frequencies
    integrated_eigenvalues(i,1) = trapz(freq_t(:,1),eigvalue(1:Nfreq/2,i));
end

percent_contri_integrated = zeros(Nblk,1);
percent_cumulative_integrated = zeros(Nblk,1);

for i = 1:Nblk
    percent_contri_integrated(i,1) = 100*(integrated_eigenvalues(i,1))/(sum(integrated_eigenvalues(1:end,1)));
end

for i = 1:Nblk
    percent_cumulative_integrated(i,1) = sum(percent_contri_integrated(1:i,1));
end

%% Plotting the cumulative energy for a particular azimuthal m integrated over frequency

close all;
figure;
h1 = plot(linspace(1,Nblk,Nblk)', percent_cumulative_integrated, 's', 'MarkerFaceColor', 'k');
hXLabel = xlabel('Number of SPOD modes','interpreter','latex','fontsize',15);
hYLabel = ylabel('Cumulative TKE','interpreter','latex','fontsize',15);
hTitle = title(strcat('Cumulative energy content integrated over $f$ for $m=$', int2str(mode)), ...
               'interpreter','latex','fontsize',15);

xlim([1 Nblk]);
ylim([40 100]);
%% Saving Plot

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, strcat('spod_cumulative_energy_m', int2str(mode),'.png'), '-dpng', '-r600');  