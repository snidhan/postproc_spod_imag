%% Written by Sheel Nidhan
% % Date - 27th November, 2018
% Comparing the POD modes for stratified case

clc; clear;
close all;

%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 5400;
stride = 100;
nstart = 2337600;
nend = nstart + (N-1)*stride;
dir2 = '/home/sheel/Work2/projects_data/spod_re5e4/fr2/spod_data/run_1.0/x_D_50/eigenspectrum/';
%% Reading the time file
time = importdata ('/home/sheel/Work/projects/spod_re5e4/post/fr2/time_stamps/time_stamp.txt');

time_spod = time(1:N,2);
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

for i = 1:Nfreq
    freq = sprintf('%04d',i);   
    filename = strcat(dir2, 'eigenvalues_freq_',freq,'.txt');
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

% Drawing a -5/3 line
count = 1;
for i = 1:size(freq,2)
    if (freq(1,i) > 0)
        xkol_plot(count,1) = freq(1,i)^(-5/3);
        freq_plot(count,1) = freq(1,i);
        count = count + 1;
        endlin
end


for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1 =  loglog(freq(i,1:end),smoothdata(eigvalue(1:end,i),'SmoothingFactor',.0001),'LineWidth',2,'Color',grey);
    %h1 = loglog(freq(i,1:end),(eigvalue(1:end,i)),'-o','MarkerSize',5,'MarkerFaceColor',grey);
    grid on;
    hold on;
end

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
hYLabel = ylabel('SPOD mode energy','interpreter','latex','fontsize',15);
hTitle = title('SPOD eigenvalue spectra for first 20 modes','interpreter','latex','fontsize',15);

xlim([0.02 10]);
ylim([0 10^-5]);

%loglog(0.01.*freq_plot(5:end),xkol_plot(5:end)/10^6,'--','LineWidth',2,'Color','r');


%% Saving images

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'spod_spectra_fr2_nfreq512_novl384_2110100_2378200_x_D_100.png','-dpng','-r300');  
% print(gcf,'spod_spectra_fr2_nfreq512_novl384_2110100_2378200_x_D_100.eps','-depsc','-r300');
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PUT IN A SEPARATE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
% Cumulative percentage contribution
%% Percentage contribution of modes

percent_contri = zeros(Nfreq,Nblk);
for i = 1:Nfreq
    s = sum(eigvalue(i,:));
    percent_contri(i,:) = 100*(eigvalue(i,:)/s); 
end

percent_cumulative = zeros(Nfreq,Nblk);
for i = 1:Nfreq
    for j = 1:Nblk
        percent_cumulative(i,j) = sum(percent_contri(i,1:j)); 
    end
end

% 
%% Summing up energy of modes from all frequencies

% Calculating mean frequency (approximation)
freq_mean = zeros(size(freq,2),1);
freq_t = freq';

for i = 1:size(freq,2)
    freq_mean(i,1) = mean(freq(:,i));
end

integrated_eigenvalues = zeros(Nblk,1);

for i = 1:Nblk              % Contribution of modes integrated over all frequencies
    integrated_eigenvalues(i,1) = trapz(freq_t(:,i),eigvalue(:,i));
end

percent_contri_integrated = zeros(Nblk,1);
percent_cumulative_integrated = zeros(Nblk,1);

for i = 1:Nblk
    percent_contri_integrated(i,1) = 100*(integrated_eigenvalues(i,1))/(sum(integrated_eigenvalues(1:end,1)));
end

for i = 1:Nblk
    percent_cumulative_integrated(i,1) = sum(percent_contri_integrated(1:i,1));
end


% figure;
% 
% h1 = plot(linspace(1,Nblk,Nblk)', percent_cumulative_integrated, 'linewidth', 2, 'color', 'k');
% 
% hXLabel = xlabel('Number of modes','interpreter','latex','fontsize',15);
% hYLabel = ylabel('Cumulative streamwise TKE','interpreter','latex','fontsize',15);
% hTitle = title('Cumulative energy content integrated over $f$ for $m=1$','interpreter','latex','fontsize',15);
% 
% 
% %% Saving Plot
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'cumulative_energy_m0.png','-dpng','-r300');  
% 
% print(gcf,'cumulative_energy_mo.eps','-depsc','-r300');
% 
