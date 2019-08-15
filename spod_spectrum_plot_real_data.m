

%% Written by Sheel Nidhan
% % Date - 27th November, 2018
% Comparing the first POD modes for different m azimuthal wavenumber

clc; clear;
close all;

%% SPOD Parameters

Nfreq = 512;
Novlp = 384;
N     = 2680;

stride = 100;
nstart = 2110100;
nend = nstart + (N-1)*stride;

dir2 = '/home/sheel/Work/projects/spod_re5e4/post/fr2/spod_data/x_D_50/nfreq512_novlp384_2110100_2378200/eigenvalues/';
%% Reading the time file
time = importdata ('/home/sheel/Work/projects/spod_re5e4/post/fr2/time_stamps/time_stamp_2110100_2378200.txt');

%Two data files missing for densp for Fr = 2 case
time([1229 1232]) = [];

%%

% begin iteration 1892600
% end   iteration 2471900

n1 =  (nstart - 2110100)/100 + 1; 
n2 = n1 + N - 1;

time = time(n1:n2,1);
dt   = time(2:end,1) - time(1:end-1,1);
plot(dt, 'bo');
pause(2);

close;
%%

Nblk = floor((N-Novlp)/(Nfreq-Novlp));

qstart = zeros(Nblk,1);
qend   = zeros(Nblk,1);

for i = 1:Nblk
    qstart(i,1) = (i-1)*(Nfreq-Novlp) + 1;
    qend(i,1)   = qstart(i) + Nfreq - 1;
end
%% Separating time data into blocks

time_blk = zeros(Nblk,Nfreq);

for i = 1:Nblk
    time_blk(i,:) = time(qstart(i,1):qend(i,1),1)';
end
%% Resolving the time series data into frequency

% Calculating the time window for each block

T = zeros(Nblk,1);
for i = 1:Nblk
    T(i,1) = time_blk(i,end) - time_blk(i,1);
end

freq  = zeros(Nblk,Nfreq);

for i = 1:Nblk
    Fs = 1/((2*pi)/Nfreq);  
    f_mapped = (Fs)*(-Nfreq/2:Nfreq/2-1)/Nfreq;
    freq(i,:) = (f_mapped.*(2*pi/T(i,1)))';

end
%% Reading SPOD eigenvalues files

eigvalue = zeros(Nfreq,Nblk);
%m = sprintf('%03d',mode);



for i = 1:Nfreq
    filename = sprintf('%04d',i);   
    filename = strcat(dir2, 'freq_eigenvalues__',filename,'_','000','.txt')
   
    A = importdata(filename);
    eigvalue(i,:) = A(:,1);
   
end

for i  = 1:Nfreq
    eigvalue(i,:) = sort(eigvalue(i,:));
end


for i = 1:Nfreq
    eigvalue(i,:) = flip(eigvalue(i,:));
end

%% Plotting SPOD eigenvalues

Nplot = 17;   % No. of modes to plot

C = repmat(linspace(1,0.1,Nplot).',1,3);

figure;

% Drawing a -5/3 line
count = 1;
for i = 1:size(freq,2)
    if (freq(1,i) > 0)
        xkol_plot(count,1) = freq(1,i)^(-5/3);
        freq_plot(count,1) = freq(1,i);
        count = count + 1;
    end
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
