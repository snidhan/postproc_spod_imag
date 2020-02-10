function [eigvalue, f] = spod_spectrum_plot_real_data(x)

%% Written by Sheel Nidhan
% % Date - February 8, 2020
% Comparing the POD modes for stratified case

close all;

%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7000;
stride = 100;
nstart = 2329600;
nend = nstart + (N-1)*stride;
dir2 = strcat('/home/sheel/Work2/projects_data/spod_re5e4/fr2/spod_data/run_3.0/x_D_', num2str(x), '/eigenspectrum/');
%% dt for calculating frequency

dt = 0.0905441280000332;

%% Fixing the frequency of SPOD spectrum

f = (0:Nfreq-1)/dt(1)/Nfreq;

if mod(Nfreq,2) == 0
    f(Nfreq/2 + 1:end) = f(Nfreq/2 + 1:end)-1/dt(1);
else
    f((Nfreq+1)/2 + 1:end) = f((Nfreq+1)/2 + 1:end) - 1/dt(1);
end

f = f';
%% Setting Nblk for the calculation

Nblk = floor((N-Novlp)/(Nfreq-Novlp));

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

for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1 =  loglog(f(1:end)',smoothdata(eigvalue(1:end,i),'SmoothingFactor',.0001),'LineWidth',2,'Color',grey);
    %h1 = loglog(freq(i,1:end),(eigvalue(1:end,i)),'-o','MarkerSize',5,'MarkerFaceColor',grey);
    grid on;
    hold on;
end

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
hYLabel = ylabel('SPOD mode energy','interpreter','latex','fontsize',15);
hTitle = title('SPOD eigenvalue spectra for first 20 modes','interpreter','latex','fontsize',15);

xlim([0.02 3]);
ylim([0 10^-3]);


%% Saving images

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,'spod_spectra_fr2_nfreq512_novl384_2110100_2378200_x_D_100.png','-dpng','-r300');  
% print(gcf,'spod_spectra_fr2_nfreq512_novl384_2110100_2378200_x_D_100.eps','-depsc','-r300'); 