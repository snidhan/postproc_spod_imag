%% Written by Sheel Nidhan
%  Plotting the SPOD modes for a specific azimuthal wavenumber for the 3dplane SPOD
%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7000;
mode  = 2;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
nr = 354;
nplanes = 22;
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
x_sampled = [20; 40; 60; 80; 100; 120];
half_width_tke = [1.9557; 2.5796; 2.9654; 3.24593; 3.4227; 3.6530];
Nrows = numvar*nr*nplanes*Nblk;
Nrows_permode = numvar*nr*nplanes;
dir2 = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/3dplane_spod/eigenmodes/';
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

eigmodes = zeros(Nrows_permode,1);
m = sprintf('%03d',mode);
freq =  sprintf('%04d',1);  
filename = strcat(dir2, 'eigenmode_freq_',freq,'_',m,'_spod26.mod');
A = importdata(filename);
eigmodes(:,1) = A(:,1) + sqrt(-1).*A(:,2);


%% Rearranging the eigmodes for visualization

eigenmodes_arranged(:,:) = reshape(eigmodes(:,1), [nr*numvar, nplanes]);
for i = 1:nplanes
    eigenmodes_separated(:,:,i) = reshape(eigenmodes_arranged(:,i), [nr, numvar]);
end

u_eigenmode = squeeze(eigenmodes_separated(:,1,:));
v_eigenmode = squeeze(eigenmodes_separated(:,2,:));
w_eigenmode = squeeze(eigenmodes_separated(:,3,:));


%% Normalizing eigenmodes


% u_eigenmode = u_eigenmode/(norm(u_eigenmode, 'fro'));
% v_eigenmode = v_eigenmode/(norm(v_eigenmode, 'fro'));
% w_eigenmode = w_eigenmode/(norm(w_eigenmode, 'fro'));

for i = 1:nplanes
    u_eigenmode(:,i) = u_eigenmode(:,i)/norm(u_eigenmode(:,i),2);
    v_eigenmode(:,i) = v_eigenmode(:,i)/norm(v_eigenmode(:,i),2);
    w_eigenmode(:,i) = w_eigenmode(:,i)/norm(w_eigenmode(:,i),2);
end


%% Plotting the eigenmodes 

figure;
hold on;
h1 = plot(rc./half_width_tke(1,1), abs(w_eigenmode(:,4)).^2, 'b-', 'Linewidth', 2);
h2 = plot(rc./half_width_tke(2,1), abs(w_eigenmode(:,8)).^2, 'r-', 'Linewidth', 2);
h3 = plot(rc./half_width_tke(3,1), abs(w_eigenmode(:,12)).^2, 'k-', 'Linewidth', 2);
h4 = plot(rc./half_width_tke(4,1), abs(w_eigenmode(:,16)).^2, 'k--', 'Linewidth', 2);
h5 = plot(rc./half_width_tke(5,1), abs(w_eigenmode(:,20)).^2, 'r--', 'Linewidth', 2);
h6 = plot(rc./half_width_tke(6,1), abs(w_eigenmode(:,22)).^2, 'b--', 'Linewidth', 2);

hXLabel = xlabel('$r/L_{Kmean}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|^{2}$','interpreter','latex','fontsize',15);
% %hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], 'x/D = 20', 'x/D = 40','x/D = 60', 'x/D = 80', 'x/D = 100', 'x/D = 120');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'phi_w_m2_st0_spodmode1_3dplane_tkebased.png','-dpng','-r600');  
print(gcf,'phi_w_m2_st0_spodmode1_3dplane_tkebased.eps','-depsc','-r600');
