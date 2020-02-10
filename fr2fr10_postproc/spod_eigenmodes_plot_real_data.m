

%% Name - Sheel Nidhan
%  Date - 8th February, 2020

%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7000;
stride = 100;
nstart = 2329600;
nend = nstart + (N-1)*stride;
nr = 333;
ntheta = 256;
numvar = 4;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nblk_sampled = 3;
Nfreq_sampled = 25;
Nrows = numvar*nr*ntheta;

dir2 = strcat('/home/sheel/Work2/projects_data/spod_re5e4/fr2/spod_data/run_3.0/x_D_', int2str(x), '/eigenmodes/');
disp(dir2);

%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  %#ok<*SAGROW> % Centered the grid faces to grid centers
end

rc = rc(1:333);

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

%% Reading SPOD eigenvalues files

eigmodes = zeros(Nrows, Nfreq_sampled, Nblk_sampled);

for i = 1:Nfreq_sampled
    for j = 1:Nblk_sampled
        
    freq = sprintf('%04d',i);   
    nblk_name = sprintf('%02d', Nblk-j+1);
    
    filename = strcat(dir2, 'eigenmode_freq_',freq,'_',nblk_name,'_real.mod');
    disp(filename);
    fid = fopen(filename,'rb');              
    h = fread(fid,1,'*int32'); % May need adjusting
    a = fread(fid, nr*ntheta*numvar, '*double');
    
    filename = strcat(dir2, 'eigenmode_freq_',freq,'_',nblk_name,'_imag.mod');
    disp(filename);
    fid = fopen(filename,'rb');              
    h = fread(fid,0,'*int32'); % May need adjusting
    b = fread(fid, nr*ntheta*numvar, '*double');
    
    eigmodes(:,i,j) = a(:,1) + sqrt(-1).*b(:,1);
    
    fclose all;
    
    end
end

%% Reshape the eigmodes array 

eigmode_reshaped = reshape(eigmodes, [nr,ntheta,numvar,Nfreq_sampled,Nblk_sampled]);


%% Plotting the eigenmodes 

theta = linspace(0,2*pi,ntheta)';
plot_mode = real(eigmodes(:,:,4,1,7));
[C,h,x,y] = polarcont(rc, theta, plot_mode, 10);

axis equal;
set(h,'Linecolor','none');
set(h,'LevelList',-0.25:0.025:0.25)
colormap jet;
colorbar;

