%% Script to run the spod_spectrum_ploot_imag_data.m for plotting the spectra


clc; clear;
close all;

x = [20; 80];
m = [0; 1; 2];

% Write out the files for plotting in paper

dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';

count = 1;
for i = 1:size(x,1)
    for j = 1:size(m,1)
        [eigvalue,f] = spod_spectrum_plot_imag_data(x(i,1), m(j,1));
        eigvalues_spectrum_plot(count,1).loc = x(i,1); %#ok<*SAGROW>
        eigvalues_spectrum_plot(count,1).azmode = m(j,1);
        eigvalues_spectrum_plot(count,1).freq = f;
        eigvalues_spectrum_plot(count,1).eigvalue = eigvalue;
        count = count + 1;
    end
end

save(strcat(dirout, 'eigvalue_spectrum_diff_loc.mat'), 'eigvalues_spectrum_plot')