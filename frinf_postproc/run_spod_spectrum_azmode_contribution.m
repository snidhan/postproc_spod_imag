%% Script to run the spod_spectrum_ploot_imag_data.m for plotting the spectra


clc; clear;
close all;

x = [20; 40; 80; 100];

dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';

for i = 1:size(x,1)
     [FREQ, MODE, eigvalue, mode, eigenvalue_fraction_contri] = spod_spectrum_contribution_of_azimuthalmode(x(i,1));
     eigvalues_contourf_plot(i,1).loc      = x(i,1);
     eigvalues_contourf_plot(i,1).FREQ     = FREQ;
     eigvalues_contourf_plot(i,1).MODE     = MODE;
     eigvalues_contourf_plot(i,1).eigvalue = eigvalue;
     
     eigvalues_bar_plot(i,1).loc      = x(i,1); %#ok<*SAGROW>
     eigvalues_bar_plot(i,1).mode      = mode;
     eigvalues_bar_plot(i,1).eigenvalue_fraction_contri = eigenvalue_fraction_contri;
end

save(strcat(dirout, 'eigvalue_contourf_diff_loc.mat'), 'eigvalues_contourf_plot')
save(strcat(dirout, 'eigvalue_barplot_diff_loc.mat'),  'eigvalues_bar_plot')