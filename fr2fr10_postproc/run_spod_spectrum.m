%% Script to run the spod_spectrum_ploot_imag_data.m for plotting the spectra


clc; clear;
close all;

x = [10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100];

% Write out the files for plotting in paper

dirout = './';

count = 1;

for i = 1:size(x,1)
    disp(i);
    [eigvalue,f] = spod_spectrum_plot_real_data(x(i,1));
    eigvalues_spectrum_plot(count,1).loc = x(i,1); %#ok<*SAGROW>
    eigvalues_spectrum_plot(count,1).freq = f;
    eigvalues_spectrum_plot(count,1).eigvalue = eigvalue;
    count = count + 1;
    fclose all;
end

save(strcat(dirout, 'fr2_eigvalue_spectrum_slice.mat'), 'eigvalues_spectrum_plot')