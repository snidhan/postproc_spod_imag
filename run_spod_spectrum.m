%% Script to run the spod_spectrum_ploot_imag_data.m for plotting the spectra


clc; clear;
close all;

x = [20];
m = [0; 1; 2];

for i = 1:size(x,1)
    for j = 1:size(m,1)
        spod_spectrum_plot_imag_data(x(i,1), m(j,1));
    end
end