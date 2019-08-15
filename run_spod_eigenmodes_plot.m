%% Script to run the spod_spectrum_ploot_imag_data.m for plotting the spectra

clc; clear;
close all;

x = [75];
m = [0; 1; 2; 3; 4; 5];
Nf_sampled = 40;
Nblk_sampled = 5;

for i = 1:size(x,1)
    for j = 1:size(m,1)
        spod_eigenmodes_plot(x(i,1), m(j,1), Nf_sampled, Nblk_sampled);
    end
end