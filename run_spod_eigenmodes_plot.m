%% Script to run the spod_spectrum_ploot_imag_data.m for plotting the spectra

clc; clear;
close all;

x = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
% x = [100];

m = [0; 1; 2; 3; 4; 5];
% m = 2;

Nf_sampled = 50;
Nblk_sampled = 5;

for i = 1:size(x,1)
    for j = 1:size(m,1)
        spod_eigenmodes_plot(x(i,1), m(j,1), Nf_sampled, Nblk_sampled);
    end
end