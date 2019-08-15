%% Script to run the spod_spectrum_ploot_imag_data.m for plotting the spectra


clc; clear;
close all;

x = [20;40;60;80;100;120];


for i = 1:size(x,1)
    spod_spectrum_azmode_contribution(x(i,1));
end