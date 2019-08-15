%% Sheel Nidhan
%  Date - December 5, 2018

%% Loading the data

clear;
clc;

case18 = load('case18_eigv.mat');
case11 = load('case11_eigv.mat');

diff = abs(case18.eigvalue(:,1:10) - case11.eigvalue(:,1:10));

perc_diff = 100*(diff ./ case18.eigvalue(:,1:10));

min_perc_diff = min(min(perc_diff));
max_perc_diff = max(max(perc_diff));