% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 3: Construction of depth-duration-frequency (DDF) curves 
% -------------------------------------------------------------------------

clear variables %clear the workspace variables
close all %close alla figures
clc %clear the command window

% -------------------------------------------------------------------------
% # 1: Calibrate the DDF curve parameters
% useful functions: linspace
% -------------------------------------------------------------------------

% import the data from previous part
load assignment1_output_part2.mat 



% you need DDF curves for each return period, so you will have a set of
% parameters for each return period

% define lists of test values for each parameter (e.g.
% c_list=linspace(0,100,3))

% implement the 'brute force' algorithm 
% -implement nested for loops over every parameter combination 
% -in the core of the loop, compute the depths estimated by the DDF curve
%   with the given return periods and parameter combination
% -compute the sum of square errors between the computed dephts and the
%   ones estimated in part 1. 
% -retain the sets that have the lowest sum of squared errors 
% -check if the best-fit errors are lower than the maximal errors given. If
%   not, repeat the loop with a larger amount of test values


% ...
% ...
% ...


