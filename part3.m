% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 2
% Part 3: Rainfall-runoff transformation for one year of rainfall data
% -------------------------------------------------------------------------

% Clear the workspace and figure
clear variables
close all

% -------------------------------------------------------------------------
% #1: Prepare rainfall event
% -------------------------------------------------------------------------

% load discretization timestep and IUHw from previous part
load Ass2_part2_result.mat

% load rainfall data (there are no empty values)
T = readtable('rainfall_data.txt',...
    'HeaderLines', 2,...
    'Format','%s%s%f'); %the format is: text string, text string and float number

% Create a vector h containing the hourly precipitation and a vector t
% containing the timestamp
h = T.rre150h0;    %rre150h0 is the MeteoSwiss code for hourly rainfall depth [mm]
t_h = datetime(T.time,'InputFormat','yyyyMMddHH'); %convert to datetime (can be slow)

% uniform downscaling of precipitation from hourly to dt time steps
t = t_h(1):hours(dt):t_h(end)+1-hours(dt); %timestamps at dt time steps
J = interp1(t_h, h, t,'previous','extrap'); %rainfall intensity in mm/h for each timestep


% -------------------------------------------------------------------------
% #2: RAINFALL SEPARATION
% -------------------------------------------------------------------------





