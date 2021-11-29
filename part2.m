% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 2
% Part 2: Rainfall-runoff transformation for one event
% -------------------------------------------------------------------------

% Clear the workspace and figure
clear variables
close all

% -------------------------------------------------------------------------
%% #1: Prepare rainfall event
% -------------------------------------------------------------------------

% Load the matrix with J, Je, I from part 1
load output_part1.mat

dt_h = 1; %hourly time steps
Je_h = Je; %rename to stress that it is Je at hourly time steps

clear Je;

%% display total volume Je_h
disp(' ')
fprintf('Volume Je = %.3f mm\n',sum(Je_h)*dt_h)

% define the time step and some time vectors of the simulation (make column vectors)
nstep=10; % [-] number of time steps per each hour (YOU CAN CHANGE THIS)
dt=1/nstep; %[h]

% create a single vector with input fluxes in dt timesteps
t_Je_h=0:1:length(Je_h)-1; %time vector for the hourly Je values
t_Je=0:dt:t_Je_h(end)+dt_h-dt;   %time vector of effective precipitation in dt timesteps
Je = interp1(t_Je_h, Je_h, t_Je ,'previous','extrap'); %rainfall intensity in mm/h for each timestep

% -------------------------------------------------------------------------
% #2: GENERATE THE INSTANTANEOUS UNIT HYDROGRAPHS (IUH)
% -------------------------------------------------------------------------
%% Watershed IUH
load IUHpars.mat % par_shape = k, par_scale = theta
x = 0:30;

fun = @(x) x.^(0.2).*exp(-x);
%Evaluate the integral from x=0 to x=Inf.
gammaK = integral(fun,0,Inf);

IUHw = (1/(gammaK * par_scale^par_shape)) .* x.^(par_shape-1) .* exp(-x/par_scale);

%IUHw = zeros(1,20);
%for k = 1:19
%    IUHw(k) = sum(gammaf(10*k: 10*k + 9))/10;
%end

figure
bar(IUHw)
title("Watershed IUH")

disp("Surface under curve of IUHw: " + sum(IUHw));

%% Channel IUH

c = 0.3*3600; %m/s * s/h ;celerity
D = 10^6; %m^2/h ;hydrodynamic dispersion
L = 7 * 1000; %km -> m

t = 1:20;
IUHc = L./(sqrt(4*pi*D).*t.^(3/2)) .* exp((-(L-c*t).^2)./(4*D.*t));
figure 
bar(IUHc)
title("Channel IUH")

disp("Surface under curve of IUHw: " + sum(IUHc));






