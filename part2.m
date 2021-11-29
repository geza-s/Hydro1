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
x = 0:dt:35;

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

disp("Surface under curve of IUHw: " + sum(IUHw)*dt);

%% Channel IUH

c = 0.3*3600; %m/s * s/h ;celerity
D = 10^6; %m^2/h ;hydrodynamic dispersion
L = 7 * 1000; %km -> m

t = 1:dt:20;
IUHc = L./(sqrt(4*pi*D).*t.^(3/2)) .* exp((-(L-c*t).^2)./(4*D.*t));
figure 
bar(IUHc)
title("Channel IUH")

disp("Surface under curve of IUHc: " + sum(IUHc)*dt);
%% Plot the two iuh in same figure, focus early part

figure
IUHcx = [IUHc, zeros(1,length(x)-length(IUHc))];
bar(x, [ IUHw; IUHcx])
title("Watershed and channel IUH")
legend("IUHw", "IUHc")

%% explicit convolution IUHw and Je

event = 1; %event number 1,2 or 3
Jei = Je(:,event);

N = length(IUHw);
M = length(Je);

Qw = zeros([1, N+M]);

%%convolution
for n = 1:N+M
    for m = 1:M
        %disp("n :" + n + " and m =" + m);
        if ((n>m) && n-m+1<N)
            Qw(n) = Qw(n) + Jei(m)*IUHw(n-m+1);
        end
    end
end

SJei = sum(Jei)*dt; %% la quantité total d'eau tombé
SJei_waited = sum(IUHw)*dt * SJei*dt; %%ce qu'on devrait attendre au vu de l'approximation du IUHw
SQw = sum(Qw)*dt*dt; %%total de pluie qu'on obtient à la sortie
disp("With a IUHw surface of " + sum(IUHw)*dt + " we should get a total water of " + SJei_waited + " but by the convolution we have : " + SQw + " !")
disp("Less than " + round(1-SQw/SJei_waited, 3)*100 + "% is lost through convulotion algorithm")

figure 
bar(Jei)
hold on
bar(Qw)
hold on
bar(IUHw)
legend("Jei", "Qw", "IUHw")
title("Convolution of IUH and Je of event " + event)

%% Channel IUH convolution

N = length(IUHc);
M = length(Qw);

Qc = zeros([1, N+M]);

%%convolution
for n = 1:N+M
    for m = 1:M
        %disp("n :" + n + " and m =" + m);
        if ((n>m) && n-m+1<N)
            Qc(n) = Qc(n) + Qw(m)*IUHc(n-m+1);
        end
    end
end

SQc = sum(Qc)*dt;
disp("Surface under Qw :" + SQw + " and under Qc :" + SQc)


figure 
bar(Qc)
hold on
bar(Qw)
hold on
bar(IUHc)
legend("Qc", "Qw", "IUHc")
title("Convolution of IUHc and Qw of event " + event)