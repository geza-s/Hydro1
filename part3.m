% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 2
% Part 3: Rainfall-runoff transformation for one year of rainfall data
% -------------------------------------------------------------------------

% Clear the workspace and figure
clear variables
close all

%% -------------------------------------------------------------------------
% #1: Prepare rainfall event
% -------------------------------------------------------------------------

% load discretization timestep and IUHw from previous part
load output_part2.mat

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


%% ------------------------------------------------------------------------
% #2: RAINFALL SEPARATION
% -------------------------------------------------------------------------

Je = 0.3*J;
I = 0.1*J;

%% ------------------------------------------------------------------------
% #3: SURFACE CONTRIBUTION
% -------------------------------------------------------------------------

N = length(IUHw); %length of IUHc
M = length(Je); %length of effective precipitation array

Qc = zeros([1, N+M]);

for n = 1:N+M
    for m = 1:M
        if ((n>m) && n-m+1<N)
            Qc(n) = Qc(n) + Je(m)*IUHw(n-m+1)*dt;
        end
    end
end
%% just to visualize 

figure
area(Qc(1:1000), 'FaceAlpha', 0.5)
hold on
area(Je(1:1000),'FaceAlpha', 0.5)
legend("Channel", "Effective precipitation")
ylabel(" magnitude [mm]")
xlabel("time in " + dt + " [hours] timestep")
%% ------------------------------------------------------------------------
% #4: SUBSURFACE CONTRIBUTION
% -------------------------------------------------------------------------

%% 4.1 create subsurface IUH 
load IUHpars.mat % par_shape = k, par_scale = theta
par_scale = par_scale*10;
x = 0:dt:35; % lets take 35 to have same length IUH's

fun = @(x) x.^(0.2).*exp(-x);
%Evaluate the integral from x=0 to x=Inf.
gammaK = integral(fun,0,Inf);

IUHsub = (1/(gammaK * par_scale^par_shape)) .* x.^(par_shape-1) .* exp(-x/par_scale);

%% check IUH

figure
area(IUHw,'FaceAlpha', 0.5)
title("Watershed IUH")

disp("Surface under curve of subsurface IUH: " + sum(IUHsub)*dt);

%% 4.2 convolve with infiltration

N = length(IUHsub); %length of IUHc
M = length(I); %length of effective precipitation array

Qsub = zeros([1, N+M]);

for n = 1:N+M
    for m = 1:M
        if ((n>m) && n-m+1<N)
            Qsub(n) = Qsub(n) + I(m)*IUHsub(n-m+1)*dt;
        end
    end
end

%% 4.3 check the result 

figure
area(Qsub(1:1000), 'FaceAlpha', 0.5)
hold on
area(I(1:1000),'FaceAlpha', 0.5)
legend("Subsurface contribution", "Infiltration rate")
ylabel(" magnitude [mm/h]")
xlabel("time in " + dt + " [hours] timestep")

%% ------------------------------------------------------------------------
% #5: TOTAL DISCHARGE
% -------------------------------------------------------------------------

totalQ = Qsub + Qc; %Since we took the same length IUH, they are same length

%% ------------------------------------------------------------------------
% #6: NOVEMBER FIGURE 
% -------------------------------------------------------------------------
i = 1:length(t); %array for simple indexing
novI = i(month(t)==11);

figure
subplot(3,1,1)

area(t(novI), totalQ(novI),'FaceColor', '#008080', 'FaceAlpha', 0.5)
hold on;
area(t(novI), Qc(novI),'FaceColor', '#0072BD', 'FaceAlpha', 0.5)
hold on;
area(t(novI), Qsub(novI),'FaceColor', '#8B4513', 'FaceAlpha', 0.9)
legend("Total runoff", "Surface runnoff", "Subsurface runoff")

title("Runoffs")

subplot(3,1,2)
a = bar(t(novI), [Qc(novI);(Qsub(novI))], 'stacked' );
a(1).FaceColor = '#0072BD';
a(2).FaceColor = '#8B4513';
legend("Surface runnoff", "Subsurface runoff")

subplot(3,1,3)
area(t(novI), Qc(novI)./totalQ(novI))
title("Proportion of surface runoff over total runoff")

%% ------------------------------------------------------------------------
% #6: NOVEMBER FIGURE 
% -------------------------------------------------------------------------

partH = length(totalQ(totalQ>0.2))/length(totalQ);
partL = length(totalQ(totalQ<0.002))/length(totalQ);

figure
pie([partL, 1-partH-partL, partH])
legend("Total discharge < 0.02", "0.02 < total discharge < 0.2", "total discharge > 0.2", 'Location', 'SouthEast') 
title("Fraction of the year of total discharge")