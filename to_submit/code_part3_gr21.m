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

N = length(IUHw); %length of watershed IUH
M = length(Je); %length of effective precipitation array

Qc = zeros([1, N+M]);

for n = 1:N+M
    for m = max(1,n-N+1):min(n,M)
        Qc(n) = Qc(n) + Je(m)*IUHw(n-m+1)*dt;
    end
end


%% just to visualize (can be removed)

figure
area(Qc(1:1000), 'FaceAlpha', 0.5)
hold on
area(Je(1:1000),'FaceAlpha', 0.5)
legend("Channel Discharge", "Effective precipitation")
ylabel(" magnitude [mm/h]")
title("Just to visualize, little time window on the result")
xlabel("time in " + dt + " [hours] timestep")

%% ------------------------------------------------------------------------
% #4: SUBSURFACE CONTRIBUTION
% -------------------------------------------------------------------------

%% 4.1 create subsurface IUH 
load IUHpars.mat % par_shape = k, par_scale = theta
par_scale = par_scale*10;
x = 0:dt:200; % lets take 35 to have same length IUH's

fun = @(x) x.^(0.2).*exp(-x);
%Evaluate the integral from x=0 to x=Inf.
gammaK = integral(fun,0,Inf);

IUHsub = (1/(gammaK * par_scale^par_shape)) .* x.^(par_shape-1) .* exp(-x/par_scale);

%% check IUH

figure
bar(x,IUHsub,'FaceAlpha', 0.5, 'Barwidth', 1)
hold on;
bar(x(1:length(IUHw)),IUHw, 'FaceAlpha', 0.5, 'Barwidth', 1)
title("Watershed IUHs")
legend("Subsurface IUH", "Surface IUH");
xlim([0,150])
xlabel("Time [hour]")
ylabel("Rate [mm/h]")
disp("Surface under curve of subsurface IUH: " + sum(IUHsub)*dt);

%% 4.2 convolve with infiltration

N = length(IUHsub); %length of IUHsub
M = length(I); %length of infiltration array

Qsub = zeros([1, N+M]);

for n = 1:N+M
    for m = max(1,n-N+1):min(n,M)
        Qsub(n) = Qsub(n) + I(m)*IUHsub(n-m+1)*dt;
    end
end

%% 4.3 check the result 

figure
area(Qsub(1:1000), 'FaceAlpha', 0.5);
hold on
area(I(1:1000),'FaceAlpha', 0.5)
legend("Subsurface contribution", "Infiltration rate")
ylabel(" magnitude [mm/h]")
xlabel("time in " + dt + " [hours] timestep")
title("Sample of the subsurface runoff response")

%% ------------------------------------------------------------------------
% #5: TOTAL DISCHARGE
% -------------------------------------------------------------------------
totalQ = zeros([1,max(length(Qc), length(Qsub))]);
totalQ(1:length(Qsub)) = Qsub; 
totalQ(1:length(Qc)) = totalQ(1:length(Qc)) + Qc;

%% ------------------------------------------------------------------------
% #6: NOVEMBER FIGURE 
% -------------------------------------------------------------------------
i = 1:length(t); %array for simple indexing
novI = i(month(t)==11);

figure

subplot(3,1,[1,2])
a = bar(t(novI), [Qsub(novI); Qc(novI)], 'stacked', 'Barwidth', 1);
a(1).FaceColor = '#8B4513';
a(2).FaceColor = '#0072BD';
legend("Subsurface portion","Surface portion")
ylim([0,0.18])
grid on;
xlabel("Time");
ylabel("Rate [mm/h]")
title("Evolution of the total runoff through November month")

subplot(3,1,3)
bar(t(novI), Qc(novI)./totalQ(novI),'Barwidth', 1)
title([" ";"Proportion of surface runoff over total runoff"])
grid on;
xlabel("Time")
ylabel("Proportion")

%% ------------------------------------------------------------------------
% #7: COMPUTE FRACTIONS 
% -------------------------------------------------------------------------

partH = length(totalQ(totalQ>0.2))/length(totalQ);
partL = length(totalQ(totalQ<0.002))/length(totalQ);

figure

map = [255,255,189;
172,255,155;
89,255,121;
45,198,127;
15,121,146;
0,75,163;
0,61,178;
23,49,191;
114,42,203;
205,35,214];
map = map/255;

pie([partL, 1-partH-partL, partH], '%.1f%%')
colormap(map);
legend("Total discharge < 0.02 mm/h", "0.02 mm/h < total discharge < 0.2 mm/h", "total discharge > 0.2 mm/h", 'Location', 'EastOutside') 
title("Fractions of the all the total discharges fluxes")


