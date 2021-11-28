% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 2
% Part 1: rainfall partitioning using the SCS Curve Number method
% -------------------------------------------------------------------------

% Clear the workspace and figure
clear variables
close all

% -------------------------------------------------------------------------
% #1: Prepare rainfall event
% -------------------------------------------------------------------------
% load rainfall depth data of events 1 and 2 (mm of rainfall every hour,
% for four hours)
load event_1
load event_2

% for event 3, use the DDF curve from previous assignment to find
% precipitation DEPTH corresponding to: a return period of 100 years and
% duration of 4 hours
duration = 4; %100 years return period but 4 hr duration
load ddf100_param.mat;
h_100 = (ddf100_param(1) .* duration)./(duration.^ddf100_param(2) + ddf100_param(3));
event_3 = round(h_100,1)*[1,1,1,1];

% NOTE: it can be useful to insert the three events into a 4x3 matrix

events = [event_1 event_2 event_3'];
leg = [ " event 1", "event 2", "event 3"];
clear event_1 event_2 event_3;
% -------------------------------------------------------------------------
%% #2-3: IMPLEMENTATION OF CURVE NUMBER METHOD
load SCSpars.mat

% average CN for the catchment.
average_CN = CN*percentages/100;

%potential maximum soil moisture retention S of the catchment
S= 25400/average_CN - 254; % units:[mm]

%% cumulative precipitation (P [mm]);
Pt = tril(ones(4))*events;
figure
plot(Pt, '-o')
title("Cumulative precipitation")
ylabel('P [mm]')
xlabel('timestep [hour]')
legend(leg)
%% The initial abstraction (Ia [mm]);
Ia = 0.2 * S;

%% The cumulative infiltration (Fa [mm]);
Fa = Pt - ((Pt - Ia).*(Pt-Ia))./(Pt + 0.8*S);
figure
plot(Fa, '-o')
title("Cumulative Infiltration Fa(t)")
ylabel('Fa [mm]')
xlabel('timestep [hour]')
legend(leg)
%% The cumulative effective precipitation (Pe [mm]);
Pet = Pt - Fa;
figure 
plot(Pet, '-o')
title("Cumulative effective precipitation")
ylabel('Pe [mm]')
xlabel('timestep [hour]')
legend(leg)
%% The infiltration intensity (I [mm/h]);

M = diag(ones([4 1])) + diag(-ones([3 1]),-1);
M = M;
I = M * Fa;
plot(I, '-o')
title("Infiltration intensity")
legend(leg)
%% The effective rainfall intensity (Je [mm/h]);

Je = M * Pet;
plot(Je, '-o')
title("Effective rainfall intensity")
ylabel("Je [mm/h]")
legend(leg)
% -------------------------------------------------------------------------
%% NB: the CN method uses cumulative fluxes
J = events;
save("output_part1.mat", "J", "Je", "I") 






