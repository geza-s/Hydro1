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
event_3 = h_100*[1,1,1,1];

% NOTE: it can be useful to insert the three events into a 4x3 matrix

events = [event_1 event_2 event_3'];
clear event_1 event_2 event_3;
% -------------------------------------------------------------------------
% #2-3: IMPLEMENTATION OF CURVE NUMBER METHOD
load SCSpars.mat
% average CN for the catchment.
average_CN = CN*(percentages/100)
%potential maximum soil moisture retention S of the catchment
S= 25400/average_CN - 254 % units:[mm]
% cumulative precipitation (P [mm]);
%The initial abstraction (Ia [mm]);
% The cumulative infiltration (Fa [mm]);
% The cumulative effective precipitation (Pe [mm]);
% The infiltration intensity (I [mm/h]);
% The effective rainfall intensity (Je [mm/h]);

% -------------------------------------------------------------------------
% NB: the CN method uses cumulative fluxes








