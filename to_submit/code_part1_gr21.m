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
h_100 = h_100/4; %because it's a total mm over 4hr
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

%% a) Cumulative precipitation (P [mm]);
Pt = tril(ones(4))*events;

figure
bar(Pt)
title("Cumulative precipitation")
ylabel('P [mm]')
xlabel('Timestep [hour]')
legend(leg, 'Location', 'NorthWest')

%% b) The initial abstraction (Ia [mm]);
Ia = 0.2 * S;

%% c) The cumulative infiltration (Fa [mm]);
Fa = S*(Pt-Ia)./(Pt-Ia+S);
Fa(Fa<0) = 0;
%We have to check that there is no more infiltration than Pt!!!!!

figure
bar(Fa)
ylabel('Fa [mm]')
xlabel('timestep [hour]')
title("Cumulative Infiltration Fa(t)")
legend(leg,'Location', 'NorthWest')


%% d) The cumulative effective precipitation (Pe [mm]);
Pet = Pt-Ia-Fa;
Pet(Pet<0)=0;

figure 
bar(Pet)
title("Cumulative effective precipitation")
ylabel('Pe [mm]')
xlabel('timestep [hour]')
legend(leg)
%% e) The infiltration intensity (I [mm/h]);

M = diag(ones([4 1])) + diag(-ones([3 1]),-1);
I = M * Fa;
bar(I)
xlabel("Timestep [hour]")
ylabel("Intensity [mm/h]")
title("Infiltration intensity")
legend(leg)
%% f) The effective rainfall intensity (Je [mm/h]);

Je = M * Pet;
%plot(Je, '-o')
bar(Je)
title("Effective rainfall intensity")
ylabel("Intensity [mm/h]")
xlabel("Timestep [hour]")
legend(leg)
% -------------------------------------------------------------------------
%% Plotting together 
ti = 1:4;

figure 

h1 = subplot(2,3,1);
a = bar(ti, [I(:,1)'; Je(:,1)'], 'stacked', 'Barwidth', 1);
a(1).FaceColor = '#8B4513';
a(2).FaceColor = '#0072BD';
hold on 
grid on
plot(events(:,1), 'o', 'Linewidth', 3)
xlabel("Timestep [hour]")
ylim([0,20])
ylabel("intensity [mm/h]")
title("Event 1")

h2 = subplot(2,3,2);
a = bar(ti, [I(:,2)'; Je(:,2)'], 'stacked', 'Barwidth', 1);
a(1).FaceColor = '#8B4513';
a(2).FaceColor = '#0072BD';
hold on 
grid on
plot(events(:,2), 'o' ,'Linewidth', 3)
xlabel("Timestep [hour]")
ylim([0,20])
ylabel("intensity [mm/h]")
%lgd = legend("Infiltration intensity portion", "Effective Precipitation intensity portion", "Rainfall Intensity", "Location", "southoutside" );
title("Event 2")

h3 = subplot(2,3,3);
a = bar(ti, [I(:,3)'; Je(:,3)'], 'stacked', 'Barwidth', 1);
a(1).FaceColor = '#8B4513';
a(2).FaceColor = '#0072BD';
hold on 
grid on
plot(events(:,3), 'o', 'Linewidth', 3)
%stairs(events(:,3))
lgd = legend("Infiltration intensity portion", "Effective Precipitation intensity portion", "Rainfall Intensity" );
xlabel("Timestep [hour]")
ylim([0,20])
ylabel("intensity [mm/h]")
title("Event 3")

hL = subplot(2,3,4.5);
poshL = get(hL,'position');     % Getting its position
posfinal = [poshL(1) + poshL(3), poshL(2) + poshL(4)/1.5, poshL(3)/2, poshL(4)/2];
set(lgd,'position', posfinal);      % Adjusting legend's position
axis(hL,'off');                 % Turning its axis off

%% Compute the total effective precipitation for each event 
%NB : events are in mm/h, therefore cumulative is in mm

figure
b = bar(Pet(4,:));
grid on;
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,1));
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

title("Total effective precipitation for each event")
xticklabels({'Event 1', 'Event 2', 'Event 3'})
ylabel("Cumulated precipitaiton [mm]")


%% Saving the output necessery for part2
J = events;
save("output_part1.mat", "J", "Je", "I") 






