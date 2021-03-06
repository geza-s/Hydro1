% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 2: Fit a Gumbel distribution and calculate critical rainfall depths
% -------------------------------------------------------------------------

clear variables %clear the workspace variables
close all %close alla figures
clc %clear the command window

% -------------------------------------------------------------------------
%% # 1: Compute the Weibull plotting position
% -------------------------------------------------------------------------

% import the data from Part1 using the function load
load assignment1_output_part1.mat
weibull_table = [];
gumbel_table = []; %table to store gumbel method parameters

%% All the computing and plotting: 

figure( 'position', [300, 50, 900, 400])
col = [    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840]; % colors for the graphes 


ranks = [1:size(annualMax, 1)]';
Fh_matrix = [];
posR = 0;
for rainfall_duration = D
    h = sort(annualMax(:,rainfall_duration));
    weibull_matrix = [ranks h];
    if (length(h)-length(unique(h)) ~= 0)
        disp("Warning: some maximal values are not unique")
        disp(rainfall_duration)
    end

    %%% Computing the empirical freqeuencies

    empirical_freq = flip(ranks)/ranks(end);    
    weibull_matrix = [weibull_matrix empirical_freq];

    %%% Computing weibul position

    weibull_position = ranks/(length(ranks)+1);
    weibull_matrix = [weibull_matrix weibull_position];

    %%% Computing the reduced variable

    yF = -log(-log(weibull_position));
    weibull_matrix = [weibull_matrix yF];

    %%% gumbel parameters: moment methods

    mh = mean(h);
    sh = std(h);

    alpha_mm= pi/(sh * sqrt(6));
    u_mm = mh-0.5772/alpha_mm;
    

    %%% Gumbel parameters : gumbel method

    %mF = mean(weibull_position); %mean of reduced variable ??????
    mYF = mean(yF);
    sYF = std(yF);

    alpha_YF = sYF/sh;
    u_YF = mh - mYF*sh/sYF;

    vect_param = [alpha_mm, u_mm, alpha_YF, u_YF]';
    gumbel_table = [gumbel_table vect_param];
    
    weibull_table = cat(3, weibull_table, weibull_matrix);
    %%% plottting comparision of gumbel functions
    x_h = [1:120];
    %gumbel_mm = exp(-exp(-gumbel_table(1,3)*(x_h-gumbel_table(2,3))));
    %gumbel_gm = exp(-exp(-gumbel_table(3,3)*(x_h-gumbel_table(4,3))));
    gumbel_mm = exp(-exp(-alpha_mm*(x_h-u_mm)));
    gumbel_gm = exp(-exp(-alpha_YF*(x_h-u_YF)));
    
    %strMM =("Moments method, duration = " + rainfall_duration + " H");
    %plot(x_h, gumbel_mm, 'k', displayName = )
    %hold('on')
    posR = posR + 1;
    strGM =("Gumbel method, duration = " + rainfall_duration + " [h]");
    plot(x_h, gumbel_gm, 'LineWidth', 1.5, 'DisplayName', strGM, 'Color', col(posR,:))
    hold('on')
    strEF =("Empirical frequency, duration = " + rainfall_duration + " [h]");
    plot(h, weibull_position, 'o', 'markersize', 3, 'Color', col(posR,:), ...
        'DisplayName', strEF);
    pause(1)
end

%Some additional elements for the graph:
title('Fitted Cumulative Gumbel Distribution', 'fontsize', 14);
xlabel('rainfall maxima h in [mm]', 'fontsize', 14);
ylabel('non-exceedance probability Fh', 'fontsize', 14);
%legend('moments method', 'gumbel method', 'empirical values',...
legend('Location','eastoutside', 'Fontsize', 12);
%%%weibull_table: [ranks h empirical_freq weibull_pos yF_reduced_variable] 
%% Compute return period for each duration

clear h Ph emp_Th;

%Th based on the weibull position
h = weibull_table(:,2,1);
Ph = weibull_table(:,4,1);
emp_Th = 1./(1-Ph); %Th associated to each empirical h vector ! 

%disp('Empirical return periods : OK') %debugging

%% 

MTT = [];
T_array = [1:1:100];
for i=[1:6]
    alpha = gumbel_table(3,i);
    u = gumbel_table(4,i);
    d = 1-(1./T_array);
    h = u - (1/alpha)*log(-log(d));
    MTT = [MTT h'];
end

MTT = [MTT T_array'];

%disp('h related to return periods T based on gumbel : Ok') %debugging

%% Create Matric H_Gum
H_Gum = [];
Ti = [10,40,100]';
for i=[1:6]
    alpha = gumbel_table(3,i);
    u = gumbel_table(4,i);
    hi = u - (1/alpha)*log(-log(1-(1./Ti)));
    H_Gum = [H_Gum hi];
end
%disp('Matric with 10,40,100 return year: ok') %debugging

%% plot this 

figure( 'position', [300, 300, 850, 500])
 % new figure for the second plot 

%col = ['#0072BD', '#EDB120', '#EDB120', '#EDB120','#EDB120']; 

for i= [1:6]
    strEmp = ("Empirical values, duration = " + D(i) + " [h]");
    h = plot(emp_Th, weibull_table(:,2,i), 'o', 'markersize', 5, 'DisplayName', strEmp);
    set(h,'Color',col(i,:));
    hold on;
    
    strG = ("Gumbel approx., duration = " + D(i)+ " [h]");
    h = plot(MTT(:,7), MTT(:,i), '-','LineWidth',1.5, 'DisplayName', strG);
    set(h,'Color',col(i,:));
    
    pHgum = plot(Ti, H_Gum(:,i), 'ok');
    if (i == 6)
        pHgum.DisplayName = "Selected values for DDF curves";
    else 
        pHgum.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    pause(1)
end

xlabel('Return Period Th [years]','fontsize',14)
ylabel('Rainfall depth h [mm]','fontsize',14)
title('Rainfall depth for each duration','fontsize',14)

legend('Location', 'eastoutside', 'Fontsize', 12)
%legend('h[mm] against empirical return periods', 'h[mm] with Gumbel distribution', 'Output Matrix: H Gum','Location','southeast', 'Fontsize', 14);

%% Saving output for next part

T = Ti;
save('assignment1_output_part2', 'H_Gum', 'T', 'D');

disp('H_Gum, T and D saved in ouput file.')

