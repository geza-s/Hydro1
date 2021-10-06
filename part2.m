% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 2: Fit a Gumbel distribution and calculate critical rainfall depths
% -------------------------------------------------------------------------

clear variables %clear the workspace variables
close all %close alla figures
clc %clear the command window

% -------------------------------------------------------------------------
% # 1: Compute the Weibull plotting position
% -------------------------------------------------------------------------

% import the data from Part1 using the function load
load assignment1_output_part1.mat

rainfall_duration = 3;
ranks = [1:size(annualMax, 1)]';
h = sort(annualMax(:,rainfall_duration));
weibull_table = [ranks h];
if (length(h)-length(unique(h)) ~= 0)
    disp("Warning: some maximal values are not unique")
end

%% Computing the empirical freqeuencies

empirical_freq = flip(ranks)/ranks(end);

weibull_table = [weibull_table empirical_freq];

%% Computing weibul position

weibull_position = ranks/(length(ranks)+1);

weibull_table = [weibull_table weibull_position];

%% Computing the reduced variable

yF = -log(-log(weibull_position));

weibull_table = [weibull_table yF];

%% Gumpel parameters: moment methods

mh = mean(h);
sh = std(h);

% WITH MATH TOOLBOX :
%syms u alpha
%eqns = [ u + 0.5772/alpha == m, pi/(sqrt(6)*alpha)  == s];
%vars = [u alpha];
%[solu, sola] = solve(eqns,vars);

alpha_mm= pi/(sh * sqrt(6));
u_mm = mh-0.5772/alpha_mm;
 
%% Gumbel parameters : gumbel method

%mF = mean(weibull_position); %mean of reduced variable ??????
mYF = mean(yF);
sYF = std(yF);

alpha_YF = sYF/sh;
u_YF = mh - mYF*sh/sYF;

%% plottting gumbel functions
x_h = [1:90];
gumpel_mm = exp(-exp(-alpha_mm*(x_h-u_mm)));
gumpel_gm = exp(-exp(-alpha_YF*(x_h-u_YF)));

plot(x_h, gumpel_mm)
hold('on')
plot(x_h, gumpel_gm, 'c')
hold('on')
plot(h, weibull_position, 'or', ...
    'markersize', 3)
legend('moments method', 'gumbel method', 'empirical values',...
    'Location','southeast', 'Fontsize', 15);
title('Gumpel distributions against rainfall depth for fixed duration');






