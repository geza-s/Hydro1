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
weibull_table = [];
gumbel_table = []; %table to store gumpel method parameters

for rainfall_duration = D
    ranks = [1:size(annualMax, 1)]';
    h = sort(annualMax(:,rainfall_duration));
    weibull_matrix = [ranks h];
    if (length(h)-length(unique(h)) ~= 0)
        disp("Warning: some maximal values are not unique")
        disp(rainfall_duration)
    end

    %% Computing the empirical freqeuencies

    empirical_freq = flip(ranks)/ranks(end);

    weibull_matrix = [weibull_matrix empirical_freq];

    %% Computing weibul position

    weibull_position = ranks/(length(ranks)+1);

    weibull_matrix = [weibull_matrix weibull_position];

    %% Computing the reduced variable

    yF = -log(-log(weibull_position));

    weibull_matrix = [weibull_matrix yF];

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

    vect_param = [alpha_mm, u_mm, alpha_YF, u_YF]';
    gumbel_table = [gumbel_table vect_param];
    
    weibull_table = cat(3, weibull_table, weibull_matrix);
    %% plottting comparision of gumbel functions
    x_h = [1:120];
    %gumpel_mm = exp(-exp(-gumbel_table(1,3)*(x_h-gumbel_table(2,3))));
    %gumpel_gm = exp(-exp(-gumbel_table(3,3)*(x_h-gumbel_table(4,3))));
    gumpel_mm = exp(-exp(-alpha_mm*(x_h-u_mm)));
    gumpel_gm = exp(-exp(-alpha_YF*(x_h-u_YF)));

    plot(x_h, gumpel_mm, 'k')
    hold('on')
    plot(x_h, gumpel_gm, 'c')
    hold('on')
    plot(h, weibull_position, 'or', ...
        'markersize', 3)
    legend('moments method', 'gumbel method', 'empirical values',...
        'Location','southeast', 'Fontsize', 15);
    title('Gumpel distributions against rainfall depth for fixed duration');

    pause(1)
end

%% Compute return period for each duration
MTh = [];
for i = [1:6]
    N = length(h);
    m = weibull_table(:,1,i);
    disp(length(m))
    Th = (N+1)./(N+1-m);
    MTh = [MTh Th];
end



