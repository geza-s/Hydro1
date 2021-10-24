% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 3: Construction of depth-duration-frequency (DDF) curves 
% -------------------------------------------------------------------------

clear variables %clear the workspace variables
close all %close alla figures
clc %clear the command window

% -------------------------------------------------------------------------
% # 1: Calibrate the DDF curve parameters
% useful functions: linspace
% -------------------------------------------------------------------------

% import the data from previous part
load assignment1_output_part2.mat 

%%
% you need DDF curves for each return period, so you will have a set of
% parameters for each return period

% define lists of test values for each parameter (e.g.
% c_list=linspace(0,100,3))

% implement the 'brute force' algorithm 
% -implement nested for loops over every parameter combination 
% -in the core of the loop, compute the depths estimated by the DDF curve
%   with the given return periods and parameter combination
% -compute the sum of square errors between the computed dephts and the
%   ones estimated in part 1. 
% -retain the sets that have the lowest sum of squared errors 
% -check if the best-fit errors are lower than the maximal errors given. If
%   not, repeat the loop with a larger amount of test values

c_list = linspace(1,100,100); %linspace(start, stop, nval)
f_list = linspace(-1,1,100); %let's begin with 5 vals in each list
e_list = linspace(0,1,100);

DDF_param = [];
Td = D;
for j = [1:3]
    err_T = 1000; %number big enough
    c = 0;
    e = 0;
    f = 0;
    hi_gum = H_Gum(j,:);
    for ci=c_list
        for fi=f_list
            for ei = e_list
                h = (ci .* Td)./(Td.^ei + fi);
                err_temp = sum((h-hi_gum).^2);
                if (err_temp < err_T)
                    c = ci;
                    e = ei;
                    f = fi;
                    err_T = err_temp;
                end
            end
        end
    end
    val = ("For T =" + T(j) + " years, the parameters c,e and f are:");
    disp(val)
    param_vect = [c, e, f, err_T];
    disp(param_vect(1:3))
    DDF_param = [DDF_param param_vect'];
end
disp('Residual sum square for each interpolation (T = 10,40,100):')
disp(DDF_param(4,:))

%% let's plot all this

col = [    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

figure
for j = [1:3]
    hi_gum = H_Gum(j,:);
    x_durations = linspace(0,100,200);
    parameters = DDF_param(:,j)';
    %disp(parameters);
    h_DDF = (parameters(1) .* x_durations)./(x_durations.^parameters(2) + parameters(3));
    str = ("DDF, Return period T = " + T(j) + " [y]");
    plot(x_durations', h_DDF, 'LineWidth', 1.5, 'color', col(j,:), 'DisplayName',str);
    hold on;
    str2 = ("Gumbel values, duration T = " + T(j) + " [y]");
    oh = plot(D, hi_gum, 'o', 'color', col(j,:), 'DisplayName',str2);
    %oh.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on;
end

grid on;
title('DDF curves for each return period given');
xlabel('Durations [hours]')
ylabel('Rainfall depth [mm]')
%legend('DDF interpolation','values from Gumbel approximation')
legend('Location', 'South East', 'Fontsize', 11)

