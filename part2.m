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
area(x,IUHw,'FaceAlpha', 0.5)
title("Watershed IUH")
xlabel("Time [h]")
ylabel("Proportion of unit discharge [mm/h]") %Important to understand that the values are per hour !!

disp("Surface under curve of IUHw: " + sum(IUHw)*dt); %Sum of all the small rectangles 

%% Channel IUH

c = 0.3*3600; %m/s * s/h ;celerity
D = 10^6; %m^2/h ;hydrodynamic dispersion
L = 7 * 1000; %km -> m

t = dt:dt:20;
IUHc = L./(sqrt(4*pi*D).*t.^(3/2)) .* exp((-(L-c*t).^2)./(4*D.*t));
IUHc = [0, IUHc]; % just because can't divide by zero 

figure 
area(IUHc,'FaceColor', '#D95319','FaceAlpha', 0.5)
title("Channel IUH")
xlabel("Time [h]")
ylabel("Proportion of unit discharge [mm/h]")

disp("Surface under curve of IUHc: " + sum(IUHc)*dt);
%% Plot the two iuh in same figure, focus early part

figure
IUHcx = [IUHc, zeros(1,length(x)-length(IUHc))];

area(x, IUHcx, 'FaceAlpha', 0.5)
hold on
area(x, IUHw, 'FaceAlpha', 0.5)
xlabel("Time [hour]")
ylabel("??")
title("Watershed and channel IUH")
legend("IUHw", "IUHc")

%% explicit convolution IUHw and Je
Qw = [];
Qc = [];
for event = [1,2,3]
    Jei = Je(:,event);

    N = length(IUHw);
    M = length(Je);

    Qwi = zeros([1, N+M]);

    %%convolution
    for n = 1:N+M
        for m = 1:M
            %disp("n :" + n + " and m =" + m);
            if ((n>m) && n-m+1<N)
                Qwi(n) = Qwi(n) + Jei(m)*dt*IUHw(n-m+1);
            end
        end
    end

    SJei = sum(Jei); %% la quantité total d'eau tombé
    SJei_waited = sum(IUHw)*dt*SJei; %%ce qu'on devrait attendre au vu de l'approximation du IUHw
    SQw = sum(Qwi); %%total de pluie qu'on obtient à la sortie
    disp("A)   With a IUHw surface of " + sum(IUHw)*dt + " we should get a total water of " + SJei_waited + " but by the convolution we have : " + SQw + " !")
    disp("     Less than " + round(1-SQw/SJei_waited, 5)*100 + "% is lost through this approximations")

    figure 
    
    subplot(2,1,1)
    
    area(Qwi,'FaceColor', '#0072BD' ,'FaceAlpha', 0.5)
    hold on
    area(Jei,'FaceColor', '#7E2F8E', 'FaceAlpha', 0.7)
    hold on
    area(IUHw,'FaceColor','#EDB120', 'FaceAlpha', 0.8)
    legend("Qw", "Jei", "IUHw")
    title("Convolution of IUH and Je of event " + event)

    %%% Channel IUH convolution

    N = length(IUHc);
    M = length(Qwi);

    Qci = zeros([1, N+M]);

    for n = 1:N+M
        for m = 1:M
            %disp("n :" + n + " and m =" + m);
            if ((n>m) && n-m+1<N)
                Qci(n) = Qci(n) + Qwi(m)*dt*IUHc(n-m+1);
            end
        end
    end

    SQc = sum(Qci);
    disp("B)   Surface under Qw :" + SQw + " and under Qc :" + SQc)
    disp("     This is also because surface under IUHc = " + sum(IUHc)*dt)
    disp("----------")
    %%%% convolution 2 plot to verify 

    subplot(2,1,2) 
    area(Qci,'FaceColor', '#D95319' ,'FaceAlpha', 0.5)
    hold on
    area(Qwi,'FaceColor', '#0072BD' ,'FaceAlpha', 0.5)
    hold on
    area(IUHc,'FaceColor','#EDB120', 'FaceAlpha', 0.8)
    legend("Qc", "Qw", "IUHc")
    title("Convolution of IUHc and Qw of event " + event)
    
    %%%saving final data we need
    Qw = [Qw Qwi'];
    Qc = [Qc Qci'];
    
end

%% test figure Fa PT reversed double plot
for event = [1,2,3]
    figure
    %plot(Fa, '-o')
    yyaxis left
    area(Je(:,event),'FaceColor', '#7E2F8E', 'FaceAlpha', 0.7)
    ylabel('Effectiv Precipitation [mm/h]')
    ylim([0,max(Je(:,event))+3])
    ax = gca;
    ax.YDir = "reverse";
    ax.YColor = '#7E2F8E';
    top = max(Je(:,event));
    ylim(ax, [0,top+top/2])



    yyaxis right
    area(Qw(:,event),'FaceColor', '#0072BD' ,'FaceAlpha', 0.5)
    hold on
    area(Qc(:,event),'FaceColor', '#D95319' ,'FaceAlpha', 0.5)
    ylabel('Q [mm/h]')
    ax.YColor = '#0072BD';

    xlabel("Time [hour]")
    legend(["Je", "Qw", "Qc"])
    title("Convolution of IUHc and Qw of event " + event)
end

%% channel maxima's

MaxC = [max(Qc); 0 0 0];
for i = 1:length(Qc)
   if Qc(i,1) == MaxC(1,1)
       MaxC(2,1) = i*dt;
       disp("Maximal Q of channel for event 1:")
       disp(MaxC(1,1) + "[mm] after " + MaxC(2,1) + " [hours] ")

   end
   if Qc(i,2) == MaxC(1,2)
       MaxC(2,2) = i*dt;
       disp("Maximal Q of channel for event 2: ")
       disp(MaxC(1,2) + "[mm] after " + MaxC(2,2) + " [hours]")
   end
   if Qc(i,3) == MaxC(1,3)
       MaxC(2,3) = i*dt;
       disp("Maximal Q of channel for event 3: ")
       disp(MaxC(1,3) + "[mm] after " + MaxC(2,3) + " [hours]")
   end
   
end

%% saving results

save("output_part2.mat", "dt", "IUHw") 
