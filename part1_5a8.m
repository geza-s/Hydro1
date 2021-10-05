
% import data into table T
T = readtable('/Users/annahalloy/Documents/EPFL/BA5/hydrology/ass1/data/data.txt',...
    'HeaderLines', 2,...
    'Format','%s%s%f',... %the format is: text string, text string and float number
    'TreatAsEmpty','-'); %this is how empty data is reported (see legend)

% Create a vector h containing the hourly precipitation and a vector t
% containing the timestamp
h = T.rre150h0;    %rre150h0 is the MeteoSwiss code for hourly rainfall depth [mm]
t = datetime(T.time,'InputFormat','yyyyMMddHH'); %convert to datetime (can be slow)
m = month(t); %gives a value in 1-12 to indicate the month of each date
y = year(t); %gives the year of each date

% fix empty values (which appear as NaN values in the Matlab)
emptyValues = isnan(h); %logical test to tell whether a value is missing or not
h(emptyValues) = 0; %give zero to those values
fprintf('%i empty values\n', sum(emptyValues)); %display how many missing values there are

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


annualMax = zeros(39,6);
maxi = zeros(39,1);
D = [1 3 6 12 24 48]; % durations
%%% boucle sur les années


for k = D
    
        for j = 1981:2019

        %initialisation de la moving window ou toutes les sommes seront
        %stockées pour apres trouver le maximum (voir maxi)
        store = zeros(600000, 1);

        % boucle 'moving window' sur toutes les valeures de h de chaque
        % année
            for i = find((y == j), 1): find((y == j), 1)+length(h(y == j))-k

                %somme sur 3 valeures (3heures) consécutives, garde les
                %valeures sur toute l'année
                store(i) =  sum(h(i:i+k-1));

            end

        maxi(j-1980) = max(store);
        % on garde uniquement le mawimum sur une année -> max de store
        annualMax(j-1980,k) = maxi(j-1980);

        end
end


save('assignment1_output_part1', 'annualMax', 'D')
