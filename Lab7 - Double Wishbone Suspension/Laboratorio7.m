
close all
clear all
clc

%% Data

teta1_min = -9;
teta1_max = -1;
teta2_min = 5;
teta2_max = 15;
teta_lbounds = [teta1_min, teta2_min]*pi*((180)^(-1)); % Lower bounds of theta1 and theta2 and converted from degrees to radians
teta_ubounds = [teta1_max, teta2_max]*pi*((180)^-1);
theta1 = linspace(teta1_min, teta1_max, 20)*pi*((180)^-1); % generates 20 linearly spaced points between teta1_min and teta1_max and converted from degrees to radians
theta2 = linspace(teta2_min, teta2_max, 20)*pi*((180)^-1);

DV = combvec(theta1, theta2); % generate all combinations of elements from the vectors theta1 and theta2
DV = DV'; % compute the matrix transpose of DV

x0 = DV;

h1 = 0.190;     % [m]
h2 = 0.270;     % [m]
g = 0.050;      % [m]
c = 0.870;      % [m]
dz = 0.05;      % [m]

% Target values that we want to reach
target_f1 = 0.100;      % [m]
target_f2 = 2*pi/180;   % [rad]



%% Genetic Algorithm:
par = 0.03;  % probability of mutation
bit = 8;

niteration = 200;

x = x0;
gen = 0;

for ii = 1:niteration

        % compute the fitness of the parents
        [F, matr_rank, obfun1, obfun2] = fitness(x, h1, h2, g, c, dz, target_f1, target_f2);
        
        gen = gen + 1
        
        if gen == 1             
            xinitial = DV;
            f1initial = obfun1;
            f2initial = obfun2;
        end
        
        % Function
        xparent = matr_rank(:, [1, 2]);
        f1parent = matr_rank(:, 3);
        f2parent = matr_rank(:, 4);
        check = find(F == 1);  % If each element of the vector F (i.e. each combination of the current generation) is equal to 1, means that we reached the optimal design variables combination
        NIR1(ii) = length(check);
        
        if F == 1       % stop condition
             break
        end
        
        % Population binary coding
        b = ga_coding(xparent, teta_lbounds, teta_ubounds, bit);
        
        % Crossover
        xbf_cross = ga_crossover(b, F, bit);
        
        % Mutation 
        xm = ga_mutation(xbf_cross, par, bit);
        
        % Get back from binary coding to decimal coding
        xson = ga_decoding(xbf_cross, teta_lbounds, teta_ubounds, bit);
        
        % Compute the fitness of the offsprings
        [Fson, matr_rank_son, obfun1, obfun2] = fitness(xson, h1, h2, g, c, dz, target_f1, target_f2);     
        
        % Discard half population (select the best individuals)
        [x_best, pop_ord] = ga_discard(xparent, F, xson, Fson);
        
        x = x_best;
              
end



%% Plotting the results

% Design Variables space:
figure
plot(x(:,1), x(:,2),'ro')
hold on
plot(x0(:,1), x0(:,2),'bo')
grid on
xlabel('X_1')
ylabel('X_2')
title('Design Variable Space')
legend('Pareto-Optimal points','Initial values')


% Objective Functions Space:
figure
plot(f1initial, f2initial,'bo')
hold on
plot(f1parent, f2parent,'ro')
grid on
ylabel('f2')
xlabel('f1')
legend('Initial values','Pareto-optimal points derived from GA')
title('Objective Functions Space')


% Rank behaviour:
figure
bar(1:gen, NIR1/length(F), 'b')
ylabel('N� of individuals with rank 1 over N� of total individuals')
xlabel('Generations')
grid on
title('Rank behaviour')



%% function definition

%% -------------- Discard (Best Half) --------------

% With the function ga_discard, the offspring (sons) set of combinations is
% appended to the parents set of combinations, generating a population of
% parents and sons, then, each person (combination) of the matrix is
% ordered in discendent way base on the fitness value and finally, only the
% best half of the population is used for the next generation, the other
% half is discarded
function [x_best, pop_ord] = ga_discard(x_dis, F_dis, xson_dis, Fson_dis)

% The matrix pop is composed by three columns, the first two are the two
% design variables, instead the third one is the correspective fitness
% value, the number of rows is obtained appending the sons combinations
% after the parents combinations:
pop = [x_dis,     F_dis;      % Parents combinations  
       xson_dis,  Fson_dis];  % Sons combinations

pop_ord = sortrows(pop, 3, 'descend');

len = size(x_dis, 1);
x_best = pop_ord(1:len, [1, 2]);  % Here, the wrost half combinations is removed from the matrix pop_ord

end




%% ----------- Fitness function --------------

% The function fitness is used to assign the fitness value to each
% combination of the determined generation, this fitness is simply computed
% as the inverse of the rank number assigned to each combination after the
% computation of the Pareto optimal set
function [F, matr_rank, obfun1, obfun2] = fitness(x, h1, h2, g, c, dz, target_f1, target_f2)

for i = 1:length(x)
    theta1 = x(i,1);
    theta2 = x(i,2);

    d(i,1) = h2/(tan(theta1 + theta2) - tan(theta1));
    h(i,1) = h1 - d(i,1)*tan(theta1);

    RCH(i,1) = c*h(i,1)/(d(i,1) + g);
    dgamma(i,1) = dz/(d(i,1) + g);

    obfun1(i,1) = ((RCH(i,1) - target_f1)^2)/target_f1^2;    % Evaluate the first objective function for each combination of design variables
    obfun2(i,1) = ((dgamma(i,1) - target_f2)^2)/target_f2^2; % Evaluate the second objective function for each combination of design variables
end

flag = zeros(length(x), 1);     % Preallocate the memory for computational efficiency
matr = [x(:,1), x(:,2), obfun1, obfun2, flag];
matr_sort = sortrows(matr, 3);  % sortrows(matr, 3) sorts matr in ascending order based on the elements of the third column
matr_rank = [];                 % Preallocate memory for matr_rank
z = 0;

while size(matr_sort(:,1)) > 0  % In this while loop, the matrix matr_sort decrease his dimension removing the optimal Pareto set computed for each loop

%% ---------- Pareto Optimal set finder ----------
%
%
    for i = 1:size(matr_sort, 1)
        for j = i:size(matr_sort, 1)
            if matr_sort(j,4) > matr_sort(i,4)
               matr_sort(j,5) = 1;
            end
        end
    end
    % Now the last column of matr_sort has many ones and zeros, the zeros
    % represent the Pareto optimal set
    index = find(matr_sort(:,5) == 0); % Indexes of the Pareto optimal set. Note that here it finds the index of all elements equal to zero in the last column of matr_sort
%
%
%% -----------------------------------------------
    rank = z + ones(length(index), 1);
    % In the following it is assigned a rank number to the curent Pareto optimal set of this loop
    matr_rank = [matr_rank; matr_sort(index,1), matr_sort(index,2), matr_sort(index,3), matr_sort(index,4), rank]; 
    z = z + 1;                % Increase the rank number
    matr_sort(index,:) = [];  % Remove from matr_sort the Pareto optimal set
    matr_sort(:,5) = 0;       % Reset to zero again the last column of matr_sort
end

% Here, for each combination of design variables it is assigned a finess
% value based on the inverse of the rank number:
F = 1./matr_rank(:,5);        % fitness function definition 

end


