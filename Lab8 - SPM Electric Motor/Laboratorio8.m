
close all
clear all
clc

%% Laboratorio 8:
%% Step 1: Dati


h_m = linspace(3.315*1e-3, 11.715*1e-3, 20);
b_m = linspace(28*1e-3, 81*1e-3, 20);
N_pc = [0, 1, 2, 3];

x_lbounds = [h_m(1), b_m(1), N_pc(1)];
x_ubounds = [h_m(end), b_m(end), N_pc(end)];


DV = combvec(h_m, b_m, N_pc);
DV = DV';

x0 = DV;


target_f1 = 1;
target_f2 = 0;   % [m^3]



%% Step2: Genetic Algorithm:

par = 0.03;         % probability of mutation
bit = [8, 8, 2];

niteration = 100;

x = x0;
gen = 0;

for ii = 1:niteration

        % compute the fitness of the parents
        [F, matr_rank, obfun1, obfun2] = fitness(x, target_f1, target_f2);
        
        gen = gen + 1
        
        if gen == 1             
            xinitial = DV;
            f1initial = obfun1;
            f2initial = obfun2;
        end
        
        % Function
        xparent = matr_rank(:, [1, 2, 3]);
        f1parent = matr_rank(:, 4);
        f2parent = matr_rank(:, 5);
        check = find(F == 1);  % If each element of the vector F (i.e. each combination of the current generation) is equal to 1, means that we reached the optimal design variables combination
        NIR1(ii) = length(check);
        
        if F == 1       % stop condition
             break
        end
        
        % population binary coding
        b = ag_coding_withDiscrete(xparent, x_lbounds, x_ubounds, bit);
        
        % crossover
        xbf_cross = ag_crossover_withDiscrete(b, F);
        
        % mutation 
        xm = ag_mutation_withDiscrete(xbf_cross, par);
        
        % get back from binary coding to decimal coding
        xson = ag_decoding_withDiscrete(xbf_cross, x_lbounds, x_ubounds, bit);
        
        % compute the fitness of the offsprings
        [Fson, matr_rank_son, obfun1, obfun2] = fitness(xson, target_f1, target_f2);     
        
        % Discard (select the best individuals)
        [x_best, pop_ord] = ga_discard(xparent, F, xson, Fson);
        
        x = x_best;
              
end



%% Step 3: Plotting the results

% Design Variables space:
figure
plot3(x0(:,1), x0(:,2), x0(:,3), 'bo')
hold on
plot3(x(:,1), x(:,2), x(:,3), 'ro')
grid on
xlabel('h_m')
ylabel('b_m')
zlabel('N_pc')
title('Design Variable Space')
legend('Initial values', 'Pareto-Optimal points')
hold off


% Objective Functions Space:
figure
plot(f1initial, f2initial, 'bo')
hold on
plot(f1parent, f2parent, 'ro')
grid on
ylabel('Volume')
xlabel('Losses')
legend('Initial values', 'Pareto-optimal points derived from GA')
title('Objective Functions Space')
xlim([0.02, 0.14])
hold off


% Rank behaviour:
figure
bar(1:gen, NIR1/length(F), 'b')
ylabel('N° of individuals with rank 1 over N° of total individuals')
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

pop_ord = sortrows(pop, 4, 'descend');

len = size(x_dis, 1);
x_best = pop_ord(1:len, [1, 2, 3]);  % Here, the wrost half combinations is removed from the matrix pop_ord

end




%% ----------- Fitness function --------------

% The function fitness is used to assign the fitness value to each
% combination of the determined generation, this fitness is simply computed
% as the inverse of the rank number assigned to each combination after the
% computation of the Pareto optimal set
function [F, matr_rank, obfun1, obfun2] = fitness(x, target_f1, target_f2)

for i = 1:length(x)
    h_m = x(i,1);
    b_m = x(i,2);
    N_pc = x(i,3);

    [eta_NP(i,1), Volume(i,1)] = function_SPM_motor(h_m, b_m, N_pc);

    obfun1(i,1) = abs(eta_NP(i,1) - target_f1);    % Evaluate the first objective function for each combination of design variables
    obfun2(i,1) = abs(Volume(i,1) - target_f2);    % Evaluate the second objective function for each combination of design variables
end

flag = zeros(1, length(x))';    % Preallocate the memory for computational efficiency
matr = [x(:,1), x(:,2), x(:,3), obfun1, obfun2, flag];
matr_sort = sortrows(matr, 4);  % sortrows(matr, 4) sorts matr in ascending order based on the elements of the fourth column
matr_rank = [];                 % Preallocate memory for matr_rank
z = 0;

while size(matr_sort(:,1), 1) > 0  % In this while loop, the matrix matr_sort decrease his dimension removing the optimal Pareto set computed for each loop

%% ---------- Pareto Optimal set finder ----------
%
%
    for i = 1:size(matr_sort, 1)
        for j = i:size(matr_sort, 1)
            if matr_sort(j,5) > matr_sort(i,5)
               matr_sort(j,6) = 1;
            end
        end
    end
    % Now the last column of matr_sort has many ones and zeros, the zeros
    % represent the Pareto optimal set
    index = find(matr_sort(:,6) == 0); % Indeces of the Pareto optimal set. Note that here it find the index of all elements equal to zero in the last column of matr_sort
%
%
%% -----------------------------------------------
    rank = z + ones(length(index), 1);
    % In the following it is assigned a rank number to the current Pareto optimal set of this loop
    matr_rank = [matr_rank; matr_sort(index,1), matr_sort(index,2), matr_sort(index,3), matr_sort(index,4), matr_sort(index,5), rank]; 
    z = z + 1;                % Increase the rank number
    matr_sort(index,:) = [];  % Remove from matr_sort the Pareto optimal set
    matr_sort(:,6) = 0;       % Reset to zero again the last column of matr_sort
end

% In the following, we want to penalize the combinations with N_pc = 0,
% because it is unfeasible, to do that, we search combinations with
% N_pc = 0 and we put to that combinations an high value of the rank:
idx_unfeasible = matr_rank(:,3)==0;
matr_rank(idx_unfeasible,6) = 100;

% Here, for each combination of design variables it is assigned a finess
% value based on the inverse of the rank number:
F = 1./matr_rank(:,6);        % fitness function definition 

end


