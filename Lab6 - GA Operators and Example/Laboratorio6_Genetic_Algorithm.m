
clear all
close all
clc

%% Crossover only

xmin = [0, 0];
xmax = [5, 5];
n = 100;
x0 = 5*rand(n, 2);   % We generate an initial population of n combinations, where each combination is characterised by two numbers between 0 and 5
bit = 8;
niteration = 200;
x_parents = x0;


for ii = 1:niteration

        % compute the fitness of the parents
        F_parents = 100 - (x_parents(:,1).^2 + x_parents(:,2).^2);

        % population binary coding
        bit_parents = ga_coding(x_parents, xmin, xmax, bit);

        % crossover
        bit_sons = ga_crossover(bit_parents, F_parents, bit);

        % get back from binary coding to decimal coding
        x_sons = ga_decoding(bit_sons, xmin, xmax, bit);

        % compute the fitness of the offsprings (prole)
        F_sons = 100 - (x_sons(:,1).^2 + x_sons(:,2).^2);

        av_F(:,ii) = mean(F_sons);
        
        x_parents = x_sons; 
        x1_animation(:,ii) = x_parents(:,1);
        x2_animation(:,ii) = x_parents(:,2);
end


figure
pause on 
for i = 1:length(x0(:,1))
    x1 = x1_animation(:,i);
    x2 = x2_animation(:,i);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k','MarkerSize', 0.01)
    hold on  
    plot(x1, x2, 'ro')
    axis_limits = [0 5 0 5];
    axis manual; 
    axis equal; 
    axis(axis_limits);
    hold off    
    drawnow;       
    pause off
end
hold on
plot(x0(:,1), x0(:,2), 'bo')
grid on
xlabel('X_1')
ylabel('X_2')
hold off
title('Design Variable Space - only Reproduction')
legend('','Last points','Initial values')


figure
plot(1:niteration, av_F, '.')
grid on
xlabel('Iterations')
ylabel('Average Fitness')
title('Convergence diagram - Only reproduction')



%% Crossover and mutation

xmin = [0, 0];
xmax = [5, 5];
n = 100;
x0 = 5*rand(n, 2);

par = 0.2;
bit = 8;

niteration = 200;
xx_parents = x0;

for ii = 1:niteration

        % compute the fitness of the parents
        FF_parents = 100 - (xx_parents(:,1).^2 + xx_parents(:,2).^2);
        
        % population binary coding
        bbit_parents = ga_coding(xx_parents, xmin, xmax, bit);
        
        % crossover
        bbit_sons = ga_crossover(bbit_parents, FF_parents, bit);
        
        % mutation 
        mutated_sons = ga_mutation(bbit_sons, par, bit);
        
        % get back from binary coding to decimal coding
        xx_sons = ga_decoding(mutated_sons, xmin, xmax, bit);
        
        % compute the fitness of the offsprings
        FF_sons = 100 - (xx_sons(:,1).^2 + xx_sons(:,2).^2);     
        
        av_FF(:,ii) = mean(FF_sons);
        xx_parents = xx_sons; 

        xx1_animation(:,ii) = xx_parents(:,1);
        xx2_animation(:,ii) = xx_parents(:,2);
end


figure
pause on 
for i = 1:length(x0(:,1))
    x1 = xx1_animation(:,i);
    x2 = xx2_animation(:,i);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k','MarkerSize', 0.01)
    hold on  
    plot(x1, x2, 'ro')
    axis_limits = [0 5 0 5];
    axis manual; 
    axis equal; 
    axis(axis_limits);
    hold off    
    drawnow;       
    pause off
end
hold on
plot(x0(:,1), x0(:,2), 'bo')
grid on
xlabel('X_1')
ylabel('X_2')
hold off
title('Design Variable Space - Reproduction and mutation')
legend('','Last points','Initial values')


figure
plot(1:niteration, av_FF, '.')
grid on
xlabel('Iterations')
ylabel('Average Fitness')
title('Convergence diagram - Reproduction and Mutation')

% In this case the mutation algorithm worsens the result simply because 
% the problem has no local minima



%% Crossover, mutation and discard

xmin = [0, 0];
xmax = [5, 5];
n = 100;
x0 = 5*rand(n, 2);

par = 0;  % We keep this parameter for the mutation null since the problem has not local minima
bit = 8;

niteration = 200;   

xxx_parents = x0;

for ii = 1:niteration

        FFF_parents = 100 - (xxx_parents(:,1).^2 + xxx_parents(:,2).^2);

        bbbit_parents = ga_coding(xxx_parents, xmin, xmax, bit);

        bbbit_sons = ga_crossover(bbbit_parents, FFF_parents, bit);
   
        mmutated_sons = ga_mutation(bbbit_sons, par, bit);
        
        xxx_sons = ga_decoding(mmutated_sons, xmin, xmax, bit);
        
        FFF_sons = 100 - (xxx_sons(:,1).^2 + xxx_sons(:,2).^2);     
        
        % Discard half population (select the best individuals)
        x_best = ga_discard(xxx_parents, FFF_parents, xxx_sons, FFF_sons);
        
        av_FFF(:,ii) = mean(FFF_sons);
        
        xxx_parents = x_best;

        xxx1_animation(:,ii) = xxx_parents(:,1);
        xxx2_animation(:,ii) = xxx_parents(:,2);
end


figure
pause on 
for i = 1:length(x0(:,1))
    x1 = xxx1_animation(:,i);
    x2 = xxx2_animation(:,i);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k','MarkerSize', 0.01)
    hold on  
    plot(x1, x2, 'ro')
    axis_limits = [0 5 0 5];
    axis manual; 
    axis equal; 
    axis(axis_limits);
    hold off    
    drawnow;       
    pause off
end
hold on
plot(x0(:,1), x0(:,2), 'bo')
grid on
xlabel('X_1')
ylabel('X_2')
hold off
title('Design Variable Space - discard BEST HALF')
legend('','Last points','Initial values')


figure
plot(1:niteration, av_FFF, '.')
grid on
xlabel('Iterations')
ylabel('Average Fitness')
title('Convergence diagram - discard BEST HALF')



%% Function definition:

% BEST HALF METHOD:
% I make a matrix with the parents individual x_parents and its sons x_sons
% and the corresponding fitness values of both parents and sons, then I will 
% come up with a matrix whose rows are double with respect to the single vector 
% x_parents (I'm joining x_parents and x_sons), then I order the individuals 
% from largest to smallest of the fitness value and I only keep the first half good.

function [x_best, pop_ord] = ga_discard(x_dis, F_dis, xson_dis, Fson_dis)

pop = [x_dis,     F_dis; 
       xson_dis,  Fson_dis];

pop_ord = sortrows(pop, 3, 'descend');

len = size(x_dis, 1);
x_best = pop_ord(1:len, [1, 2]);

end


