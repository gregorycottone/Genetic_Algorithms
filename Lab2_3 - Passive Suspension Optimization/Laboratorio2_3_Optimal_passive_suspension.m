
clear all
close all
clc


%% Lab 2

Ab = 1.4e-5;                                  % m
v = 30;                                       % m/s vehicle velocity
m1 = 30;                                      % kg wheel, unsprung mass
m2 = 230;                                     % kg quarter car vehicle model, sprung mass
k1 = 120000;                                  % N/m tyre radial stiffness
A = sqrt(Ab*v*0.5);
k2min = 1;                                 % N/m suspension minimum stiffness
r2min = 1;                                  % Ns/m suspension minimum dumping
k2max = 50000;                                % N/m suspension maximum stiffness
r2max = 5000;                                 % Ns/m suspension maximum dumping
n = 100;                                    % number of elements
k2 = linspace(k2min, k2max, n);             % design variable vector NB:k2 and r2 must have the same length!!
r2 = linspace(r2min, r2max, n);             % design variable vector
DV = combvec(k2,r2);                          % all possible combinations of k2 and r2
DV = DV';
N = length(DV);                               % number of combinations



%% Discomfort and Road holding

% For each combination of k2 and r2 we compute the corrispective values of
% the Discomfort and Road holding

tic

for i = 1:length(DV(:,1))
    x = DV(i,:);
    Dis(i,1) = Discomfort(x, A, m1, m2, k1);
    Rh(i,1) = Road_holding(x, A, m1, m2, k1);
end

all_obj_values = toc;



%% Big matrix generation

tic

tabel = [DV, Dis, Rh, zeros(N,1)];  % the last column of zeros is to preallocate the memory for computational efficiency
bigmat = sortrows(tabel, 3);        % sort all the rows in an ascendent way with respect to Dis value


%% Pareto optimal set finder

for ii = 1:N
    for jj = ii:N
        if bigmat(jj, 4) > bigmat(ii, 4)
            bigmat(jj, 5) = 1;
        end
    end
end

flag = bigmat(:,5);
ind = find(flag == 0);  % find the index of all elements equal to zero in the last column of bigmat

pareto_toc = toc;
first_case_toc = all_obj_values + pareto_toc


%% Figures

figure
plot(bigmat(:,3), bigmat(:,4), 'b.')
hold on
grid on
plot(bigmat(ind,3), bigmat(ind,4), 'ro')
ylabel('Road Holding')
xlabel('Discomfort')
legend('All the possible configurations', 'Pareto-optimal set')
title('Objective functions space')

figure
plot(bigmat(:,1), bigmat(:,2), 'b.')
hold on
grid on
plot(bigmat(ind,1), bigmat(ind,2), 'ro')
ylabel('Damping r2')
xlabel('Stiffness k2')
legend('All the possible configurations', 'Pareto-optimal set')
title('Design variables space')



%% Lab 3 - wheight method

tic

Disbigmat = bigmat(ind, 3);  % all the Pareto points of the Discomfort objective function 
Rhbigmat = bigmat(ind, 4);   % all the Pareto points of the Road holding objective function

% ATTENZIONE !!! 
% Chiarire se nella normalizzazione ci vuole il massimo
% valore delle funzioni obbiettivo o il massimo valore delle stesse ma nel
% set di Optimal Pareto

% maxDis = max(Disbigmat);  % finding the maximum of discomfort Pareto set optimal points
% maxRh = max(Rhbigmat);    % finding the maximum of road holding Pareto set optimal points
maxDis = max(Dis);          % finding the maximum of discomfort function 
maxRh = max(Rh);            % finding the maximum of road holding function
x0 = [2000, 500];           % defining the initial guess for the iterative cycle - IT MUST BE REASONABLE
x = x0;
lambdas = linspace(0, 1, 1000); % definition of the lambda vecotor --> equals to the number of points that weight method will outcome 


for iii = 1:length(lambdas)
    lambda1 = lambdas(iii);
    lambda2 = 1 - lambda1;    % Because the sum of lambads should be always 1
    x = fminsearch(@(x) funweight(x, lambda1, lambda2, maxDis, maxRh, A, m1, m2, k1), x0);  
    x_weight(iii,:) = x;       % matrix that accumulates the Pareto optimal points computed cycle by cycle
    x0 = x;                    % every cycle i choose to assign the initial guess of the next cycle - IMPORTANT

    % Here, the correspective values of Discomfort and Road holding is
    % computed for the weight method solution in each iteration 
    Disf(iii) = Discomfort(x, A, m1, m2, k1);
    Rhf(iii) = Road_holding(x, A, m1, m2, k1);
end

weight_toc = toc;
second_case_toc = all_obj_values + weight_toc


figure
plot(x_weight(:,1), x_weight(:,2), 'b.', bigmat(ind,1), bigmat(ind,2), 'ro')
grid on
ylabel('Damping r2')
xlabel('Stiffness k2')
title('Design variables space')
legend('Weighted method set','Pareto-optimal set')


figure
plot(Disf, Rhf, 'b.', bigmat(ind,3), bigmat(ind,4), 'ro')
ylabel('Road Holding')
xlabel('Discomfort')
legend('Weighted method set','Pareto-optimal set')
title('Objective functions space')
grid on


%% Lab 3 - constraint method

tic 

% Also here it is not clear if I should use the Pareto set or the set of
% Road holding solutions:
% minRh = min(Rhbigmat);                % finding the minimum of Road holding Pareto set optimal points
minRh = min(Rh);                        % finding the minimum of road holding function

epsilon = linspace(minRh, maxRh, 1000);  % creation of a vector of the epsilon that the iterative cycle must slide --> equals to the number of points that constraints method will outcome 
xx0 =  x0;

for jjj = 1:length(epsilon)
    options = optimset('Algorithm','interior-point','Display','off');
%     options = optimset('Algorithm','sqp','Display','off');
    xx = fmincon(@(xx) Discomfort(xx,A,m1,m2,k1), xx0, [], [], [], [], [k2min,r2min], [k2max,r2max], @(xx) Road_holding_eps(xx,epsilon(jjj),A,m1,m2,k1), options);
    constraintmet(jjj,:) = xx; % matrix that accumulates the pareto optimal points computed cycle by cycle

    % Here, the correspective values of Discomfort and Road holding is
    % computed for the constraint method solution in each iteration 
    Discons(jjj) = Discomfort(xx, A, m1, m2, k1);
    Rhcons(jjj) = Road_holding(xx, A, m1, m2, k1);
end

constr_toc = toc;
third_case_toc = all_obj_values + constr_toc


figure
plot(constraintmet(:,1), constraintmet(:,2), 'b.', bigmat(ind,1), bigmat(ind,2), 'ro')
grid on
ylabel('Damping r2')
xlabel('Stiffness k2')
title('Design variables space')
legend('Constraints method set','Pareto-optimal set')

figure
plot(Discons, Rhcons, 'b.', bigmat(ind,3), bigmat(ind,4), 'ro')
ylabel('Road Holding')
xlabel('Discomfort')
legend('Constraints method set','Pareto-optimal set')
title('Objective functions space')
grid on



%% Analytical solution

% defining the symbolic variables

syms k2 r2 real

Dis = A*sqrt(((m1+m2)/((m2^2)*r2))*(k2^2) + (k1*r2)/(m2^2));
Rh =  A*sqrt((((m1+m2)^3)/((m2^2)*r2))*(k2^2) - ((2*m1*k1*(m1+m2))/(m2*r2))*k2 + (k1*r2*(m1+m2)^2)/(m2^2) + ((k1^2) * m1)/(r2));

dRh_dk2 = diff(Rh,k2);    % partial derivative of Road holding with respect to k2
dRh_dr2 = diff(Rh,r2);    % partial derivative of Road holding with respect to r2
dDis_dk2 = diff(Dis,k2);  % partial derivative of Discomfort with respect to k2
dDis_dr2 = diff(Dis,r2);  % partial derivative of Discomfort with respect to r2

pareto = simplify(dDis_dk2*dRh_dr2 - dDis_dr2*dRh_dk2);  % This is a function of k2 and r2
DVspacek2 = solve(pareto, k2);           % Solve the previous equation with respect to k2, note that in output there are two possible solutions, the correct one is the first


figure
ezplot(DVspacek2(1), [0, 3000])  % limiting r2 up to 3000 for a better plotting proportionality
title('Desing variables - Analytical solution')
xlabel('r2 [N/m]')
ylabel('k2 [Ns/m]')




Disr2 = simplify(subs(Dis, k2, DVspacek2(1)));  % Dis as a function of only r2 --> Dis = f(r2)
Rhr2 = simplify(subs(Rh, k2, DVspacek2(1)));    % Rh as a function of only r2 --> Rh = f(r2)
r2inv = simplify(finverse(Disr2, r2));          % make the inverse of the discomfort --> r2 = f^-1(Dis)
RhfDis = simplify(subs(Rhr2, r2, r2inv));       % substitute r2 = f^-1(Dis) inside Rh = f(r2) in order to obtain Rh = f(Dis)


figure
ezplot(RhfDis, [0.0005, 1.4])                   % limiting the discomfort in particular the lower bond must be greater than 0 beacause it would tend to infinite
title('Objective functions - Analytical solutions')
xlabel('Discomfort')
ylabel('Road holding')


%% Minimum of the Road holding

% minimum for Discomfort (when k2 = r2 = 0) implies a road holdind that tends
% to infinite, it has no sense compute the minimum for discomfort, it is for 0 values of DVs 
minRhdfk2 = solve(simplify(dRh_dk2), r2);  % partial derivative of Road holding with respect to k2 = 0 --> expressed as r2 = f(k2)
minRhdfr2 = solve(simplify(dRh_dr2), r2);  % partial derivative of Road holding with respect to r2 = 0 --> expressed as r2 = f(k2)
% we find the places where the minimum (or maximum) are located of the two
% partial derivatives 
% the partial derivative of Road holding with respect to k2 has no minimum(or maximum), instead the other shows a place of minimum points where the Road holding fucntion must be re-evaluated 
Rh1 = simplify(subs(Rh, r2, minRhdfr2(2)));  % substitution of minimum points found before -- > Rh = f(k2)


figure
ezplot(Rh1, [0, k2max])  % here are plotted the minimum points of Road holding: here we can estimate the minimum value of Rh with respect to k2
title('Road holding - local minimum point')
xlabel('k2')
ylabel('Road holding')


k2 = 10:100:k2max;
xxx = eval(Rh1);  % evaluation of all the points once defined k2
[minRh, indk2] = min(xxx);
maxk2 = k2(indk2);

r2 = 10:10:r2max;
yyy = eval(simplify(DVspacek2(1)));      % computation of the k2 vector as a function of r2: pareto optimal set in DV space
jj = find(abs((yyy/maxk2) - 1)<= 0.01);  % finding the value of r2 that gives back k2max that minimizes Rh

figure
hold on
ezplot(DVspacek2(1), [0, 3000])
plot(r2(jj(end)), yyy(jj(end)), 'ro', 0, 0, 'ro')
title('Desing variables - Analytical solution')
xlabel('r2 [N/m]')
ylabel('k2 [Ns/m]')
hold off

zzz = eval(Rhr2);   % evaluation of the Rh expressed as Rh = f(r2)
eee = eval(Disr2);  % evaluation of the Dis expressed as Dis = f(r2)
ii = find(abs((zzz/minRh) - 1) <= 0.00001);  % finding the value of r2 that gives back minRh that minimizes Rh

figure
hold on
ezplot(RhfDis, [0.0005, 1.4]) 
plot(eee(ii(end)), minRh, 'ro')
title('Objective functions - Analytical solutions')
xlabel('Discomfort')
ylabel('Road holding')
hold off


%% Functions

function phi = funweight(x, lambda1, lambda2, maxDis, maxRh, A, m1, m2, k1)
Disf = Discomfort(x, A, m1, m2, k1);
Rhf = Road_holding(x, A, m1, m2, k1);
phi = lambda1*(Disf/maxDis) + lambda2*(Rhf/maxRh); % Weighted sum of Discomfort and Road holding functions
end


function [Rh, ceq] = Road_holding_eps(x, epsilon, A, m1, m2, k1)
Rh = Road_holding(x, A, m1, m2, k1) - epsilon;
ceq = [];
end


function Dis = Discomfort(x, A, m1, m2, k1)
Dis = A*sqrt(((m1+m2)/((m2^2)*x(2)))*(x(1)^2) + (k1*x(2))/(m2^2));
end


function Rh = Road_holding(x, A, m1, m2, k1)
Rh =  A*sqrt((((m1+m2)^3)/((m2^2)*x(2)))*(x(1)^2) - ((2*m1*k1*(m1+m2))/(m2*x(2)))*x(1) + (k1*x(2)*(m1+m2)^2)/(m2^2) + ((k1^2) * m1)/(x(2)));
end


