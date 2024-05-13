
clear all
close all
clc


%% Lab 4

Ab = 1.4e-5;                      % m
v = 30;                           % m/s vehicle velocity
m1 = 30;                          % kg wheel, unsprung mass
m2 = 230;                         % kg quarter car vehicle model, sprung mass
k1 = 120000;                      % N/m tyre radial stiffness
A = sqrt(Ab*v*0.5);
g1min = 0;                        % N/m suspension maximum stiffness
g2min = 0;                        % Ns/m suspension maximum dumping
g1max = 400000;                   % N/m suspension maximum stiffness
g2max = 100000;                   % Ns/m suspension maximum dumping
n = 100;                          % number of elements
g1 = linspace(g1min, g1max, n);
g2 = linspace(g2min, g2max, n);
DV = combvec(g1, g2);             % all possible combinations of k2 and r2
DV = DV';
N = length(DV);                   % number of combinations



%% Discomfort and Road holding

% For each combination of k2 and r2 we compute the corrispective values of
% the Discomfort and Road holding

tic

for i = 1:length(DV(:,1))
    x = DV(i,:);
    Dis(i,1) = Discomfort(x, A, m2, k1);
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
ylabel('g1')
xlabel('g2')
legend('All the possible configurations', 'Pareto-optimal set')
title('Design variables space')



%% Functions


function Dis = Discomfort(x, A, m2, k1)
Dis = A*sqrt((k1*x(1))/(m2*x(2)));
end


function Rh = Road_holding(x, A, m1, m2, k1)
Rh = A*sqrt(((x(1)*k1*((m1+m2)^2))/(m2*x(2))) - (((k1^2)*(2*m1+m2))/(x(2))) + ((k1^3)*m2)/(x(1)*x(2)) + (((k1^2)*m1*x(2))/(x(1)*m2)));
end


