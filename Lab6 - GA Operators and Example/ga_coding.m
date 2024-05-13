function b = ga_coding(x, xmin, xmax, bit)

% This function performs the binary coding and concatenates the variables (Levi)

% INPUT:
%   x: matrix of individuals of population (design variables stored in columns) [x1, x2, ...]
%   xmin: vector of lower bounds [1, ndv]
%   xmax: vector of upper bounds [1, ndv]
%   bit: number of bits

% OUTPUT:
%   b: matrix of individuals of population in binary code [b1; b2; ...] 


[n_pop, n_var] = size(x);
% Note:
% n_pop = number of population
% n_var = number of variables

% normalization
X = round((x - repmat(xmin, n_pop, 1))./repmat(xmax - xmin, n_pop, 1)*(2^(bit) - 1)); 
% Note:
% Y = round(X) rounds each element of X to the nearest integer
% B = repmat(xmax - xmin, n_pop, 1) returns an array containing n_pop copies of [xmax - xmin] in the row and column dimensions. The size of B is size([xmax - xmin])*n_pop when [xmax - xmin] is a matrix.

% converte in binario
B = dec2bin(X, bit);

% concatenates
b = [];

for i = 1:n_var
    b = [b, B((i - 1)*n_pop + 1:i*n_pop, :)];
end

