function x = ga_decoding(b, xmin, xmax, bit)


% This function performs the decodification from binary code and separates the variables (Levi) 

% INPUT:
%   b: matrix of individuals in binary code [b1; b2; ...] 
%   xmin: vector of lower bounds [1, ndv]
%   xmax: vector of upper bounds [1, ndv]
%   bit: number of bits

% OUTPUT:
%   x: matrix of individuals (design variables are stored in columns) [x1, x2, ...]


n_var = length(xmin);
n_pop = size(b, 1);

% Separates the variables and performs the decoding
for i = 1:n_var
    X(:,i) = bin2dec(b(:,(i-1)*bit+1:i*bit));
end


% denormalization
x = X./(2^bit - 1).*repmat(xmax - xmin, n_pop, 1) + repmat(xmin, n_pop, 1);



