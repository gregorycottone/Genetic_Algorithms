function b=ga_coding(x,xmin,xmax,bit);

% b=ga_coding(x,F,bit)
%
% this function performs the binary coding and concatenates the variables (Levi)
%
% input:
%   x: matrix of individuals of population (design variables stored in columns) [x1 x2 ...]
%   xmin: vector of lower bounds [1,ndv]
%   xmax: vector of upper bounds [1,ndv]
%   bit: number of bits
% output:
%   b: matrix of individuals of population in binary code [b1; b2; ...] 

[n_pop,n_var]=size(x);

% normalization
X=round((x-repmat(xmin,n_pop,1))./repmat(xmax-xmin,n_pop,1)*(2^(bit)-1));
% converte in binario
B=dec2bin(X,bit);
% concatenates
b=[];
for i=1:n_var
    b=[b B((i-1)*n_pop+1:i*n_pop,:)];
end