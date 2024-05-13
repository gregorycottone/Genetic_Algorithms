function b=ag_coding_withDiscrete(x,xmin,xmax,bit)

% b=ag_coding_withDiscrete(x,xmin,xmax,bit)
%
% Ver 2.1 04/2022
%
% Matlab function that implements binary coding and concatenates variables
% [modified from Levi]
%
% input:
%   x: matrix of individuals (design variables arranged columnwise) [x1 x2 ...]
%   xmin: row vector with lower bounds of each design variable [1,ndv]
%   xmax: row vector with upper bounds of each design variable [1,ndv]
%   bit: row vector with number of bits of each design variable [1,ndv]
% output:
%   b: matrix of individuals in binary coding (each row is an individual) [b1; b2; ...] 

[n_pop,n_var]=size(x);

% normalization
X=[];
for jj=1:n_var
    Xi=round((x(:,jj)-repmat(xmin(jj),n_pop,1))./repmat(xmax(jj)-xmin(jj),n_pop,1)*(2^(bit(jj))-1));
    X=[X,Xi];
end
b=[];

% binary conversion and concatenation
for i=1:n_var
    Bi=dec2bin(X(:,i),bit(i));
    b=[b Bi];
end


