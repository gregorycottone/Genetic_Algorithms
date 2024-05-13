function x=ag_decoding_withDiscrete(b,xmin,xmax,bit)

% b=ag_coding_withDiscrete(x,xmin,xmax,bit)
%
% Ver 2.1 04/2022
%
% Matlab function that implements decimal decoding from binary string [modified from Levi]
%
% input:
%   b: matrix with binary population [b1; b2; ...] 
%   xmin: row vector with lower bounds of each design variable [1,ndv]
%   xmax: row vector with upper bounds of each design variable [1,ndv]
%   bit: row vector with number of bits of each design variable [1,ndv]% output:
% Output:
%   x: matrix with DVs (columnwise) [x1 x2 ...]

n_var=length(xmin);
n_pop=size(b,1);

% divide the string and decodes from binary
istart=1;
for i=1:n_var
    deltai=bit(i);
    X(:,i)=bin2dec(b(:,istart:istart+deltai-1));
    istart=istart+deltai;
end
% denormalization
x=[];
for jj=1:n_var
    xi=X(:,jj)./(2^bit(jj)-1).*repmat(xmax(jj)-xmin(jj),n_pop,1)+repmat(xmin(jj),n_pop,1);
    x=[x,xi];
end
end