function xb=ga_mutation(xb,par,bit)

% xbm=ga_mutation(xb,par,bit)
%
% this function performs the mutation of the genes
%
% input:
%   xb: matrix of the individuals [xb1; xb2; ...]
%   par: mutation parameter pmut=par*pmin+(1-par)*pmax
%        pmin=min([1/nind; 1/len(chromosome)]) and the same for pmax
%   bit: number of bits
% output:
%   xb: matrix of mutated individuals [xb1; xb2; ...] 

[nind,ngeni]=size(xb);
pmin=min([1/nind 1/ngeni]);
pmax=max([1/nind 1/ngeni]);
p=par*pmin+(1-par)*pmax; % probability of mutation
mut=rand(size(xb)); % casuality of mutation of each gene
imut=mut<p; % indicates which are the genes to be mutated
xb(imut)=dec2bin(abs(1-bin2dec(xb(imut)))); % mutation of the genes