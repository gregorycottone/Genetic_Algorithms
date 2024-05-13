function [xb,p]=ag_mutation_withDiscrete(xb,par)

% [xb,p]=ag_mutation_withDiscrete(xb,par)
%
% Matlab function that implements random mutation of chromosomes 
% Ver 2.1 04/2022
%
% input:
%   xb: population individuals [xb1; xb2; ...]
%   par: parameter to set mutation probability-> pmut=par*pmin+(1-par)*pmax
%        pmin=min([1/populationSize; 1/len(chromosome)]),
%        pmin=max([1/populationSize; 1/len(chromosome)])
% output:
%   xb: matrice with mutated individuals [xb1; xb2; ...]
%   p: mutation probability

[nind,ngeni]=size(xb);
pmin=min([1/nind 1/ngeni]);
pmax=max([1/nind 1/ngeni]);
p=par*pmin+(1-par)*pmax; % probabilità che il gene muti
mut=rand(size(xb)); % casualità di mutazione di ogni gene
imut=mut<p; % indica quali geni devono mutare
xb(imut)=dec2bin(abs(1-bin2dec(xb(imut)))); % muta i geni