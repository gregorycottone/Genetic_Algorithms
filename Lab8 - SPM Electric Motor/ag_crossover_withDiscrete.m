function xbf=ag_crossover_withDiscrete(xb,F)

% xbf=ag_crossover_withDiscrete(x,F)
%
% Matlab function that implements single point crossover with Roulette Wheel 
% Ver 2.1 04/2022
%
% input:
%   xb: matrix of individuals in binary coding [xb1; xb2; ...]
%   F: fitness of each individual (vettore colonna)
% output:
%   xbf: matrix of offsprings [xbf1; xbf2; ...] 

[nind,ngeni]=size(xb);
p=F/sum(F); % probabilità di essere scelto per generare figli
ruota=cumsum(p); % costruisce la "ruota" dei genitori
gen=rand(size(F)); % sceglie punti sulla "ruota" dei genitori
for i=1:length(F)
    igen(i)=min(find(ruota>=gen(i))); % indica quali sono i genitori
end

cross=round(ngeni*rand(1,nind/2));% sceglie dove tagliare i cromosomi
for i=1:length(igen)/2
    xbf(2*i-1,:)=[xb(igen(2*i-1),1:cross(i)) xb(igen(2*i),(cross(i)+1):end)];
    xbf(2*i,:)=[xb(igen(2*i),1:cross(i)) xb(igen(2*i-1),(cross(i)+1):end) ];
end