function xbf=ga_crossover(xb,F,bit)

% xbf=ga_crossover(x,F,bit)
%
% this function generates the offspring population with a single point crossover
%
% input:
%   xb: matrix of individuals [xb1; xb2; ...]
%   F: fitness of each individual (column vector)
%   bit: number of bits
% output:
%   xbf: matrix of offsprings [xbf1; xbf2; ...] 

[nind,ngeni]=size(xb);
p=F/sum(F); % probability of being selected for crossover
ruota=cumsum(p); % builds the roulette wheel
gen=rand(size(F)); % chooses the points on the roulette wheel of the parents
for i=1:length(F)
    igen(i)=min(find(ruota>=gen(i))); % indicates the parents
end
cross=round(ngeni*rand(1,nind/2));% selection of the point in which the chromosomes are cut 
for i=1:length(igen)/2
    xbf(2*i-1,:)=[xb(igen(2*i-1),1:cross(i)) xb(igen(2*i),(cross(i)+1):end)];
    xbf(2*i,:)=[xb(igen(2*i),1:cross(i)) xb(igen(2*i-1),(cross(i)+1):end) ];
end