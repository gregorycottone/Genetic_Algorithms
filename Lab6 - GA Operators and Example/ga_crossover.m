function xbf = ga_crossover(xb, F, bit)


% This function generates the offspring (prole) population with a single
% point crossover (incrocio)

% INPUT:
%   xb: matrix of individuals in binary code [xb1; xb2; ...], note that xb is the output of the function ga_coding
%   F: fitness of each individual (column vector of dimension length(xb) x 1)
%   bit: number of bits

% OUTPUT:
%   xbf: matrix of offsprings in binary code [xbf1; xbf2; ...]

[nind, ngeni] = size(xb); % nind = length(xb) while ngeni = bit*2
p = F/sum(F);             % for each individual it findes the corresponding probability of being selected for crossover (vector of dimension nind x 1)
ruota = cumsum(p);        % builds the roulette wheel (vector of dimension nind x 1)
% Note:
% cumsum(p) returns a vector containing the cumulative sum of the elements of p

gen = rand(size(F));      % chooses the points on the roulette wheel of the parents (vector of dimension nind x 1)

% ??????????????????????????????????
for i = 1:length(F)
    igen(i) = min(find(ruota >= gen(i)));  % indicates the parents
end

cross = round(ngeni*rand(1, nind/2));  % selection of the point in which the chromosomes are cut 

for i = 1:length(igen)/2
    xbf(2*i-1,:) = [xb(igen(2*i-1), 1:cross(i)), xb(igen(2*i), (cross(i)+1):end)];
    xbf(2*i,:) = [xb(igen(2*i), 1:cross(i)), xb(igen(2*i-1), (cross(i)+1):end)];
end
% ??????????????????????????????????

end

