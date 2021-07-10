function [iCa, mCa] = ICa(V, mt, Co, Cs)
%ICa computes the calcium current
%   V : membrane voltage at time t
%   mt : value of activation variable at time t
%   Co : concentration of Ca outside the cell
%   Cs : concentration of Ca just below the membrane

% constants
gCa = 1.2;

aCa = 300 * (80 - V) / (exp((80 - V) / 25) - 1);
bCa = 1000 / (1 + exp((V + 38)/ 7));

ECa = 12.9 * log(Co / Cs); 

mCah = exp((40-V)/18)/(1+exp((40-V)/18));

iCa = gCa * mt ^ 4 * mCah * (V - ECa);
mCa = aCa * (1 - mt) - bCa * mt + mt;
end

