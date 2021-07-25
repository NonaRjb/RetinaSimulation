function [iK_Ca, mKc] = IK_Ca(V, mt, Cs)
%IK_Ca computes Ca-dependent K current
%   V : membrane voltage at time t
%   mt : value of activation variable at time t
%   Cs : concentration of Ca in sub-membrane area
% constants 
gKc = 0.5;
Ek = -80;

aKc = 15 * (80 - V) / (exp((80 - V)/40) - 1);
bKc = 20 * exp(-V/35);

mKc1 = Cs / (Cs + 0.3);

iK_Ca = gKc * mt ^ 2 * mKc1 * (V - Ek);
mKc = aKc * (1 - mt) - bKc * mt + mt;
end

