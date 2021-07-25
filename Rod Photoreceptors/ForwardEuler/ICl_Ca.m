function [iCl_Ca] = ICl_Ca(V, Cs)
%ICl_Ca computes Ca-dependent Cl current
%   V : membrane voltage at time t
%   Cs : concentration of Ca just below the membrane
% constants
gCl = 6.5;
ECl = -45;
Clh = 0.37;

mCl = 1 / (1 + exp((Clh - Cs)/0.09));
iCl_Ca = gCl * mCl * (V - ECl);
end

