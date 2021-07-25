function [ih, M] = Ih(V, M0)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% constants
gh = 0.305;
Eh = -32;

ah = 3 / (exp((V + 88) / 14) + 1);
bh = 20 / (exp(-(V + 18) / 20) + 1);

K = [-4*ah, bh, 0, 0, 0;...
    4*ah, -(bh+3*ah), 2*bh, 0, 0;...
    0, 3*ah, -(2*bh+2*ah), 3*bh, 0;...
    0, 0, 2*ah, -(3*bh+ah), 4*bh;...
    0, 0, 0, ah, -4*bh];

M = M0 + K * M0;
mh = M0(3) + M0(4) + M0(5);
ih = gh * mh * (V - Eh);
end

