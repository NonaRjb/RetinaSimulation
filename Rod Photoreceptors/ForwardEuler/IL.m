function [iL] = IL(V)
%IL computes leakage current
%   V : membrane voltage at time t
% constants
gl = 0.5;
El = -55;

iL = gl * (V - El);
end

