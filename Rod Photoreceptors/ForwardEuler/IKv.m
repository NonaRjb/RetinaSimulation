function [iKv, mK, hK] = IKv(V, mt, ht)
%IKv computes the delayed rectifying potassium current
%   V : membrane voltage at time t
%   mt : value of activation variable at time t
%   ht : value of inactivation variable at time t
%   iKv : delayed rectifying potassium current at time t
%   mK : value of activation variable at time t+1
%   hK : value of inactivation variable at time t+1

% constants
gKv = 2.0;
EK = -80;

amK = 5 * (100 - V)/ (exp((100 - V) / 42) - 1);
bmK = 9 * exp(-(V - 20) / 40);

ahK = 0.15 * exp(-V / 22);
bhK = 0.4125 / (exp((10 - V) / 7) + 1);

iKv = gKv * mt ^ 3 * ht * (V - EK);
mK = amK * (1 - mt) - bmK * mt + mt;
hK = ahK * (1 - ht) - bhK * ht + ht;

end

