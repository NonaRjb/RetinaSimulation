function out = firing_rate(s, th, L, f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    error('Insufficient input')
elseif nargin == 3
    f = gausswin(L);
end

spikes = double(s > th);
out = filter(f,1,spikes);
end

