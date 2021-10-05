function out = firing_rate(s, th, L, f)
%FIRING_RATE calculates the spike rate of the input signal
%   s : input signal
%   th: threshold
%   L : length of gaussian window
%   f : filtering window
if nargin < 3
    error('Insufficient input')
elseif nargin == 3
    f = gausswin(L);
end

[~, locs] = findpeaks(s);
raw_spks = s;
raw_spks(setdiff(1:length(s),locs)) = th-1;
spikes = double(raw_spks > th);
out = filter(f,1,spikes);
end

