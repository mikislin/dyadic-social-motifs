function y = segmentwise_lowpass(x, passFreq, fps, nPad)
%SEGMENTWISE_LOWPASS Fast lowpass using persistent butter filter

persistent b a lastFreq lastFps

% Recompute filter only if params changed
if isempty(b) || isempty(a) || passFreq ~= lastFreq || fps ~= lastFps
    Wn = passFreq / (fps/2); % normalized cutoff
    [b,a] = butter(2, Wn);   % 2nd order Butterworth
    lastFreq = passFreq;
    lastFps = fps;
end

% Apply zero-phase filtering
y = filtfilt(b, a, x);

% Restore edges (same as before)
if nPad > 0
    nPad = min([nPad, numel(x)]);
    y(1:nPad) = x(1:nPad);
    y(end-nPad+1:end) = x(end-nPad+1:end);
end
end