function y = lowpass_filter(x, fs, fc)
% LOWPASS_FILTER - Zero-phase Butterworth lowpass filter
%
% Usage:
%   y = lowpass_filter(x, fs, fc)
%
% Applies a 4th-order zero-phase Butterworth lowpass filter.
%
% Inputs:
%   x  - Input signal vector
%   fs - Sampling frequency (Hz)
%   fc - Cutoff frequency (Hz)
%
% Output:
%   y - Filtered signal

    [b, a] = butter(4, fc / (fs/2), 'low');
    y = filtfilt(b, a, x);
end
