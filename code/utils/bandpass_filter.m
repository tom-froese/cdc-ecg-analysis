function y = bandpass_filter(x, fs, f_low, f_high)
% BANDPASS_FILTER - Zero-phase Butterworth bandpass filter
%
% Usage:
%   y = bandpass_filter(x, fs, f_low, f_high)
%
% Applies a 4th-order zero-phase Butterworth bandpass filter.
% Implemented as cascaded high-pass and low-pass for numerical stability.
%
% Inputs:
%   x      - Input signal vector
%   fs     - Sampling frequency (Hz)
%   f_low  - Lower cutoff frequency (Hz)
%   f_high - Upper cutoff frequency (Hz)
%
% Output:
%   y - Filtered signal

    [b, a] = butter(4, f_low / (fs/2), 'high');
    y = filtfilt(b, a, x);
    [b, a] = butter(4, f_high / (fs/2), 'low');
    y = filtfilt(b, a, y);
end
