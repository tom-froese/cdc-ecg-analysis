function r_peaks = detect_r_peaks(ecg, fs)
% DETECT_R_PEAKS - Pan-Tompkins-inspired R-peak detection
%
% Usage:
%   r_peaks = detect_r_peaks(ecg, fs)
%
% Pan-Tompkins-inspired algorithm with adaptive two-level thresholding:
%   1. Bandpass 5-15 Hz (QRS energy enhancement)
%   2. Differentiate
%   3. Square
%   4. Moving window integration (150 ms)
%   5. Adaptive thresholding (signal/noise level tracking)
%   6. Refine to actual R-peak maxima in original ECG
%
% Inputs:
%   ecg - ECG signal vector (in physical units, e.g. mV)
%   fs  - Sampling frequency (Hz)
%
% Output:
%   r_peaks - Column vector of R-peak sample indices
%
% Reference:
%   Pan J, Tompkins WJ. A real-time QRS detection algorithm.
%   IEEE Trans Biomed Eng. 1985;32(3):230-236.

    %% 1. Bandpass 5-15 Hz
    [b, a] = butter(2, [5 / (fs/2), 15 / (fs/2)], 'bandpass');
    filtered = filtfilt(b, a, ecg);

    %% 2. Differentiate
    diff_ecg = diff(filtered);

    %% 3. Square
    squared = diff_ecg .^ 2;

    %% 4. Moving window integration (150 ms)
    win_len = round(0.150 * fs);
    kernel = ones(win_len, 1) / win_len;
    integrated = conv(squared, kernel, 'same');

    %% 5. Find peaks with initial threshold
    % Initial threshold from first 10 seconds
    init_samp = min(length(integrated), 10 * fs);
    threshold = 0.3 * max(integrated(1:init_samp));

    % Minimum distance between peaks: 200 ms (300 bpm ceiling)
    min_dist = round(0.200 * fs);

    % Find all candidate peaks above a low threshold
    [pks, locs] = findpeaks(integrated, 'MinPeakHeight', threshold * 0.2, ...
                            'MinPeakDistance', min_dist);

    if isempty(locs)
        r_peaks = [];
        return;
    end

    %% 6. Adaptive thresholding (two-level, as in Pan-Tompkins)
    spki = threshold;       % Signal peak running estimate
    npki = threshold * 0.1; % Noise peak running estimate
    threshold1 = npki + 0.25 * (spki - npki);

    keep = false(size(locs));
    for i = 1:length(locs)
        if pks(i) > threshold1
            keep(i) = true;
            spki = 0.125 * pks(i) + 0.875 * spki;
        else
            npki = 0.125 * pks(i) + 0.875 * npki;
        end
        threshold1 = npki + 0.25 * (spki - npki);
    end

    locs = locs(keep);

    %% 7. Refine to actual R-peak maxima in ECG (±75 ms)
    search_hw = round(0.075 * fs);
    r_peaks = zeros(size(locs));
    for i = 1:length(locs)
        s0 = max(1, locs(i) - search_hw);
        s1 = min(length(ecg), locs(i) + search_hw);
        [~, idx] = max(ecg(s0:s1));
        r_peaks(i) = s0 + idx - 1;
    end

    % Remove duplicates
    r_peaks = unique(r_peaks);

    % Final pass: remove any peaks closer than 200 ms
    if length(r_peaks) > 1
        keep = true(size(r_peaks));
        for i = 2:length(r_peaks)
            if (r_peaks(i) - r_peaks(i-1)) < min_dist
                keep(i) = false;
            end
        end
        r_peaks = r_peaks(keep);
    end

    r_peaks = r_peaks(:);
end
