function [t_end, t_peak, details, reason] = detect_t_end(ecg, r_peaks, beat_idx, fs)
% DETECT_T_END - T-wave end detection using the tangent method
%
% Usage:
%   [t_end, t_peak, details, reason] = detect_t_end(ecg, r_peaks, beat_idx, fs)
%
% Detects the T-wave end point by:
%   1. Defining a search window from 100 ms after the R-peak to 70% of RR
%   2. Finding the T-wave peak (largest deflection from baseline)
%   3. Computing the maximum slope on the trailing edge
%   4. Extrapolating the tangent line to the baseline
%
% Inputs:
%   ecg      - Filtered ECG signal vector
%   r_peaks  - Vector of all R-peak sample indices
%   beat_idx - Index into r_peaks for the current beat (1-based)
%   fs       - Sampling frequency (Hz)
%
% Outputs:
%   t_end   - T-end sample index (NaN on failure)
%   t_peak  - T-peak sample index (NaN on failure)
%   details - Struct with intermediate values (baseline, slope point, etc.)
%   reason  - String describing outcome ('OK' or failure reason)
%
% Reference:
%   Laguna P, Thakor NV, et al. New algorithm for QT interval analysis
%   in 24-hour Holter ECG. Med Biol Eng Comput. 1990;28(1):67-73.

    t_end  = NaN;
    t_peak = NaN;
    details = struct();
    reason = '';

    r_samp = r_peaks(beat_idx);
    n = length(ecg);

    % Search window: 100 ms after R to 70% of RR
    qrs_end = r_samp + round(0.10 * fs);
    if beat_idx < length(r_peaks)
        rr_samp = r_peaks(beat_idx + 1) - r_samp;
        search_end = r_samp + round(0.70 * rr_samp);
        search_end = min(search_end, r_peaks(beat_idx + 1) - round(0.05 * fs));
    else
        search_end = r_samp + round(0.6 * fs);
    end
    qrs_end    = max(qrs_end, 1);
    search_end = min(search_end, n);

    if search_end - qrs_end < round(0.08 * fs)
        reason = 'window too short';
        return;
    end

    window = ecg(qrs_end:search_end);
    w_idx  = qrs_end:search_end;

    % Baseline from last 30 ms of search window
    bl_n = max(1, round(0.03 * fs));
    baseline = median(window(end - bl_n + 1 : end));

    % T-peak: dominant peak relative to baseline
    [pos_val, pos_idx] = max(window - baseline);
    [neg_val, neg_idx] = min(window - baseline);

    if abs(pos_val) >= abs(neg_val)
        tp_local = pos_idx;
        polarity = 1;
    else
        tp_local = neg_idx;
        polarity = -1;
    end
    t_peak = w_idx(tp_local);

    % Check trailing edge has enough room
    if tp_local >= length(window) - round(0.04 * fs)
        reason = 'T-peak at window end';
        return;
    end

    % Derivative of trailing edge
    trailing = window(tp_local:end);
    dt = 1 / fs;
    deriv = diff(trailing) / dt;

    if polarity == 1
        [max_slope_val, msi] = min(deriv);
    else
        [max_slope_val, msi] = max(deriv);
    end

    if abs(max_slope_val) < 0.01
        reason = 'slope too flat';
        return;
    end

    sp_local = tp_local + msi - 1;
    sp_global = w_idx(sp_local);
    sp_val = window(sp_local);

    % Tangent extrapolation to baseline
    delta_samples = (baseline - sp_val) / (max_slope_val * dt);
    te_local = sp_local + round(delta_samples);

    if te_local <= tp_local
        reason = 'T-end before T-peak';
        return;
    end

    if te_local <= length(window)
        t_end = w_idx(te_local);
    else
        t_end = w_idx(end) + (te_local - length(window));
    end
    t_end = min(t_end, n);

    % Validate RT duration
    rt_s = (t_end - r_samp) / fs;
    if rt_s < 0.20 || rt_s > 0.60
        reason = sprintf('RT=%.0fms', rt_s * 1000);
        t_end = NaN;
        return;
    end

    % Must not pass next R
    if beat_idx < length(r_peaks) && t_end >= r_peaks(beat_idx + 1)
        reason = 'T-end past next R';
        t_end = NaN;
        return;
    end

    details.baseline  = baseline;
    details.slope_pt  = sp_global;
    details.slope_val = sp_val;
    details.max_slope = max_slope_val;
    reason = 'OK';
end
