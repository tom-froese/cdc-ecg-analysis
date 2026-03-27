function [r_peaks, t_ends] = extract_rt_pairs(ann_samples, ann_symbols)
% EXTRACT_RT_PAIRS - Extract matched R-peak and T-end sample pairs
%
% Usage:
%   [r_peaks, t_ends] = extract_rt_pairs(ann_samples, ann_symbols)
%
% Given arrays of annotation samples and symbols (from read_wfdb_annotations),
% extracts R-peak locations (symbol 'N') and matches each to its
% corresponding T-wave offset (symbol ')' following a 't' symbol).
%
% Only returns pairs where the T-end falls between the current R-peak
% and the next R-peak (i.e., physiologically plausible matches).
%
% Inputs:
%   ann_samples - Column vector of annotation sample numbers
%   ann_symbols - Cell array of annotation symbols
%
% Outputs:
%   r_peaks - Column vector of matched R-peak sample numbers
%   t_ends  - Column vector of matched T-end sample numbers

    % Extract all R-peaks
    r_peaks_all = [];
    for i = 1:length(ann_symbols)
        if strcmp(ann_symbols{i}, 'N')
            r_peaks_all(end+1) = ann_samples(i);
        end
    end

    % Extract all T-end markers: ')' following 't'
    t_ends_all = [];
    for i = 1:length(ann_symbols)-1
        if strcmp(ann_symbols{i}, 't') && strcmp(ann_symbols{i+1}, ')')
            t_ends_all(end+1) = ann_samples(i+1);
        end
    end

    r_peaks_all = r_peaks_all(:);
    t_ends_all  = t_ends_all(:);

    % Match each R-peak to the nearest following T-end
    % (must fall before the next R-peak)
    matched_r = [];
    matched_t = [];

    for i = 1:length(r_peaks_all)
        r = r_peaks_all(i);
        t_idx = find(t_ends_all > r, 1, 'first');
        if ~isempty(t_idx)
            if i < length(r_peaks_all)
                if t_ends_all(t_idx) < r_peaks_all(i+1)
                    matched_r(end+1) = r;
                    matched_t(end+1) = t_ends_all(t_idx);
                end
            else
                matched_r(end+1) = r;
                matched_t(end+1) = t_ends_all(t_idx);
            end
        end
    end

    r_peaks = matched_r(:);
    t_ends  = matched_t(:);
end
