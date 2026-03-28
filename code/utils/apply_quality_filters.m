function [valid_mask, filter_report] = apply_quality_filters(beats, varargin)
% APPLY_QUALITY_FILTERS - Uniform quality filters for CDC analysis
%
% Usage:
%   valid_mask = apply_quality_filters(beats)
%   [valid_mask, report] = apply_quality_filters(beats, 'Name', Value, ...)
%
% Applies physiologically motivated quality filters uniformly across all
% datasets. These filters are applied at analysis time (not during export)
% to ensure transparency and reproducibility.
%
% Two stages of filtering:
%
%   Stage 1 — Beat-level (returned in valid_mask):
%     RR interval:  0.3 s < RR < 3.0 s   (20-200 bpm range)
%     RT interval:  0.15 s < RT < 0.6 s   (physiological QT bounds)
%     CDC ratio:    0.2 < RT/RR < 0.6     (physiological duty cycle range)
%     RT < RR:      T-end must precede next R-peak
%
%   Stage 2 — Subject-level (returned in filter_report):
%     MinBeatsPerSubject: Minimum number of valid beats required per
%     subject for a meaningful subject-level median. Subjects with
%     fewer valid beats are flagged in filter_report.excluded_subjects.
%     The caller is responsible for excluding these subjects after
%     aggregation; the beat-level valid_mask is unaffected.
%
% Inputs:
%   beats - Table with columns: r_sample, t_end_sample, next_r_sample, fs,
%           unique_subject_id
%
% Optional Name-Value pairs:
%   'RRMin'              - Minimum RR interval in seconds (default: 0.3)
%   'RRMax'              - Maximum RR interval in seconds (default: 3.0)
%   'RTMin'              - Minimum RT interval in seconds (default: 0.15)
%   'RTMax'              - Maximum RT interval in seconds (default: 0.6)
%   'CDCMin'             - Minimum CDC ratio (default: 0.2)
%   'CDCMax'             - Maximum CDC ratio (default: 0.6)
%   'MinBeatsPerSubject' - Minimum valid beats per subject (default: 2)
%   'Verbose'            - Print filter report (default: true)
%
% Outputs:
%   valid_mask    - Logical vector (true = beat passes all beat-level filters)
%   filter_report - Struct with rejection counts per criterion and:
%       .min_beats_per_subject  - The threshold used
%       .excluded_subjects      - Cell array of unique_subject_id values
%                                 with fewer valid beats than the threshold
%       .n_subjects_excluded    - Count of excluded subjects

    %% Parse optional arguments
    p = inputParser;
    addParameter(p, 'RRMin', 0.3);
    addParameter(p, 'RRMax', 3.0);
    addParameter(p, 'RTMin', 0.15);
    addParameter(p, 'RTMax', 0.6);
    addParameter(p, 'CDCMin', 0.2);
    addParameter(p, 'CDCMax', 0.6);
    addParameter(p, 'MinBeatsPerSubject', 2);
    addParameter(p, 'Verbose', true);
    parse(p, varargin{:});
    opts = p.Results;

    %% Compute intervals
    n_beats = height(beats);
    fs = beats.fs;

    rr_s = (beats.next_r_sample - beats.r_sample) ./ fs;
    rt_s = (beats.t_end_sample - beats.r_sample) ./ fs;
    cdc  = rt_s ./ rr_s;

    %% Apply filters sequentially (for reporting)
    f_rr_pos   = rr_s > 0;
    f_rt_pos   = rt_s > 0;
    f_rt_lt_rr = rt_s < rr_s;
    f_rr_min   = rr_s > opts.RRMin;
    f_rr_max   = rr_s < opts.RRMax;
    f_rt_min   = rt_s > opts.RTMin;
    f_rt_max   = rt_s < opts.RTMax;
    f_cdc_min  = cdc > opts.CDCMin;
    f_cdc_max  = cdc < opts.CDCMax;

    valid_mask = f_rr_pos & f_rt_pos & f_rt_lt_rr & ...
                 f_rr_min & f_rr_max & f_rt_min & f_rt_max & ...
                 f_cdc_min & f_cdc_max;

    %% Build report
    filter_report.n_total     = n_beats;
    filter_report.n_valid     = sum(valid_mask);
    filter_report.n_rejected  = n_beats - sum(valid_mask);
    filter_report.pct_valid   = 100 * sum(valid_mask) / n_beats;

    filter_report.rejected_rr_nonpositive = sum(~f_rr_pos);
    filter_report.rejected_rt_nonpositive = sum(~f_rt_pos);
    filter_report.rejected_rt_exceeds_rr  = sum(f_rr_pos & f_rt_pos & ~f_rt_lt_rr);
    filter_report.rejected_rr_too_short   = sum(f_rr_pos & ~f_rr_min);
    filter_report.rejected_rr_too_long    = sum(f_rr_pos & f_rr_min & ~f_rr_max);
    filter_report.rejected_rt_too_short   = sum(f_rt_pos & ~f_rt_min);
    filter_report.rejected_rt_too_long    = sum(f_rt_pos & f_rt_min & ~f_rt_max);
    filter_report.rejected_cdc_too_low    = sum(f_rr_pos & f_rt_pos & f_rt_lt_rr & ~f_cdc_min);
    filter_report.rejected_cdc_too_high   = sum(f_rr_pos & f_rt_pos & f_rt_lt_rr & f_cdc_min & ~f_cdc_max);

    %% Subject-level filter: minimum valid beats per subject
    min_beats = opts.MinBeatsPerSubject;
    filter_report.min_beats_per_subject = min_beats;

    if min_beats > 0 && ismember('unique_subject_id', beats.Properties.VariableNames)
        [unique_subjects, ~, subj_idx] = unique(beats.unique_subject_id);
        n_valid_per_subject = accumarray(subj_idx, valid_mask, [], @sum);
        below_threshold = n_valid_per_subject < min_beats;
        filter_report.excluded_subjects    = unique_subjects(below_threshold);
        filter_report.n_subjects_excluded  = sum(below_threshold);
        filter_report.n_subjects_total     = length(unique_subjects);
    else
        filter_report.excluded_subjects    = {};
        filter_report.n_subjects_excluded  = 0;
        filter_report.n_subjects_total     = 0;
    end

    %% Print report
    if opts.Verbose
        fprintf('  Quality filter: %d/%d beats pass (%.1f%%)\n', ...
                filter_report.n_valid, filter_report.n_total, filter_report.pct_valid);
        if filter_report.n_rejected > 0
            fprintf('    Rejected: RR<=0: %d, RT<=0: %d, RT>=RR: %d\n', ...
                    filter_report.rejected_rr_nonpositive, ...
                    filter_report.rejected_rt_nonpositive, ...
                    filter_report.rejected_rt_exceeds_rr);
            fprintf('    Rejected: RR<%.1fs: %d, RR>%.1fs: %d\n', ...
                    opts.RRMin, filter_report.rejected_rr_too_short, ...
                    opts.RRMax, filter_report.rejected_rr_too_long);
            fprintf('    Rejected: RT<%.2fs: %d, RT>%.1fs: %d\n', ...
                    opts.RTMin, filter_report.rejected_rt_too_short, ...
                    opts.RTMax, filter_report.rejected_rt_too_long);
            fprintf('    Rejected: CDC<%.1f: %d, CDC>%.1f: %d\n', ...
                    opts.CDCMin, filter_report.rejected_cdc_too_low, ...
                    opts.CDCMax, filter_report.rejected_cdc_too_high);
        end
        if filter_report.n_subjects_excluded > 0
            fprintf('    Subjects with <%d valid beats: %d/%d excluded\n', ...
                    min_beats, filter_report.n_subjects_excluded, ...
                    filter_report.n_subjects_total);
        end
    end
end