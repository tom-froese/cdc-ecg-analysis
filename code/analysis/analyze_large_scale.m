function results = analyze_large_scale()
% ANALYZE_LARGE_SCALE - CDC analysis using large algorithmically annotated databases
%
% This analysis uses datasets where fiducial points were obtained through
% validated algorithms or hybrid approaches:
%
%   1. Fantasia        - HC: Database R-peaks, tangent T-end (healthy volunteers)
%   2. Autonomic Aging - HC: Fully automatic (healthy volunteers)
%   3. PTB             - HC vs Path: Manual T-end, algorithmic R-peak
%   4. PTB-XL          - CN vs Path: ECGDeli validated algorithm
%
% Group classification follows the hierarchical model (analyze_hierarchical_model.m):
%   Healthy Control   (HC)  = Fantasia, Autonomic Aging, PTB 'healthy'
%                             (all verified volunteers per source documentation)
%   Clinically Normal (CN)  = PTB-XL 'healthy' (hospital patients, normal ECG)
%   Pathological      (Path) = PTB, PTB-XL 'pathological'
%
% Purpose: Confirm the gold-standard findings at scale. These datasets
% provide statistical power that the small manually annotated databases
% cannot. Uniform quality filters are applied to all datasets, including
% a minimum of 2 valid beats per subject for a meaningful median.
%
% Subject aggregation uses unique_subject_id (format: Database_RecordID)
% for consistency with the hierarchical model pipeline.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    n_bootstrap = 5000;
    inv_e = 1 / exp(1);

    fprintf('================================================================\n');
    fprintf('LARGE-SCALE CDC CONFIRMATION (Algorithmic Annotation)\n');
    fprintf('================================================================\n');
    fprintf('Reference value: 1/e = %.4f\n', inv_e);
    fprintf('Uniform quality filters applied to all datasets.\n');
    fprintf('================================================================\n\n');

    %% Load datasets
    ptb_beats = readtable(paths.csv_ptb);
    ptbxl_beats = readtable(paths.csv_ptbxl);
    fantasia_beats = readtable(paths.csv_fantasia);
    aa_beats = readtable(paths.csv_autonomic_aging);

    fprintf('PTB:             %d beats loaded\n', height(ptb_beats));
    fprintf('PTB-XL:          %d beats loaded\n', height(ptbxl_beats));
    fprintf('Fantasia:        %d beats loaded\n', height(fantasia_beats));
    fprintf('Autonomic Aging: %d beats loaded\n\n', height(aa_beats));

    %% Apply uniform quality filters
    fprintf('Applying quality filters...\n\n');

    fprintf('  PTB:\n');
    [ptb_mask, ptb_filt] = apply_quality_filters(ptb_beats);
    ptb_beats = ptb_beats(ptb_mask, :);
    ptb_excluded = ptb_filt.excluded_subjects;

    fprintf('  PTB-XL:\n');
    [ptbxl_mask, ptbxl_filt] = apply_quality_filters(ptbxl_beats);
    ptbxl_beats = ptbxl_beats(ptbxl_mask, :);
    ptbxl_excluded = ptbxl_filt.excluded_subjects;

    fprintf('  Fantasia:\n');
    [fan_mask, fan_filt] = apply_quality_filters(fantasia_beats);
    fantasia_beats = fantasia_beats(fan_mask, :);
    fan_excluded = fan_filt.excluded_subjects;

    fprintf('  Autonomic Aging:\n');
    [aa_mask, aa_filt] = apply_quality_filters(aa_beats);
    aa_beats = aa_beats(aa_mask, :);
    aa_excluded = aa_filt.excluded_subjects;

    fprintf('\nAfter filtering:\n');
    fprintf('  PTB:             %d beats\n', height(ptb_beats));
    fprintf('  PTB-XL:          %d beats\n', height(ptbxl_beats));
    fprintf('  Fantasia:        %d beats\n', height(fantasia_beats));
    fprintf('  Autonomic Aging: %d beats\n\n', height(aa_beats));

    %% Analysis 1: Fantasia - Healthy volunteers
    fprintf('================================================================\n');
    fprintf('FANTASIA: HEALTHY VOLUNTEERS\n');
    fprintf('(Database R-peaks, tangent T-end, 2-hour recordings)\n');
    fprintf('================================================================\n');

    fantasia_results = analyze_healthy_cohort(fantasia_beats, 'Fantasia', n_bootstrap, inv_e, fan_excluded);

    %% Analysis 2: Autonomic Aging - Healthy volunteers
    fprintf('\n================================================================\n');
    fprintf('AUTONOMIC AGING: HEALTHY VOLUNTEERS\n');
    fprintf('(Fully automatic, ages 18-92)\n');
    fprintf('================================================================\n');

    aa_results = analyze_healthy_cohort(aa_beats, 'Autonomic Aging', n_bootstrap, inv_e, aa_excluded);

    %% Analysis 3: PTB - Healthy Control vs Pathological
    fprintf('\n================================================================\n');
    fprintf('PTB: HEALTHY CONTROL vs PATHOLOGICAL\n');
    fprintf('(Manual T-end, algorithmic R-peak)\n');
    fprintf('================================================================\n');

    ptb_results = analyze_hc_vs_pathological(ptb_beats, 'PTB', n_bootstrap, inv_e, ptb_excluded);

    %% Analysis 4: PTB-XL - Clinically normal vs Pathological
    fprintf('\n================================================================\n');
    fprintf('PTB-XL: CLINICALLY NORMAL vs PATHOLOGICAL\n');
    fprintf('(ECGDeli automatic annotation, N ~ 18,000 patients)\n');
    fprintf('================================================================\n');

    ptbxl_results = analyze_cn_vs_pathological(ptbxl_beats, 'PTB-XL', n_bootstrap, inv_e, ptbxl_excluded);

    %% Integrated summary
    fprintf('\n================================================================\n');
    fprintf('LARGE-SCALE SUMMARY\n');
    fprintf('================================================================\n\n');

    fprintf('                              Mode      95%% CI          n      %% > 1/e\n');
    fprintf('------------------------------------------------------------------------\n');
    fprintf('Fantasia HC                   %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            fantasia_results.mode_all, fantasia_results.ci_all(1), ...
            fantasia_results.ci_all(2), fantasia_results.n_all, ...
            fantasia_results.pct_above);
    fprintf('Autonomic Aging HC            %.4f    [%.3f, %.3f]  %5d    %.1f%%\n', ...
            aa_results.mode_all, aa_results.ci_all(1), ...
            aa_results.ci_all(2), aa_results.n_all, aa_results.pct_above);
    fprintf('------------------------------------------------------------------------\n');
    if ~isnan(ptb_results.mode_hc)
        fprintf('PTB Healthy control           %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
                ptb_results.mode_hc, ptb_results.ci_hc(1), ...
                ptb_results.ci_hc(2), ptb_results.n_hc, ...
                ptb_results.pct_hc_above);
    end
    fprintf('PTB Pathological              %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            ptb_results.mode_path, ptb_results.ci_path(1), ...
            ptb_results.ci_path(2), ptb_results.n_path, ...
            ptb_results.pct_path_above);
    fprintf('------------------------------------------------------------------------\n');
    fprintf('PTB-XL Clinically normal      %.4f    [%.3f, %.3f]  %5d    %.1f%%\n', ...
            ptbxl_results.mode_cn, ptbxl_results.ci_cn(1), ...
            ptbxl_results.ci_cn(2), ptbxl_results.n_cn, ...
            ptbxl_results.pct_cn_above);
    fprintf('PTB-XL Pathological           %.4f    [%.3f, %.3f]  %5d    %.1f%%\n', ...
            ptbxl_results.mode_path, ptbxl_results.ci_path(1), ...
            ptbxl_results.ci_path(2), ptbxl_results.n_path, ...
            ptbxl_results.pct_path_above);
    fprintf('------------------------------------------------------------------------\n');
    fprintf('Reference: 1/e = %.4f\n', inv_e);
    fprintf('================================================================\n');

    %% Package results
    results.fantasia = fantasia_results;
    results.autonomic_aging = aa_results;
    results.ptb = ptb_results;
    results.ptbxl = ptbxl_results;
    results.inv_e = inv_e;

    save(fullfile(paths.results, 'large_scale_results.mat'), 'results');
    fprintf('\nResults saved to: %s\n', fullfile(paths.results, 'large_scale_results.mat'));
end

%% ========================================================================
%  DATASET-SPECIFIC ANALYSIS FUNCTIONS
%  ========================================================================

function results = analyze_healthy_cohort(beats, dataset_name, n_bootstrap, inv_e, excluded_subjects)
% ANALYZE_HEALTHY_COHORT - Single-group analysis for healthy populations
%
% Used for Fantasia and Autonomic Aging: all subjects are healthy
% volunteers, so we compute a single distribution and test proximity
% to 1/e. Age-related effects are handled by the hierarchical model.

    fprintf('Processing %s...\n', dataset_name);

    [unique_subjects, ~, subject_idx] = unique(beats.unique_subject_id);
    n_subjects = length(unique_subjects);

    % Identify subjects excluded by minimum-beats filter
    if nargin >= 5 && ~isempty(excluded_subjects)
        is_excluded = ismember(unique_subjects, excluded_subjects);
    else
        is_excluded = false(n_subjects, 1);
    end

    median_ratio = zeros(n_subjects, 1);
    n_beats_per = zeros(n_subjects, 1);

    for i = 1:n_subjects
        mask = subject_idx == i;
        rec = beats(mask, :);
        rt_ms = (rec.t_end_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        rr_ms = (rec.next_r_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        valid = rr_ms > 0;
        if sum(valid) >= 1
            median_ratio(i) = median(rt_ms(valid) ./ rr_ms(valid));
            n_beats_per(i) = sum(valid);
        else
            median_ratio(i) = NaN;
        end
    end

    valid = ~isnan(median_ratio) & n_beats_per > 0 & ~is_excluded;
    median_ratio = median_ratio(valid);

    n_all = length(median_ratio);
    fprintf('  Valid subjects: %d\n', n_all);

    [mode_all, ci_all] = bootstrap_mode(median_ratio, n_bootstrap);
    pct_above = 100 * sum(median_ratio > inv_e) / n_all;

    fprintf('  Mode: %.4f [%.4f, %.4f]\n', mode_all, ci_all(1), ci_all(2));
    fprintf('  Mean: %.4f, Median: %.4f\n', mean(median_ratio), median(median_ratio));
    fprintf('  1/e in CI: %s\n', yesno(inv_e >= ci_all(1) && inv_e <= ci_all(2)));

    results = struct('n_all', n_all, 'mode_all', mode_all, 'ci_all', ci_all, ...
        'pct_above', pct_above, 'all_ratios', median_ratio);
end

function results = analyze_cn_vs_pathological(beats, dataset_name, n_bootstrap, inv_e, excluded_subjects)
% ANALYZE_CN_VS_PATHOLOGICAL - Two-group analysis: Clinically Normal vs Pathological
%
% Used for PTB-XL. The 'healthy' group label denotes hospital patients
% with normal ECG findings — Clinically Normal, not Healthy Control.

    fprintf('Processing %s...\n', dataset_name);

    [unique_subjects, ~, subject_idx] = unique(beats.unique_subject_id);
    n_subjects = length(unique_subjects);

    % Identify subjects excluded by minimum-beats filter
    if nargin >= 5 && ~isempty(excluded_subjects)
        is_excluded = ismember(unique_subjects, excluded_subjects);
    else
        is_excluded = false(n_subjects, 1);
    end

    median_ratio = zeros(n_subjects, 1);
    n_beats_per = zeros(n_subjects, 1);
    group = cell(n_subjects, 1);

    for i = 1:n_subjects
        mask = subject_idx == i;
        rec = beats(mask, :);
        rt_ms = (rec.t_end_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        rr_ms = (rec.next_r_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        valid = (rr_ms > 0) & (rt_ms > 0) & (rt_ms < rr_ms);
        if sum(valid) >= 1
            median_ratio(i) = median(rt_ms(valid) ./ rr_ms(valid));
            n_beats_per(i) = sum(valid);
        else
            median_ratio(i) = NaN;
        end
        group{i} = get_string_field(rec, 'group');
    end

    valid = ~isnan(median_ratio) & n_beats_per > 0 & ~is_excluded;
    median_ratio = median_ratio(valid);
    group = group(valid);

    cn_ratios   = median_ratio(strcmp(group, 'healthy'));
    path_ratios = median_ratio(~strcmp(group, 'healthy'));
    n_cn   = length(cn_ratios);
    n_path = length(path_ratios);

    fprintf('  Clinically normal: n = %d\n', n_cn);
    fprintf('  Pathological:      n = %d\n', n_path);

    if n_cn >= 3
        [mode_cn, ci_cn] = bootstrap_mode(cn_ratios, n_bootstrap);
        pct_cn_above = 100 * sum(cn_ratios > inv_e) / n_cn;
    else
        mode_cn = NaN; ci_cn = [NaN NaN]; pct_cn_above = NaN;
    end

    [mode_path, ci_path] = bootstrap_mode(path_ratios, n_bootstrap);
    pct_path_above = 100 * sum(path_ratios > inv_e) / n_path;

    fprintf('  CN mode:   %.4f [%.4f, %.4f]\n', mode_cn, ci_cn(1), ci_cn(2));
    fprintf('  Path mode: %.4f [%.4f, %.4f]\n', mode_path, ci_path(1), ci_path(2));

    if n_cn >= 3 && n_path >= 3
        [p_val, ~, stats] = ranksum(cn_ratios, path_ratios);
        fprintf('  Wilcoxon: p = %.2e%s\n', p_val, stars(p_val));
        % Rank-biserial effect size
        if isfield(stats, 'ranksum')
            U = stats.ranksum - n_cn*(n_cn+1)/2;
            effect_size = 1 - 2*U/(n_cn*n_path);
        else
            effect_size = NaN;
        end
        fprintf('  Effect size (rank-biserial r): %.3f\n', effect_size);
    else
        p_val = NaN; effect_size = NaN;
    end

    results = struct('n_cn', n_cn, 'n_path', n_path, ...
        'mode_cn', mode_cn, 'ci_cn', ci_cn, ...
        'mode_path', mode_path, 'ci_path', ci_path, ...
        'pct_cn_above', pct_cn_above, 'pct_path_above', pct_path_above, ...
        'p_value', p_val, 'effect_size', effect_size, ...
        'cn_ratios', cn_ratios, 'path_ratios', path_ratios);
end

function results = analyze_hc_vs_pathological(beats, dataset_name, n_bootstrap, inv_e, excluded_subjects)
% ANALYZE_HC_VS_PATHOLOGICAL - Two-group analysis: Healthy Control vs Pathological
%
% Used for PTB. The 'healthy' group label denotes verified healthy
% volunteers — Healthy Control — per the PhysioNet source documentation.

    fprintf('Processing %s...\n', dataset_name);

    [unique_subjects, ~, subject_idx] = unique(beats.unique_subject_id);
    n_subjects = length(unique_subjects);

    % Identify subjects excluded by minimum-beats filter
    if nargin >= 5 && ~isempty(excluded_subjects)
        is_excluded = ismember(unique_subjects, excluded_subjects);
    else
        is_excluded = false(n_subjects, 1);
    end

    median_ratio = zeros(n_subjects, 1);
    n_beats_per = zeros(n_subjects, 1);
    group = cell(n_subjects, 1);

    for i = 1:n_subjects
        mask = subject_idx == i;
        rec = beats(mask, :);
        rt_ms = (rec.t_end_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        rr_ms = (rec.next_r_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        valid = (rr_ms > 0) & (rt_ms > 0) & (rt_ms < rr_ms);
        if sum(valid) >= 1
            median_ratio(i) = median(rt_ms(valid) ./ rr_ms(valid));
            n_beats_per(i) = sum(valid);
        else
            median_ratio(i) = NaN;
        end
        group{i} = get_string_field(rec, 'group');
    end

    valid = ~isnan(median_ratio) & n_beats_per > 0 & ~is_excluded;
    median_ratio = median_ratio(valid);
    group = group(valid);

    hc_ratios   = median_ratio(strcmp(group, 'healthy'));
    path_ratios = median_ratio(~strcmp(group, 'healthy'));
    n_hc   = length(hc_ratios);
    n_path = length(path_ratios);

    fprintf('  Healthy control: n = %d\n', n_hc);
    fprintf('  Pathological:    n = %d\n', n_path);

    if n_hc >= 3
        [mode_hc, ci_hc] = bootstrap_mode(hc_ratios, n_bootstrap);
        pct_hc_above = 100 * sum(hc_ratios > inv_e) / n_hc;
    else
        mode_hc = NaN; ci_hc = [NaN NaN]; pct_hc_above = NaN;
    end

    [mode_path, ci_path] = bootstrap_mode(path_ratios, n_bootstrap);
    pct_path_above = 100 * sum(path_ratios > inv_e) / n_path;

    fprintf('  HC mode:   %.4f [%.4f, %.4f]\n', mode_hc, ci_hc(1), ci_hc(2));
    fprintf('  Path mode: %.4f [%.4f, %.4f]\n', mode_path, ci_path(1), ci_path(2));

    if n_hc >= 3 && n_path >= 3
        [p_val, ~, stats] = ranksum(hc_ratios, path_ratios);
        fprintf('  Wilcoxon: p = %.2e%s\n', p_val, stars(p_val));
        if isfield(stats, 'ranksum')
            U = stats.ranksum - n_hc*(n_hc+1)/2;
            effect_size = 1 - 2*U/(n_hc*n_path);
        else
            effect_size = NaN;
        end
        fprintf('  Effect size (rank-biserial r): %.3f\n', effect_size);
    else
        p_val = NaN; effect_size = NaN;
    end

    results = struct('n_hc', n_hc, 'n_path', n_path, ...
        'mode_hc', mode_hc, 'ci_hc', ci_hc, ...
        'mode_path', mode_path, 'ci_path', ci_path, ...
        'pct_hc_above', pct_hc_above, 'pct_path_above', pct_path_above, ...
        'p_value', p_val, 'effect_size', effect_size, ...
        'hc_ratios', hc_ratios, 'path_ratios', path_ratios);
end

%% ========================================================================
%  HELPERS
%  ========================================================================

function s = get_string_field(rec, field_name)
% GET_STRING_FIELD - Robustly extract a char value from the first row of a
% table column that may be stored as cell, string, or categorical.
    val = rec.(field_name)(1);
    if iscell(val)
        s = val{1};
    elseif isstring(val)
        s = char(val);
    elseif iscategorical(val)
        s = char(val);
    else
        s = char(val);
    end
end

function s = yesno(b)
    if b, s = 'YES'; else, s = 'NO'; end
end

function s = stars(p)
    if p < 0.001, s = ' ***';
    elseif p < 0.01, s = ' **';
    elseif p < 0.05, s = ' *';
    else, s = ' (ns)'; end
end