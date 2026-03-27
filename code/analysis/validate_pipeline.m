function results = validate_pipeline()
% VALIDATE_PIPELINE - Validate automatic annotation pipeline against gold standards
%
% Compares the automatic R-peak (Pan-Tompkins) and T-end (tangent method)
% detectors against expert manual annotations in two gold-standard databases:
%
%   1. LUDB  (Tier 1: full manual, fs = 500 Hz, lead II, n = 200)
%   2. QTDB  (Tier 1: full manual, fs = 250 Hz, q1c annotator, n ~ 105)
%
% Two levels of analysis are reported:
%   - Beat-level:    per-beat agreement (full transparency)
%   - Subject-level: per-subject median CDC agreement (analysis-relevant)
%
% The critical outcome is subject-level CDC agreement, since subject-level
% median CDC is what enters all downstream analyses.
%
% QTDB note: The QT Database contains 15-minute recordings but only
% 30-100 beats per record were manually annotated, concentrated in the
% final 5 minutes (Laguna et al., Comput. Cardiol. 1997). The automatic
% detector runs on the full signal but comparison is restricted to the
% annotated segment (+/- 2 seconds buffer).
%
% Outputs:
%   results - Struct with per-database and pooled agreement statistics
%   Console output suitable for Supplementary Table
%   CSV:    results/validation_results.csv
%   Figure: results/figures/SI_Fig_validation_bland_altman.{pdf,png,fig}
%
% Usage:
%   results = validate_pipeline();
%
% Tom Froese, OIST Embodied Cognitive Science Unit, March 2026

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('PIPELINE VALIDATION: Automatic vs Manual Annotations\n');
    fprintf('============================================================\n\n');

    paths = config();

    %% ================================================================
    %  VALIDATE LUDB
    %  ================================================================
    fprintf('────────────────────────────────────────────────────────────\n');
    fprintf('VALIDATING LUDB (Tier 1: full manual, fs=500 Hz)\n');
    fprintf('────────────────────────────────────────────────────────────\n');

    ludb = validate_ludb(paths);

    %% ================================================================
    %  VALIDATE QTDB
    %  ================================================================
    fprintf('\n────────────────────────────────────────────────────────────\n');
    fprintf('VALIDATING QTDB (Tier 1: full manual, fs=250 Hz)\n');
    fprintf('  Note: only 30-100 beats per record were manually annotated\n');
    fprintf('  (final 5 min of 15-min recordings; Laguna et al. 1997).\n');
    fprintf('  Auto detector restricted to annotated segment.\n');
    fprintf('────────────────────────────────────────────────────────────\n');

    qtdb = validate_qtdb(paths);

    %% ================================================================
    %  POOLED RESULTS
    %  ================================================================
    fprintf('\n────────────────────────────────────────────────────────────\n');
    fprintf('POOLED RESULTS\n');
    fprintf('────────────────────────────────────────────────────────────\n');

    pooled = pool_results(ludb, qtdb);

    %% ================================================================
    %  SUMMARY TABLES
    %  ================================================================
    fprintf('\n============================================================\n');
    fprintf('TABLE 1: Beat-level agreement (Auto vs Manual)\n');
    fprintf('============================================================\n\n');
    print_beat_level_table(ludb, qtdb, pooled);

    fprintf('\n============================================================\n');
    fprintf('TABLE 2: Subject-level agreement (median CDC per subject)\n');
    fprintf('============================================================\n\n');
    print_subject_level_table(ludb, qtdb);

    %% ================================================================
    %  BLAND-ALTMAN PLOTS
    %  ================================================================
    fprintf('\nGenerating Bland-Altman plots...\n');
    plot_bland_altman(ludb, qtdb, paths);

    %% ================================================================
    %  SAVE RESULTS
    %  ================================================================
    results.ludb   = ludb;
    results.qtdb   = qtdb;
    results.pooled = pooled;

    save(fullfile(paths.results, 'validation_results.mat'), 'results');
    save_results_csv(ludb, qtdb, pooled, paths);

    fprintf('\nValidation complete.\n');
    fprintf('  MAT: %s\n', fullfile(paths.results, 'validation_results.mat'));
    fprintf('  CSV: %s\n', fullfile(paths.results, 'validation_results.csv'));
end


%% ========================================================================
%  LUDB VALIDATION
%  ========================================================================
function res = validate_ludb(paths)

    base_path = paths.raw_ludb;
    fs = 500;
    lead_name = 'ii';
    n_subjects = 200;

    r_tol = round(0.075 * fs);
    t_tol = round(0.100 * fs);

    all_r_err = []; all_t_err = [];
    all_cdc_manual = []; all_cdc_auto = [];
    n_manual_r = 0; n_matched_r = 0; n_auto_r = 0;
    n_matched_t = 0; n_records = 0; n_auto_t_fail = 0;

    subj_cdc_manual = []; subj_cdc_auto = [];

    for subj = 1:n_subjects
        hea_file = fullfile(base_path, 'data', [num2str(subj) '.hea']);
        dat_file = fullfile(base_path, 'data', [num2str(subj) '.dat']);
        ann_file = fullfile(base_path, 'data', [num2str(subj) '.' lead_name]);

        if ~exist(hea_file, 'file') || ~exist(dat_file, 'file') || ...
           ~exist(ann_file, 'file')
            continue;
        end

        try
            [~, n_samp, n_sig, gains, baselines, sig_names, fmt] = ...
                read_wfdb_header(hea_file);
            lead_idx = find_lead_index(sig_names, lead_name);
            if isempty(lead_idx), continue; end

            raw = read_wfdb_signal(dat_file, n_samp, n_sig, fmt);
            ecg_raw = (raw(:, lead_idx) - baselines(lead_idx)) / gains(lead_idx);

            [ann_samples, ann_symbols] = read_wfdb_annotations(ann_file, 'ludb');
            [man_r, man_t] = extract_rt_pairs(ann_samples, ann_symbols);
            if length(man_r) < 3, continue; end

            ecg_filt = bandpass_filter(ecg_raw, fs, 0.5, 40);
            auto_r = detect_r_peaks(ecg_filt, fs);
            if isempty(auto_r), continue; end

            auto_t = nan(size(auto_r));
            for b = 1:length(auto_r)
                [te, ~, ~, ~] = detect_t_end(ecg_filt, auto_r, b, fs);
                if ~isnan(te)
                    auto_t(b) = te;
                else
                    n_auto_t_fail = n_auto_t_fail + 1;
                end
            end

            n_manual_r = n_manual_r + length(man_r);
            n_auto_r   = n_auto_r + length(auto_r);

            [r_err, t_err, cdc_man, cdc_aut, n_mr, n_mt] = ...
                match_and_compare(man_r, man_t, auto_r, auto_t, r_tol, t_tol, fs);

            all_r_err      = [all_r_err; r_err];
            all_t_err      = [all_t_err; t_err];
            all_cdc_manual = [all_cdc_manual; cdc_man];
            all_cdc_auto   = [all_cdc_auto;   cdc_aut];
            n_matched_r = n_matched_r + n_mr;
            n_matched_t = n_matched_t + n_mt;
            n_records   = n_records + 1;

            % --- Subject-level median CDC ---
            man_cdc_subj = compute_median_cdc(man_r, man_t, fs);
            aut_cdc_subj = compute_median_cdc_auto(auto_r, auto_t, fs);

            if ~isnan(man_cdc_subj) && ~isnan(aut_cdc_subj)
                subj_cdc_manual(end+1) = man_cdc_subj;
                subj_cdc_auto(end+1)   = aut_cdc_subj;
            end

        catch ME
            fprintf('  LUDB subj %d: %s\n', subj, ME.message);
        end
    end

    fprintf('  Records processed: %d\n', n_records);
    fprintf('  Manual R-peaks: %d, Auto R-peaks: %d\n', n_manual_r, n_auto_r);
    fprintf('  R-peak matches: %d (sensitivity: %.1f%%)\n', ...
        n_matched_r, 100 * n_matched_r / max(n_manual_r, 1));
    fprintf('  T-end matches (beat-level): %d, Auto T-end failures: %d\n', ...
        n_matched_t, n_auto_t_fail);

    res = compute_agreement(all_r_err, all_t_err, all_cdc_manual, all_cdc_auto);
    res.database = 'LUDB';
    res.n_records = n_records;
    res.n_manual_r = n_manual_r;
    res.n_auto_r = n_auto_r;
    res.n_matched_r = n_matched_r;
    res.n_matched_t = n_matched_t;
    res.r_sensitivity = 100 * n_matched_r / max(n_manual_r, 1);
    res.subj = compute_subject_agreement(subj_cdc_manual(:), subj_cdc_auto(:));
    res.subj_cdc_manual = subj_cdc_manual(:);
    res.subj_cdc_auto   = subj_cdc_auto(:);

    print_agreement(res);
end


%% ========================================================================
%  QTDB VALIDATION
%  ========================================================================
function res = validate_qtdb(paths)

    base_path = paths.raw_qtdb;
    fs = 250;
    records = get_qtdb_records();

    r_tol = round(0.075 * fs);
    t_tol = round(0.100 * fs);

    all_r_err = []; all_t_err = [];
    all_cdc_manual = []; all_cdc_auto = [];
    n_manual_r = 0; n_matched_r = 0; n_auto_r = 0;
    n_matched_t = 0; n_records = 0; n_auto_t_fail = 0;

    subj_cdc_manual = []; subj_cdc_auto = [];

    for i = 1:length(records)
        record = records{i};
        hea_file = fullfile(base_path, [record '.hea']);
        dat_file = fullfile(base_path, [record '.dat']);

        if ~exist(hea_file, 'file') || ~exist(dat_file, 'file')
            continue;
        end

        ann_file = '';
        for ann_ext = {'q1c', 'pu0', 'qt1', 'man'}
            candidate = fullfile(base_path, [record '.' ann_ext{1}]);
            if exist(candidate, 'file')
                ann_file = candidate;
                break;
            end
        end
        if isempty(ann_file), continue; end

        try
            [fs_rec, n_samp, n_sig, gains, baselines, sig_names, fmt] = ...
                read_wfdb_header(hea_file);
            lead_idx = 1;
            raw = read_wfdb_signal(dat_file, n_samp, n_sig, fmt);
            ecg_raw = (raw(:, lead_idx) - baselines(lead_idx)) / gains(lead_idx);

            if fs_rec > 0 && ~isnan(fs_rec)
                fs_actual = fs_rec;
            else
                fs_actual = fs;
            end

            if any(~isfinite(ecg_raw))
                fprintf('  QTDB %s: non-finite signal values, skipping\n', record);
                continue;
            end

            [ann_samples, ann_symbols] = read_wfdb_annotations(ann_file, 'qtdb');
            [man_r, man_t] = extract_rt_pairs(ann_samples, ann_symbols);
            if length(man_r) < 3, continue; end

            % --- RESTRICT to annotated segment (+/- 2 s buffer) ---
            seg_start = max(1, man_r(1) - round(2 * fs_actual));
            seg_end   = min(length(ecg_raw), man_r(end) + round(2 * fs_actual));

            % Run auto pipeline on full signal
            ecg_filt = bandpass_filter(ecg_raw, fs_actual, 0.5, 40);
            auto_r_full = detect_r_peaks(ecg_filt, fs_actual);
            if isempty(auto_r_full), continue; end

            % Keep only auto R-peaks within annotated segment
            in_segment = auto_r_full >= seg_start & auto_r_full <= seg_end;
            auto_r = auto_r_full(in_segment);
            if length(auto_r) < 3, continue; end

            % Detect T-ends using full context for RR intervals
            auto_t = nan(size(auto_r));
            seg_indices = find(in_segment);
            for k = 1:length(seg_indices)
                bi = seg_indices(k);
                [te, ~, ~, ~] = detect_t_end(ecg_filt, auto_r_full, bi, fs_actual);
                if ~isnan(te)
                    auto_t(k) = te;
                else
                    n_auto_t_fail = n_auto_t_fail + 1;
                end
            end

            n_manual_r = n_manual_r + length(man_r);
            n_auto_r   = n_auto_r + length(auto_r);

            [r_err, t_err, cdc_man, cdc_aut, n_mr, n_mt] = ...
                match_and_compare(man_r, man_t, auto_r, auto_t, ...
                                  r_tol, t_tol, fs_actual);

            all_r_err      = [all_r_err; r_err];
            all_t_err      = [all_t_err; t_err];
            all_cdc_manual = [all_cdc_manual; cdc_man];
            all_cdc_auto   = [all_cdc_auto;   cdc_aut];
            n_matched_r = n_matched_r + n_mr;
            n_matched_t = n_matched_t + n_mt;
            n_records   = n_records + 1;

            % Record-level median CDC (matched beats only)
            if length(cdc_man) >= 3 && length(cdc_aut) >= 3
                subj_cdc_manual(end+1) = median(cdc_man);
                subj_cdc_auto(end+1)   = median(cdc_aut);
            end

        catch ME
            fprintf('  QTDB %s: %s\n', record, ME.message);
        end
    end

    fprintf('  Records processed: %d\n', n_records);
    fprintf('  Manual R-peaks: %d, Auto R-peaks (in segment): %d\n', ...
        n_manual_r, n_auto_r);
    fprintf('  R-peak matches: %d (sensitivity: %.1f%%)\n', ...
        n_matched_r, 100 * n_matched_r / max(n_manual_r, 1));
    fprintf('  T-end matches (beat-level): %d, Auto T-end failures: %d\n', ...
        n_matched_t, n_auto_t_fail);

    res = compute_agreement(all_r_err, all_t_err, all_cdc_manual, all_cdc_auto);
    res.database = 'QTDB';
    res.n_records = n_records;
    res.n_manual_r = n_manual_r;
    res.n_auto_r = n_auto_r;
    res.n_matched_r = n_matched_r;
    res.n_matched_t = n_matched_t;
    res.r_sensitivity = 100 * n_matched_r / max(n_manual_r, 1);
    res.subj = compute_subject_agreement(subj_cdc_manual(:), subj_cdc_auto(:));
    res.subj_cdc_manual = subj_cdc_manual(:);
    res.subj_cdc_auto   = subj_cdc_auto(:);

    print_agreement(res);
end


%% ========================================================================
%  MATCHING
%  ========================================================================
function [r_err, t_err, cdc_man, cdc_aut, n_matched_r, n_matched_t] = ...
    match_and_compare(man_r, man_t, auto_r, auto_t, r_tol, ~, fs)

    r_err = []; t_err = [];
    cdc_man = []; cdc_aut = [];
    n_matched_r = 0; n_matched_t = 0;

    matched_auto_idx = nan(length(man_r), 1);

    for i = 1:length(man_r)
        diffs = abs(auto_r - man_r(i));
        [min_diff, best_idx] = min(diffs);
        if min_diff <= r_tol
            if ~ismember(best_idx, matched_auto_idx(1:i-1))
                matched_auto_idx(i) = best_idx;
                n_matched_r = n_matched_r + 1;
                r_err(end+1, 1) = (auto_r(best_idx) - man_r(i)) / fs * 1000;
            end
        end
    end

    for i = 1:length(man_r)
        ai = matched_auto_idx(i);
        if isnan(ai), continue; end
        if i > length(man_t) || isnan(man_t(i)), continue; end
        if isnan(auto_t(ai)), continue; end
        if i >= length(man_r), continue; end

        next_ai = matched_auto_idx(i + 1);
        if isnan(next_ai), continue; end

        man_rt = (man_t(i) - man_r(i)) / fs;
        man_rr = (man_r(i+1) - man_r(i)) / fs;
        aut_rt = (auto_t(ai) - auto_r(ai)) / fs;
        aut_rr = (auto_r(next_ai) - auto_r(ai)) / fs;

        if man_rr <= 0 || aut_rr <= 0, continue; end
        if man_rt <= 0 || aut_rt <= 0, continue; end
        if man_rt >= man_rr || aut_rt >= aut_rr, continue; end

        t_err(end+1, 1) = (auto_t(ai) - man_t(i)) / fs * 1000;
        cdc_man(end+1, 1) = man_rt / man_rr;
        cdc_aut(end+1, 1) = aut_rt / aut_rr;
        n_matched_t = n_matched_t + 1;
    end
end


%% ========================================================================
%  SUBJECT-LEVEL CDC HELPERS
%  ========================================================================
function med_cdc = compute_median_cdc(r_peaks, t_ends, fs)
% Median CDC from annotation pairs with quality filters
    n = min(length(r_peaks)-1, length(t_ends));
    if n < 2, med_cdc = NaN; return; end

    cdc_vals = [];
    for i = 1:n
        rt = (t_ends(i) - r_peaks(i)) / fs;
        rr = (r_peaks(i+1) - r_peaks(i)) / fs;
        if rr > 0 && rt > 0 && rt < rr
            c = rt / rr;
            if c > 0.2 && c < 0.6
                cdc_vals(end+1) = c;
            end
        end
    end
    if isempty(cdc_vals)
        med_cdc = NaN;
    else
        med_cdc = median(cdc_vals);
    end
end


function med_cdc = compute_median_cdc_auto(auto_r, auto_t, fs)
% Median CDC from automatic detections (all valid beats)
    valid = find(~isnan(auto_t));
    if length(valid) < 2, med_cdc = NaN; return; end

    cdc_vals = [];
    for k = 1:length(valid)-1
        bi = valid(k);
        % Find next auto R-peak (not necessarily in valid list)
        if bi < length(auto_r)
            rt = (auto_t(bi) - auto_r(bi)) / fs;
            rr = (auto_r(bi+1) - auto_r(bi)) / fs;
            if rr > 0 && rt > 0 && rt < rr
                c = rt / rr;
                if c > 0.2 && c < 0.6
                    cdc_vals(end+1) = c;
                end
            end
        end
    end
    if isempty(cdc_vals)
        med_cdc = NaN;
    else
        med_cdc = median(cdc_vals);
    end
end


%% ========================================================================
%  AGREEMENT STATISTICS
%  ========================================================================
function res = compute_agreement(r_err_ms, t_err_ms, cdc_man, cdc_aut)
    res.r_n = length(r_err_ms);
    res.r_mae = mean(abs(r_err_ms));
    res.r_bias = mean(r_err_ms);
    res.r_std = std(r_err_ms);
    res.r_loa_lo = res.r_bias - 1.96 * res.r_std;
    res.r_loa_hi = res.r_bias + 1.96 * res.r_std;

    res.t_n = length(t_err_ms);
    res.t_mae = mean(abs(t_err_ms));
    res.t_bias = mean(t_err_ms);
    res.t_std = std(t_err_ms);
    res.t_loa_lo = res.t_bias - 1.96 * res.t_std;
    res.t_loa_hi = res.t_bias + 1.96 * res.t_std;

    cdc_err = cdc_aut - cdc_man;
    res.cdc_n = length(cdc_err);
    res.cdc_mae = mean(abs(cdc_err));
    res.cdc_bias = mean(cdc_err);
    res.cdc_std = std(cdc_err);
    res.cdc_loa_lo = res.cdc_bias - 1.96 * res.cdc_std;
    res.cdc_loa_hi = res.cdc_bias + 1.96 * res.cdc_std;
    if res.cdc_n >= 3
        res.cdc_corr = corr(cdc_man, cdc_aut);
    else
        res.cdc_corr = NaN;
    end

    res.r_err_ms = r_err_ms;
    res.t_err_ms = t_err_ms;
    res.cdc_man = cdc_man;
    res.cdc_aut = cdc_aut;
end


function s = compute_subject_agreement(man_med, aut_med)
    s.n = length(man_med);
    err = aut_med - man_med;
    s.mae = mean(abs(err));
    s.bias = mean(err);
    s.std = std(err);
    s.loa_lo = s.bias - 1.96 * s.std;
    s.loa_hi = s.bias + 1.96 * s.std;
    s.median_ae = median(abs(err));
    if s.n >= 3
        s.corr = corr(man_med, aut_med);
    else
        s.corr = NaN;
    end
    fprintf('\n  --- Subject-level CDC (median per subject, N=%d) ---\n', s.n);
    fprintf('  MAE=%.4f, Bias=%+.4f, LoA=[%.4f, %.4f], r=%.3f\n', ...
        s.mae, s.bias, s.loa_lo, s.loa_hi, s.corr);
    fprintf('  Median AE=%.4f\n', s.median_ae);
end


%% ========================================================================
%  POOL RESULTS
%  ========================================================================
function pooled = pool_results(ludb, qtdb)
    r_err_all = [ludb.r_err_ms; qtdb.r_err_ms];
    t_err_all = [ludb.t_err_ms; qtdb.t_err_ms];
    cdc_man_all = [ludb.cdc_man; qtdb.cdc_man];
    cdc_aut_all = [ludb.cdc_aut; qtdb.cdc_aut];

    pooled.r_n = length(r_err_all);
    pooled.r_mae = mean(abs(r_err_all));
    pooled.r_bias = mean(r_err_all);
    pooled.r_std = std(r_err_all);
    pooled.r_loa_lo = pooled.r_bias - 1.96 * pooled.r_std;
    pooled.r_loa_hi = pooled.r_bias + 1.96 * pooled.r_std;

    pooled.t_n = length(t_err_all);
    pooled.t_mae = mean(abs(t_err_all));
    pooled.t_bias = mean(t_err_all);
    pooled.t_std = std(t_err_all);
    pooled.t_loa_lo = pooled.t_bias - 1.96 * pooled.t_std;
    pooled.t_loa_hi = pooled.t_bias + 1.96 * pooled.t_std;

    cdc_err_all = cdc_aut_all - cdc_man_all;
    pooled.cdc_n = length(cdc_err_all);
    pooled.cdc_mae = mean(abs(cdc_err_all));
    pooled.cdc_bias = mean(cdc_err_all);
    pooled.cdc_std = std(cdc_err_all);
    pooled.cdc_loa_lo = pooled.cdc_bias - 1.96 * pooled.cdc_std;
    pooled.cdc_loa_hi = pooled.cdc_bias + 1.96 * pooled.cdc_std;
    if pooled.cdc_n >= 3
        pooled.cdc_corr = corr(cdc_man_all, cdc_aut_all);
    else
        pooled.cdc_corr = NaN;
    end

    pooled.n_matched_r = ludb.n_matched_r + qtdb.n_matched_r;
    pooled.n_manual_r = ludb.n_manual_r + qtdb.n_manual_r;
    pooled.r_sensitivity = 100 * pooled.n_matched_r / max(pooled.n_manual_r, 1);
    pooled.n_matched_t = ludb.n_matched_t + qtdb.n_matched_t;

    pooled.r_err_ms = r_err_all;
    pooled.t_err_ms = t_err_all;
    pooled.cdc_man = cdc_man_all;
    pooled.cdc_aut = cdc_aut_all;

    subj_man_all = [ludb.subj_cdc_manual; qtdb.subj_cdc_manual];
    subj_aut_all = [ludb.subj_cdc_auto;   qtdb.subj_cdc_auto];
    pooled.subj = compute_subject_agreement(subj_man_all, subj_aut_all);
    pooled.subj_cdc_manual = subj_man_all;
    pooled.subj_cdc_auto   = subj_aut_all;
end


%% ========================================================================
%  BLAND-ALTMAN PLOTS (4 panels: 3 beat-level + 1 subject-level)
%  ========================================================================
function plot_bland_altman(ludb, qtdb, paths)

    fig = figure('Position', [50 100 1600 900], 'Color', 'w');

    col_ludb = [0.20, 0.47, 0.76];
    col_qtdb = [0.85, 0.40, 0.25];
    ms = 6;

    % ---- Panel a: R-peak ----
    subplot(2, 2, 1); hold on;
    n_l = length(ludb.r_err_ms);
    if n_l > 0
        scatter(1:n_l, ludb.r_err_ms, ms, col_ludb, 'filled', 'MarkerFaceAlpha', 0.4);
    end
    if ~isempty(qtdb.r_err_ms)
        scatter(n_l+(1:length(qtdb.r_err_ms)), qtdb.r_err_ms, ms, col_qtdb, ...
            'filled', 'MarkerFaceAlpha', 0.4);
    end
    all_r = [ludb.r_err_ms; qtdb.r_err_ms];
    add_ba_lines(all_r);
    ylabel('R-peak error (ms): Auto − Manual');
    xlabel('Beat index');
    title('a  R-peak timing (beat-level)', 'FontWeight', 'bold');
    set(gca, 'FontSize', 10, 'Box', 'on');
    annotate_panel(mean(all_r), mean(abs(all_r)), []);

    % ---- Panel b: T-end ----
    subplot(2, 2, 2); hold on;
    if ~isempty(ludb.t_err_ms)
        scatter(ludb.cdc_man, ludb.t_err_ms, ms, col_ludb, 'filled', 'MarkerFaceAlpha', 0.4);
    end
    if ~isempty(qtdb.t_err_ms)
        scatter(qtdb.cdc_man, qtdb.t_err_ms, ms, col_qtdb, 'filled', 'MarkerFaceAlpha', 0.4);
    end
    all_t = [ludb.t_err_ms; qtdb.t_err_ms];
    add_ba_lines(all_t);
    xlabel('Manual CDC');
    ylabel('T-end error (ms): Auto − Manual');
    title('b  T-end timing (beat-level)', 'FontWeight', 'bold');
    set(gca, 'FontSize', 10, 'Box', 'on');
    annotate_panel(mean(all_t), mean(abs(all_t)), []);

    % ---- Panel c: CDC beat-level ----
    subplot(2, 2, 3); hold on;
    cdc_man_all = [ludb.cdc_man; qtdb.cdc_man];
    cdc_aut_all = [ludb.cdc_aut; qtdb.cdc_aut];
    cdc_mean = (cdc_man_all + cdc_aut_all) / 2;
    cdc_diff = cdc_aut_all - cdc_man_all;
    n_l_c = length(ludb.cdc_man);
    if n_l_c > 0
        scatter(cdc_mean(1:n_l_c), cdc_diff(1:n_l_c), ms, col_ludb, ...
            'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'LUDB');
    end
    if ~isempty(qtdb.cdc_man)
        scatter(cdc_mean(n_l_c+1:end), cdc_diff(n_l_c+1:end), ms, col_qtdb, ...
            'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'QTDB');
    end
    add_ba_lines(cdc_diff);
    xlabel('Mean CDC (Manual + Auto) / 2');
    ylabel('\DeltaCDC: Auto − Manual');
    title('c  CDC (beat-level)', 'FontWeight', 'bold');
    legend('Location', 'southeast', 'FontSize', 8);
    set(gca, 'FontSize', 10, 'Box', 'on');
    annotate_panel(mean(cdc_diff), mean(abs(cdc_diff)), corr(cdc_man_all, cdc_aut_all));

    % ---- Panel d: CDC subject-level ----
    subplot(2, 2, 4); hold on;
    subj_man_all = [ludb.subj_cdc_manual; qtdb.subj_cdc_manual];
    subj_aut_all = [ludb.subj_cdc_auto;   qtdb.subj_cdc_auto];
    subj_mean = (subj_man_all + subj_aut_all) / 2;
    subj_diff = subj_aut_all - subj_man_all;
    n_l_s = length(ludb.subj_cdc_manual);
    if n_l_s > 0
        scatter(subj_mean(1:n_l_s), subj_diff(1:n_l_s), 25, col_ludb, ...
            'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'LUDB');
    end
    if ~isempty(qtdb.subj_cdc_manual)
        scatter(subj_mean(n_l_s+1:end), subj_diff(n_l_s+1:end), 25, col_qtdb, ...
            'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'QTDB');
    end
    add_ba_lines(subj_diff);
    xlabel('Mean median CDC (Manual + Auto) / 2');
    ylabel('\DeltaCDC: Auto − Manual');
    title('d  CDC (subject-level median)', 'FontWeight', 'bold');
    legend('Location', 'southeast', 'FontSize', 8);
    set(gca, 'FontSize', 10, 'Box', 'on');
    annotate_panel(mean(subj_diff), mean(abs(subj_diff)), corr(subj_man_all, subj_aut_all));

    % --- Save with proper page sizing ---
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0 0 40 22]);
    set(fig, 'PaperSize', [40 22]);

    out_pdf = fullfile(paths.figures, 'SI_Fig_validation_bland_altman.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig_validation_bland_altman.png');
    out_fig = fullfile(paths.figures, 'SI_Fig_validation_bland_altman.fig');

    print(fig, out_pdf, '-dpdf', '-painters');
    print(fig, out_png, '-dpng', '-r300');
    savefig(fig, out_fig);

    fprintf('  Saved: %s  (vector)\n', out_pdf);
    fprintf('  Saved: %s  (raster, 300 dpi)\n', out_png);
    fprintf('  Saved: %s  (editable)\n', out_fig);
end


function add_ba_lines(data)
    yline(mean(data), 'k-', 'LineWidth', 1.2);
    yline(mean(data) + 1.96*std(data), 'k--', 'LineWidth', 0.8);
    yline(mean(data) - 1.96*std(data), 'k--', 'LineWidth', 0.8);
    yline(0, 'Color', [0.6 0.6 0.6], 'LineStyle', ':');
end


function annotate_panel(bias_val, mae_val, r_val)
    if abs(bias_val) < 1
        text(0.03, 0.95, sprintf('Bias: %+.4f', bias_val), ...
            'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top');
        text(0.03, 0.87, sprintf('MAE: %.4f', mae_val), ...
            'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top');
    else
        text(0.03, 0.95, sprintf('Bias: %+.1f ms', bias_val), ...
            'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top');
        text(0.03, 0.87, sprintf('MAE: %.1f ms', mae_val), ...
            'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top');
    end
    if ~isempty(r_val) && ~isnan(r_val)
        text(0.03, 0.79, sprintf('r = %.3f', r_val), ...
            'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top');
    end
end


%% ========================================================================
%  PRINT FUNCTIONS
%  ========================================================================
function print_agreement(res)
    fprintf('\n  --- %s Beat-level Agreement ---\n', res.database);
    fprintf('  R-peak:  N=%d, MAE=%.1f ms, Bias=%+.1f ms, LoA=[%.1f, %.1f] ms\n', ...
        res.r_n, res.r_mae, res.r_bias, res.r_loa_lo, res.r_loa_hi);
    fprintf('  T-end:   N=%d, MAE=%.1f ms, Bias=%+.1f ms, LoA=[%.1f, %.1f] ms\n', ...
        res.t_n, res.t_mae, res.t_bias, res.t_loa_lo, res.t_loa_hi);
    fprintf('  CDC:     N=%d, MAE=%.4f, Bias=%+.4f, LoA=[%.4f, %.4f], r=%.3f\n', ...
        res.cdc_n, res.cdc_mae, res.cdc_bias, res.cdc_loa_lo, res.cdc_loa_hi, res.cdc_corr);
    fprintf('  R-peak sensitivity: %.1f%%\n', res.r_sensitivity);
end


function print_beat_level_table(ludb, qtdb, pooled)
    fprintf('%-12s | %6s | %8s | %8s | %20s | %6s\n', ...
        'Metric', 'N', 'MAE', 'Bias', '95%% LoA', 'r');
    fprintf('%s\n', repmat('-', 1, 72));
    fprintf('%-12s\n', 'R-PEAK');
    pr('LUDB',   ludb.r_n,   ludb.r_mae,   ludb.r_bias,   ludb.r_loa_lo,   ludb.r_loa_hi,   NaN, 'ms');
    pr('QTDB',   qtdb.r_n,   qtdb.r_mae,   qtdb.r_bias,   qtdb.r_loa_lo,   qtdb.r_loa_hi,   NaN, 'ms');
    pr('Pooled', pooled.r_n, pooled.r_mae, pooled.r_bias, pooled.r_loa_lo, pooled.r_loa_hi, NaN, 'ms');
    fprintf('%s\n', repmat('-', 1, 72));
    fprintf('%-12s\n', 'T-END');
    pr('LUDB',   ludb.t_n,   ludb.t_mae,   ludb.t_bias,   ludb.t_loa_lo,   ludb.t_loa_hi,   NaN, 'ms');
    pr('QTDB',   qtdb.t_n,   qtdb.t_mae,   qtdb.t_bias,   qtdb.t_loa_lo,   qtdb.t_loa_hi,   NaN, 'ms');
    pr('Pooled', pooled.t_n, pooled.t_mae, pooled.t_bias, pooled.t_loa_lo, pooled.t_loa_hi, NaN, 'ms');
    fprintf('%s\n', repmat('-', 1, 72));
    fprintf('%-12s\n', 'CDC');
    pr('LUDB',   ludb.cdc_n,   ludb.cdc_mae,   ludb.cdc_bias,   ludb.cdc_loa_lo,   ludb.cdc_loa_hi,   ludb.cdc_corr,   'cdc');
    pr('QTDB',   qtdb.cdc_n,   qtdb.cdc_mae,   qtdb.cdc_bias,   qtdb.cdc_loa_lo,   qtdb.cdc_loa_hi,   qtdb.cdc_corr,   'cdc');
    pr('Pooled', pooled.cdc_n, pooled.cdc_mae, pooled.cdc_bias, pooled.cdc_loa_lo, pooled.cdc_loa_hi, pooled.cdc_corr, 'cdc');
    fprintf('%s\n', repmat('-', 1, 72));
    fprintf('\nR-peak sensitivity: LUDB=%.1f%%, QTDB=%.1f%%, Pooled=%.1f%%\n', ...
        ludb.r_sensitivity, qtdb.r_sensitivity, pooled.r_sensitivity);
end

function pr(label, n, mae, bias, loa_lo, loa_hi, r_val, unit)
    if strcmp(unit, 'ms')
        fprintf('  %-10s | %6d | %6.1f ms | %+6.1f ms | [%+6.1f, %+6.1f] ms |      -\n', ...
            label, n, mae, bias, loa_lo, loa_hi);
    else
        fprintf('  %-10s | %6d | %8.4f | %+8.4f | [%+7.4f, %+7.4f] | %6.3f\n', ...
            label, n, mae, bias, loa_lo, loa_hi, r_val);
    end
end

function print_subject_level_table(ludb, qtdb)
    fprintf('%-12s | %6s | %8s | %8s | %20s | %6s | %8s\n', ...
        'Database', 'N', 'MAE', 'Bias', '95%% LoA', 'r', 'Med AE');
    fprintf('%s\n', repmat('-', 1, 80));
    s = ludb.subj;
    fprintf('  %-10s | %6d | %8.4f | %+8.4f | [%+7.4f, %+7.4f] | %6.3f | %8.4f\n', ...
        'LUDB', s.n, s.mae, s.bias, s.loa_lo, s.loa_hi, s.corr, s.median_ae);
    s = qtdb.subj;
    fprintf('  %-10s | %6d | %8.4f | %+8.4f | [%+7.4f, %+7.4f] | %6.3f | %8.4f\n', ...
        'QTDB', s.n, s.mae, s.bias, s.loa_lo, s.loa_hi, s.corr, s.median_ae);
    fprintf('%s\n', repmat('-', 1, 80));
    fprintf('\nContext: between-group CDC differences from hierarchical model:\n');
    fprintf('  HC vs CN:   Delta_CDC ~ 0.040  (subject-level MAE should be << this)\n');
    fprintf('  HC vs Path: Delta_CDC ~ 0.088  (subject-level MAE should be << this)\n');
end


%% ========================================================================
%  SAVE CSV
%  ========================================================================
function save_results_csv(ludb, qtdb, pooled, paths)
    fid = fopen(fullfile(paths.results, 'validation_results.csv'), 'w');
    fprintf(fid, 'level,database,metric,n,mae,bias,loa_lo,loa_hi,corr,sensitivity_pct\n');

    % Beat-level
    fprintf(fid, 'beat,LUDB,R-peak,%d,%.2f,%.2f,%.2f,%.2f,,%.1f\n', ludb.r_n, ludb.r_mae, ludb.r_bias, ludb.r_loa_lo, ludb.r_loa_hi, ludb.r_sensitivity);
    fprintf(fid, 'beat,LUDB,T-end,%d,%.2f,%.2f,%.2f,%.2f,,\n', ludb.t_n, ludb.t_mae, ludb.t_bias, ludb.t_loa_lo, ludb.t_loa_hi);
    fprintf(fid, 'beat,LUDB,CDC,%d,%.5f,%.5f,%.5f,%.5f,%.4f,\n', ludb.cdc_n, ludb.cdc_mae, ludb.cdc_bias, ludb.cdc_loa_lo, ludb.cdc_loa_hi, ludb.cdc_corr);
    fprintf(fid, 'beat,QTDB,R-peak,%d,%.2f,%.2f,%.2f,%.2f,,%.1f\n', qtdb.r_n, qtdb.r_mae, qtdb.r_bias, qtdb.r_loa_lo, qtdb.r_loa_hi, qtdb.r_sensitivity);
    fprintf(fid, 'beat,QTDB,T-end,%d,%.2f,%.2f,%.2f,%.2f,,\n', qtdb.t_n, qtdb.t_mae, qtdb.t_bias, qtdb.t_loa_lo, qtdb.t_loa_hi);
    fprintf(fid, 'beat,QTDB,CDC,%d,%.5f,%.5f,%.5f,%.5f,%.4f,\n', qtdb.cdc_n, qtdb.cdc_mae, qtdb.cdc_bias, qtdb.cdc_loa_lo, qtdb.cdc_loa_hi, qtdb.cdc_corr);
    fprintf(fid, 'beat,Pooled,R-peak,%d,%.2f,%.2f,%.2f,%.2f,,%.1f\n', pooled.r_n, pooled.r_mae, pooled.r_bias, pooled.r_loa_lo, pooled.r_loa_hi, pooled.r_sensitivity);
    fprintf(fid, 'beat,Pooled,T-end,%d,%.2f,%.2f,%.2f,%.2f,,\n', pooled.t_n, pooled.t_mae, pooled.t_bias, pooled.t_loa_lo, pooled.t_loa_hi);
    fprintf(fid, 'beat,Pooled,CDC,%d,%.5f,%.5f,%.5f,%.5f,%.4f,\n', pooled.cdc_n, pooled.cdc_mae, pooled.cdc_bias, pooled.cdc_loa_lo, pooled.cdc_loa_hi, pooled.cdc_corr);

    % Subject-level
    s = ludb.subj;
    fprintf(fid, 'subject,LUDB,CDC,%d,%.5f,%.5f,%.5f,%.5f,%.4f,\n', s.n, s.mae, s.bias, s.loa_lo, s.loa_hi, s.corr);
    s = qtdb.subj;
    fprintf(fid, 'subject,QTDB,CDC,%d,%.5f,%.5f,%.5f,%.5f,%.4f,\n', s.n, s.mae, s.bias, s.loa_lo, s.loa_hi, s.corr);
    s = pooled.subj;
    fprintf(fid, 'subject,Pooled,CDC,%d,%.5f,%.5f,%.5f,%.5f,%.4f,\n', s.n, s.mae, s.bias, s.loa_lo, s.loa_hi, s.corr);

    fclose(fid);
end


%% ========================================================================
%  HELPERS
%  ========================================================================
function idx = find_lead_index(sig_names, lead_name)
    idx = [];
    for i = 1:length(sig_names)
        if strcmpi(strtrim(sig_names{i}), lower(lead_name))
            idx = i; return;
        end
    end
    variants = {lead_name, upper(lead_name), 'II', 'ii'};
    for v = 1:length(variants)
        for i = 1:length(sig_names)
            if strcmpi(strtrim(sig_names{i}), variants{v})
                idx = i; return;
            end
        end
    end
    if isempty(idx) && ~isempty(sig_names), idx = 1; end
end

function records = get_qtdb_records()
    records = [{'sel100','sel102','sel103','sel104','sel114','sel116', ...
                'sel117','sel123','sel213','sel221','sel223','sel230', ...
                'sel231','sel232','sel233'}, ...  % arrhythmia
               {'sel301','sel302','sel306','sel307','sel308','sel310'}, ...  % ST change
               {'sel803','sel808','sel811','sel820','sel821','sel840', ...
                'sel847','sel853','sel871','sel872','sel873','sel883','sel891'}, ...  % SVA
               {'sel16265','sel16272','sel16273','sel16420','sel16483', ...
                'sel16539','sel16773','sel16786','sel16795','sel17453'}, ...  % healthy
               {'sele0104','sele0106','sele0107','sele0110','sele0111', ...
                'sele0112','sele0114','sele0116','sele0121','sele0122', ...
                'sele0124','sele0126','sele0129','sele0133','sele0136', ...
                'sele0166','sele0170','sele0203','sele0210','sele0211', ...
                'sele0303','sele0405','sele0406','sele0409','sele0411', ...
                'sele0509','sele0603','sele0604','sele0606','sele0607', ...
                'sele0609','sele0612','sele0704'}, ...  % European ST-T
               {'sel30','sel31','sel32','sel33','sel34','sel35', ...
                'sel36','sel37','sel38','sel39','sel40','sel41', ...
                'sel42','sel43','sel44','sel45','sel46','sel47', ...
                'sel48','sel49','sel50','sel51','sel52','sel17152'}, ...  % sudden death
               {'sel14046','sel14157','sel14172','sel15814'}];  % long-term
end
