function results = analyze_gold_standard_validation()
% ANALYZE_GOLD_STANDARD_VALIDATION - Validate automatic annotation pipeline
%   against expert manual annotations in gold-standard databases
%
% Compares the automatic R-peak (Pan-Tompkins) and T-end (tangent method)
% detectors against expert manual annotations in two gold-standard databases:
%
%   1. LUDB  (Tier 1: full manual, fs = 500 Hz, lead II, n = 200)
%   2. QTDB  (Tier 1: full manual, fs = 250 Hz, q1c annotator, n ~ 105)
%
% Two levels of analysis:
%   - Beat-level:    per-beat agreement (full transparency)
%   - Subject-level: per-subject median CDC agreement (analysis-relevant)
%
% QTDB note: Each 15-min recording has only 30-100 manually annotated beats,
% concentrated in the final 5 minutes (Laguna et al., Comput. Cardiol. 1997).
% The auto detector runs on the full signal but comparison is restricted to
% the annotated segment (+/- 30 s buffer).
%
% Outputs:
%   results/validation_results.mat   - Full results struct
%   results/validation_results.csv   - Summary table for reproducibility
%
% Companion figure: plot_SI_Fig3.m (Bland-Altman plots)
%
% Tom Froese, OIST Embodied Cognitive Science Unit, March 2026

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('PIPELINE VALIDATION: Automatic vs Manual Annotations\n');
    fprintf('============================================================\n\n');

    paths = config();

    % --- LUDB ---
    fprintf('────────────────────────────────────────────────────────────\n');
    fprintf('VALIDATING LUDB (Tier 1: full manual, fs=500 Hz)\n');
    fprintf('────────────────────────────────────────────────────────────\n');
    ludb = validate_ludb(paths);

    % --- QTDB ---
    fprintf('\n────────────────────────────────────────────────────────────\n');
    fprintf('VALIDATING QTDB (Tier 1: full manual, fs=250 Hz)\n');
    fprintf('  Note: 30-100 beats/record manually annotated (final 5 min).\n');
    fprintf('  Auto detector restricted to annotated segment.\n');
    fprintf('────────────────────────────────────────────────────────────\n');
    qtdb = validate_qtdb(paths);

    % --- POOLED ---
    fprintf('\n────────────────────────────────────────────────────────────\n');
    fprintf('POOLED RESULTS\n');
    fprintf('────────────────────────────────────────────────────────────\n');
    pooled = pool_results(ludb, qtdb);

    % --- TABLES ---
    fprintf('\n============================================================\n');
    fprintf('TABLE 1: Beat-level agreement\n');
    fprintf('============================================================\n\n');
    print_beat_level_table(ludb, qtdb, pooled);

    fprintf('\n============================================================\n');
    fprintf('TABLE 2: Subject-level agreement (median CDC per subject)\n');
    fprintf('============================================================\n\n');
    print_subject_level_table(ludb, qtdb, pooled);

    % --- SAVE ---
    results.ludb   = ludb;
    results.qtdb   = qtdb;
    results.pooled = pooled;

    save(fullfile(paths.results, 'validation_results.mat'), 'results');
    save_results_csv(ludb, qtdb, pooled, paths);

    fprintf('\nValidation complete.\n');
    fprintf('  MAT: %s\n', fullfile(paths.results, 'validation_results.mat'));
    fprintf('  CSV: %s\n', fullfile(paths.results, 'validation_results.csv'));
    fprintf('\nTo generate the Bland-Altman figure, run: plot_SI_Fig3()\n');
end


%% ========================================================================
%  LUDB VALIDATION
%  ========================================================================
function res = validate_ludb(paths)

    base_path = paths.raw_ludb;
    fs = 500;
    lead_name = 'ii';
    n_subjects = 200;

    r_tol = round(0.075 * fs);   % +/-75 ms

    % Accumulators
    all_r_err = []; all_t_err = [];
    all_cdc_man = []; all_cdc_aut = [];
    n_man_r = 0; n_aut_r = 0; n_match_r = 0;
    n_match_t = 0; n_recs = 0; n_tfail = 0;
    subj_man = []; subj_aut = [];

    for subj = 1:n_subjects
        hea_file = fullfile(base_path, 'data', [num2str(subj) '.hea']);
        dat_file = fullfile(base_path, 'data', [num2str(subj) '.dat']);
        ann_file = fullfile(base_path, 'data', [num2str(subj) '.' lead_name]);

        if ~exist(hea_file, 'file') || ~exist(dat_file, 'file') || ...
           ~exist(ann_file, 'file')
            continue;
        end

        try
            % Load signal
            [~, n_samp, n_sig, gains, baselines, sig_names, fmt] = ...
                read_wfdb_header(hea_file);
            li = find_lead(sig_names, lead_name);
            if isempty(li), continue; end
            raw = read_wfdb_signal(dat_file, n_samp, n_sig, fmt);
            ecg = (raw(:, li) - baselines(li)) / gains(li);

            % Load manual annotations
            [as, ay] = read_wfdb_annotations(ann_file, 'ludb');
            [mr, mt] = extract_rt_pairs(as, ay);
            if length(mr) < 3, continue; end

            % Run auto pipeline
            ef = bandpass_filter(ecg, fs, 0.5, 40);
            ar = detect_r_peaks(ef, fs);
            if isempty(ar), continue; end

            at = nan(size(ar));
            for b = 1:length(ar)
                [te, ~, ~, ~] = detect_t_end(ef, ar, b, fs);
                if ~isnan(te), at(b) = te; else, n_tfail = n_tfail + 1; end
            end

            % Beat-level comparison
            n_man_r = n_man_r + length(mr);
            n_aut_r = n_aut_r + length(ar);
            [re, te2, cm, ca, nmr, nmt] = ...
                do_match(mr, mt, ar, at, r_tol, fs);
            all_r_err = [all_r_err; re];
            all_t_err = [all_t_err; te2];
            all_cdc_man = [all_cdc_man; cm];
            all_cdc_aut = [all_cdc_aut; ca];
            n_match_r = n_match_r + nmr;
            n_match_t = n_match_t + nmt;
            n_recs = n_recs + 1;

            % Subject-level median CDC
            sm = med_cdc(mr, mt, fs);
            sa = med_cdc_auto(ar, at, fs);
            if ~isnan(sm) && ~isnan(sa)
                subj_man(end+1) = sm;
                subj_aut(end+1) = sa;
            end

        catch ME
            fprintf('  LUDB subj %d: %s\n', subj, ME.message);
        end
    end

    fprintf('  Records: %d\n', n_recs);
    fprintf('  Manual R: %d, Auto R: %d, Matched: %d (%.1f%%)\n', ...
        n_man_r, n_aut_r, n_match_r, 100*n_match_r/max(n_man_r,1));
    fprintf('  T-end matches: %d, T-end failures: %d\n', n_match_t, n_tfail);

    res = make_res('LUDB', all_r_err, all_t_err, all_cdc_man, all_cdc_aut, ...
        n_recs, n_man_r, n_aut_r, n_match_r, n_match_t, subj_man(:), subj_aut(:));
    print_res(res);
end


%% ========================================================================
%  QTDB VALIDATION
%  ========================================================================
function res = validate_qtdb(paths)

    base_path = paths.raw_qtdb;
    fs = 250;
    records = qtdb_records();

    r_tol = round(0.075 * fs);

    all_r_err = []; all_t_err = [];
    all_cdc_man = []; all_cdc_aut = [];
    n_man_r = 0; n_aut_r = 0; n_match_r = 0;
    n_match_t = 0; n_recs = 0; n_tfail = 0;
    subj_man = []; subj_aut = [];

    fprintf('  Processing %d QTDB records (annotated segment only + 30 s buffer)...\n', length(records));

    for i = 1:length(records)
        rec = records{i};
        fprintf('    [%3d/%d] %s ... ', i, length(records), rec);

        hea_file = fullfile(base_path, [rec '.hea']);
        dat_file = fullfile(base_path, [rec '.dat']);
        if ~exist(hea_file, 'file') || ~exist(dat_file, 'file')
            fprintf('missing files → skip\n');
            continue;
        end

        % Find annotation file
        af = '';
        for ext = {'q1c', 'pu0', 'qt1', 'man'}
            c = fullfile(base_path, [rec '.' ext{1}]);
            if exist(c, 'file'), af = c; break; end
        end
        if isempty(af)
            fprintf('no annotations → skip\n');
            continue;
        end

        try
            [fsr, n_samp, n_sig, gains, baselines, ~, fmt] = read_wfdb_header(hea_file);
            if fsr > 0 && ~isnan(fsr)
                fsa = fsr;
            else
                fsa = fs;
            end

            raw = read_wfdb_signal(dat_file, n_samp, n_sig, fmt);
            ecg = (raw(:,1) - baselines(1)) / gains(1);

            if any(~isfinite(ecg))
                fprintf('non-finite → skip\n');
                continue;
            end

            [as, ay] = read_wfdb_annotations(af, 'qtdb');
            [mr, mt] = extract_rt_pairs(as, ay);
            if length(mr) < 3
                fprintf('too few manual beats → skip\n');
                continue;
            end

            % Crop to annotated segment + buffer
            buffer = round(30 * fsa);
            seg_s = max(1, mr(1) - buffer);
            seg_e = min(length(ecg), mr(end) + buffer);
            ecg_seg = ecg(seg_s:seg_e);

            ef = bandpass_filter(ecg_seg, fsa, 0.5, 40);
            ar_seg = detect_r_peaks(ef, fsa);   % LOCAL indices in cropped ef

            if isempty(ar_seg)
                fprintf('no R-peaks → skip\n');
                continue;
            end

            ar_full = ar_seg + (seg_s - 1);     % global indices

            in_seg = ar_full >= mr(1) & ar_full <= mr(end);
            ar = ar_full(in_seg);
            if length(ar) < 3
                fprintf('too few auto beats → skip\n');
                continue;
            end

            % T-end detection using correct local indices
            at = nan(size(ar));
            seg_idx = find(in_seg);
            for k = 1:length(ar)
                local_bi = seg_idx(k);
                [te_local, ~, ~, ~] = detect_t_end(ef, ar_seg, local_bi, fsa);
                if ~isnan(te_local)
                    at(k) = te_local + (seg_s - 1);       % store global index
                else
                    n_tfail = n_tfail + 1;
                end
            end

            % Beat-level comparison
            n_man_r = n_man_r + length(mr);
            n_aut_r = n_aut_r + length(ar);
            [re, te2, cm, ca, nmr, nmt] = do_match(mr, mt, ar, at, r_tol, fsa);

            all_r_err = [all_r_err; re];
            all_t_err = [all_t_err; te2];
            all_cdc_man = [all_cdc_man; cm];
            all_cdc_aut = [all_cdc_aut; ca];
            n_match_r = n_match_r + nmr;
            n_match_t = n_match_t + nmt;
            n_recs = n_recs + 1;

            if length(cm) >= 3 && length(ca) >= 3
                subj_man(end+1) = median(cm);
                subj_aut(end+1) = median(ca);
            end

            fprintf('OK (%d matched beats, %d T-ends)\n', nmr, nmt);

        catch ME
            fprintf('ERROR: %s\n', ME.message);
        end
    end

    fprintf('  Records: %d\n', n_recs);
    fprintf('  Manual R: %d, Auto R (segment): %d, Matched: %d (%.1f%%)\n', ...
        n_man_r, n_aut_r, n_match_r, 100*n_match_r/max(n_man_r,1));
    fprintf('  T-end matches: %d, T-end failures: %d\n', n_match_t, n_tfail);

    res = make_res('QTDB', all_r_err, all_t_err, all_cdc_man, all_cdc_aut, ...
        n_recs, n_man_r, n_aut_r, n_match_r, n_match_t, subj_man(:), subj_aut(:));
    print_res(res);
end


%% ========================================================================
%  BEAT MATCHING
%  ========================================================================
function [r_err, t_err, cdc_m, cdc_a, nm_r, nm_t] = ...
    do_match(man_r, man_t, aut_r, aut_t, r_tol, fs)

    r_err = []; t_err = []; cdc_m = []; cdc_a = [];
    nm_r = 0; nm_t = 0;
    matched = nan(length(man_r), 1);

    % Match each manual R to nearest auto R within tolerance
    for i = 1:length(man_r)
        d = abs(aut_r - man_r(i));
        [md, bi] = min(d);
        if md <= r_tol && ~ismember(bi, matched(1:i-1))
            matched(i) = bi;
            nm_r = nm_r + 1;
            r_err(end+1, 1) = (aut_r(bi) - man_r(i)) / fs * 1000;
        end
    end

    % Compare T-ends and CDC for consecutive matched pairs
    for i = 1:length(man_r)-1
        ai = matched(i);
        ai_next = matched(i+1);
        if isnan(ai) || isnan(ai_next), continue; end
        if i > length(man_t) || isnan(man_t(i)), continue; end
        if isnan(aut_t(ai)), continue; end

        m_rt = (man_t(i) - man_r(i)) / fs;
        m_rr = (man_r(i+1) - man_r(i)) / fs;
        a_rt = (aut_t(ai) - aut_r(ai)) / fs;
        a_rr = (aut_r(ai_next) - aut_r(ai)) / fs;

        if m_rr <= 0 || a_rr <= 0 || m_rt <= 0 || a_rt <= 0, continue; end
        if m_rt >= m_rr || a_rt >= a_rr, continue; end

        t_err(end+1, 1) = (aut_t(ai) - man_t(i)) / fs * 1000;
        cdc_m(end+1, 1) = m_rt / m_rr;
        cdc_a(end+1, 1) = a_rt / a_rr;
        nm_t = nm_t + 1;
    end
end


%% ========================================================================
%  SUBJECT-LEVEL CDC HELPERS
%  ========================================================================
function mc = med_cdc(rp, te, fs)
    n = min(length(rp)-1, length(te));
    if n < 2, mc = NaN; return; end
    v = [];
    for i = 1:n
        rt = (te(i) - rp(i)) / fs;
        rr = (rp(i+1) - rp(i)) / fs;
        if rr > 0 && rt > 0 && rt < rr
            c = rt / rr;
            if c > 0.2 && c < 0.6, v(end+1) = c; end
        end
    end
    if isempty(v), mc = NaN; else, mc = median(v); end
end

function mc = med_cdc_auto(ar, at, fs)
    vi = find(~isnan(at));
    if length(vi) < 2, mc = NaN; return; end
    v = [];
    for k = 1:length(vi)-1
        bi = vi(k);
        if bi < length(ar)
            rt = (at(bi) - ar(bi)) / fs;
            rr = (ar(bi+1) - ar(bi)) / fs;
            if rr > 0 && rt > 0 && rt < rr
                c = rt / rr;
                if c > 0.2 && c < 0.6, v(end+1) = c; end
            end
        end
    end
    if isempty(v), mc = NaN; else, mc = median(v); end
end


%% ========================================================================
%  AGREEMENT STATS
%  ========================================================================
function res = make_res(db, r_err, t_err, cm, ca, nr, nmr, nar, nxr, nxt, sm, sa)
    res.database = db;
    res.n_records = nr;
    res.n_manual_r = nmr; res.n_auto_r = nar;
    res.n_matched_r = nxr; res.n_matched_t = nxt;
    res.r_sensitivity = 100 * nxr / max(nmr, 1);

    res.r_n = length(r_err); res.r_mae = mn(abs(r_err)); res.r_bias = mn(r_err);
    res.r_std = sd(r_err);
    res.r_loa_lo = res.r_bias - 1.96*res.r_std;
    res.r_loa_hi = res.r_bias + 1.96*res.r_std;

    res.t_n = length(t_err); res.t_mae = mn(abs(t_err)); res.t_bias = mn(t_err);
    res.t_std = sd(t_err);
    res.t_loa_lo = res.t_bias - 1.96*res.t_std;
    res.t_loa_hi = res.t_bias + 1.96*res.t_std;

    ce = ca - cm;
    res.cdc_n = length(ce); res.cdc_mae = mn(abs(ce)); res.cdc_bias = mn(ce);
    res.cdc_std = sd(ce);
    res.cdc_loa_lo = res.cdc_bias - 1.96*res.cdc_std;
    res.cdc_loa_hi = res.cdc_bias + 1.96*res.cdc_std;
    if res.cdc_n >= 3, res.cdc_corr = corr(cm, ca); else, res.cdc_corr = NaN; end

    res.r_err_ms = r_err; res.t_err_ms = t_err;
    res.cdc_man = cm; res.cdc_aut = ca;

    % Subject-level
    res.subj_cdc_manual = sm; res.subj_cdc_auto = sa;
    se = sa - sm;
    res.subj_n = length(se); res.subj_mae = mn(abs(se)); res.subj_bias = mn(se);
    res.subj_std = sd(se); res.subj_median_ae = mdn(abs(se));
    res.subj_loa_lo = res.subj_bias - 1.96*res.subj_std;
    res.subj_loa_hi = res.subj_bias + 1.96*res.subj_std;
    if res.subj_n >= 3, res.subj_corr = corr(sm, sa); else, res.subj_corr = NaN; end
end

function v = mn(x)
    if isempty(x), v = NaN; else, v = mean(x); end
end

function v = sd(x)
    if isempty(x), v = NaN; else, v = std(x); end
end

function v = mdn(x)
    if isempty(x), v = NaN; else, v = median(x); end
end


%% ========================================================================
%  POOL
%  ========================================================================
function p = pool_results(a, b)
    re = [a.r_err_ms; b.r_err_ms];
    te = [a.t_err_ms; b.t_err_ms];
    cm = [a.cdc_man; b.cdc_man];
    ca = [a.cdc_aut; b.cdc_aut];
    sm = [a.subj_cdc_manual; b.subj_cdc_manual];
    sa = [a.subj_cdc_auto;   b.subj_cdc_auto];
    p = make_res('Pooled', re, te, cm, ca, ...
        a.n_records+b.n_records, a.n_manual_r+b.n_manual_r, ...
        a.n_auto_r+b.n_auto_r, a.n_matched_r+b.n_matched_r, ...
        a.n_matched_t+b.n_matched_t, sm, sa);

    fprintf('\n  --- Pooled subject-level (N=%d) ---\n', p.subj_n);
    fprintf('  MAE=%.4f, Bias=%+.4f, LoA=[%.4f, %.4f], r=%.3f\n', ...
        p.subj_mae, p.subj_bias, p.subj_loa_lo, p.subj_loa_hi, p.subj_corr);
end


%% ========================================================================
%  PRINT FUNCTIONS
%  ========================================================================
function print_res(r)
    fprintf('\n  --- %s Beat-level ---\n', r.database);
    fprintf('  R:   N=%d, MAE=%.1f ms, Bias=%+.1f ms, LoA=[%.1f, %.1f]\n', ...
        r.r_n, r.r_mae, r.r_bias, r.r_loa_lo, r.r_loa_hi);
    fprintf('  T:   N=%d, MAE=%.1f ms, Bias=%+.1f ms, LoA=[%.1f, %.1f]\n', ...
        r.t_n, r.t_mae, r.t_bias, r.t_loa_lo, r.t_loa_hi);
    fprintf('  CDC: N=%d, MAE=%.4f, Bias=%+.4f, r=%.3f\n', ...
        r.cdc_n, r.cdc_mae, r.cdc_bias, r.cdc_corr);
    fprintf('  Sensitivity: %.1f%%\n', r.r_sensitivity);
    fprintf('\n  --- %s Subject-level (N=%d) ---\n', r.database, r.subj_n);
    fprintf('  MAE=%.4f, Bias=%+.4f, LoA=[%.4f, %.4f], r=%.3f, MedAE=%.4f\n', ...
        r.subj_mae, r.subj_bias, r.subj_loa_lo, r.subj_loa_hi, r.subj_corr, r.subj_median_ae);
end

function print_beat_level_table(a, b, p)
    hdr = '%-10s | %6s | %8s | %9s | %22s | %6s';
    fmt_ms  = '  %-8s | %6d | %6.1f ms | %+7.1f ms | [%+7.1f, %+7.1f] ms |      -';
    fmt_cdc = '  %-8s | %6d | %8.4f | %+9.4f | [%+8.4f, %+8.4f] | %6.3f';
    sep = repmat('-', 1, 72);

    fprintf([hdr '\n'], 'Metric', 'N', 'MAE', 'Bias', '95% LoA', 'r');
    fprintf('%s\nR-PEAK\n', sep);
    fprintf([fmt_ms '\n'], 'LUDB', a.r_n, a.r_mae, a.r_bias, a.r_loa_lo, a.r_loa_hi);
    fprintf([fmt_ms '\n'], 'QTDB', b.r_n, b.r_mae, b.r_bias, b.r_loa_lo, b.r_loa_hi);
    fprintf([fmt_ms '\n'], 'Pooled', p.r_n, p.r_mae, p.r_bias, p.r_loa_lo, p.r_loa_hi);
    fprintf('%s\nT-END\n', sep);
    fprintf([fmt_ms '\n'], 'LUDB', a.t_n, a.t_mae, a.t_bias, a.t_loa_lo, a.t_loa_hi);
    fprintf([fmt_ms '\n'], 'QTDB', b.t_n, b.t_mae, b.t_bias, b.t_loa_lo, b.t_loa_hi);
    fprintf([fmt_ms '\n'], 'Pooled', p.t_n, p.t_mae, p.t_bias, p.t_loa_lo, p.t_loa_hi);
    fprintf('%s\nCDC\n', sep);
    fprintf([fmt_cdc '\n'], 'LUDB', a.cdc_n, a.cdc_mae, a.cdc_bias, a.cdc_loa_lo, a.cdc_loa_hi, a.cdc_corr);
    fprintf([fmt_cdc '\n'], 'QTDB', b.cdc_n, b.cdc_mae, b.cdc_bias, b.cdc_loa_lo, b.cdc_loa_hi, b.cdc_corr);
    fprintf([fmt_cdc '\n'], 'Pooled', p.cdc_n, p.cdc_mae, p.cdc_bias, p.cdc_loa_lo, p.cdc_loa_hi, p.cdc_corr);
    fprintf('%s\n', sep);
    fprintf('Sensitivity: LUDB=%.1f%%, QTDB=%.1f%%, Pooled=%.1f%%\n', ...
        a.r_sensitivity, b.r_sensitivity, p.r_sensitivity);
end

function print_subject_level_table(a, b, p)
    fprintf('%-10s | %5s | %8s | %9s | %22s | %6s | %8s\n', ...
        'Database', 'N', 'MAE', 'Bias', '95% LoA', 'r', 'MedAE');
    sep = repmat('-', 1, 78);
    fprintf('%s\n', sep);
    fmt = '  %-8s | %5d | %8.4f | %+9.4f | [%+8.4f, %+8.4f] | %6.3f | %8.4f';
    fprintf([fmt '\n'], 'LUDB', a.subj_n, a.subj_mae, a.subj_bias, a.subj_loa_lo, a.subj_loa_hi, a.subj_corr, a.subj_median_ae);
    fprintf([fmt '\n'], 'QTDB', b.subj_n, b.subj_mae, b.subj_bias, b.subj_loa_lo, b.subj_loa_hi, b.subj_corr, b.subj_median_ae);
    fprintf([fmt '\n'], 'Pooled', p.subj_n, p.subj_mae, p.subj_bias, p.subj_loa_lo, p.subj_loa_hi, p.subj_corr, p.subj_median_ae);
    fprintf('%s\n', sep);
    fprintf('\nContext (between-group CDC differences):\n');
    fprintf('  HC vs CN:   ~0.040   (subject MAE should be << this)\n');
    fprintf('  HC vs Path: ~0.088   (subject MAE should be << this)\n');
end


%% ========================================================================
%  SAVE CSV
%  ========================================================================
function save_results_csv(a, b, p, paths)
    f = fopen(fullfile(paths.results, 'validation_results.csv'), 'w');
    fprintf(f, 'level,database,metric,n,mae,median_ae,bias,loa_lo,loa_hi,corr,sensitivity\n');
    wr_beat(f, a); wr_beat(f, b); wr_beat(f, p);
    wr_subj(f, 'LUDB', a); wr_subj(f, 'QTDB', b); wr_subj(f, 'Pooled', p);
    fclose(f);
end

function wr_beat(f, r)
    fprintf(f, 'beat,%s,R-peak,%d,%.2f,,%.2f,%.2f,%.2f,,%.1f\n', r.database, r.r_n, r.r_mae, r.r_bias, r.r_loa_lo, r.r_loa_hi, r.r_sensitivity);
    fprintf(f, 'beat,%s,T-end,%d,%.2f,,%.2f,%.2f,%.2f,,\n', r.database, r.t_n, r.t_mae, r.t_bias, r.t_loa_lo, r.t_loa_hi);
    fprintf(f, 'beat,%s,CDC,%d,%.5f,,%.5f,%.5f,%.5f,%.4f,\n', r.database, r.cdc_n, r.cdc_mae, r.cdc_bias, r.cdc_loa_lo, r.cdc_loa_hi, r.cdc_corr);
end

function wr_subj(f, lab, r)
    fprintf(f, 'subject,%s,CDC,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.4f,\n', lab, r.subj_n, r.subj_mae, r.subj_median_ae, r.subj_bias, r.subj_loa_lo, r.subj_loa_hi, r.subj_corr);
end

%% ========================================================================
%  HELPERS
%  ========================================================================
function idx = find_lead(names, target)
    idx = [];
    for i = 1:length(names)
        if strcmpi(strtrim(names{i}), target), idx = i; return; end
    end
    for i = 1:length(names)
        if strcmpi(strtrim(names{i}), 'II'), idx = i; return; end
    end
    if isempty(idx) && ~isempty(names), idx = 1; end
end

function r = qtdb_records()
    r = {'sel100','sel102','sel103','sel104','sel114','sel116','sel117','sel123', ...
         'sel213','sel221','sel223','sel230','sel231','sel232','sel233', ...
         'sel301','sel302','sel306','sel307','sel308','sel310', ...
         'sel803','sel808','sel811','sel820','sel821','sel840','sel847', ...
         'sel853','sel871','sel872','sel873','sel883','sel891', ...
         'sel16265','sel16272','sel16273','sel16420','sel16483', ...
         'sel16539','sel16773','sel16786','sel16795','sel17453', ...
         'sele0104','sele0106','sele0107','sele0110','sele0111','sele0112', ...
         'sele0114','sele0116','sele0121','sele0122','sele0124','sele0126', ...
         'sele0129','sele0133','sele0136','sele0166','sele0170','sele0203', ...
         'sele0210','sele0211','sele0303','sele0405','sele0406','sele0409', ...
         'sele0411','sele0509','sele0603','sele0604','sele0606','sele0607', ...
         'sele0609','sele0612','sele0704', ...
         'sel30','sel31','sel32','sel33','sel34','sel35','sel36','sel37', ...
         'sel38','sel39','sel40','sel41','sel42','sel43','sel44','sel45', ...
         'sel46','sel47','sel48','sel49','sel50','sel51','sel52','sel17152', ...
         'sel14046','sel14157','sel14172','sel15814'};
end