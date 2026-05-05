function results = analyze_gold_standard()
% ANALYZE_GOLD_STANDARD - CDC analysis using fully manually annotated databases
%
% This analysis uses ONLY datasets where both R-peaks and T-wave endpoints
% were placed by expert cardiologists (annotation_method = 'manual'):
%
%   1. LUDB  - Healthy controls vs Patients (pathological ECG) (n ~ 200)
%   2. QTDB  - Patients (non-pathological ECG) vs Patients (pathological ECG)
%             vs Sudden death (n ~ 100)
%
% LUDB "Healthy" are genuine healthy volunteers (Kalyakulina et al. 2020,
% IEEE Access: "healthy volunteers and patients of the Nizhny Novgorod
% City Hospital No 5").
%
% QTDB three-group classification:
%   Patients (non-pathological ECG) — internal categorical 'NonPathECG'
%       - MIT-BIH Normal Sinus Rhythm records (source group='healthy').
%         Hospital referrals to the Beth Israel Arrhythmia Laboratory
%         with no significant arrhythmias — NOT healthy volunteers.
%   Patients (pathological ECG) — internal categorical 'PathECG'
%       - MIT-BIH Arrhythmia, SVA, ST Change, European ST-T, Long-Term
%         records (source group='pathological').
%   Sudden death      - BIH Sudden Death records (group='sudden_death')
%
% Purpose: Establish the core CDC pattern using gold-standard fiducial
% points that cannot be attributed to algorithmic bias. These datasets
% are small but the annotation quality is beyond reproach.
%
% No beat-level quality filters are applied to manually annotated data,
% since expert annotation implicitly excludes artefactual beats.
%
% Subject aggregation uses unique_subject_id (format: Database_RecordID)
% for consistency with the hierarchical model pipeline.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    n_bootstrap = 5000;
    inv_e = 1 / exp(1);

    fprintf('================================================================\n');
    fprintf('GOLD-STANDARD CDC ANALYSIS (Manual Annotation Only)\n');
    fprintf('================================================================\n');
    fprintf('Reference value: 1/e = %.4f\n', inv_e);
    fprintf('Annotation tier: MANUAL (expert-placed fiducial points)\n');
    fprintf('================================================================\n\n');

    %% Load datasets
    ludb_beats = readtable(paths.csv_ludb);
    qtdb_beats = readtable(paths.csv_qtdb);

    fprintf('LUDB: %d beats loaded\n', height(ludb_beats));
    fprintf('QTDB: %d beats loaded\n\n', height(qtdb_beats));

    %% Verify annotation method
    assert(all(strcmp(get_col_as_char(ludb_beats, 'annotation_method'), 'manual')), ...
           'LUDB contains non-manual annotations');
    assert(all(strcmp(get_col_as_char(qtdb_beats, 'annotation_method'), 'manual')), ...
           'QTDB contains non-manual annotations');

    %% Analysis 1: LUDB - Healthy controls vs Patients (pathological ECG)
    fprintf('================================================================\n');
    fprintf('LUDB: HEALTHY CONTROLS vs PATHOLOGICAL\n');
    fprintf('================================================================\n');

    ludb_results = analyze_healthy_vs_pathological(ludb_beats, 'LUDB', ...
                   'Healthy controls', n_bootstrap, inv_e);

    %% Analysis 2: QTDB - Patients (non-pathological ECG) vs Patients (pathological ECG) vs Sudden death
    fprintf('\n================================================================\n');
    fprintf('QTDB: PATIENTS (NON-PATHOLOGICAL ECG) vs PATIENTS (PATHOLOGICAL ECG) vs SUDDEN DEATH\n');
    fprintf('================================================================\n');

    qtdb_results = analyze_qtdb_three_groups(qtdb_beats, n_bootstrap, inv_e);

    %% Integrated summary
    fprintf('\n================================================================\n');
    fprintf('GOLD-STANDARD SUMMARY\n');
    fprintf('================================================================\n\n');

    fprintf('                                      Mode      95%% CI          n      %% > 1/e\n');
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('LUDB Healthy controls                 %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            ludb_results.mode_healthy, ludb_results.ci_healthy(1), ...
            ludb_results.ci_healthy(2), ludb_results.n_healthy, ...
            ludb_results.pct_healthy_above);
    fprintf('LUDB Patients (pathological ECG)      %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            ludb_results.mode_pe, ludb_results.ci_pe(1), ...
            ludb_results.ci_pe(2), ludb_results.n_pe, ...
            ludb_results.pct_pe_above);
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('QTDB Patients (non-pathological ECG)  %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            qtdb_results.mode_npe, qtdb_results.ci_npe(1), ...
            qtdb_results.ci_npe(2), qtdb_results.n_npe, ...
            qtdb_results.pct_npe_above);
    fprintf('QTDB Patients (pathological ECG)      %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            qtdb_results.mode_pe, qtdb_results.ci_pe(1), ...
            qtdb_results.ci_pe(2), qtdb_results.n_pe, ...
            qtdb_results.pct_pe_above);
    fprintf('QTDB Sudden death                     %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            qtdb_results.mode_fatal, qtdb_results.ci_fatal(1), ...
            qtdb_results.ci_fatal(2), qtdb_results.n_fatal, ...
            qtdb_results.pct_fatal_above);
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('Reference: 1/e = %.4f\n\n', inv_e);

    fprintf('Statistical comparisons:\n');
    fprintf('  LUDB: p = %.2e (healthy vs pathological)\n', ludb_results.p_value);
    fprintf('  QTDB: p = %.2e (Patients (non-pathological ECG) vs Patients (pathological ECG))\n', qtdb_results.p_npe_vs_pe);
    fprintf('  QTDB: p = %.2e (Patients (pathological ECG) vs sudden death)\n', qtdb_results.p_pe_vs_fatal);
    fprintf('  QTDB: p = %.2e (Kruskal-Wallis, all three groups)\n', qtdb_results.p_kruskal);

    %% Package results
    results.ludb = ludb_results;
    results.qtdb = qtdb_results;
    results.inv_e = inv_e;

    save(fullfile(paths.results, 'gold_standard_results.mat'), 'results');
    fprintf('\nResults saved to: %s\n', fullfile(paths.results, 'gold_standard_results.mat'));
end

%% ========================================================================
%  HEALTHY vs PATHOLOGICAL (used for LUDB)
%  ========================================================================
function results = analyze_healthy_vs_pathological(beats, dataset_name, healthy_label, n_bootstrap, inv_e)

    fprintf('Processing %s...\n', dataset_name);

    %% Compute subject-level medians using unique_subject_id
    [unique_subjects, ~, subject_idx] = unique(beats.unique_subject_id);
    n_subjects = length(unique_subjects);

    median_ratio = zeros(n_subjects, 1);
    n_beats = zeros(n_subjects, 1);
    group = cell(n_subjects, 1);

    for i = 1:n_subjects
        mask = subject_idx == i;
        rec = beats(mask, :);

        rt_ms = (rec.t_end_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        rr_ms = (rec.next_r_sample - rec.r_sample) ./ rec.fs(1) * 1000;

        valid = rr_ms > 0;
        if sum(valid) >= 1
            ratios = rt_ms(valid) ./ rr_ms(valid);
            median_ratio(i) = median(ratios);
            n_beats(i) = sum(valid);
        else
            median_ratio(i) = NaN;
            n_beats(i) = 0;
        end

        group{i} = get_string_field(rec, 'group');
    end

    valid = ~isnan(median_ratio) & n_beats > 0;
    median_ratio = median_ratio(valid);
    group = group(valid);

    %% Split
    healthy_idx = strcmp(group, 'healthy');
    pe_idx = ~healthy_idx;

    healthy_ratios = median_ratio(healthy_idx);
    pe_ratios = median_ratio(pe_idx);

    n_healthy = length(healthy_ratios);
    n_pe = length(pe_ratios);

    fprintf('  %s: n = %d\n', healthy_label, n_healthy);
    fprintf('  Patients (pathological ECG): n = %d\n', n_pe);

    %% Statistics
    [mode_healthy, ci_healthy] = bootstrap_mode(healthy_ratios, n_bootstrap);
    pct_healthy_above = 100 * sum(healthy_ratios > inv_e) / n_healthy;

    fprintf('\n  %s:\n', upper(healthy_label));
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_healthy, ci_healthy(1), ci_healthy(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(healthy_ratios), median(healthy_ratios));
    fprintf('    1/e in CI: %s\n', yesno(inv_e >= ci_healthy(1) && inv_e <= ci_healthy(2)));

    [mode_pe, ci_pe] = bootstrap_mode(pe_ratios, n_bootstrap);
    pct_pe_above = 100 * sum(pe_ratios > inv_e) / n_pe;

    fprintf('\n  PATIENTS (PATHOLOGICAL ECG):\n');
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_pe, ci_pe(1), ci_pe(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(pe_ratios), median(pe_ratios));
    fprintf('    1/e in CI: %s\n', yesno(inv_e >= ci_pe(1) && inv_e <= ci_pe(2)));

    %% Comparison
    if n_healthy >= 3 && n_pe >= 3
        [p_wilcoxon, ~] = ranksum(healthy_ratios, pe_ratios);
        fprintf('\n  Wilcoxon: p = %.2e%s\n', p_wilcoxon, stars(p_wilcoxon));
    else
        p_wilcoxon = NaN;
    end

    %% Package
    results.n_healthy = n_healthy;
    results.n_pe = n_pe;
    results.mode_healthy = mode_healthy;
    results.ci_healthy = ci_healthy;
    results.mode_pe = mode_pe;
    results.ci_pe = ci_pe;
    results.pct_healthy_above = pct_healthy_above;
    results.pct_pe_above = pct_pe_above;
    results.p_value = p_wilcoxon;
    results.healthy_ratios = healthy_ratios;
    results.pe_ratios = pe_ratios;
end

%% ========================================================================
%  QTDB THREE-GROUP ANALYSIS
%  ========================================================================
function results = analyze_qtdb_three_groups(beats, n_bootstrap, inv_e)
% Three-group split (source CSV value → internal categorical → display label):
%   'healthy'      → 'NonPathECG'  → "Patients (non-pathological ECG)"
%                    MIT-BIH Normal Sinus Rhythm (n~10), hospital
%                    referrals with no arrhythmias, NOT healthy volunteers.
%   'pathological' → 'PathECG'     → "Patients (pathological ECG)"
%                    Arrhythmia, SVA, ST Change, European ST-T,
%                    Long-Term records.
%   'sudden_death' → 'SuddenDeath' → "Sudden death"
%                    BIH Sudden Death records.

    fprintf('Processing QTDB (three-group split)...\n');

    [unique_subjects, ~, subject_idx] = unique(beats.unique_subject_id);
    n_subjects = length(unique_subjects);

    median_ratio = zeros(n_subjects, 1);
    n_beats = zeros(n_subjects, 1);
    group = cell(n_subjects, 1);

    for i = 1:n_subjects
        mask = subject_idx == i;
        rec = beats(mask, :);

        rt_ms = (rec.t_end_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        rr_ms = (rec.next_r_sample - rec.r_sample) ./ rec.fs(1) * 1000;

        valid = rr_ms > 0;
        if sum(valid) >= 1
            ratios = rt_ms(valid) ./ rr_ms(valid);
            median_ratio(i) = median(ratios);
            n_beats(i) = sum(valid);
        else
            median_ratio(i) = NaN;
            n_beats(i) = 0;
        end
        group{i} = get_string_field(rec, 'group');
    end

    valid = ~isnan(median_ratio) & n_beats > 0;
    median_ratio = median_ratio(valid);
    group = group(valid);

    %% Three-way split
    normal_idx = strcmp(group, 'healthy');
    fatal_idx  = strcmp(group, 'sudden_death');
    patho_idx  = ~normal_idx & ~fatal_idx;

    npe_ratios = median_ratio(normal_idx);
    pe_ratios  = median_ratio(patho_idx);
    fatal_ratios  = median_ratio(fatal_idx);

    n_npe = length(npe_ratios);
    n_pe = length(pe_ratios);
    n_fatal = length(fatal_ratios);

    fprintf('  Patients (non-pathological ECG): n = %d (MIT-BIH Normal Sinus Rhythm)\n', n_npe);
    fprintf('  Patients (pathological ECG):     n = %d (Arrhythmia, SVA, ST, European, Long-Term)\n', n_pe);
    fprintf('  Sudden death:                    n = %d (BIH Sudden Death)\n', n_fatal);

    %% Statistics — Patients (non-pathological ECG)
    if n_npe >= 3
        [mode_npe, ci_npe] = bootstrap_mode(npe_ratios, n_bootstrap);
    else
        mode_npe = median(npe_ratios); ci_npe = [NaN NaN];
    end
    pct_npe_above = 100 * sum(npe_ratios > inv_e) / n_npe;

    fprintf('\n  PATIENTS (NON-PATHOLOGICAL ECG):\n');
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_npe, ci_npe(1), ci_npe(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(npe_ratios), median(npe_ratios));
    fprintf('    %% > 1/e: %.1f%%\n', pct_npe_above);

    %% Statistics — Patients (pathological ECG)
    [mode_pe, ci_pe] = bootstrap_mode(pe_ratios, n_bootstrap);
    pct_pe_above = 100 * sum(pe_ratios > inv_e) / n_pe;

    fprintf('\n  PATIENTS (PATHOLOGICAL ECG):\n');
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_pe, ci_pe(1), ci_pe(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(pe_ratios), median(pe_ratios));
    fprintf('    %% > 1/e: %.1f%%\n', pct_pe_above);

    %% Statistics — Sudden death
    [mode_fatal, ci_fatal] = bootstrap_mode(fatal_ratios, n_bootstrap);
    pct_fatal_above = 100 * sum(fatal_ratios > inv_e) / n_fatal;

    fprintf('\n  SUDDEN DEATH:\n');
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_fatal, ci_fatal(1), ci_fatal(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(fatal_ratios), median(fatal_ratios));
    fprintf('    %% > 1/e: %.1f%%\n', pct_fatal_above);

    %% Pairwise comparisons
    if n_npe >= 3 && n_pe >= 3
        p_npe_vs_pe = ranksum(npe_ratios, pe_ratios);
        fprintf('\n  Patients (non-pathological ECG) vs Patients (pathological ECG): p = %.2e%s\n', ...
                p_npe_vs_pe, stars(p_npe_vs_pe));
    else
        p_npe_vs_pe = NaN;
        fprintf('\n  Patients (non-pathological ECG) vs Patients (pathological ECG): insufficient n\n');
    end

    p_pe_vs_fatal = ranksum(pe_ratios, fatal_ratios);
    fprintf('  Patients (pathological ECG) vs Sudden death:                    p = %.2e%s\n', ...
            p_pe_vs_fatal, stars(p_pe_vs_fatal));

    %% Kruskal-Wallis (overall three-group test)
    all_ratios = [npe_ratios; pe_ratios; fatal_ratios];
    group_labels = [repmat({'NonPathECG'}, n_npe, 1); ...
                    repmat({'PathECG'}, n_pe, 1); ...
                    repmat({'SuddenDeath'}, n_fatal, 1)];
    p_kruskal = kruskalwallis(all_ratios, group_labels, 'off');
    fprintf('  Kruskal-Wallis (3 groups):         p = %.2e%s\n', ...
            p_kruskal, stars(p_kruskal));

    %% Package
    results.n_npe = n_npe;
    results.n_pe = n_pe;
    results.n_fatal = n_fatal;
    results.mode_npe = mode_npe;
    results.ci_npe = ci_npe;
    results.mode_pe = mode_pe;
    results.ci_pe = ci_pe;
    results.mode_fatal = mode_fatal;
    results.ci_fatal = ci_fatal;
    results.pct_npe_above = pct_npe_above;
    results.pct_pe_above = pct_pe_above;
    results.pct_fatal_above = pct_fatal_above;
    results.p_npe_vs_pe = p_npe_vs_pe;
    results.p_pe_vs_fatal = p_pe_vs_fatal;
    results.p_kruskal = p_kruskal;
    results.npe_ratios = npe_ratios;
    results.pe_ratios = pe_ratios;
    results.fatal_ratios = fatal_ratios;
end

%% ========================================================================
%  HELPERS
%  ========================================================================

function s = get_string_field(rec, field_name)
% GET_STRING_FIELD - Robustly extract a char value from the first row of a
% table column that may be stored as cell, string, or categorical.
% Matches the pattern used in analyze_hierarchical_model.m.
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

function col = get_col_as_char(tbl, field_name)
% GET_COL_AS_CHAR - Return an entire table column as a cell array of char,
% regardless of whether it was read as cell, string, or categorical.
    raw = tbl.(field_name);
    if iscell(raw)
        col = raw;
    elseif isstring(raw)
        col = cellstr(raw);
    elseif iscategorical(raw)
        col = cellstr(raw);
    else
        col = cellstr(string(raw));
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
