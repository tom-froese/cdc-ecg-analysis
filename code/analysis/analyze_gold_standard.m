function results = analyze_gold_standard()
% ANALYZE_GOLD_STANDARD - CDC analysis using fully manually annotated databases
%
% This analysis uses ONLY datasets where both R-peaks and T-wave endpoints
% were placed by expert cardiologists (annotation_method = 'manual'):
%
%   1. LUDB  - Healthy controls vs Pathological (n ~ 200)
%   2. QTDB  - Clinically normal vs Pathological vs Sudden death (n ~ 100)
%
% LUDB "Healthy" are genuine healthy volunteers (Kalyakulina et al. 2020,
% IEEE Access: "healthy volunteers and patients of the Nizhny Novgorod
% City Hospital No 5").
%
% QTDB three-group classification:
%   Clinically normal - MIT-BIH Normal Sinus Rhythm records (group='healthy')
%                       Hospital referrals to the Beth Israel Arrhythmia
%                       Laboratory with no significant arrhythmias — NOT
%                       healthy volunteers, hence "Clinically Normal"
%   Pathological      - MIT-BIH Arrhythmia, SVA, ST Change, European ST-T,
%                       Long-Term records (group='pathological')
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

    %% Analysis 1: LUDB - Healthy controls vs Pathological
    fprintf('================================================================\n');
    fprintf('LUDB: HEALTHY CONTROLS vs PATHOLOGICAL\n');
    fprintf('================================================================\n');

    ludb_results = analyze_healthy_vs_pathological(ludb_beats, 'LUDB', ...
                   'Healthy controls', n_bootstrap, inv_e);

    %% Analysis 2: QTDB - Clinically normal vs Pathological vs Sudden death
    fprintf('\n================================================================\n');
    fprintf('QTDB: CLINICALLY NORMAL vs PATHOLOGICAL vs SUDDEN DEATH\n');
    fprintf('================================================================\n');

    qtdb_results = analyze_qtdb_three_groups(qtdb_beats, n_bootstrap, inv_e);

    %% Integrated summary
    fprintf('\n================================================================\n');
    fprintf('GOLD-STANDARD SUMMARY\n');
    fprintf('================================================================\n\n');

    fprintf('                              Mode      95%% CI          n      %% > 1/e\n');
    fprintf('------------------------------------------------------------------------\n');
    fprintf('LUDB Healthy controls         %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            ludb_results.mode_healthy, ludb_results.ci_healthy(1), ...
            ludb_results.ci_healthy(2), ludb_results.n_healthy, ...
            ludb_results.pct_healthy_above);
    fprintf('LUDB Pathological             %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            ludb_results.mode_patient, ludb_results.ci_patient(1), ...
            ludb_results.ci_patient(2), ludb_results.n_patient, ...
            ludb_results.pct_patient_above);
    fprintf('------------------------------------------------------------------------\n');
    fprintf('QTDB Clinically normal        %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            qtdb_results.mode_normal, qtdb_results.ci_normal(1), ...
            qtdb_results.ci_normal(2), qtdb_results.n_normal, ...
            qtdb_results.pct_normal_above);
    fprintf('QTDB Pathological             %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            qtdb_results.mode_pathological, qtdb_results.ci_pathological(1), ...
            qtdb_results.ci_pathological(2), qtdb_results.n_pathological, ...
            qtdb_results.pct_pathological_above);
    fprintf('QTDB Sudden death             %.4f    [%.3f, %.3f]    %3d    %.1f%%\n', ...
            qtdb_results.mode_fatal, qtdb_results.ci_fatal(1), ...
            qtdb_results.ci_fatal(2), qtdb_results.n_fatal, ...
            qtdb_results.pct_fatal_above);
    fprintf('------------------------------------------------------------------------\n');
    fprintf('Reference: 1/e = %.4f\n\n', inv_e);

    fprintf('Statistical comparisons:\n');
    fprintf('  LUDB: p = %.2e (healthy vs pathological)\n', ludb_results.p_value);
    fprintf('  QTDB: p = %.2e (clinically normal vs pathological)\n', qtdb_results.p_normal_vs_patho);
    fprintf('  QTDB: p = %.2e (pathological vs sudden death)\n', qtdb_results.p_patho_vs_fatal);
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
    patient_idx = ~healthy_idx;

    healthy_ratios = median_ratio(healthy_idx);
    patient_ratios = median_ratio(patient_idx);

    n_healthy = length(healthy_ratios);
    n_patient = length(patient_ratios);

    fprintf('  %s: n = %d\n', healthy_label, n_healthy);
    fprintf('  Pathological: n = %d\n', n_patient);

    %% Statistics
    [mode_healthy, ci_healthy] = bootstrap_mode(healthy_ratios, n_bootstrap);
    pct_healthy_above = 100 * sum(healthy_ratios > inv_e) / n_healthy;

    fprintf('\n  %s:\n', upper(healthy_label));
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_healthy, ci_healthy(1), ci_healthy(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(healthy_ratios), median(healthy_ratios));
    fprintf('    1/e in CI: %s\n', yesno(inv_e >= ci_healthy(1) && inv_e <= ci_healthy(2)));

    [mode_patient, ci_patient] = bootstrap_mode(patient_ratios, n_bootstrap);
    pct_patient_above = 100 * sum(patient_ratios > inv_e) / n_patient;

    fprintf('\n  PATHOLOGICAL:\n');
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_patient, ci_patient(1), ci_patient(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(patient_ratios), median(patient_ratios));
    fprintf('    1/e in CI: %s\n', yesno(inv_e >= ci_patient(1) && inv_e <= ci_patient(2)));

    %% Comparison
    if n_healthy >= 3 && n_patient >= 3
        [p_wilcoxon, ~] = ranksum(healthy_ratios, patient_ratios);
        fprintf('\n  Wilcoxon: p = %.2e%s\n', p_wilcoxon, stars(p_wilcoxon));
    else
        p_wilcoxon = NaN;
    end

    %% Package
    results.n_healthy = n_healthy;
    results.n_patient = n_patient;
    results.mode_healthy = mode_healthy;
    results.ci_healthy = ci_healthy;
    results.mode_patient = mode_patient;
    results.ci_patient = ci_patient;
    results.pct_healthy_above = pct_healthy_above;
    results.pct_patient_above = pct_patient_above;
    results.p_value = p_wilcoxon;
    results.healthy_ratios = healthy_ratios;
    results.patient_ratios = patient_ratios;
end

%% ========================================================================
%  QTDB THREE-GROUP ANALYSIS
%  ========================================================================
function results = analyze_qtdb_three_groups(beats, n_bootstrap, inv_e)
% Three-group split:
%   'healthy'       → Clinically normal (MIT-BIH Normal Sinus Rhythm, n~10)
%                     Hospital referrals with no arrhythmias, NOT healthy volunteers
%   'pathological'  → Pathological (Arrhythmia, SVA, ST Change, European ST-T, Long-Term)
%   'sudden_death'  → Sudden death (BIH Sudden Death)

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

    normal_ratios = median_ratio(normal_idx);
    patho_ratios  = median_ratio(patho_idx);
    fatal_ratios  = median_ratio(fatal_idx);

    n_normal = length(normal_ratios);
    n_pathological = length(patho_ratios);
    n_fatal = length(fatal_ratios);

    fprintf('  Clinically normal: n = %d (MIT-BIH Normal Sinus Rhythm)\n', n_normal);
    fprintf('  Pathological:      n = %d (Arrhythmia, SVA, ST, European, Long-Term)\n', n_pathological);
    fprintf('  Sudden death:      n = %d (BIH Sudden Death)\n', n_fatal);

    %% Statistics — Clinically normal
    if n_normal >= 3
        [mode_normal, ci_normal] = bootstrap_mode(normal_ratios, n_bootstrap);
    else
        mode_normal = median(normal_ratios); ci_normal = [NaN NaN];
    end
    pct_normal_above = 100 * sum(normal_ratios > inv_e) / n_normal;

    fprintf('\n  CLINICALLY NORMAL:\n');
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_normal, ci_normal(1), ci_normal(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(normal_ratios), median(normal_ratios));
    fprintf('    %% > 1/e: %.1f%%\n', pct_normal_above);

    %% Statistics — Pathological
    [mode_pathological, ci_pathological] = bootstrap_mode(patho_ratios, n_bootstrap);
    pct_pathological_above = 100 * sum(patho_ratios > inv_e) / n_pathological;

    fprintf('\n  PATHOLOGICAL:\n');
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_pathological, ci_pathological(1), ci_pathological(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(patho_ratios), median(patho_ratios));
    fprintf('    %% > 1/e: %.1f%%\n', pct_pathological_above);

    %% Statistics — Sudden death
    [mode_fatal, ci_fatal] = bootstrap_mode(fatal_ratios, n_bootstrap);
    pct_fatal_above = 100 * sum(fatal_ratios > inv_e) / n_fatal;

    fprintf('\n  SUDDEN DEATH:\n');
    fprintf('    Mode: %.4f [%.4f, %.4f]\n', mode_fatal, ci_fatal(1), ci_fatal(2));
    fprintf('    Mean: %.4f, Median: %.4f\n', mean(fatal_ratios), median(fatal_ratios));
    fprintf('    %% > 1/e: %.1f%%\n', pct_fatal_above);

    %% Pairwise comparisons
    if n_normal >= 3 && n_pathological >= 3
        p_normal_vs_patho = ranksum(normal_ratios, patho_ratios);
        fprintf('\n  Clinically normal vs Pathological: p = %.2e%s\n', ...
                p_normal_vs_patho, stars(p_normal_vs_patho));
    else
        p_normal_vs_patho = NaN;
        fprintf('\n  Clinically normal vs Pathological: insufficient n\n');
    end

    p_patho_vs_fatal = ranksum(patho_ratios, fatal_ratios);
    fprintf('  Pathological vs Sudden death:      p = %.2e%s\n', ...
            p_patho_vs_fatal, stars(p_patho_vs_fatal));

    %% Kruskal-Wallis (overall three-group test)
    all_ratios = [normal_ratios; patho_ratios; fatal_ratios];
    group_labels = [repmat({'ClinicallyNormal'}, n_normal, 1); ...
                    repmat({'Pathological'}, n_pathological, 1); ...
                    repmat({'SuddenDeath'}, n_fatal, 1)];
    p_kruskal = kruskalwallis(all_ratios, group_labels, 'off');
    fprintf('  Kruskal-Wallis (3 groups):         p = %.2e%s\n', ...
            p_kruskal, stars(p_kruskal));

    %% Package
    results.n_normal = n_normal;
    results.n_pathological = n_pathological;
    results.n_fatal = n_fatal;
    results.mode_normal = mode_normal;
    results.ci_normal = ci_normal;
    results.mode_pathological = mode_pathological;
    results.ci_pathological = ci_pathological;
    results.mode_fatal = mode_fatal;
    results.ci_fatal = ci_fatal;
    results.pct_normal_above = pct_normal_above;
    results.pct_pathological_above = pct_pathological_above;
    results.pct_fatal_above = pct_fatal_above;
    results.p_normal_vs_patho = p_normal_vs_patho;
    results.p_patho_vs_fatal = p_patho_vs_fatal;
    results.p_kruskal = p_kruskal;
    results.normal_ratios = normal_ratios;
    results.pathological_ratios = patho_ratios;
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
