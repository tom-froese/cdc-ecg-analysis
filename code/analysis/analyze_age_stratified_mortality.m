function results = analyze_age_stratified_mortality()
% ANALYZE_AGE_STRATIFIED_MORTALITY - Proximity to 1/e predicts survival
%   at every age decade
%
% Tests whether patients whose CDC is closer to the thermodynamic optimum
% (1/e) have lower all-cause mortality within each age decade, using the
% full CODE-15% cohort with pathological status as a covariate.
%
% Primary analysis: within-stratum tertiles of |CDC - 1/e|
%   For each age decade, patients are split into tertiles of deviation
%   from 1/e. This controls for age-related drift and asks: among
%   people of the same age, do those closer to 1/e fare better?
%
% Robustness: fixed cutoffs (|CDC - 1/e| < 0.03 / 0.03-0.06 / > 0.06)
%   Confirms results are not an artefact of tertile construction.
%
% Statistical models:
%   - Cochran-Mantel-Haenszel test (T1 vs T3 across strata)
%   - Cox PH with CDC_dev × Age interaction
%   - Cox PH main effects
%   - Within-stratum proportion tests
%
% Subject aggregation uses unique_subject_id for consistency with the
% unified beat table format (cf. analyze_cdc_vs_hr.m).
%
% Dependencies: config.m, apply_quality_filters.m
% Data:         CODE-15% beat extraction + exams.csv mortality labels
% Output:       age_stratified_mortality_results.mat
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    fprintf('================================================================\n');
    fprintf('AGE-STRATIFIED MORTALITY: PROXIMITY TO 1/e AS PROTECTION\n');
    fprintf('================================================================\n');
    fprintf('Reference optimum: 1/e = %.4f\n\n', inv_e);

    %% ================================================================
    %  1. LOAD AND AGGREGATE TO SUBJECT LEVEL
    %  ================================================================

    fprintf('Loading beat-level data...\n');
    beats = readtable(paths.csv_code15, 'TextType', 'string');
    fprintf('  %d beats loaded\n', height(beats));

    fprintf('Applying quality filters (no CDC ratio cutoffs)...\n');
    [valid_mask, ~] = apply_quality_filters(beats, ...
        'CDCMin', 0, 'CDCMax', 1);
    beats = beats(valid_mask, :);
    fprintf('  %d beats after filtering\n', height(beats));

    % --- Subject-level aggregation via unique_subject_id ---
    % unique_subject_id format: CODE15_<patient_id>
    % This is the project-wide convention; see analyze_cdc_vs_hr.m
    fs = beats.fs(1);
    rt_ms = (beats.t_end_sample - beats.r_sample) / fs * 1000;
    rr_ms = (beats.next_r_sample - beats.r_sample) / fs * 1000;
    beat_cdc = rt_ms ./ rr_ms;
    valid_beats = (rr_ms > 0) & (rt_ms > 0) & (rt_ms < rr_ms);

    [G, subject_ids] = findgroups(beats.unique_subject_id);

    CDC = splitapply(@(x, m) median(x(m)), beat_cdc, valid_beats, G);
    HR  = splitapply(@(x, m) 60000 ./ median(x(m)), rr_ms, valid_beats, G);
    Age = splitapply(@(x) x(1), beats.age, G);
    Sex_str = splitapply(@(x) x(1), beats.sex, G);

    % Clinical group from source_subset
    is_normal_beat = (beats.source_subset == "normal");
    total_per  = splitapply(@numel, beats.record_id, G);
    normal_per = splitapply(@sum, is_normal_beat, G);
    is_cn = (total_per == normal_per);

    Group_str = strings(length(subject_ids), 1);
    Group_str(is_cn)  = "ClinicallyNormal";
    Group_str(~is_cn) = "Pathological";

    % Extract numeric patient ID for merge with exams.csv
    % unique_subject_id = "CODE15_<patient_id>"
    if isstring(subject_ids)
        pid_str = extractAfter(subject_ids, "CODE15_");
        merge_id = str2double(pid_str);
    elseif iscell(subject_ids)
        pid_str = cellfun(@(x) extractAfter(string(x), "CODE15_"), ...
            subject_ids, 'UniformOutput', false);
        merge_id = str2double(pid_str);
    else
        merge_id = subject_ids;
    end

    beat_data = table(subject_ids, merge_id, Age, CDC, HR, Sex_str, Group_str, ...
        'VariableNames', {'unique_subject_id', 'merge_id', ...
                          'Age', 'CDC', 'HR', 'Sex', 'Group'});

    fprintf('  %d unique subjects after aggregation\n', height(beat_data));

    %% ================================================================
    %  2. MERGE WITH MORTALITY ANNOTATIONS
    %  ================================================================

    if isfield(paths, 'csv_code15_exams')
        exams_csv = paths.csv_code15_exams;
    else
        exams_csv = fullfile(paths.raw_code15, 'exams.csv');
    end
    fprintf('\nLoading exams metadata: %s\n', exams_csv);

    exams = readtable(exams_csv, 'TextType', 'string', ...
        'TreatAsMissing', '', 'VariableNamingRule', 'preserve');
    fprintf('  %d exam records loaded\n', height(exams));

    % Keep first exam per patient
    exams = sortrows(exams, {'patient_id', 'exam_id'}, {'ascend', 'ascend'});
    [~, first_idx] = unique(exams.patient_id, 'first');
    exams_first = exams(first_idx, :);

    if isnumeric(exams_first.patient_id)
        exams_first.merge_id = exams_first.patient_id;
    else
        exams_first.merge_id = str2double(exams_first.patient_id);
    end

    mortality_vars = intersect(exams_first.Properties.VariableNames, ...
        {'death', 'timey', 'normal_ecg', 'is_male'});

    merged = innerjoin(beat_data, exams_first, ...
        'LeftKeys', 'merge_id', 'RightKeys', 'merge_id', ...
        'RightVariables', mortality_vars);

    fprintf('  Merged records: %d\n', height(merged));

    % Parse mortality outcome
    if isstring(merged.death) || iscell(merged.death)
        merged.Deceased = lower(string(merged.death)) == "true" | ...
                          merged.death == "1";
    elseif islogical(merged.death)
        merged.Deceased = merged.death;
    else
        merged.Deceased = (merged.death == 1);
    end

    has_timey = ismember('timey', merged.Properties.VariableNames);
    if has_timey
        if isstring(merged.timey) || iscell(merged.timey)
            merged.FollowUp_yrs = str2double(merged.timey);
        else
            merged.FollowUp_yrs = merged.timey;
        end
    end

    valid_mort = ~isnan(double(merged.Deceased));
    if has_timey
        valid_mort = valid_mort & ~isnan(merged.FollowUp_yrs) & ...
                     (merged.FollowUp_yrs > 0);
    end
    valid_mort = valid_mort & ~isnan(merged.CDC) & ~isnan(merged.Age) & ...
                 (merged.Age > 0);

    data = merged(valid_mort, :);
    data.Group = categorical(data.Group);
    data.Sex = categorical(data.Sex);
    data.CDC_dev = abs(data.CDC - inv_e);
    data.Pathological = double(data.Group == 'Pathological');

    fprintf('  Valid records: %d (Deaths: %d, %.2f%%)\n', ...
        height(data), sum(data.Deceased), 100 * mean(data.Deceased));

    %% ================================================================
    %  3. DEFINE AGE STRATA AND DEVIATION GROUPS
    %  ================================================================

    age_edges  = [18, 40, 50, 60, 70, 80, Inf];
    age_labels = {'18-39', '40-49', '50-59', '60-69', '70-79', '80+'};
    n_age = length(age_labels);

    data.AgeBin = discretize(data.Age, age_edges);

    % --- (A) Primary: within-stratum tertiles ---
    data.DevTertile = NaN(height(data), 1);
    tertile_bounds = NaN(n_age, 2);

    for ai = 1:n_age
        mask = data.AgeBin == ai;
        dev = data.CDC_dev(mask);
        q33 = quantile(dev, 1/3);
        q67 = quantile(dev, 2/3);
        tertile_bounds(ai, :) = [q33, q67];

        t = ones(sum(mask), 1);
        t(dev >= q33) = 2;
        t(dev >= q67) = 3;
        data.DevTertile(mask) = t;
    end

    tert_labels = {'Near 1/e', 'Moderate', 'Far from 1/e'};

    fprintf('\nTertile boundaries (|CDC - 1/e|):\n');
    fprintf('  %-8s  %10s  %10s\n', 'Age', 'T1/T2', 'T2/T3');
    fprintf('  %s\n', repmat('-', 1, 32));
    for ai = 1:n_age
        fprintf('  %-8s  %10.4f  %10.4f\n', age_labels{ai}, ...
            tertile_bounds(ai, 1), tertile_bounds(ai, 2));
    end

    % --- (B) Robustness: fixed cutoffs ---
    dev_edges = [0, 0.03, 0.06, Inf];
    data.DevFixed = discretize(data.CDC_dev, dev_edges);
    fix_labels = {'Near 1/e (<0.03)', 'Moderate (0.03-0.06)', 'Far from 1/e (>0.06)'};

    %% ================================================================
    %  4. TABULATE MORTALITY
    %  ================================================================

    fprintf('\n================================================================\n');
    fprintf('PRIMARY: WITHIN-STRATUM TERTILES\n');
    fprintf('================================================================\n\n');

    [mort_rate_tert, mort_n_tert, mort_d_tert] = tabulate_mortality( ...
        data, n_age, 3, age_labels, tert_labels, 'AgeBin', 'DevTertile');

    fprintf('\n================================================================\n');
    fprintf('ROBUSTNESS: FIXED CUTOFFS\n');
    fprintf('================================================================\n\n');

    [mort_rate_fix, mort_n_fix, mort_d_fix] = tabulate_mortality( ...
        data, n_age, 3, age_labels, fix_labels, 'AgeBin', 'DevFixed');

    %% ================================================================
    %  5. CMH TEST (T1 vs T3 across age strata)
    %  ================================================================

    fprintf('\n--- Cochran-Mantel-Haenszel test: T1 vs T3 across age strata ---\n');
    cmh_num = 0;
    cmh_denom = 0;
    for ai = 1:n_age
        a = mort_d_tert(ai, 1);   % deaths in T1
        b = mort_n_tert(ai, 1) - a;  % survivors in T1
        c = mort_d_tert(ai, 3);   % deaths in T3
        d = mort_n_tert(ai, 3) - c;  % survivors in T3
        n_str = a + b + c + d;
        if n_str > 0 && mort_n_tert(ai, 1) >= 10 && mort_n_tert(ai, 3) >= 10
            cmh_num   = cmh_num   + (a - (a+b)*(a+c)/n_str);
            cmh_denom = cmh_denom + (a+b)*(c+d)*(a+c)*(b+d) / (n_str^2*(n_str-1));
        end
    end
    cmh_chi2 = cmh_num^2 / cmh_denom;
    cmh_p = 1 - chi2cdf(cmh_chi2, 1);
    fprintf('  CMH chi2 = %.2f, df = 1, p = %.2e\n', cmh_chi2, cmh_p);

    %% ================================================================
    %  6. COX PH MODELS
    %  ================================================================

    if has_timey
        fprintf('\n================================================================\n');
        fprintf('COX PROPORTIONAL HAZARDS MODELS\n');
        fprintf('================================================================\n\n');

        data.is_male = double(data.Sex == 'Male' | data.Sex == 'M');
        sd_dev = std(data.CDC_dev);
        censoring = ~data.Deceased;
        time = data.FollowUp_yrs;

        % --- 6a. Interaction model ---
        fprintf('--- Cox PH: CDC_dev × Age interaction ---\n\n');
        X_cox = [data.CDC_dev, data.Age, data.CDC_dev .* data.Age, ...
                 data.is_male, data.Pathological];
        cox_names = {'CDC_dev', 'Age', 'CDC_dev x Age', 'Male', 'Pathological'};

        [b_cox, ~, ~, stats_cox] = coxphfit(X_cox, time, ...
            'Censoring', censoring);
        print_cox(cox_names, b_cox, stats_cox);
        fprintf('  Per-SD HR for CDC_dev: %.3f\n\n', exp(b_cox(1) * sd_dev));

        % --- 6b. Residualised interaction ---
        fprintf('--- Cox PH: residualised CDC_dev × Age ---\n\n');
        lm_dev_age = fitlm(data.Age, data.CDC_dev);
        data.CDC_dev_resid = lm_dev_age.Residuals.Raw;
        fprintf('  Age explains %.1f%% of variance in |CDC - 1/e|\n\n', ...
            100 * lm_dev_age.Rsquared.Ordinary);

        X_resid = [data.CDC_dev_resid, data.Age, ...
                   data.CDC_dev_resid .* data.Age, ...
                   data.is_male, data.Pathological];
        resid_names = {'CDC_dev_resid', 'Age', 'Resid x Age', ...
                       'Male', 'Pathological'};

        [b_resid, ~, ~, stats_resid] = coxphfit(X_resid, time, ...
            'Censoring', censoring);
        print_cox(resid_names, b_resid, stats_resid);

        fprintf('\n  Residualised interaction p: %.2e (raw: %.2e)\n', ...
            stats_resid.p(3), stats_cox.p(3));

        % --- 6c. Main effects only ---
        fprintf('\n--- Cox PH: main effects only ---\n\n');
        X_main = [data.CDC_dev, data.Age, data.is_male, data.Pathological];
        main_names = {'CDC_dev', 'Age', 'Male', 'Pathological'};

        [b_main, ~, ~, stats_main] = coxphfit(X_main, time, ...
            'Censoring', censoring);
        print_cox(main_names, b_main, stats_main);
        fprintf('  Per-SD HR for CDC_dev: %.3f\n', exp(b_main(1) * sd_dev));
    end

    %% ================================================================
    %  7. WITHIN-STRATUM PROPORTION TESTS
    %  ================================================================

    fprintf('\n--- Within-stratum tests: T1 vs T3 (tertiles) ---\n');
    pvals_tert = run_stratum_tests(data, n_age, age_labels, ...
        'DevTertile', 1, 3, 'T1', 'T3');

    fprintf('\n--- Within-stratum tests: Near vs Far (fixed cutoffs) ---\n');
    pvals_fix = run_stratum_tests(data, n_age, age_labels, ...
        'DevFixed', 1, 3, 'Near', 'Far');

    %% ================================================================
    %  8. SAVE RESULTS
    %  ================================================================

    results.data = data;
    results.inv_e = inv_e;
    results.age_edges = age_edges;
    results.age_labels = age_labels;
    results.n_total = height(data);
    results.n_deceased = sum(data.Deceased);

    % Primary: tertile results
    results.tertile.bounds = tertile_bounds;
    results.tertile.labels = tert_labels;
    results.tertile.mort_rate = mort_rate_tert;
    results.tertile.mort_n = mort_n_tert;
    results.tertile.mort_d = mort_d_tert;
    results.tertile.pvals = pvals_tert;

    % CMH test
    results.tertile.cmh_chi2 = cmh_chi2;
    results.tertile.cmh_p = cmh_p;

    % Robustness: fixed cutoffs
    results.fixed.dev_edges = dev_edges;
    results.fixed.labels = {'Near 1/e', 'Moderate', 'Far from 1/e'};
    results.fixed.mort_rate = mort_rate_fix;
    results.fixed.mort_n = mort_n_fix;
    results.fixed.mort_d = mort_d_fix;
    results.fixed.pvals = pvals_fix;

    % Cox models
    if exist('b_cox', 'var')
        results.cox_interaction.beta = b_cox;
        results.cox_interaction.stats = stats_cox;
    end
    if exist('b_resid', 'var')
        results.cox_resid.beta = b_resid;
        results.cox_resid.stats = stats_resid;
        results.lm_dev_on_age = lm_dev_age;
    end
    if exist('b_main', 'var')
        results.cox_main.beta = b_main;
        results.cox_main.stats = stats_main;
    end

    mat_file = fullfile(paths.results, 'age_stratified_mortality_results.mat');
    save(mat_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', mat_file);
    fprintf('================================================================\n');
end


%% ========================================================================
%  HELPERS
%  ========================================================================

function [mort_rate, mort_n, mort_d] = tabulate_mortality( ...
        data, n_age, n_dev, age_labels, dev_labels, age_col, dev_col)
% TABULATE_MORTALITY - Cross-tabulate mortality by age × deviation group

    mort_rate = NaN(n_age, n_dev);
    mort_n    = zeros(n_age, n_dev);
    mort_d    = zeros(n_age, n_dev);

    fprintf('%-8s', 'Age');
    for di = 1:n_dev
        fprintf('  %-28s', dev_labels{di});
    end
    fprintf('  %-10s\n', 'Fold');
    fprintf('%s\n', repmat('-', 1, 110));

    for ai = 1:n_age
        fprintf('%-8s', age_labels{ai});
        for di = 1:n_dev
            mask = data.(age_col) == ai & data.(dev_col) == di;
            n_cell = sum(mask);
            d_cell = sum(data.Deceased(mask));
            mort_n(ai, di) = n_cell;
            mort_d(ai, di) = d_cell;
            if n_cell > 0
                mort_rate(ai, di) = d_cell / n_cell;
                fprintf('  %5.2f%% (%5d/%6d)         ', ...
                    100 * mort_rate(ai, di), d_cell, n_cell);
            else
                fprintf('  %28s', '—');
            end
        end
        if mort_rate(ai, 1) > 0 && ~isnan(mort_rate(ai, n_dev))
            fprintf('  %.1f\times', mort_rate(ai, n_dev) / mort_rate(ai, 1));
        end
        fprintf('\n');
    end
end

function print_cox(names, beta, stats)
    fprintf('  %-20s %10s %10s %10s %10s\n', 'Covariate', 'beta', 'HR', 'z', 'p');
    fprintf('  %s\n', repmat('-', 1, 60));
    for i = 1:length(beta)
        z = beta(i) / stats.se(i);
        fprintf('  %-20s %10.4f %10.3f %10.3f %10.2e\n', ...
            names{i}, beta(i), exp(beta(i)), z, stats.p(i));
    end
end

function pvals = run_stratum_tests(data, n_age, age_labels, ...
        dev_col, bin_lo, bin_hi, label_lo, label_hi)
    pvals = NaN(n_age, 1);
    for ai = 1:n_age
        mask_lo = data.AgeBin == ai & data.(dev_col) == bin_lo;
        mask_hi = data.AgeBin == ai & data.(dev_col) == bin_hi;
        n_lo = sum(mask_lo);  d_lo = sum(data.Deceased(mask_lo));
        n_hi = sum(mask_hi);  d_hi = sum(data.Deceased(mask_hi));

        if n_lo > 10 && n_hi > 10 && d_lo > 0 && d_hi > 0
            [~, p_z] = prop_test([d_lo, d_hi], [n_lo, n_hi]);
            pvals(ai) = p_z;
            fprintf('  %s: %s=%.2f%% (%d/%d), %s=%.2f%% (%d/%d), p=%.2e\n', ...
                age_labels{ai}, ...
                label_lo, 100*d_lo/n_lo, d_lo, n_lo, ...
                label_hi, 100*d_hi/n_hi, d_hi, n_hi, p_z);
        else
            fprintf('  %s: Insufficient data for test\n', age_labels{ai});
        end
    end
end

function [z, p] = prop_test(successes, totals)
    p1 = successes(1) / totals(1);
    p2 = successes(2) / totals(2);
    pp = sum(successes) / sum(totals);
    se = sqrt(pp * (1 - pp) * (1/totals(1) + 1/totals(2)));
    z = (p1 - p2) / se;
    p = 2 * (1 - normcdf(abs(z)));
end
