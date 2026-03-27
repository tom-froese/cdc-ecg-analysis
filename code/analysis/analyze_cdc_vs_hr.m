function results = analyze_cdc_vs_hr()
% ANALYZE_CDC_VS_HR - CDC carries independent information beyond heart rate
%
%   PART 1 — PREDICTING THE HEALTHY HR FROM FIRST PRINCIPLES
%     Healthy controls (LUDB, Fantasia, Autonomic Aging; N ~ 1,165).
%     The intersection of the empirical CDC-HR regression with
%     CDC = 1/e predicts the optimal resting HR from first principles.
%     Bootstrap CI quantifies uncertainty.
%
%   PART 2 — HIDDEN RISK IN CLINICALLY NORMAL PATIENTS
%     CODE-15% patients with normal ECG only. Tests whether CDC
%     deviation predicts mortality independently of HR, even among
%     patients who passed standard clinical screening.
%
%     Subject aggregation uses unique_subject_id (CODE15_<patient_id>)
%     via vectorised splitapply for performance on large datasets.
%
% Dependencies:
%   config.m, apply_quality_filters.m
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    fprintf('================================================================\n');
    fprintf('CDC vs HEART RATE: INDEPENDENCE ANALYSIS\n');
    fprintf('================================================================\n');
    fprintf('Reference: 1/e = %.4f\n', inv_e);
    fprintf('================================================================\n\n');

    %% ================================================================
    %  PART 1: PREDICTING THE HEALTHY HR FROM FIRST PRINCIPLES
    %  ================================================================

    fprintf('================================================================\n');
    fprintf('PART 1: PREDICTING THE HEALTHY HR FROM FIRST PRINCIPLES\n');
    fprintf('================================================================\n\n');

    S_hier = load(fullfile(paths.results, 'hierarchical_results.mat'), 'all_data');
    all_data = S_hier.all_data;

    healthy = all_data(all_data.Group == 'HealthyControl', :);
    fprintf('Healthy control subjects: N = %d\n', height(healthy));
    fprintf('  Datasets: %s\n', strjoin(categories(removecats(healthy.Dataset)), ', '));
    fprintf('  Age range: %.0f - %.0f years (mean %.1f +/- %.1f)\n', ...
        min(healthy.Age), max(healthy.Age), mean(healthy.Age), std(healthy.Age));
    fprintf('  HR range:  %.1f - %.1f bpm (mean %.1f +/- %.1f)\n', ...
        min(healthy.HR), max(healthy.HR), mean(healthy.HR), std(healthy.HR));
    fprintf('  CDC range: %.4f - %.4f (mean %.4f +/- %.4f)\n\n', ...
        min(healthy.CDC), max(healthy.CDC), mean(healthy.CDC), std(healthy.CDC));

    % --- 1a. Correlations ---
    fprintf('--- 1a. CDC-HR correlations ---\n');
    [r_pearson, p_pearson] = corr(healthy.CDC, healthy.HR);
    [r_spearman, p_spearman] = corr(healthy.CDC, healthy.HR, 'Type', 'Spearman');
    fprintf('  Pearson:  r = %.4f, p = %.2e\n', r_pearson, p_pearson);
    fprintf('  Spearman: rho = %.4f, p = %.2e\n\n', r_spearman, p_spearman);

    fprintf('  Per-dataset:\n');
    datasets = categories(removecats(healthy.Dataset));
    for d = 1:length(datasets)
        mask_d = healthy.Dataset == datasets{d};
        if sum(mask_d) >= 10
            [r_d, p_d] = corr(healthy.CDC(mask_d), healthy.HR(mask_d));
            fprintf('    %-16s  r = %.4f, p = %.2e  (N = %d)\n', ...
                datasets{d}, r_d, p_d, sum(mask_d));
        end
    end

    % --- 1b. Linear regression and intersection ---
    fprintf('\n--- 1b. Linear regression: CDC ~ HR ---\n');
    mdl_cdc_hr = fitlm(healthy.HR, healthy.CDC);
    b0 = mdl_cdc_hr.Coefficients.Estimate(1);
    b1 = mdl_cdc_hr.Coefficients.Estimate(2);
    fprintf('  Slope:     %.6f per bpm (SE = %.6f, p = %.2e)\n', ...
        b1, mdl_cdc_hr.Coefficients.SE(2), mdl_cdc_hr.Coefficients.pValue(2));
    fprintf('  Intercept: %.4f\n', b0);
    fprintf('  R-squared: %.4f\n', mdl_cdc_hr.Rsquared.Ordinary);
    fprintf('  => HR explains %.1f%% of CDC variance in healthy controls.\n\n', ...
        100 * mdl_cdc_hr.Rsquared.Ordinary);

    hr_opt = (inv_e - b0) / b1;
    fprintf('  Intersection with CDC = 1/e:\n');
    fprintf('    HR_opt = %.2f bpm\n', hr_opt);

    % --- 1c. Bootstrap CI for intersection HR ---
    fprintf('\n--- 1c. Bootstrap CI for intersection HR ---\n');
    n_boot = 5000;
    n_healthy = height(healthy);
    rng(42);
    boot_hr_opt = zeros(n_boot, 1);
    hr_vec = healthy.HR;
    cdc_vec = healthy.CDC;

    for bi = 1:n_boot
        idx = randi(n_healthy, n_healthy, 1);
        X = [ones(n_healthy,1), hr_vec(idx)];
        coeffs = X \ cdc_vec(idx);      % fast OLS via backslash
        if coeffs(2) ~= 0
            boot_hr_opt(bi) = (inv_e - coeffs(1)) / coeffs(2);
        else
            boot_hr_opt(bi) = NaN;
        end
    end
    boot_hr_opt = boot_hr_opt(~isnan(boot_hr_opt));
    hr_opt_ci = [prctile(boot_hr_opt, 2.5), prctile(boot_hr_opt, 97.5)];

    fprintf('    HR_opt = %.2f bpm [%.2f, %.2f] (95%% bootstrap CI, B=%d)\n', ...
        hr_opt, hr_opt_ci(1), hr_opt_ci(2), n_boot);

    if hr_opt >= 60 && hr_opt <= 80
        fprintf('    => Falls within the optimal healthy range (60-80 bpm).\n');
    elseif hr_opt >= 50 && hr_opt <= 100
        fprintf('    => Falls within the normal resting range (50-100 bpm).\n');
    end

    % --- 1d. CDC across HR bins ---
    fprintf('\n--- 1d. CDC across HR bins ---\n');
    hr_min_bin = 5 * floor(min(healthy.HR) / 5);
    hr_max_bin = 5 * ceil(max(healthy.HR) / 5);
    hc_hr_edges = hr_min_bin:5:hr_max_bin;
    n_hc_bins = length(hc_hr_edges) - 1;

    hc_bin_stats = struct('label', {}, 'center', {}, 'n', {}, ...
        'mean_cdc', {}, 'median_cdc', {}, 'std_cdc', {}, ...
        'ci_lo', {}, 'ci_hi', {}, 'mean_hr', {});

    fprintf('  %-12s  %6s  %10s  %10s  %10s\n', 'HR bin', 'N', 'Mean CDC', 'Median CDC', 'SD');
    fprintf('  %s\n', repmat('-', 1, 55));

    for hi = 1:n_hc_bins
        mask = healthy.HR >= hc_hr_edges(hi) & healthy.HR < hc_hr_edges(hi+1);
        n_bin = sum(mask);
        lbl = sprintf('%d-%d', hc_hr_edges(hi), hc_hr_edges(hi+1));
        hc_bin_stats(hi).label = lbl;
        hc_bin_stats(hi).center = (hc_hr_edges(hi) + hc_hr_edges(hi+1)) / 2;
        hc_bin_stats(hi).n = n_bin;
        if n_bin >= 5
            hc_bin_stats(hi).mean_cdc   = mean(healthy.CDC(mask));
            hc_bin_stats(hi).median_cdc = median(healthy.CDC(mask));
            hc_bin_stats(hi).std_cdc    = std(healthy.CDC(mask));
            se = hc_bin_stats(hi).std_cdc / sqrt(n_bin);
            hc_bin_stats(hi).ci_lo = hc_bin_stats(hi).mean_cdc - 1.96 * se;
            hc_bin_stats(hi).ci_hi = hc_bin_stats(hi).mean_cdc + 1.96 * se;
            hc_bin_stats(hi).mean_hr = mean(healthy.HR(mask));
            fprintf('  %-12s  %6d  %10.4f  %10.4f  %10.4f\n', ...
                lbl, n_bin, hc_bin_stats(hi).mean_cdc, ...
                hc_bin_stats(hi).median_cdc, hc_bin_stats(hi).std_cdc);
        else
            hc_bin_stats(hi).mean_cdc = NaN; hc_bin_stats(hi).median_cdc = NaN;
            hc_bin_stats(hi).std_cdc = NaN;  hc_bin_stats(hi).ci_lo = NaN;
            hc_bin_stats(hi).ci_hi = NaN;    hc_bin_stats(hi).mean_hr = NaN;
            fprintf('  %-12s  %6d          -           -           -\n', lbl, n_bin);
        end
    end

    % --- 1e. One-sample t-tests ---
    fprintf('\n--- 1e. One-sample t-tests: mean CDC vs 1/e per HR bin ---\n');
    for hi = 1:n_hc_bins
        mask = healthy.HR >= hc_hr_edges(hi) & healthy.HR < hc_hr_edges(hi+1);
        if sum(mask) >= 10
            [~, p_t, ~, t_stats] = ttest(healthy.CDC(mask), inv_e);
            fprintf('    HR %s:  mean = %.4f,  t(%d) = %+.2f,  p = %.2e\n', ...
                hc_bin_stats(hi).label, hc_bin_stats(hi).mean_cdc, ...
                t_stats.df, t_stats.tstat, p_t);
        end
    end

    %% ================================================================
    %  PART 2: HIDDEN RISK IN CLINICALLY NORMAL PATIENTS (CODE-15%)
    %  ================================================================

    fprintf('\n\n================================================================\n');
    fprintf('PART 2: HIDDEN RISK IN CLINICALLY NORMAL PATIENTS\n');
    fprintf('================================================================\n\n');

    fprintf('Loading CODE-15%% beat-level data...\n');
    beats = readtable(paths.csv_code15, 'TextType', 'string');
    fprintf('  %d beats loaded\n', height(beats));

    fprintf('Applying quality filters (no CDC ratio cutoffs)...\n');
    [valid_mask, ~] = apply_quality_filters(beats, 'CDCMin', 0, 'CDCMax', 1);
    beats = beats(valid_mask, :);
    fprintf('  %d beats after filtering\n', height(beats));

    % --- Vectorised subject-level aggregation via splitapply ---
    %  This replaces the per-subject for-loop, reducing runtime from
    %  ~30 min to ~30 sec on typical hardware.
    fprintf('  Aggregating to subject level (vectorised)...\n');

    fs = beats.fs(1);
    rt_ms = (beats.t_end_sample - beats.r_sample) / fs * 1000;
    rr_ms = (beats.next_r_sample - beats.r_sample) / fs * 1000;
    beat_ratios = rt_ms ./ rr_ms;
    valid_beats = (rr_ms > 0) & (rt_ms > 0) & (rt_ms < rr_ms);

    [G, subject_ids] = findgroups(beats.unique_subject_id);

    CDC = splitapply(@(x, v) median(x(v)), beat_ratios, valid_beats, G);
    HR  = splitapply(@(x, v) 60000 ./ median(x(v)), rr_ms, valid_beats, G);
    Age = splitapply(@(x) x(1), beats.age, G);

    % Sex: first value per subject
    if isstring(beats.sex) || iscell(beats.sex)
        Sex_str = splitapply(@(x) x(1), string(beats.sex), G);
    else
        Sex_str = splitapply(@(x) string(x(1)), beats.sex, G);
    end

    % Record ID (patient_id for mortality merge)
    Record_id = splitapply(@(x) x(1), string(beats.record_id), G);

    % Group: ClinicallyNormal if ALL beats have source_subset == "normal"
    is_normal_beat = (beats.source_subset == "normal");
    total_per   = splitapply(@numel, beats.record_id, G);
    normal_per  = splitapply(@sum, is_normal_beat, G);
    is_cn = (total_per == normal_per);

    Group_str = strings(length(subject_ids), 1);
    Group_str(is_cn)  = "ClinicallyNormal";
    Group_str(~is_cn) = "Pathological";

    fprintf('  %d unique subjects\n', length(subject_ids));
    fprintf('    ClinicallyNormal: %d\n', sum(is_cn));
    fprintf('    Pathological:     %d\n', sum(~is_cn));

    beat_data = table(Record_id, subject_ids, Age, CDC, HR, Sex_str, Group_str, ...
        'VariableNames', {'record_id', 'unique_subject_id', 'Age', 'CDC', ...
                          'HR', 'Sex', 'Group'});

    % --- Merge with mortality annotations ---
    if isfield(paths, 'csv_code15_exams')
        exams_csv_path = paths.csv_code15_exams;
    else
        exams_csv_path = fullfile(paths.raw_code15, 'exams.csv');
    end
    fprintf('\nLoading exams metadata from:\n  %s\n', exams_csv_path);

    exams = readtable(exams_csv_path, 'TextType', 'string', ...
        'TreatAsMissing', '', 'VariableNamingRule', 'preserve');
    has_timey = ismember('timey', exams.Properties.VariableNames);

    exams = sortrows(exams, {'patient_id', 'exam_id'}, {'ascend', 'ascend'});
    [~, first_idx] = unique(exams.patient_id, 'first');
    exams_first = exams(first_idx, :);

    exams_first.merge_id = string(exams_first.patient_id);
    beat_data.merge_id   = string(beat_data.record_id);

    mortality_vars = intersect(exams_first.Properties.VariableNames, {'death', 'timey'});
    merged = innerjoin(beat_data, exams_first, ...
        'LeftKeys', 'merge_id', 'RightKeys', 'merge_id', ...
        'RightVariables', mortality_vars);
    fprintf('  Merged records: %d\n', height(merged));

    % Parse mortality
    if isstring(merged.death) || iscell(merged.death)
        merged.Deceased = lower(string(merged.death)) == "true" | merged.death == "1";
    elseif islogical(merged.death)
        merged.Deceased = merged.death;
    else
        merged.Deceased = (merged.death == 1);
    end
    if has_timey && ismember('timey', merged.Properties.VariableNames)
        if isstring(merged.timey) || iscell(merged.timey)
            merged.FollowUp_yrs = str2double(merged.timey);
        else
            merged.FollowUp_yrs = merged.timey;
        end
    end

    valid_mort = ~isnan(double(merged.Deceased));
    if has_timey && ismember('FollowUp_yrs', merged.Properties.VariableNames)
        valid_mort = valid_mort & ~isnan(merged.FollowUp_yrs) & (merged.FollowUp_yrs > 0);
    end
    valid_mort = valid_mort & ~isnan(merged.CDC) & ~isnan(merged.Age) & ...
                 ~isnan(merged.HR) & (merged.Age > 0);

    data_all = merged(valid_mort, :);
    data_all.Group = categorical(data_all.Group);
    data_all.Sex = categorical(data_all.Sex);

    fprintf('  Full cohort: N = %d  (Deceased: %d, %.2f%%)\n', ...
        height(data_all), sum(data_all.Deceased), 100*mean(data_all.Deceased));

    % --- RESTRICT TO CLINICALLY NORMAL ---
    data = data_all(data_all.Group == 'ClinicallyNormal', :);
    data.Age_c   = data.Age - mean(data.Age);
    data.CDC_dev = abs(data.CDC - inv_e);
    data.is_male = double(data.Sex == 'Male' | data.Sex == 'M');
    n_cn = height(data);
    n_cn_dead = sum(data.Deceased);

    fprintf('\n  *** Restricting to ClinicallyNormal ***\n');
    fprintf('  N = %d  (Deceased: %d, %.2f%%)\n', n_cn, n_cn_dead, 100*mean(data.Deceased));
    fprintf('  HR range: %.1f - %.1f bpm (mean %.1f +/- %.1f)\n', ...
        min(data.HR), max(data.HR), mean(data.HR), std(data.HR));
    fprintf('  CDC range: %.4f - %.4f (mean %.4f +/- %.4f)\n', ...
        min(data.CDC), max(data.CDC), mean(data.CDC), std(data.CDC));

    % --- 2a. VIF and partial correlations ---
    fprintf('\n--- 2a. VIF and partial correlations (ClinicallyNormal) ---\n');
    [r_cn, p_cn] = corr(data.CDC, data.HR);
    fprintf('  CDC-HR correlation: r = %.4f, p = %.2e\n', r_cn, p_cn);
    r_cdc_hr_dev = corr(data.CDC_dev, data.HR);
    vif = 1 / (1 - r_cdc_hr_dev^2);
    fprintf('  CDC_dev-HR correlation: r = %.4f,  VIF = %.2f\n', r_cdc_hr_dev, vif);

    covariates_cn = [data.Age_c, data.is_male];
    r_cdc_mort_partial = partial_corr(data.CDC_dev, double(data.Deceased), ...
        [data.HR, covariates_cn]);
    r_hr_mort_partial = partial_corr(data.HR, double(data.Deceased), ...
        [data.CDC_dev, covariates_cn]);
    fprintf('  Partial r (CDC_dev-Mortality | HR, Age, Sex): %.4f\n', r_cdc_mort_partial);
    fprintf('  Partial r (HR-Mortality | CDC_dev, Age, Sex): %.4f\n', r_hr_mort_partial);

    % --- 2b. Within-HR-stratum mortality by CDC tertile ---
    fprintf('\n--- 2b. Within-HR-stratum mortality (ClinicallyNormal) ---\n\n');

    hr_edges = [50 55 60 65 70 75 80 85 90 95 100];
    hr_labels = arrayfun(@(i) sprintf('%d-%d', hr_edges(i), hr_edges(i+1)), ...
        1:length(hr_edges)-1, 'UniformOutput', false)';
    n_hr_bins = length(hr_labels);

    cdc_tert_edges = quantile(data.CDC_dev, [1/3, 2/3]);
    tert_labels = {'Near 1/e', 'Moderate', 'Far from 1/e'};
    fprintf('  CDC_dev tertile thresholds: [%.4f, %.4f]\n\n', ...
        cdc_tert_edges(1), cdc_tert_edges(2));

    mort_rate = NaN(n_hr_bins, 3);
    mort_n    = zeros(n_hr_bins, 3);
    mort_dead = zeros(n_hr_bins, 3);

    fprintf('  %-10s  %-15s  %6s  %6s  %10s\n', 'HR bin', 'CDC tertile', 'N', 'Deaths', 'Rate (%)');
    fprintf('  %s\n', repmat('-', 1, 55));

    for hi = 1:n_hr_bins
        mask_hr = data.HR >= hr_edges(hi) & data.HR < hr_edges(hi+1);
        for ti = 1:3
            if ti == 1,     mask_t = data.CDC_dev <= cdc_tert_edges(1);
            elseif ti == 2, mask_t = data.CDC_dev > cdc_tert_edges(1) & data.CDC_dev <= cdc_tert_edges(2);
            else,           mask_t = data.CDC_dev > cdc_tert_edges(2);
            end
            mask = mask_hr & mask_t; n_cell = sum(mask);
            mort_n(hi, ti) = n_cell;
            if n_cell >= 10
                mort_dead(hi, ti) = sum(data.Deceased(mask));
                mort_rate(hi, ti) = mean(data.Deceased(mask));
                fprintf('  %-10s  %-15s  %6d  %6d  %10.2f\n', ...
                    hr_labels{hi}, tert_labels{ti}, n_cell, mort_dead(hi,ti), 100*mort_rate(hi,ti));
            else
                fprintf('  %-10s  %-15s  %6d      -         -\n', hr_labels{hi}, tert_labels{ti}, n_cell);
            end
        end
        fprintf('\n');
    end

    % Fold-change summary
    fprintf('  --- Fold-change (Far vs Near 1/e) per HR stratum ---\n');
    for hi = 1:n_hr_bins
        if mort_n(hi,1) >= 30 && mort_n(hi,3) >= 30 && mort_rate(hi,1) > 0 && ~isnan(mort_rate(hi,3))
            fold = mort_rate(hi,3) / mort_rate(hi,1);
            [z, p] = prop_test([mort_dead(hi,3), mort_dead(hi,1)], [mort_n(hi,3), mort_n(hi,1)]);
            fprintf('    HR %s bpm:  %.2f%% vs %.2f%%  (%.1fx, z=%.2f, p=%.2e)\n', ...
                hr_labels{hi}, 100*mort_rate(hi,3), 100*mort_rate(hi,1), fold, z, p);
        end
    end

    % CMH test
    fprintf('\n  --- Cochran-Mantel-Haenszel test ---\n');
    cmh_num = 0; cmh_denom = 0;
    for hi = 1:n_hr_bins
        if mort_n(hi,1) >= 10 && mort_n(hi,3) >= 10
            a = mort_dead(hi,1); b = mort_n(hi,1)-a;
            c = mort_dead(hi,3); d = mort_n(hi,3)-c;
            n_str = a+b+c+d;
            if n_str > 0
                cmh_num   = cmh_num + (a - (a+b)*(a+c)/n_str);
                cmh_denom = cmh_denom + (a+b)*(c+d)*(a+c)*(b+d)/(n_str^2*(n_str-1));
            end
        end
    end
    cmh_chi2 = cmh_num^2 / cmh_denom;
    cmh_p = 1 - chi2cdf(cmh_chi2, 1);
    fprintf('    CMH chi2 = %.2f, df = 1, p = %.2e\n', cmh_chi2, cmh_p);

    % --- 2c. Nested Cox models ---
    fprintf('\n--- 2c. Nested Cox models (ClinicallyNormal only) ---\n\n');
    has_fu = ismember('FollowUp_yrs', data.Properties.VariableNames);
    cox_results = struct();

    if has_fu
        censoring = ~data.Deceased; time = data.FollowUp_yrs;

        X_A = [data.HR, data.Age_c, data.is_male];
        [b_A, logL_A, ~, stats_A] = coxphfit(X_A, time, 'Censoring', censoring);
        fprintf('Model A: survival ~ HR + Age + Sex\n');
        print_cox_table({'HR','Age_c','Male'}, b_A, stats_A);
        fprintf('  Log-L: %.2f\n\n', logL_A);

        X_B = [data.HR, data.CDC_dev, data.Age_c, data.is_male];
        [b_B, logL_B, ~, stats_B] = coxphfit(X_B, time, 'Censoring', censoring);
        fprintf('Model B: survival ~ HR + CDC_dev + Age + Sex\n');
        print_cox_table({'HR','CDC_dev','Age_c','Male'}, b_B, stats_B);
        fprintf('  Log-L: %.2f\n\n', logL_B);

        X_C = [data.CDC_dev, data.Age_c, data.is_male];
        [b_C, logL_C, ~, stats_C] = coxphfit(X_C, time, 'Censoring', censoring);
        fprintf('Model C: survival ~ CDC_dev + Age + Sex\n');
        print_cox_table({'CDC_dev','Age_c','Male'}, b_C, stats_C);
        fprintf('  Log-L: %.2f\n\n', logL_C);

        LR_BA = 2*(logL_B - logL_A); p_BA = 1 - chi2cdf(LR_BA, 1);
        LR_BC = 2*(logL_B - logL_C); p_BC = 1 - chi2cdf(LR_BC, 1);
        fprintf('--- Nested model comparisons ---\n');
        fprintf('  B vs A (CDC_dev added to HR):  LR chi2=%.2f, p=%.2e\n', LR_BA, p_BA);
        fprintf('  B vs C (HR added to CDC_dev):  LR chi2=%.2f, p=%.2e\n\n', LR_BC, p_BC);

        sd_hr = std(data.HR); sd_cdc = std(data.CDC_dev);
        fprintf('  Per-SD HRs in Model B: HR=%.3f (SD=%.1f), CDC_dev=%.3f (SD=%.4f)\n', ...
            exp(b_B(1)*sd_hr), sd_hr, exp(b_B(2)*sd_cdc), sd_cdc);

        cox_results.logL_A=logL_A; cox_results.logL_B=logL_B; cox_results.logL_C=logL_C;
        cox_results.LR_BA=LR_BA; cox_results.p_BA=p_BA;
        cox_results.LR_BC=LR_BC; cox_results.p_BC=p_BC;
        cox_results.b_A=b_A; cox_results.stats_A=stats_A;
        cox_results.b_B=b_B; cox_results.stats_B=stats_B;
        cox_results.b_C=b_C; cox_results.stats_C=stats_C;
    else
        fprintf('  Follow-up time not available.\n');
    end

    %% ================================================================
    %  SAVE RESULTS
    %  ================================================================

    results.inv_e = inv_e;

    results.healthy.n = height(healthy);
    results.healthy.r_pearson = r_pearson; results.healthy.p_pearson = p_pearson;
    results.healthy.r_spearman = r_spearman; results.healthy.p_spearman = p_spearman;
    results.healthy.r2_cdc_on_hr = mdl_cdc_hr.Rsquared.Ordinary;
    results.healthy.slope = b1; results.healthy.intercept = b0;
    results.healthy.slope_p = mdl_cdc_hr.Coefficients.pValue(2);
    results.healthy.hr_opt = hr_opt;
    results.healthy.hr_opt_ci = hr_opt_ci;
    results.healthy.hr_edges = hc_hr_edges;
    results.healthy.bin_stats = hc_bin_stats;
    results.healthy.data = healthy;

    results.code15_cn.data = data;
    results.code15_cn.n = n_cn;
    results.code15_cn.n_deceased = n_cn_dead;
    results.code15_cn.coupling.r_cdc_hr = r_cn;
    results.code15_cn.coupling.r_cdc_hr_dev = r_cdc_hr_dev;
    results.code15_cn.coupling.vif = vif;
    results.code15_cn.coupling.r_cdc_mort_partial = r_cdc_mort_partial;
    results.code15_cn.coupling.r_hr_mort_partial = r_hr_mort_partial;
    if has_fu, results.code15_cn.cox = cox_results; end
    results.code15_cn.hr_strata.hr_edges = hr_edges;
    results.code15_cn.hr_strata.hr_labels = hr_labels;
    results.code15_cn.hr_strata.cdc_tert_edges = cdc_tert_edges;
    results.code15_cn.hr_strata.tert_labels = tert_labels;
    results.code15_cn.hr_strata.mort_rate = mort_rate;
    results.code15_cn.hr_strata.mort_n = mort_n;
    results.code15_cn.hr_strata.mort_dead = mort_dead;
    results.code15_cn.hr_strata.cmh_chi2 = cmh_chi2;
    results.code15_cn.hr_strata.cmh_p = cmh_p;

    mat_file = fullfile(paths.results, 'cdc_vs_hr_results.mat');
    save(mat_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', mat_file);
    fprintf('================================================================\n');
end


%% ========================================================================
%  HELPERS
%  ========================================================================

function r = partial_corr(x, y, Z)
    x = x(:); y = y(:); Z = [ones(size(Z,1),1), Z];
    r = corr(x - Z*(Z\x), y - Z*(Z\y));
end

function print_cox_table(names, b, stats)
    fprintf('  %-20s %10s %10s %10s %10s\n', 'Covariate', 'beta', 'HR', 'z', 'p');
    fprintf('  %s\n', repmat('-', 1, 60));
    for i = 1:length(b)
        fprintf('  %-20s %10.4f %10.3f %10.3f %10.2e\n', ...
            names{i}, b(i), exp(b(i)), b(i)/stats.se(i), stats.p(i));
    end
end

function [z, p] = prop_test(successes, totals)
    p1 = successes(1)/totals(1); p2 = successes(2)/totals(2);
    pp = sum(successes)/sum(totals);
    se = sqrt(pp*(1-pp)*(1/totals(1)+1/totals(2)));
    z = (p1-p2)/se; p = 2*(1-normcdf(abs(z)));
end
