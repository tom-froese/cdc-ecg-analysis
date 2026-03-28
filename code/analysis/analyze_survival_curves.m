function results = analyze_survival_curves()
% ANALYZE_SURVIVAL_CURVES - Kaplan-Meier survival analysis by CDC deviation
%
% Computes Kaplan-Meier survival curves for the CODE-15% cohort,
% stratified by CDC-deviation tertile and (optionally) by sex.
%
% Analyses:
%   1. Overall KM curves by within-subject CDC-deviation tertile
%   2. Log-rank tests (pairwise and omnibus) for tertile separation
%   3. Sex-stratified KM curves by CDC-deviation tertile
%   4. Sex-stratified log-rank tests
%   5. Median survival time estimates per stratum
%   6. Restricted mean survival time (RMST) differences at fixed horizon
%
% The deviation tertiles are computed globally (not within age strata)
% to provide a single, interpretable stratification for the KM curves.
% This complements the age-stratified tertile analysis in Figure 2
% (analyze_age_stratified_mortality.m), which uses within-stratum
% tertiles to control for confounding.
%
% Subject aggregation uses unique_subject_id for consistency with the
% unified beat table format (cf. analyze_age_stratified_mortality.m).
%
% Dependencies: config.m, apply_quality_filters.m
% Data:         CODE-15% beat extraction + exams.csv mortality labels
% Output:       survival_curve_results.mat
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    fprintf('================================================================\n');
    fprintf('KAPLAN-MEIER SURVIVAL ANALYSIS: CDC DEVIATION TERTILES\n');
    fprintf('================================================================\n');
    fprintf('Reference optimum: 1/e = %.4f\n\n', inv_e);

    %% ================================================================
    %  1. LOAD AND AGGREGATE TO SUBJECT LEVEL
    %  ================================================================
    %  Identical to analyze_age_stratified_mortality.m — ensures that the
    %  same subject-level dataset is used for both analyses.

    fprintf('Loading beat-level data...\n');
    beats = readtable(paths.csv_code15, 'TextType', 'string');
    fprintf('  %d beats loaded\n', height(beats));

    fprintf('Applying quality filters (no CDC ratio cutoffs)...\n');
    [valid_mask, ~] = apply_quality_filters(beats, ...
        'CDCMin', 0, 'CDCMax', 1);
    beats = beats(valid_mask, :);
    fprintf('  %d beats after filtering\n', height(beats));

    % --- Subject-level aggregation via unique_subject_id ---
    fs = beats.fs(1);
    rt_ms = (beats.t_end_sample - beats.r_sample) / fs * 1000;
    rr_ms = (beats.next_r_sample - beats.r_sample) / fs * 1000;
    beat_cdc = rt_ms ./ rr_ms;
    valid_beats = (rr_ms > 0) & (rt_ms > 0) & (rt_ms < rr_ms);

    [G, subject_ids] = findgroups(beats.unique_subject_id);

    CDC = splitapply(@(x, m) median(x(m)), beat_cdc, valid_beats, G);
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

    beat_data = table(subject_ids, merge_id, Age, CDC, Sex_str, Group_str, ...
        'VariableNames', {'unique_subject_id', 'merge_id', ...
                          'Age', 'CDC', 'Sex', 'Group'});

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

    % Parse follow-up time (required for survival analysis)
    has_timey = ismember('timey', merged.Properties.VariableNames);
    if ~has_timey
        error(['exams.csv must contain a "timey" column (follow-up years) ', ...
               'for survival analysis.']);
    end
    if isstring(merged.timey) || iscell(merged.timey)
        merged.FollowUp_yrs = str2double(merged.timey);
    else
        merged.FollowUp_yrs = merged.timey;
    end

    % Apply validity filters
    valid = ~isnan(double(merged.Deceased)) & ...
            ~isnan(merged.FollowUp_yrs) & (merged.FollowUp_yrs > 0) & ...
            ~isnan(merged.CDC) & ~isnan(merged.Age) & (merged.Age > 0);

    data = merged(valid, :);
    data.Group = categorical(data.Group);
    data.Sex = categorical(data.Sex);
    data.CDC_dev = abs(data.CDC - inv_e);
    data.Pathological = double(data.Group == 'Pathological');
    data.Censored = ~data.Deceased;  % ecmnfish convention: 1 = censored

    fprintf('  Valid records: %d (Deaths: %d, %.2f%%)\n', ...
        height(data), sum(data.Deceased), 100 * mean(data.Deceased));
    fprintf('  Median follow-up: %.2f years\n', median(data.FollowUp_yrs));

    %% ================================================================
    %  3. DEFINE CDC-DEVIATION TERTILES (GLOBAL)
    %  ================================================================

    q33 = quantile(data.CDC_dev, 1/3);
    q67 = quantile(data.CDC_dev, 2/3);

    data.DevTertile = ones(height(data), 1);
    data.DevTertile(data.CDC_dev >= q33) = 2;
    data.DevTertile(data.CDC_dev >= q67) = 3;

    tert_labels = {'Near 1/e (T1)', 'Moderate (T2)', 'Far from 1/e (T3)'};

    fprintf('\nGlobal CDC-deviation tertile boundaries:\n');
    fprintf('  T1/T2 boundary: |CDC - 1/e| = %.4f\n', q33);
    fprintf('  T2/T3 boundary: |CDC - 1/e| = %.4f\n', q67);

    for ti = 1:3
        mask = data.DevTertile == ti;
        fprintf('  %s: N = %d, Deaths = %d (%.2f%%)\n', ...
            tert_labels{ti}, sum(mask), sum(data.Deceased(mask)), ...
            100 * mean(data.Deceased(mask)));
    end

    %% ================================================================
    %  4. KAPLAN-MEIER ESTIMATION — OVERALL
    %  ================================================================

    fprintf('\n================================================================\n');
    fprintf('KAPLAN-MEIER ESTIMATION — BY CDC-DEVIATION TERTILE\n');
    fprintf('================================================================\n\n');

    km_overall = struct();
    for ti = 1:3
        mask = data.DevTertile == ti;
        [km_t, km_s, km_lo, km_hi, km_at_risk] = kaplan_meier( ...
            data.FollowUp_yrs(mask), data.Deceased(mask));

        km_overall(ti).t = km_t;
        km_overall(ti).s = km_s;
        km_overall(ti).lo = km_lo;
        km_overall(ti).hi = km_hi;
        km_overall(ti).at_risk = km_at_risk;
        km_overall(ti).n = sum(mask);
        km_overall(ti).d = sum(data.Deceased(mask));
        km_overall(ti).label = tert_labels{ti};

        med_surv = median_survival(km_t, km_s);
        km_overall(ti).median_surv = med_surv;

        if isnan(med_surv)
            fprintf('  %s: median survival not reached\n', tert_labels{ti});
        else
            fprintf('  %s: median survival = %.2f years\n', tert_labels{ti}, med_surv);
        end
    end

    %% ================================================================
    %  5. LOG-RANK TESTS — OVERALL
    %  ================================================================

    fprintf('\n--- Log-rank tests (overall) ---\n');

    % Omnibus: all three groups
    [chi2_omni, p_omni] = logrank_test( ...
        data.FollowUp_yrs, data.Deceased, data.DevTertile);
    fprintf('  Omnibus (3-group): chi2 = %.2f, df = 2, p = %.2e\n', ...
        chi2_omni, p_omni);

    % Pairwise: T1 vs T3
    mask_13 = (data.DevTertile == 1) | (data.DevTertile == 3);
    grp_13 = data.DevTertile(mask_13);
    [chi2_13, p_13] = logrank_test( ...
        data.FollowUp_yrs(mask_13), data.Deceased(mask_13), grp_13);
    fprintf('  T1 vs T3:          chi2 = %.2f, df = 1, p = %.2e\n', ...
        chi2_13, p_13);

    % Pairwise: T1 vs T2
    mask_12 = (data.DevTertile == 1) | (data.DevTertile == 2);
    grp_12 = data.DevTertile(mask_12);
    [chi2_12, p_12] = logrank_test( ...
        data.FollowUp_yrs(mask_12), data.Deceased(mask_12), grp_12);
    fprintf('  T1 vs T2:          chi2 = %.2f, df = 1, p = %.2e\n', ...
        chi2_12, p_12);

    % Pairwise: T2 vs T3
    mask_23 = (data.DevTertile == 2) | (data.DevTertile == 3);
    grp_23 = data.DevTertile(mask_23);
    [chi2_23, p_23] = logrank_test( ...
        data.FollowUp_yrs(mask_23), data.Deceased(mask_23), grp_23);
    fprintf('  T2 vs T3:          chi2 = %.2f, df = 1, p = %.2e\n', ...
        chi2_23, p_23);

    %% ================================================================
    %  6. RESTRICTED MEAN SURVIVAL TIME (RMST)
    %  ================================================================
    %  RMST is the area under the KM curve up to a fixed time horizon.
    %  Unlike median survival (which may not be reached), RMST always
    %  yields a finite estimand. We use the minimum of the maximum
    %  observed times across groups as the truncation point.

    fprintf('\n--- Restricted Mean Survival Time (RMST) ---\n');

    tau_max = Inf;
    for ti = 1:3
        tau_max = min(tau_max, max(km_overall(ti).t));
    end
    % Use a round number close to the minimum max follow-up
    tau = floor(tau_max);
    fprintf('  Truncation horizon: tau = %d years\n', tau);

    for ti = 1:3
        rmst_val = compute_rmst(km_overall(ti).t, km_overall(ti).s, tau);
        km_overall(ti).rmst = rmst_val;
        fprintf('  %s: RMST(%.0f yr) = %.3f years\n', ...
            tert_labels{ti}, tau, rmst_val);
    end

    rmst_diff_13 = km_overall(1).rmst - km_overall(3).rmst;
    fprintf('  RMST difference (T1 - T3): %.3f years\n', rmst_diff_13);

    %% ================================================================
    %  7. SEX-STRATIFIED KAPLAN-MEIER
    %  ================================================================

    fprintf('\n================================================================\n');
    fprintf('SEX-STRATIFIED KAPLAN-MEIER\n');
    fprintf('================================================================\n\n');

    % Determine sex encoding
    data.is_male = (data.Sex == 'Male') | (data.Sex == 'M') | ...
                   (data.Sex == 'male') | (data.Sex == 'm');

    sex_labels = {'Female', 'Male'};
    sex_masks  = {~data.is_male, data.is_male};

    km_by_sex = struct();

    for si = 1:2
        fprintf('--- %s (N = %d, Deaths = %d) ---\n', ...
            sex_labels{si}, sum(sex_masks{si}), ...
            sum(data.Deceased(sex_masks{si})));

        km_sex_group = struct();

        for ti = 1:3
            mask = sex_masks{si} & (data.DevTertile == ti);
            [km_t, km_s, km_lo, km_hi, km_ar] = kaplan_meier( ...
                data.FollowUp_yrs(mask), data.Deceased(mask));

            km_sex_group(ti).t = km_t;
            km_sex_group(ti).s = km_s;
            km_sex_group(ti).lo = km_lo;
            km_sex_group(ti).hi = km_hi;
            km_sex_group(ti).at_risk = km_ar;
            km_sex_group(ti).n = sum(mask);
            km_sex_group(ti).d = sum(data.Deceased(mask));
            km_sex_group(ti).label = tert_labels{ti};

            med_surv = median_survival(km_t, km_s);
            km_sex_group(ti).median_surv = med_surv;

            rmst_val = compute_rmst(km_t, km_s, tau);
            km_sex_group(ti).rmst = rmst_val;

            fprintf('  %s: N=%d, Deaths=%d (%.2f%%), RMST=%.3f\n', ...
                tert_labels{ti}, sum(mask), sum(data.Deceased(mask)), ...
                100 * mean(data.Deceased(mask)), rmst_val);
        end

        % Log-rank test within this sex
        sex_data = data(sex_masks{si}, :);
        [chi2_sex, p_sex] = logrank_test( ...
            sex_data.FollowUp_yrs, sex_data.Deceased, sex_data.DevTertile);
        fprintf('  Log-rank (omnibus): chi2 = %.2f, p = %.2e\n\n', ...
            chi2_sex, p_sex);

        km_by_sex(si).label = sex_labels{si};
        km_by_sex(si).km = km_sex_group;
        km_by_sex(si).logrank_chi2 = chi2_sex;
        km_by_sex(si).logrank_p = p_sex;
        km_by_sex(si).n = sum(sex_masks{si});
        km_by_sex(si).d = sum(data.Deceased(sex_masks{si}));
    end

    %% ================================================================
    %  8. NUMBER-AT-RISK TABLE
    %  ================================================================
    %  Pre-compute number-at-risk at fixed time points for the plotting
    %  script to display beneath the KM curves.

    risk_times = 0:1:tau;

    fprintf('\n--- Number at risk (overall) ---\n');
    fprintf('%-20s', 'Time (yr)');
    fprintf('  %8d', risk_times);
    fprintf('\n');

    at_risk_overall = NaN(3, length(risk_times));
    for ti = 1:3
        mask = data.DevTertile == ti;
        t_i = data.FollowUp_yrs(mask);
        for ri = 1:length(risk_times)
            at_risk_overall(ti, ri) = sum(t_i >= risk_times(ri));
        end
        fprintf('%-20s', tert_labels{ti});
        fprintf('  %8d', at_risk_overall(ti, :));
        fprintf('\n');
    end

    % Sex-stratified at-risk tables
    at_risk_sex = cell(2, 1);
    for si = 1:2
        at_risk_sex{si} = NaN(3, length(risk_times));
        for ti = 1:3
            mask = sex_masks{si} & (data.DevTertile == ti);
            t_i = data.FollowUp_yrs(mask);
            for ri = 1:length(risk_times)
                at_risk_sex{si}(ti, ri) = sum(t_i >= risk_times(ri));
            end
        end
    end

    %% ================================================================
    %  9. SAVE RESULTS
    %  ================================================================

    results.data = data;
    results.inv_e = inv_e;
    results.tertile_bounds = [q33, q67];
    results.tert_labels = tert_labels;
    results.n_total = height(data);
    results.n_deceased = sum(data.Deceased);
    results.median_followup = median(data.FollowUp_yrs);

    % Overall KM
    results.km_overall = km_overall;
    results.logrank_omnibus = struct('chi2', chi2_omni, 'p', p_omni);
    results.logrank_t1_vs_t3 = struct('chi2', chi2_13, 'p', p_13);
    results.logrank_t1_vs_t2 = struct('chi2', chi2_12, 'p', p_12);
    results.logrank_t2_vs_t3 = struct('chi2', chi2_23, 'p', p_23);
    results.rmst_tau = tau;
    results.rmst_diff_t1_t3 = rmst_diff_13;

    % Sex-stratified KM
    results.km_by_sex = km_by_sex;
    results.sex_labels = sex_labels;

    % At-risk tables
    results.risk_times = risk_times;
    results.at_risk_overall = at_risk_overall;
    results.at_risk_sex = at_risk_sex;

    mat_file = fullfile(paths.results, 'survival_curve_results.mat');
    save(mat_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', mat_file);
    fprintf('================================================================\n');
end


%% ========================================================================
%  KAPLAN-MEIER ESTIMATOR
%  ========================================================================

function [t_km, s_km, ci_lo, ci_hi, n_at_risk] = kaplan_meier(time, event)
% KAPLAN_MEIER - Non-parametric survival function with Greenwood CI
%
%   time  - follow-up time (positive real)
%   event - logical: true = event (death), false = censored
%
%   Returns:
%     t_km     - sorted unique event times (with leading 0)
%     s_km     - survival probability at each time
%     ci_lo    - lower 95% Greenwood CI
%     ci_hi    - upper 95% Greenwood CI (capped at 1.0)
%     n_at_risk - number at risk just before each event time

    time = double(time);
    event = logical(event);

    % Sort by time
    [time_sorted, ord] = sort(time);
    event_sorted = event(ord);

    % Unique event times (exclude censoring-only times)
    event_times = unique(time_sorted(event_sorted));

    n_total = length(time);

    % Pre-allocate with t = 0, S = 1
    n_events = length(event_times);
    t_km = [0; event_times];
    s_km = ones(n_events + 1, 1);
    n_at_risk = [n_total; zeros(n_events, 1)];

    % Greenwood variance accumulator
    greenwood_sum = 0;
    var_km = zeros(n_events + 1, 1);

    s_prev = 1;
    for i = 1:n_events
        ti = event_times(i);

        % Number at risk: subjects with time >= ti
        n_i = sum(time_sorted >= ti);
        % Number of events at this time
        d_i = sum(time_sorted == ti & event_sorted);

        n_at_risk(i + 1) = n_i;

        % KM step
        if n_i > 0
            s_prev = s_prev * (1 - d_i / n_i);
        end
        s_km(i + 1) = s_prev;

        % Greenwood variance (cumulative)
        if n_i > d_i && n_i > 0
            greenwood_sum = greenwood_sum + d_i / (n_i * (n_i - d_i));
        end
        var_km(i + 1) = s_km(i + 1)^2 * greenwood_sum;
    end

    % 95% CI using log-log transform (more stable at tails)
    se_km = sqrt(var_km);
    z = 1.96;

    ci_lo = zeros(size(s_km));
    ci_hi = ones(size(s_km));

    for i = 2:length(s_km)
        if s_km(i) > 0 && s_km(i) < 1 && se_km(i) > 0
            log_log_s = log(-log(s_km(i)));
            se_log_log = se_km(i) / (s_km(i) * abs(log(s_km(i))));
            ci_lo(i) = exp(-exp(log_log_s + z * se_log_log));
            ci_hi(i) = exp(-exp(log_log_s - z * se_log_log));
        elseif s_km(i) == 1
            ci_lo(i) = 1;
            ci_hi(i) = 1;
        else
            ci_lo(i) = 0;
            ci_hi(i) = 0;
        end
    end

    ci_lo(1) = 1;
    ci_hi(1) = 1;

    % Ensure CI bounds are valid
    ci_lo = max(ci_lo, 0);
    ci_hi = min(ci_hi, 1);
end


%% ========================================================================
%  LOG-RANK TEST
%  ========================================================================

function [chi2, p] = logrank_test(time, event, group)
% LOGRANK_TEST - Log-rank (Mantel-Cox) test for survival differences
%
%   time  - follow-up time
%   event - logical: true = event
%   group - integer group labels (1, 2, ... k)
%
%   Returns chi-squared statistic and p-value (df = k-1)

    time = double(time);
    event = logical(event);
    group = double(group);

    groups = unique(group);
    k = length(groups);

    % All unique event times
    event_times = sort(unique(time(event)));

    % Accumulate observed - expected for each group
    O = zeros(k, 1);  % observed events
    E = zeros(k, 1);  % expected events under H0
    V = zeros(k, k);  % variance-covariance matrix

    for j = 1:length(event_times)
        tj = event_times(j);

        d_total = sum(time == tj & event);      % total events at tj
        n_total = sum(time >= tj);               % total at risk at tj

        if n_total <= 1 || d_total == 0
            continue;
        end

        for gi = 1:k
            mask_g = (group == groups(gi));
            d_gi = sum(time == tj & event & mask_g);
            n_gi = sum(time >= tj & mask_g);

            O(gi) = O(gi) + d_gi;
            E(gi) = E(gi) + n_gi * d_total / n_total;
        end

        % Variance contribution (hypergeometric)
        if n_total > 1
            for gi = 1:k
                n_gi = sum(time >= tj & (group == groups(gi)));
                for gj = gi:k
                    n_gj = sum(time >= tj & (group == groups(gj)));
                    if gi == gj
                        v_ij = n_gi * (n_total - n_gi) * d_total * ...
                               (n_total - d_total) / (n_total^2 * (n_total - 1));
                    else
                        v_ij = -n_gi * n_gj * d_total * ...
                               (n_total - d_total) / (n_total^2 * (n_total - 1));
                    end
                    V(gi, gj) = V(gi, gj) + v_ij;
                    if gi ~= gj
                        V(gj, gi) = V(gj, gi) + v_ij;
                    end
                end
            end
        end
    end

    % Use first k-1 groups for the test statistic
    U = O(1:k-1) - E(1:k-1);
    V_sub = V(1:k-1, 1:k-1);

    % Regularise if near-singular
    if rcond(V_sub) < 1e-12
        V_sub = V_sub + 1e-10 * eye(k-1);
    end

    chi2 = U' / V_sub * U;
    p = 1 - chi2cdf(chi2, k - 1);
end


%% ========================================================================
%  MEDIAN SURVIVAL
%  ========================================================================

function med = median_survival(t_km, s_km)
% MEDIAN_SURVIVAL - Time at which KM curve first drops to or below 0.5
    idx = find(s_km <= 0.5, 1, 'first');
    if isempty(idx)
        med = NaN;  % not reached
    else
        med = t_km(idx);
    end
end


%% ========================================================================
%  RESTRICTED MEAN SURVIVAL TIME (RMST)
%  ========================================================================

function rmst = compute_rmst(t_km, s_km, tau)
% COMPUTE_RMST - Area under the KM curve up to truncation time tau

    % Truncate KM curve at tau
    t_trunc = t_km(t_km <= tau);
    s_trunc = s_km(1:length(t_trunc));

    % Append the value at tau if needed
    if t_trunc(end) < tau
        % S(tau) is the last known value (step function)
        t_trunc = [t_trunc; tau];
        s_trunc = [s_trunc; s_trunc(end)];
    end

    % Area under step function: sum of S(t_i) * (t_{i+1} - t_i)
    rmst = 0;
    for i = 1:length(t_trunc) - 1
        rmst = rmst + s_trunc(i) * (t_trunc(i + 1) - t_trunc(i));
    end
end
