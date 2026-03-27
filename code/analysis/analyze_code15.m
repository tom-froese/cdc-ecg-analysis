function results = analyze_code15()
% ANALYZE_CODE15 - CDC analysis for the CODE-15% dataset
%
% Standalone analysis of the CODE-15% database (annotation_method =
% 'tangent_automatic'). Fits a linear model:
%   CDC ~ Age_c * Group * Sex
%
% where Group is ClinicallyNormal vs Pathological, determined by
% whether all beats for a patient have source_subset = 'normal'.
%
% Subject aggregation uses unique_subject_id (format: CODE15_PatientID)
% for consistency with the hierarchical model pipeline.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);
    n_bootstrap = 5000;

    fprintf('================================================================\n');
    fprintf('CODE-15%%: CLINICALLY NORMAL vs PATHOLOGICAL\n');
    fprintf('================================================================\n');
    fprintf('Reference: 1/e = %.4f\n', inv_e);
    fprintf('================================================================\n\n');

    %% Load data
    fprintf('Loading CODE-15%% dataset...\n');
    beats = readtable(paths.csv_code15, 'TextType', 'string');
    fprintf('  %d beats loaded\n', height(beats));

    %% Apply quality filters
    fprintf('Applying quality filters...\n');
    [valid_mask, ~] = apply_quality_filters(beats);
    beats = beats(valid_mask, :);
    fprintf('  %d beats after filtering\n\n', height(beats));

    %% Compute patient-level CDC using unique_subject_id
    fs = beats.fs(1);
    rt = (beats.t_end_sample - beats.r_sample) / fs * 1000;
    rr = (beats.next_r_sample - beats.r_sample) / fs * 1000;
    beat_ratios = rt ./ rr;

    valid_beats = (rr > 0) & (rt > 0) & (rt < rr);

    [G, patient_ids] = findgroups(beats.unique_subject_id);

    CDC = splitapply(@(x, m) median(x(m)), beat_ratios, valid_beats, G);
    HR  = splitapply(@(x, m) 60000 ./ median(x(m)), rr, valid_beats, G);
    Age = splitapply(@(x) x(1), beats.age, G);
    Sex_str = splitapply(@(x) x(1), beats.sex, G);
    Sex = categorical(Sex_str);

    % Determine group
    is_normal_beat = (beats.source_subset == "normal");
    total_beats_per = splitapply(@numel, beats.unique_subject_id, G);
    normal_beats_per = splitapply(@sum, is_normal_beat, G);
    is_clinically_normal = (total_beats_per == normal_beats_per);

    Group_str = strings(length(patient_ids), 1);
    Group_str(is_clinically_normal) = "ClinicallyNormal";
    Group_str(~is_clinically_normal) = "Pathological";
    Group = categorical(Group_str);

    all_data = table(Age, CDC, HR, Sex, Group);

    % Clean
    valid_subjects = ~isnan(all_data.CDC) & ~isnan(all_data.Age) & (all_data.Age > 0);
    all_data = all_data(valid_subjects, :);
    all_data.Age_c = all_data.Age - mean(all_data.Age);

    %% Summary
    n_cn   = sum(all_data.Group == 'ClinicallyNormal');
    n_path = sum(all_data.Group == 'Pathological');

    fprintf('Valid subjects: %d\n', height(all_data));
    fprintf('  ClinicallyNormal: %d\n', n_cn);
    fprintf('  Pathological:     %d\n', n_path);
    fprintf('  Age range: %.0f-%.0f (Mean: %.1f +/- %.1f)\n', ...
        min(all_data.Age), max(all_data.Age), mean(all_data.Age), std(all_data.Age));

    %% Distribution statistics (for SI figure)
    fprintf('\nDistribution statistics:\n');

    cn_ratios   = all_data.CDC(all_data.Group == 'ClinicallyNormal');
    path_ratios = all_data.CDC(all_data.Group == 'Pathological');

    [mode_cn, ci_cn]     = bootstrap_mode(cn_ratios, n_bootstrap);
    [mode_path, ci_path] = bootstrap_mode(path_ratios, n_bootstrap);

    [p_ranksum, ~] = ranksum(cn_ratios, path_ratios);

    fprintf('  CN mode:   %.4f [%.4f, %.4f]  dCDC=%+.4f\n', ...
            mode_cn, ci_cn(1), ci_cn(2), mode_cn - inv_e);
    fprintf('  Path mode: %.4f [%.4f, %.4f]  dCDC=%+.4f\n', ...
            mode_path, ci_path(1), ci_path(2), mode_path - inv_e);
    fprintf('  Wilcoxon rank-sum: p = %.2e\n', p_ranksum);

    %% Fit linear model
    fprintf('\nFitting: CDC ~ Age_c * Group * Sex\n');
    lm = fitlm(all_data, 'CDC ~ Age_c * Group * Sex');
    disp(lm);

    fprintf('\n=== ANOVA ===\n');
    disp(anova(lm));

    %% Key results
    coeffs = lm.Coefficients;

    fprintf('\n================================================================\n');
    fprintf('KEY RESULTS\n');
    fprintf('================================================================\n\n');

    idx = strcmp(coeffs.Properties.RowNames, '(Intercept)');
    fprintf('Intercept (ClinicallyNormal, Female, mean age): %.4f (SE=%.4f, p=%.2e)\n\n', ...
        coeffs.Estimate(idx), coeffs.SE(idx), coeffs.pValue(idx));

    idx = strcmp(coeffs.Properties.RowNames, 'Group_Pathological');
    fprintf('Pathological vs ClinicallyNormal: beta=%.4f, SE=%.4f, p=%.2e\n\n', ...
        coeffs.Estimate(idx), coeffs.SE(idx), coeffs.pValue(idx));

    idx_age = strcmp(coeffs.Properties.RowNames, 'Age_c');
    fprintf('Age slope (ClinicallyNormal): %.5f/yr (SE=%.5f, p=%.2e)\n', ...
        coeffs.Estimate(idx_age), coeffs.SE(idx_age), coeffs.pValue(idx_age));

    idx_age_p = strcmp(coeffs.Properties.RowNames, 'Age_c:Group_Pathological');
    if any(idx_age_p)
        fprintf('Additional slope (Pathological): %.5f/yr (SE=%.5f, p=%.2e)\n', ...
            coeffs.Estimate(idx_age_p), coeffs.SE(idx_age_p), coeffs.pValue(idx_age_p));
    end

    %% Marginal means
    fprintf('\nObserved means:\n');
    for g = {'ClinicallyNormal', 'Pathological'}
        mask = all_data.Group == g{1};
        fprintf('  %-20s  %.4f +/- %.4f  (N=%d)\n', ...
            g{1}, mean(all_data.CDC(mask)), std(all_data.CDC(mask)), sum(mask));
    end
    fprintf('  1/e reference:       %.4f\n', inv_e);

    %% One-sample t-tests against 1/e
    fprintf('\nOne-sample t-tests vs 1/e:\n');
    for g = {'ClinicallyNormal', 'Pathological'}
        mask = all_data.Group == g{1};
        cdc_vals = all_data.CDC(mask);
        [~, p, ~, tstats] = ttest(cdc_vals, inv_e);
        fprintf('  %-20s  t(%d)=%.3f, p=%.2e\n', g{1}, tstats.df, tstats.tstat, p);
    end

    %% Save
    results.lm = lm;
    results.all_data = all_data;
    results.inv_e = inv_e;

    % Distribution statistics (for plotting)
    results.mode_cn = mode_cn;
    results.ci_cn = ci_cn;
    results.mode_path = mode_path;
    results.ci_path = ci_path;
    results.p_ranksum = p_ranksum;
    results.n_cn = n_cn;
    results.n_path = n_path;

    save(fullfile(paths.results, 'code15_results.mat'), '-struct', 'results');
    fprintf('\nResults saved.\n');
end
