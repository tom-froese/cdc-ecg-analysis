function results = analyze_code15()
% ANALYZE_CODE15 - CDC analysis for the CODE-15% dataset
%
% Standalone analysis of the CODE-15% database (annotation_method =
% 'tangent_automatic'). Fits a linear model:
%   CDC ~ Age_c * Group * Sex
%
% where Group has two levels (internal categorical → display label):
%   NonPathECG → "Patients (non-pathological ECG)"
%   PathECG    → "Patients (pathological ECG)"
% Group assignment: a patient is NonPathECG iff all their beats have
% source_subset = 'normal'; otherwise PathECG.
%
% Subject aggregation uses unique_subject_id (format: CODE15_PatientID)
% for consistency with the hierarchical model pipeline.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);
    n_bootstrap = 5000;

    fprintf('================================================================\n');
    fprintf('CODE-15%%: PATIENTS (NON-PATHOLOGICAL ECG) vs PATIENTS (PATHOLOGICAL ECG)\n');
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
    is_npe = (total_beats_per == normal_beats_per);

    Group_str = strings(length(patient_ids), 1);
    Group_str(is_npe) = "NonPathECG";
    Group_str(~is_npe) = "PathECG";
    Group = categorical(Group_str);

    all_data = table(Age, CDC, HR, Sex, Group);

    % Clean
    valid_subjects = ~isnan(all_data.CDC) & ~isnan(all_data.Age) & (all_data.Age > 0);
    all_data = all_data(valid_subjects, :);
    all_data.Age_c = all_data.Age - mean(all_data.Age);

    %% Summary
    n_npe   = sum(all_data.Group == 'NonPathECG');
    n_pe = sum(all_data.Group == 'PathECG');

    fprintf('Valid subjects: %d\n', height(all_data));
    fprintf('  Patients (non-pathological ECG): %d\n', n_npe);
    fprintf('  Patients (pathological ECG):     %d\n', n_pe);
    fprintf('  Age range: %.0f-%.0f (Mean: %.1f +/- %.1f)\n', ...
        min(all_data.Age), max(all_data.Age), mean(all_data.Age), std(all_data.Age));

    %% Distribution statistics (for SI figure)
    fprintf('\nDistribution statistics:\n');

    npe_ratios   = all_data.CDC(all_data.Group == 'NonPathECG');
    pe_ratios = all_data.CDC(all_data.Group == 'PathECG');

    [mode_npe, ci_npe]     = bootstrap_mode(npe_ratios, n_bootstrap);
    [mode_pe, ci_pe] = bootstrap_mode(pe_ratios, n_bootstrap);

    [p_ranksum, ~] = ranksum(npe_ratios, pe_ratios);

    fprintf('  NonPathECG mode: %.4f [%.4f, %.4f]  dCDC=%+.4f\n', ...
            mode_npe, ci_npe(1), ci_npe(2), mode_npe - inv_e);
    fprintf('  PathECG mode:    %.4f [%.4f, %.4f]  dCDC=%+.4f\n', ...
            mode_pe, ci_pe(1), ci_pe(2), mode_pe - inv_e);
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
    fprintf('Intercept (Patients (non-pathological ECG), Female, mean age): %.4f (SE=%.4f, p=%.2e)\n\n', ...
        coeffs.Estimate(idx), coeffs.SE(idx), coeffs.pValue(idx));

    idx = strcmp(coeffs.Properties.RowNames, 'Group_PathECG');
    fprintf('Patients (pathological ECG) vs Patients (non-pathological ECG): beta=%.4f, SE=%.4f, p=%.2e\n\n', ...
        coeffs.Estimate(idx), coeffs.SE(idx), coeffs.pValue(idx));

    idx_age = strcmp(coeffs.Properties.RowNames, 'Age_c');
    fprintf('Age slope (Patients (non-pathological ECG)): %.5f/yr (SE=%.5f, p=%.2e)\n', ...
        coeffs.Estimate(idx_age), coeffs.SE(idx_age), coeffs.pValue(idx_age));

    idx_age_pe = strcmp(coeffs.Properties.RowNames, 'Age_c:Group_PathECG');
    if any(idx_age_pe)
        fprintf('Additional slope (Patients (pathological ECG)): %.5f/yr (SE=%.5f, p=%.2e)\n', ...
            coeffs.Estimate(idx_age_pe), coeffs.SE(idx_age_pe), coeffs.pValue(idx_age_pe));
    end

    %% Marginal means
    group_ids    = {'NonPathECG', 'PathECG'};
    group_labels = {'Patients (non-pathological ECG)', 'Patients (pathological ECG)'};
    fprintf('\nObserved means:\n');
    for gi = 1:numel(group_ids)
        mask = all_data.Group == group_ids{gi};
        fprintf('  %-32s  %.4f +/- %.4f  (N=%d)\n', ...
            group_labels{gi}, mean(all_data.CDC(mask)), std(all_data.CDC(mask)), sum(mask));
    end
    fprintf('  1/e reference:                    %.4f\n', inv_e);

    %% One-sample t-tests against 1/e
    fprintf('\nOne-sample t-tests vs 1/e:\n');
    for gi = 1:numel(group_ids)
        mask = all_data.Group == group_ids{gi};
        cdc_vals = all_data.CDC(mask);
        [~, p, ~, tstats] = ttest(cdc_vals, inv_e);
        fprintf('  %-32s  t(%d)=%.3f, p=%.2e\n', group_labels{gi}, tstats.df, tstats.tstat, p);
    end

    %% Save
    results.lm = lm;
    results.all_data = all_data;
    results.inv_e = inv_e;

    % Distribution statistics (for plotting)
    results.mode_npe = mode_npe;
    results.ci_npe = ci_npe;
    results.mode_pe = mode_pe;
    results.ci_pe = ci_pe;
    results.p_ranksum = p_ranksum;
    results.n_npe = n_npe;
    results.n_pe = n_pe;

    save(fullfile(paths.results, 'code15_results.mat'), '-struct', 'results');
    fprintf('\nResults saved.\n');
end
