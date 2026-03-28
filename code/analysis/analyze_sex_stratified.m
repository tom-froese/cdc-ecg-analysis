function results = analyze_sex_stratified()
% ANALYZE_SEX_STRATIFIED - Sex-stratified CDC aging trajectories
%
% Tests whether the core Figure 1 finding — healthy hearts converge on
% 1/e while pathological hearts deviate with age — holds separately in
% females and males.
%
% Analyses:
%   1. Descriptive statistics (N, mean CDC, mean ΔCDC) by Group × Sex
%   2. OLS regression of ΔCDC on Age within each Group × Sex cell
%   3. OLS regression of diastole on Age within each Group × Sex cell
%   4. Sex × Age interaction within each clinical group (tests whether
%      the age slope differs between sexes)
%   5. One-sample t-tests: mean CDC vs 1/e by Group × Sex
%   6. Cohen's d for sex differences in CDC and HR (quantifies the
%      "CDC is sex-invariant, HR is not" finding)
%
% This analysis uses the five-database hierarchical dataset (same as
% Figure 1), not the CODE-15% cohort. The sex-stratified KM survival
% curves (SI Fig 7) already demonstrate sex consistency in the mortality
% analysis; this figure demonstrates it in the aging trajectory analysis.
%
% Data source: hierarchical_results.mat (via analyze_hierarchical_model.m)
% Output:      sex_stratified_results.mat
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    fprintf('================================================================\n');
    fprintf('SEX-STRATIFIED CDC AGING TRAJECTORIES\n');
    fprintf('================================================================\n');
    fprintf('Reference optimum: 1/e = %.4f\n\n', inv_e);

    %% ================================================================
    %  1. LOAD DATA
    %  ================================================================

    S = load(fullfile(paths.results, 'hierarchical_results.mat'), ...
             'all_data', 'lme');
    D = S.all_data;

    D.delta_CDC = D.CDC - inv_e;
    D.RR_ms     = 60000 ./ D.HR;
    D.Dias_ms   = D.RR_ms .* (1 - D.CDC);

    fprintf('Loaded %d subjects from hierarchical dataset\n', height(D));

    %% ================================================================
    %  2. GROUP × SEX DESCRIPTIVE STATISTICS
    %  ================================================================

    groups       = {'HealthyControl', 'ClinicallyNormal', 'Pathological'};
    group_labels = {'Healthy Control', 'Clinically Normal', 'Pathological'};
    sexes        = {'F', 'M'};
    sex_labels   = {'Female', 'Male'};

    fprintf('\n--- Descriptive statistics ---\n');
    fprintf('%-22s %-8s  %6s  %8s  %8s  %8s  %8s\n', ...
        'Group', 'Sex', 'N', 'CDC', 'dCDC', 'HR', 'Dias(ms)');
    fprintf('%s\n', repmat('-', 1, 78));

    desc = struct();

    for gi = 1:3
        for si = 1:2
            mask = (D.Group == groups{gi}) & (D.Sex == sexes{si});
            n = sum(mask);

            desc(gi, si).group     = group_labels{gi};
            desc(gi, si).sex       = sex_labels{si};
            desc(gi, si).n         = n;
            desc(gi, si).mean_cdc  = mean(D.CDC(mask));
            desc(gi, si).std_cdc   = std(D.CDC(mask));
            desc(gi, si).mean_dcdc = mean(D.delta_CDC(mask));
            desc(gi, si).std_dcdc  = std(D.delta_CDC(mask));
            desc(gi, si).mean_hr   = mean(D.HR(mask));
            desc(gi, si).std_hr    = std(D.HR(mask));
            desc(gi, si).mean_dias = mean(D.Dias_ms(mask));

            fprintf('%-22s %-8s  %6d  %8.4f  %+8.4f  %8.1f  %8.1f\n', ...
                group_labels{gi}, sex_labels{si}, n, ...
                desc(gi, si).mean_cdc, desc(gi, si).mean_dcdc, ...
                desc(gi, si).mean_hr, desc(gi, si).mean_dias);
        end
    end

    %% ================================================================
    %  3. OLS REGRESSIONS: ΔCDC ~ Age WITHIN GROUP × SEX
    %  ================================================================

    fprintf('\n--- OLS: dCDC ~ Age within each Group × Sex ---\n');
    fprintf('%-22s %-8s  %12s  %12s  %10s  %10s\n', ...
        'Group', 'Sex', 'Slope(/yr)', 'Intercept', 'R^2', 'p(slope)');
    fprintf('%s\n', repmat('-', 1, 78));

    ols_dcdc = struct();

    for gi = 1:3
        for si = 1:2
            mask = (D.Group == groups{gi}) & (D.Sex == sexes{si});
            mdl = fitlm(D.Age(mask), D.delta_CDC(mask));

            ols_dcdc(gi, si).slope     = mdl.Coefficients.Estimate(2);
            ols_dcdc(gi, si).intercept = mdl.Coefficients.Estimate(1);
            ols_dcdc(gi, si).p_slope   = mdl.Coefficients.pValue(2);
            ols_dcdc(gi, si).r2        = mdl.Rsquared.Ordinary;

            fprintf('%-22s %-8s  %+12.6f  %+12.4f  %10.4f  %10.2e\n', ...
                group_labels{gi}, sex_labels{si}, ...
                ols_dcdc(gi, si).slope, ols_dcdc(gi, si).intercept, ...
                ols_dcdc(gi, si).r2, ols_dcdc(gi, si).p_slope);
        end
    end

    %% ================================================================
    %  4. OLS REGRESSIONS: DIASTOLE ~ Age WITHIN GROUP × SEX
    %  ================================================================

    fprintf('\n--- OLS: Diastole(ms) ~ Age within each Group × Sex ---\n');
    fprintf('%-22s %-8s  %12s  %10s  %10s\n', ...
        'Group', 'Sex', 'Slope(ms/yr)', 'R^2', 'p(slope)');
    fprintf('%s\n', repmat('-', 1, 66));

    ols_dias = struct();

    for gi = 1:3
        for si = 1:2
            mask = (D.Group == groups{gi}) & (D.Sex == sexes{si});
            mdl = fitlm(D.Age(mask), D.Dias_ms(mask));

            ols_dias(gi, si).slope   = mdl.Coefficients.Estimate(2);
            ols_dias(gi, si).p_slope = mdl.Coefficients.pValue(2);
            ols_dias(gi, si).r2      = mdl.Rsquared.Ordinary;

            fprintf('%-22s %-8s  %+12.4f  %10.4f  %10.2e\n', ...
                group_labels{gi}, sex_labels{si}, ...
                ols_dias(gi, si).slope, ols_dias(gi, si).r2, ...
                ols_dias(gi, si).p_slope);
        end
    end

    %% ================================================================
    %  5. SEX × AGE INTERACTION WITHIN EACH GROUP
    %  ================================================================
    %  Tests whether the age slope differs between F and M within each
    %  clinical group. A significant interaction means the sexes age at
    %  different rates for that clinical tier.

    fprintf('\n--- Sex x Age interaction tests (within group) ---\n');

    D.is_male = double(D.Sex == 'M');
    D.Age_c = D.Age - mean(D.Age);

    interaction = struct();

    for gi = 1:3
        d_sub = D(D.Group == groups{gi}, :);

        mdl_int = fitlm(d_sub, 'delta_CDC ~ Age_c * is_male');
        coeff_names = mdl_int.CoefficientNames;
        int_idx = find(contains(coeff_names, ':'));

        if ~isempty(int_idx)
            int_beta = mdl_int.Coefficients.Estimate(int_idx);
            int_p    = mdl_int.Coefficients.pValue(int_idx);
        else
            int_beta = NaN;
            int_p = NaN;
        end

        interaction(gi).group = group_labels{gi};
        interaction(gi).beta  = int_beta;
        interaction(gi).p     = int_p;

        fprintf('  %-22s  Sex x Age beta = %+.6f, p = %.2e', ...
            group_labels{gi}, int_beta, int_p);
        if int_p < 0.05
            fprintf('  *significant*\n');
        else
            fprintf('  (n.s.)\n');
        end
    end

    %% ================================================================
    %  6. ONE-SAMPLE T-TESTS: MEAN CDC VS 1/e BY GROUP × SEX
    %  ================================================================

    fprintf('\n--- One-sample t-tests: CDC vs 1/e ---\n');

    ttest_results = struct();

    for gi = 1:3
        for si = 1:2
            mask = (D.Group == groups{gi}) & (D.Sex == sexes{si});
            [~, p, ~, stats] = ttest(D.CDC(mask), inv_e);

            ttest_results(gi, si).t    = stats.tstat;
            ttest_results(gi, si).df   = stats.df;
            ttest_results(gi, si).p    = p;
            ttest_results(gi, si).mean = mean(D.CDC(mask));

            fprintf('  %-22s %-8s  mean=%.4f, t(%d)=%+.2f, p=%.2e\n', ...
                group_labels{gi}, sex_labels{si}, ...
                ttest_results(gi, si).mean, stats.df, stats.tstat, p);
        end
    end

    %% ================================================================
    %  7. COHEN'S d: SEX DIFFERENCES IN CDC vs HR
    %  ================================================================
    %  Demonstrates "CDC is sex-invariant, HR is not": the ratio is
    %  preserved across sexes while the underlying rate differs.

    fprintf('\n--- Cohen''s d for sex differences (within group) ---\n');
    fprintf('%-22s  %10s  %10s\n', 'Group', 'd(CDC)', 'd(HR)');
    fprintf('%s\n', repmat('-', 1, 46));

    cohens_d_results = struct();

    for gi = 1:3
        mask_f = (D.Group == groups{gi}) & (D.Sex == 'F');
        mask_m = (D.Group == groups{gi}) & (D.Sex == 'M');

        d_cdc = cohens_d_pooled(D.CDC(mask_m), D.CDC(mask_f));
        d_hr  = cohens_d_pooled(D.HR(mask_m),  D.HR(mask_f));

        cohens_d_results(gi).group = group_labels{gi};
        cohens_d_results(gi).d_cdc = d_cdc;
        cohens_d_results(gi).d_hr  = d_hr;

        fprintf('%-22s  %+10.3f  %+10.3f\n', group_labels{gi}, d_cdc, d_hr);
    end

    fprintf('\n  Interpretation: |d(CDC)| << |d(HR)| confirms that the\n');
    fprintf('  ratio is preserved across sexes while the underlying\n');
    fprintf('  heart rate differs — the components shift together.\n');

    %% ================================================================
    %  8. HIERARCHICAL MODEL: SEX-RELATED TERMS
    %  ================================================================
    %  Report the relevant coefficients from the full LME for the
    %  SI legend (already fitted in analyze_hierarchical_model.m).

    fprintf('\n--- From hierarchical LME: sex-related terms ---\n');

    has_lme = isfield(S, 'lme') && ~isempty(S.lme) && isobject(S.lme);
    if has_lme
        lme_coeffs = S.lme.Coefficients;
        fprintf('%-45s %10s %10s %10s\n', 'Term', 'Estimate', 'SE', 'p');
        fprintf('%s\n', repmat('-', 1, 80));

        for i = 1:height(lme_coeffs)
            name = lme_coeffs.Name{i};
            if contains(name, 'Sex', 'IgnoreCase', true) || ...
               contains(name, 'Male', 'IgnoreCase', true)
                fprintf('%-45s %+10.5f %10.5f %10.2e\n', ...
                    name, lme_coeffs.Estimate(i), ...
                    lme_coeffs.SE(i), lme_coeffs.pValue(i));
            end
        end
    else
        fprintf('  LME object not available in hierarchical_results.mat\n');
        fprintf('  (Re-run analyze_hierarchical_model.m with -v7.3 save)\n');
    end

    %% ================================================================
    %  9. SAVE RESULTS
    %  ================================================================

    results.all_data       = D;
    results.inv_e          = inv_e;
    results.groups         = groups;
    results.group_labels   = group_labels;
    results.sexes          = sexes;
    results.sex_labels     = sex_labels;
    results.desc           = desc;
    results.ols_dcdc       = ols_dcdc;
    results.ols_dias       = ols_dias;
    results.interaction    = interaction;
    results.ttest_vs_1e    = ttest_results;
    results.cohens_d       = cohens_d_results;

    mat_file = fullfile(paths.results, 'sex_stratified_results.mat');
    save(mat_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', mat_file);
    fprintf('================================================================\n');
end


%% ========================================================================
%  HELPERS
%  ========================================================================

function d = cohens_d_pooled(x1, x2)
% COHENS_D_POOLED - Pooled-variance Cohen's d
    n1 = length(x1);  n2 = length(x2);
    sp = sqrt(((n1-1)*var(x1) + (n2-1)*var(x2)) / (n1 + n2 - 2));
    d = (mean(x1) - mean(x2)) / sp;
end
