function results = analyze_hierarchical_model()
% ANALYZE_HIERARCHICAL_MODEL - Pooled hierarchical LMM across all datasets
%
% Fits a linear mixed-effects model to subject-level median CDC values:
%
%   CDC ~ Age_c * Group * Sex + (1|Dataset)
%
% where Group has three levels:
%   HealthyControl    - Verified healthy volunteers (LUDB, Fantasia, Autonomic Aging)
%   ClinicallyNormal  - Hospital patients with normal ECG findings (PTB, PTB-XL)
%   Pathological      - Hospital patients with cardiac pathology
%
% Dataset is a random intercept absorbing annotation-method variance and
% other between-database differences.
%
% Subject aggregation uses unique_subject_id (format: Database_RecordID)
% to ensure unambiguous subject-level grouping across datasets with
% heterogeneous ID formats.
%
% Also extracts subject-level median heart rate (HR) for exploratory
% analysis of the "gait adjustment" hypothesis: the healthy aging heart
% may preserve its optimal CDC ratio by slowing its rate.
%
% Quality filters are applied uniformly to algorithmically annotated
% datasets; manually annotated datasets (LUDB, QTDB) are used as-is.
% QTDB is excluded from the hierarchical model because it lacks the
% 'healthy' group needed for the three-way comparison.
%
% A supplementary analysis adds annotation_method as a fixed effect
% to verify that the core findings are not driven by detection method.
%
% Dependencies:
%   config.m, apply_quality_filters.m, build_beats_table.m
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    fprintf('================================================================\n');
    fprintf('HIERARCHICAL MODEL: POOLED ACROSS ALL DATASETS\n');
    fprintf('================================================================\n');
    fprintf('Reference: 1/e = %.4f\n', inv_e);
    fprintf('================================================================\n\n');

    %% ====================================================================
    %  Step 1: Load and process each dataset to subject level
    %  ====================================================================
    fprintf('Loading and processing datasets...\n');

    ludb_data = process_dataset(paths.csv_ludb, 'LUDB', ...
        'ApplyFilters', false, 'GroupMap', @map_ludb_group);

    ptb_data = process_dataset(paths.csv_ptb, 'PTB', ...
        'ApplyFilters', true, 'GroupMap', @map_clinical_group);

    ptbxl_data = process_dataset(paths.csv_ptbxl, 'PTB-XL', ...
        'ApplyFilters', true, 'GroupMap', @map_clinical_group, ...
        'ClampAge', {300, 90});

    fantasia_data = process_dataset(paths.csv_fantasia, 'Fantasia', ...
        'ApplyFilters', true, 'GroupMap', @map_healthy_only);

    aa_data = process_dataset(paths.csv_autonomic_aging, 'AutonomicAging', ...
        'ApplyFilters', true, 'GroupMap', @map_healthy_only);

    %% ====================================================================
    %  Step 2: Combine into single table and prepare model variables
    %  ====================================================================
    all_data = [ludb_data; ptb_data; ptbxl_data; fantasia_data; aa_data];

    % Standardise Group to 3 levels (merge Fantasia age subgroups)
    all_data.Group = categorical(all_data.Group);
    all_data.Group = mergecats(all_data.Group, ...
        intersect(categories(all_data.Group), ...
                  {'HealthyControl', 'Young Healthy', 'Elderly Healthy'}), ...
        'HealthyControl');

    % Mean-centre age for interpretable intercept
    all_data.Age_c = all_data.Age - mean(all_data.Age);

    fprintf('\nTotal subjects: %d\n', height(all_data));
    fprintf('   HealthyControl:   %d\n', sum(all_data.Group == 'HealthyControl'));
    fprintf('   ClinicallyNormal: %d\n', sum(all_data.Group == 'ClinicallyNormal'));
    fprintf('   Pathological:     %d\n\n', sum(all_data.Group == 'Pathological'));

    % Per-dataset breakdown
    fprintf('Per-dataset breakdown:\n');
    datasets = categories(all_data.Dataset);
    for d = 1:length(datasets)
        ds = datasets{d};
        mask = all_data.Dataset == ds;
        fprintf('  %-16s  N=%5d  (HC=%d, CN=%d, Path=%d)  Age: %.1f+/-%.1f  Male: %.1f%%\n', ...
            ds, sum(mask), ...
            sum(mask & all_data.Group == 'HealthyControl'), ...
            sum(mask & all_data.Group == 'ClinicallyNormal'), ...
            sum(mask & all_data.Group == 'Pathological'), ...
            mean(all_data.Age(mask)), std(all_data.Age(mask)), ...
            100 * mean(all_data.Sex(mask) == 'M'));
    end

    %% ====================================================================
    %  Step 3: Primary model
    %  ====================================================================
    fprintf('\n================================================================\n');
    fprintf('PRIMARY MODEL: CDC ~ Age_c * Group * Sex + (1|Dataset)\n');
    fprintf('================================================================\n');

    lme = fitlme(all_data, 'CDC ~ Age_c * Group * Sex + (1|Dataset)');
    disp(lme);

    fprintf('\n=== TYPE III ANOVA ===\n');
    anova_tbl = anova(lme); %#ok<NASGU>
    disp(anova(lme));

    %% Extract and report key effects
    [~, ~, stats] = fixedEffects(lme);

    idx_intercept = find(strcmp(stats.Name, '(Intercept)'));
    idx_cn   = find(strcmp(stats.Name, 'Group_ClinicallyNormal'));
    idx_path = find(strcmp(stats.Name, 'Group_Pathological'));
    idx_age  = find(strcmp(stats.Name, 'Age_c'));
    idx_sex  = find(strcmp(stats.Name, 'Sex_M'));

    fprintf('\n================================================================\n');
    fprintf('KEY RESULTS\n');
    fprintf('================================================================\n\n');

    fprintf('Intercept (HealthyControl, Female, mean age): %.4f (SE=%.4f, p=%.2e)\n\n', ...
        stats.Estimate(idx_intercept), stats.SE(idx_intercept), stats.pValue(idx_intercept));

    fprintf('Group effects (relative to HealthyControl):\n');
    fprintf('  ClinicallyNormal: beta=%.4f, SE=%.4f, p=%.2e\n', ...
        stats.Estimate(idx_cn), stats.SE(idx_cn), stats.pValue(idx_cn));
    fprintf('  Pathological:     beta=%.4f, SE=%.4f, p=%.2e\n\n', ...
        stats.Estimate(idx_path), stats.SE(idx_path), stats.pValue(idx_path));

    fprintf('Age slope (HealthyControl): %.5f/yr (SE=%.5f, p=%.2e)\n', ...
        stats.Estimate(idx_age), stats.SE(idx_age), stats.pValue(idx_age));
    fprintf('Sex (Male vs Female):       %.4f   (SE=%.4f, p=%.2e)\n\n', ...
        stats.Estimate(idx_sex), stats.SE(idx_sex), stats.pValue(idx_sex));

    %% ====================================================================
    %  Step 4: Observed means and one-sample t-tests against 1/e
    %  ====================================================================
    fprintf('================================================================\n');
    fprintf('OBSERVED MEANS PER GROUP (at mean age)\n');
    fprintf('================================================================\n');

    groups = {'HealthyControl', 'ClinicallyNormal', 'Pathological'};
    for g = 1:length(groups)
        mask = all_data.Group == groups{g};
        fprintf('  %-20s  %.4f +/- %.4f  (N=%d)\n', ...
            groups{g}, mean(all_data.CDC(mask)), std(all_data.CDC(mask)), sum(mask));
    end
    fprintf('  1/e reference:       %.4f\n\n', inv_e);

    fprintf('One-sample t-tests: CDC vs 1/e\n');
    for g = 1:length(groups)
        mask = all_data.Group == groups{g};
        cdc_vals = all_data.CDC(mask);
        [~, p, ci, tstats] = ttest(cdc_vals, inv_e);
        fprintf('  %-20s  t(%d)=%.3f, p=%.2e, diff=%.4f [%.4f, %.4f]\n', ...
            groups{g}, tstats.df, tstats.tstat, p, ...
            mean(cdc_vals) - inv_e, ci(1) - inv_e, ci(2) - inv_e);
    end

    %% Random effects
    fprintf('\nRandom effects (Dataset intercept):\n');
    [~, ~, re_stats] = randomEffects(lme);
    disp(re_stats);

    %% Model fit
    fprintf('Model fit: AIC=%.1f, BIC=%.1f, LogL=%.1f, Residual SD=%.5f\n', ...
        lme.ModelCriterion.AIC, lme.ModelCriterion.BIC, ...
        lme.LogLikelihood, sqrt(lme.MSE));

    %% ====================================================================
    %  Step 5: Supplementary — annotation method robustness check
    %  ====================================================================
    fprintf('\n================================================================\n');
    fprintf('SUPPLEMENTARY: ANNOTATION METHOD ROBUSTNESS CHECK\n');
    fprintf('================================================================\n');

    all_data.AnnotTier = categorical(all_data.AnnotTier);
    lme_annot = fitlme(all_data, 'CDC ~ Age_c * Group * Sex + AnnotTier + (1|Dataset)');

    fprintf('Model with AnnotTier as fixed effect:\n');
    [~, ~, stats_annot] = fixedEffects(lme_annot);

    annot_rows = contains(stats_annot.Name, 'AnnotTier');
    if any(annot_rows)
        fprintf('  Annotation tier effects:\n');
        for j = find(annot_rows)'
            fprintf('    %s: beta=%.4f, SE=%.4f, p=%.2e\n', ...
                stats_annot.Name{j}, stats_annot.Estimate(j), ...
                stats_annot.SE(j), stats_annot.pValue(j));
        end
    end

    fprintf('\n  AIC without AnnotTier: %.1f\n', lme.ModelCriterion.AIC);
    fprintf('  AIC with AnnotTier:    %.1f\n', lme_annot.ModelCriterion.AIC);
    fprintf('  Delta AIC: %.1f (positive = primary model preferred)\n', ...
            lme_annot.ModelCriterion.AIC - lme.ModelCriterion.AIC);

    %% ====================================================================
    %  Step 6: Exploratory — heart rate vs age by clinical group
    %  ====================================================================
    fprintf('\n================================================================\n');
    fprintf('EXPLORATORY: HEART RATE vs AGE BY CLINICAL GROUP\n');
    fprintf('================================================================\n');

    fprintf('\nHR summary per group:\n');
    for g = 1:length(groups)
        mask = all_data.Group == groups{g};
        hr = all_data.HR(mask);
        fprintf('  %-20s  %.1f +/- %.1f bpm  [%.1f - %.1f]  (N=%d)\n', ...
            groups{g}, mean(hr), std(hr), min(hr), max(hr), sum(mask));
    end

    fprintf('\nAge-HR correlations per group:\n');
    for g = 1:length(groups)
        mask = all_data.Group == groups{g};
        [r, p] = corr(all_data.Age(mask), all_data.HR(mask));
        fprintf('  %-20s  r = %.3f, p = %.2e\n', groups{g}, r, p);
    end

    %% HR hierarchical model (supplementary)
    fprintf('\n================================================================\n');
    fprintf('SUPPLEMENTARY: HR ~ Age_c * Group * Sex + (1|Dataset)\n');
    fprintf('================================================================\n');

    lme_hr = fitlme(all_data, 'HR ~ Age_c * Group * Sex + (1|Dataset)');
    disp(lme_hr);

    fprintf('\n=== TYPE III ANOVA (HR model) ===\n');
    anova_hr = anova(lme_hr); %#ok<NASGU>
    disp(anova(lme_hr));

    %% ====================================================================
    %  Step 7: Save results
    %  ====================================================================
    results.lme       = lme;
    results.all_data  = all_data;
    results.inv_e     = inv_e;
    results.lme_annot = lme_annot;
    results.lme_hr    = lme_hr;

    % Save as v7.3 (HDF5-based) to preserve MATLAB objects and tables.
    % The default v5 format silently drops LinearMixedModel objects and
    % tables with categorical columns.
    mat_file = fullfile(paths.results, 'hierarchical_results.mat');
    save(mat_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', mat_file);

    % Also export the subject-level data as CSV for portability.
    % This allows collaborators using R or Python to inspect the pooled
    % dataset without requiring MATLAB.
    csv_file = fullfile(paths.results, 'hierarchical_subject_data.csv');
    export_data = all_data;
    export_data.Group   = string(export_data.Group);
    export_data.Sex     = string(export_data.Sex);
    export_data.Dataset = string(export_data.Dataset);
    export_data.AnnotTier = string(export_data.AnnotTier);
    writetable(export_data, csv_file);
    fprintf('Subject-level data saved to: %s\n', csv_file);
end


%% ========================================================================
%  UNIFIED DATASET PROCESSING
%  ========================================================================

function data = process_dataset(filename, dataset_label, varargin)
% PROCESS_DATASET - Load a beat-level CSV and aggregate to subject level
%
% Reads the unified CDC beat table, optionally applies quality filters,
% computes per-subject median CDC and HR, and maps diagnostic groups via
% a caller-supplied function handle.
%
% Subject-level aggregation uses the unique_subject_id column, which
% disambiguates subjects across databases with overlapping record_id
% formats (e.g., LUDB's numeric '1', '2' vs PTB-XL's 'patient15709').
%
% Inputs:
%   filename      - Path to beat-level CSV file
%   dataset_label - Label for the Dataset column (e.g., 'PTB-XL')
%
% Optional Name-Value pairs:
%   'ApplyFilters' - Logical, whether to apply quality filters (default: false)
%   'GroupMap'      - Function handle: group_label = f(group_string)
%                     Maps the raw 'group' column value to a model-level
%                     group label (e.g., 'HealthyControl', 'Pathological')
%   'ClampAge'      - {sentinel_value, replacement_value} for age clamping
%                     (e.g., {300, 90} for PTB-XL's coded age cap)
%
% Output:
%   data - Table with columns: Age, CDC, HR, Sex, Group, Dataset, AnnotTier

    p = inputParser;
    addParameter(p, 'ApplyFilters', false, @islogical);
    addParameter(p, 'GroupMap', @map_passthrough, @(x) isa(x, 'function_handle'));
    addParameter(p, 'ClampAge', {}, @iscell);
    parse(p, varargin{:});

    apply_filt = p.Results.ApplyFilters;
    group_map  = p.Results.GroupMap;
    clamp_age  = p.Results.ClampAge;

    % Load beats
    beats = readtable(filename);

    % Age clamping (e.g., PTB-XL encodes 90+ as 300)
    if ~isempty(clamp_age)
        sentinel    = clamp_age{1};
        replacement = clamp_age{2};
        beats.age(beats.age == sentinel) = replacement;
    end

    % Quality filters
    if apply_filt
        [mask, ~] = apply_quality_filters(beats, 'Verbose', false);
        beats = beats(mask, :);
    end

    % Aggregate to subject level using unique_subject_id
    [unique_subjects, ~, subject_idx] = unique(beats.unique_subject_id);
    n = length(unique_subjects);

    Age       = zeros(n, 1);
    CDC       = zeros(n, 1);
    HR        = zeros(n, 1);
    Sex       = cell(n, 1);
    Group     = cell(n, 1);
    AnnotTier = cell(n, 1);

    for i = 1:n
        mask = (subject_idx == i);
        rec  = beats(mask, :);

        % Compute RT and RR intervals in milliseconds
        rt_ms = (rec.t_end_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        rr_ms = (rec.next_r_sample - rec.r_sample) ./ rec.fs(1) * 1000;
        valid = (rr_ms > 0) & (rt_ms > 0) & (rt_ms < rr_ms);

        if sum(valid) >= 1
            CDC(i) = median(rt_ms(valid) ./ rr_ms(valid));
            HR(i)  = 60000 / median(rr_ms(valid));
        else
            CDC(i) = NaN;
            HR(i)  = NaN;
        end

        Age(i)       = rec.age(1);
        Sex{i}       = get_sex(rec);
        AnnotTier{i} = get_string_field(rec, 'annotation_method');
        Group{i}     = group_map(get_string_field(rec, 'group'));
    end

    % Exclude subjects with missing data
    valid = ~isnan(CDC) & ~isnan(HR) & ~isnan(Age) & (Age > 0);

    data = table(Age(valid), CDC(valid), HR(valid), ...
                 categorical(Sex(valid)), categorical(Group(valid)), ...
                 categorical(repmat({dataset_label}, sum(valid), 1)), ...
                 categorical(AnnotTier(valid)), ...
                 'VariableNames', {'Age','CDC','HR','Sex','Group','Dataset','AnnotTier'});

    fprintf('  %-16s  %5d subjects (%d excluded)\n', ...
            dataset_label, height(data), n - sum(valid));
end


%% ========================================================================
%  GROUP MAPPING FUNCTIONS
%  ========================================================================

function label = map_ludb_group(raw_group)
% LUDB healthy controls are verified volunteers
    if strcmpi(raw_group, 'healthy')
        label = 'HealthyControl';
    else
        label = 'Pathological';
    end
end

function label = map_clinical_group(raw_group)
% PTB and PTB-XL 'healthy' means clinically normal ECG, not volunteer status
    if strcmpi(raw_group, 'healthy')
        label = 'ClinicallyNormal';
    else
        label = 'Pathological';
    end
end

function label = map_healthy_only(~)
% Fantasia and Autonomic Aging: all subjects are healthy volunteers
    label = 'HealthyControl';
end

function label = map_passthrough(raw_group)
% Default: pass through the raw group label unchanged
    label = raw_group;
end


%% ========================================================================
%  FIELD EXTRACTION HELPERS
%  ========================================================================

function s = get_sex(rec)
% GET_SEX - Extract sex string from first row, handling cell/string types
    if iscell(rec.sex) && ~isempty(rec.sex)
        s = rec.sex{1};
    elseif isstring(rec.sex)
        s = char(rec.sex(1));
    else
        s = '';
    end
end

function s = get_string_field(rec, field_name)
% GET_STRING_FIELD - Robustly extract a char value from the first row of a
% table column that may be stored as cell, string, or categorical
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