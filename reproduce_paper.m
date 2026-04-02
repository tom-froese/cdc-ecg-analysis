% REPRODUCE_PAPER - One-command reproduction of all analyses, figures, and tables
%
% Reproduces every result in:
%   Froese, T. et al. (in prep.). Nature Aging.
%
% Prerequisites:
%   1. Run setup.m once to add all project paths.
%   2. Place the preprocessed CSV files (*_beats.csv, code15_exams.csv)
%      in data/preprocessed/. Download from Zenodo:
%      https://doi.org/10.5281/zenodo.19246123
%
% Usage:
%   >> setup
%   >> reproduce_paper
%
% Output:
%   results/*.mat              — analysis results
%   results/figures/*.pdf/png  — all main and supplementary figures
%
% Runtime: ~10–30 minutes depending on hardware (the CODE-15% analyses
% are the bottleneck due to n = 214,176 patients).
%
% Tom Froese, OIST Embodied Cognitive Science Unit

fprintf('\n');
fprintf('============================================================\n');
fprintf('CDC-ECG ANALYSIS — REPRODUCING NATURE AGING PAPER\n');
fprintf('============================================================\n\n');

paths = config();

%% ---- Verify preprocessed data is available ----

required_csvs = { ...
    paths.csv_ludb, paths.csv_qtdb, paths.csv_ptb, ...
    paths.csv_ptbxl, paths.csv_fantasia, ...
    paths.csv_autonomic_aging, paths.csv_code15, ...
    paths.csv_code15_exams};

missing = {};
for i = 1:length(required_csvs)
    if ~isfile(required_csvs{i})
        missing{end+1} = required_csvs{i}; %#ok<AGROW>
    end
end

if ~isempty(missing)
    fprintf('ERROR: Missing preprocessed CSV files:\n');
    for i = 1:length(missing)
        fprintf('  %s\n', missing{i});
    end
    fprintf('\nDownload from Zenodo and place in: %s\n', paths.preprocessed);
    fprintf('https://doi.org/10.5281/zenodo.19246123\n');
    error('Cannot reproduce paper without preprocessed data.');
end

fprintf('All preprocessed CSV files found.\n\n');

%% ---- STEP 1: Run all analyses ----
% Each analysis loads its required CSVs and saves a .mat results file.
% Ordering respects data dependencies (plot scripts load these .mat files).

fprintf('============================================================\n');
fprintf('STEP 1: STATISTICAL ANALYSES\n');
fprintf('============================================================\n\n');

fprintf('[1/7] Gold-standard analysis (LUDB, QTDB — manual annotations)...\n');
analyze_gold_standard();

fprintf('\n[2/7] Large-scale analysis (Fantasia, Autonomic Aging, PTB, PTB-XL)...\n');
analyze_large_scale();

fprintf('\n[3/7] Hierarchical linear model (all five databases pooled)...\n');
analyze_hierarchical_model();

fprintf('\n[4/7] CDC vs heart rate independence analysis...\n');
analyze_cdc_vs_hr();

fprintf('\n[5/7] CODE-15%% standalone analysis...\n');
analyze_code15();

fprintf('\n[6/7] Age-stratified mortality analysis (CODE-15%%)...\n');
analyze_age_stratified_mortality();

fprintf('\n[7/7] Survival curves and sex-stratified analyses (CODE-15%%)...\n');
analyze_survival_curves();
analyze_sex_stratified();

%% ---- STEP 2: Generate all figures ----

fprintf('\n============================================================\n');
fprintf('STEP 2: FIGURES\n');
fprintf('============================================================\n\n');

% Main figures
fprintf('Plotting Fig. 1  — CDC aging trajectories...\n');
plot_Fig1();

fprintf('Plotting Fig. 2  — Mortality by CDC-deviation tertile...\n');
plot_Fig2();

% Supplementary figures
fprintf('Plotting SI Fig. 1 — CDC vs heart rate...\n');
plot_SI_Fig1();

fprintf('Plotting SI Fig. 2 — Systole ceiling and RR compensation...\n');
plot_SI_Fig2();

% SI Fig. 3 requires raw data — see note below
fprintf('Skipping SI Fig. 3 — pipeline validation (requires raw databases).\n');
fprintf('  To generate, download LUDB and QTDB from PhysioNet into data/raw/,\n');
fprintf('  then run: analyze_gold_standard_validation(); plot_SI_Fig3();\n');

fprintf('Plotting SI Fig. 4 — Gold-standard CDC distributions...\n');
plot_SI_Fig4();

fprintf('Plotting SI Fig. 5 — Large-scale CDC distributions...\n');
plot_SI_Fig5();

fprintf('Plotting SI Fig. 6 — CODE-15%% CDC distributions...\n');
plot_SI_Fig6();

fprintf('Plotting SI Fig. 7 — Kaplan-Meier survival curves...\n');
plot_SI_Fig7();

fprintf('Plotting SI Fig. 8 — Sex-stratified CDC aging...\n');
plot_SI_Fig8();

%% ---- Done ----

fprintf('\n============================================================\n');
fprintf('REPRODUCTION COMPLETE\n');
fprintf('============================================================\n');
fprintf('Results : %s\n', paths.results);
fprintf('Figures : %s\n', paths.figures);
fprintf('\nNote: SI Fig. 3 (pipeline validation) requires raw PhysioNet\n');
fprintf('databases. See instructions above to generate it separately.\n\n');
