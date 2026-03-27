% REPRODUCE_FIGURES - One-command reproduction of ALL figures and tables

fprintf('============================================================\n');
fprintf('CDC-ECG ANALYSIS — REPRODUCING NATURE AGING PAPER\n');
fprintf('============================================================\n\n');

paths = config();

% 1. Export beat data (creates all *_beats.csv in data/preprocessed/)
fprintf('1. Exporting beat data from all databases...\n');
export_all_beats();

% 2. Run all analyses (these generate the figures and tables)
fprintf('2. Running gold-standard analysis (Act 1)...\n');
analyze_gold_standard();

fprintf('3. Running large-scale analysis (Act 2)...\n');
analyze_large_scale();

fprintf('4. Running hierarchical model (Act 3)...\n');
analyze_hierarchical_model();

fprintf('5. Running CODE-15%% standalone analysis...\n');
analyze_code15();

fprintf('\n============================================================\n');
fprintf('REPRODUCTION COMPLETE!\n');
fprintf('============================================================\n');
fprintf('→ All figures are in: results/figures/\n');
fprintf('→ All tables   are in: results/tables/\n\n');
fprintf('You can now open any figure directly in MATLAB.\n');