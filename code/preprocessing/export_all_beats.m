function export_all_beats()
% EXPORT_ALL_BEATS - Master export script for all CDC datasets
%
% Usage:
%   export_all_beats()
%
% Runs the beat extraction pipeline for each dataset and produces
% unified CSV files in data/preprocessed/. Each CSV contains one row
% per beat with the following columns:
%
%   record_id         - Subject/record identifier
%   source_database   - Database name (e.g. 'LUDB', 'QTDB', 'Fantasia')
%   source_subset     - Diagnostic category or sub-population
%   group             - 'healthy', 'pathological', or 'sudden_death'
%   age               - Age in years (NaN if unavailable)
%   sex               - 'M', 'F', or '' if unavailable
%   fs                - Sampling frequency (Hz)
%   recording_id      - Recording identifier within subject
%   beat_idx          - Beat index within recording (1, 2, 3, ...)
%   r_sample          - R-peak sample number
%   t_end_sample      - T-wave end sample number
%   next_r_sample     - Next R-peak sample number
%   annotation_method - Provenance of fiducial points (see below)
%
% Annotation method values:
%   'manual'            - Both R-peak and T-end from expert annotations
%   'manual_t_auto_r'   - T-end from expert annotations, R-peak algorithmic
%   'manual_r_auto_t'   - R-peak from database annotations, T-end automatic
%   'ecgdeli_automatic' - Both from ECGDeli validated algorithm
%   'tangent_automatic' - R-peak detected, T-end via tangent method
%
% NO beat-level quality filters are applied during export. All filtering
% is deferred to the analysis stage for transparency and reproducibility.
%
% Tom Froese, OIST Embodied Cognitive Science Unit
% Created: January-February 2026

    paths = config();

    fprintf('============================================================\n');
    fprintf('CDC BEAT EXPORT: ALL DATASETS\n');
    fprintf('============================================================\n');
    fprintf('Output directory: %s\n\n', paths.preprocessed);

    total_beats = 0;
    total_records = 0;

    %% Dataset 1: LUDB (Tier 1 - fully manual)
    if exist(paths.raw_ludb, 'dir')
        fprintf('────────────────────────────────────────────────────────────\n');
        fprintf('EXPORTING LUDB\n');
        fprintf('────────────────────────────────────────────────────────────\n');
        [nb, nr] = export_ludb(paths.raw_ludb, paths.csv_ludb);
        total_beats = total_beats + nb;
        total_records = total_records + nr;
    else
        warning('LUDB not found: %s', paths.raw_ludb);
    end

    %% Dataset 2: QTDB (Tier 1 - fully manual)
    if exist(paths.raw_qtdb, 'dir')
        fprintf('\n────────────────────────────────────────────────────────────\n');
        fprintf('EXPORTING QTDB\n');
        fprintf('────────────────────────────────────────────────────────────\n');
        [nb, nr] = export_qtdb(paths.raw_qtdb, paths.csv_qtdb);
        total_beats = total_beats + nb;
        total_records = total_records + nr;
    else
        warning('QTDB not found: %s', paths.raw_qtdb);
    end

    %% Dataset 3: PTB Diagnostic (Tier 2 - T-end manual, R-peak algorithmic)
    if exist(paths.raw_ptb, 'dir')
        fprintf('\n────────────────────────────────────────────────────────────\n');
        fprintf('EXPORTING PTB\n');
        fprintf('────────────────────────────────────────────────────────────\n');
        [nb, nr] = export_ptb(paths.raw_ptb, paths.csv_ptb);
        total_beats = total_beats + nb;
        total_records = total_records + nr;
    else
        warning('PTB not found: %s', paths.raw_ptb);
    end

    %% Dataset 4: PTB-XL (Tier 2 - ECGDeli automatic)
    if exist(paths.raw_ptbxl, 'dir') && exist(paths.raw_ptbxl_features, 'dir')
        fprintf('\n────────────────────────────────────────────────────────────\n');
        fprintf('EXPORTING PTB-XL\n');
        fprintf('────────────────────────────────────────────────────────────\n');
        [nb, nr] = export_ptbxl(paths.raw_ptbxl, paths.raw_ptbxl_features, paths.csv_ptbxl);
        total_beats = total_beats + nb;
        total_records = total_records + nr;
    else
        warning('PTB-XL or features dataset not found.');
    end

    %% Dataset 5: Fantasia (Tier 3 - R from database, T-end tangent)
    if exist(paths.raw_fantasia, 'dir')
        fprintf('\n────────────────────────────────────────────────────────────\n');
        fprintf('EXPORTING FANTASIA\n');
        fprintf('────────────────────────────────────────────────────────────\n');
        [nb, nr] = export_fantasia(paths.raw_fantasia, paths.csv_fantasia);
        total_beats = total_beats + nb;
        total_records = total_records + nr;
    else
        warning('Fantasia not found: %s', paths.raw_fantasia);
    end

    %% Dataset 6: Autonomic Aging (Tier 3 - fully automatic)
    if exist(paths.raw_autonomic_aging, 'dir')
        fprintf('\n────────────────────────────────────────────────────────────\n');
        fprintf('EXPORTING AUTONOMIC AGING\n');
        fprintf('────────────────────────────────────────────────────────────\n');
        [nb, nr] = export_autonomic_aging(paths.raw_autonomic_aging, paths.csv_autonomic_aging);
        total_beats = total_beats + nb;
        total_records = total_records + nr;
    else
        warning('Autonomic Aging not found: %s', paths.raw_autonomic_aging);
    end

    %% Dataset 7: CODE-15% (Tier 3 - fully automatic)
    if exist(paths.raw_code15, 'dir')
        fprintf('\n────────────────────────────────────────────────────────────\n');
        fprintf('EXPORTING CODE-15%%\n');
        fprintf('────────────────────────────────────────────────────────────\n');
        [nb, nr] = export_code15(paths.raw_code15, paths.csv_code15);
        total_beats = total_beats + nb;
        total_records = total_records + nr;
    else
        warning('CODE-15%% not found: %s', paths.raw_code15);
    end

    %% Summary
    fprintf('\n============================================================\n');
    fprintf('EXPORT COMPLETE\n');
    fprintf('============================================================\n');
    fprintf('Total: %d beats from %d records across all datasets\n', ...
            total_beats, total_records);
    fprintf('============================================================\n');
end