function paths = config()
% CONFIG - Central path configuration for the CDC Analysis project
%
% Usage:
%   paths = config();
%
% Returns a struct with all project-relative paths.
% Set paths.root to the absolute path of your project folder,
% then all other paths are derived automatically.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    %% Project root (EDIT THIS to match your local setup)
    paths.root = fileparts(mfilename('fullpath'));

    %% Raw data directories
    paths.raw = fullfile(paths.root, 'data', 'raw');

    paths.raw_ludb     = fullfile(paths.raw, 'lobachevsky-university-electrocardiography-database-1.0.1');
    paths.raw_qtdb     = fullfile(paths.raw, 'qt-database-1.0.0');
    paths.raw_ptb      = fullfile(paths.raw, 'ptb-diagnostic-ecg-database-1.0.0');
    paths.raw_ptbxl    = fullfile(paths.raw, 'ptb-xl-a-large-publicly-available-electrocardiography-dataset-1.0.3');
    paths.raw_ptbxl_features = fullfile(paths.raw, 'ptb-xl-a-comprehensive-electrocardiographic-feature-dataset-1.0.1');
    paths.raw_fantasia = fullfile(paths.raw, 'fantasia-database-1.0.0');
    paths.raw_autonomic_aging = fullfile(paths.raw, ...
        'autonomic-aging-a-dataset-to-quantify-changes-of-cardiovascular-autonomic-function-during-healthy-aging-1.0.0');
    paths.raw_code15   = fullfile(paths.raw, 'code-15-1.0.0');

    %% Preprocessed data (output CSVs)
    paths.preprocessed = fullfile(paths.root, 'data', 'preprocessed');

    paths.csv_ludb     = fullfile(paths.preprocessed, 'ludb_beats.csv');
    paths.csv_qtdb     = fullfile(paths.preprocessed, 'qtdb_beats.csv');
    paths.csv_ptb      = fullfile(paths.preprocessed, 'ptb_beats.csv');
    paths.csv_ptbxl    = fullfile(paths.preprocessed, 'ptbxl_beats.csv');
    paths.csv_fantasia = fullfile(paths.preprocessed, 'fantasia_beats.csv');
    paths.csv_autonomic_aging = fullfile(paths.preprocessed, 'autonomic_aging_beats.csv');
    paths.csv_code15   = fullfile(paths.preprocessed, 'code15_beats.csv');
    paths.csv_code15_exams = fullfile(paths.preprocessed, 'code15_exams.csv');

    %% Results
    paths.results = fullfile(paths.root, 'results');
    paths.figures = fullfile(paths.results, 'figures');

    %% Ensure output directories exist
    dirs_to_create = {paths.preprocessed, paths.figures};
    for i = 1:length(dirs_to_create)
        if ~exist(dirs_to_create{i}, 'dir')
            mkdir(dirs_to_create{i});
        end
    end

    %% Add utility functions to path
    addpath(fullfile(paths.root, 'code', 'utils'));
end
