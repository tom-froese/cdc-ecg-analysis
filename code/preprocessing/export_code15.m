function [n_beats_out, n_records_out] = export_code15(base_path, output_file)
% EXPORT_CODE15 - Export CODE-15% beat data in unified CDC format
%
% Annotation tier: TANGENT_AUTOMATIC
%   R-peaks are detected using the Pan-Tompkins algorithm.
%   T-wave endpoints are detected using the tangent method.
%   (Same pipeline as Autonomic Aging export.)
%
% The CODE-15% dataset is a large-scale clinical ECG database from the
% Telehealth Network of Minas Gerais (TNMG), Brazil, containing 345,779
% 12-lead ECG exams from 233,770 patients sampled at 400 Hz.
%
% Data format:
%   exams.csv          - Metadata (exam_id, patient_id, age, is_male,
%                        normal_ecg, diagnostic columns, trace_file)
%   exams_part{i}.hdf5 - HDF5 files with datasets 'tracings' (N x 4096 x 12)
%                        and 'exam_id' (N x 1). Signals are in 1e-4 V units.
%                        Leads: DI, DII, DIII, AVR, AVL, AVF, V1-V6.
%                        Short recordings are zero-padded symmetrically.
%
% Patient-level inclusion criteria:
%   Many patients have multiple exams across the 2010-2016 collection
%   period. To ensure one independent data point per patient:
%
%   1. EXAM SELECTION: Only the first exam per patient (lowest exam_id)
%      is processed. This provides a baseline cross-sectional snapshot
%      and preserves the mortality follow-up data (death, timey), which
%      is only available for each patient's first exam.
%
%   2. GROUP ASSIGNMENT: The normal_ecg flag of the selected (first)
%      exam determines group membership at the time of measurement:
%        - normal_ecg = true  -> source_subset = 'normal',  group = 'healthy'
%        - normal_ecg = false -> source_subset = 'abnormal', group = 'pathological'
%      Some patients classified as normal at first exam may develop
%      pathology in later exams; this is standard for cross-sectional
%      cohort designs and makes the group comparison conservative.
%
%   These criteria are applied before any HDF5 files are opened, so
%   only ~233k exams are processed instead of ~346k.
%
% Lead selection: DII (index 2) is used for R-peak and T-end detection,
% as it typically provides the clearest P-QRS-T morphology.
%
% Zero-padding handling: The dataset pads short recordings (7 s = 2800
% samples) symmetrically to 4096. We detect and strip padding so that
% the Pan-Tompkins detector does not produce spurious detections in the
% flat zero regions.
%
% IMPORTANT: No beat-level quality filters are applied during export.
% All beats where detect_t_end returns a valid (non-NaN) result are
% exported. Quality filtering is deferred to the analysis stage.
%
% Memory strategy: HDF5 tracings are read in batches (configurable) to
% avoid loading the entire ~4 GB file at once.
%
% Usage:
%   [n_beats, n_records] = export_code15(base_path, output_file)
%
%   Or via the project pipeline:
%   paths = config();
%   export_code15(paths.raw_code15, paths.csv_code15);
%
% Tom Froese, OIST Embodied Cognitive Science Unit
% Created: March 2026

    if ~exist(base_path, 'dir')
        error('CODE-15%% database not found: %s', base_path);
    end

    %% ====================================================================
    %  CONFIGURATION
    %  ====================================================================
    FS            = 400;    % Sampling frequency (fixed for CODE-15%)
    N_SAMPLES     = 4096;   % Samples per tracing (after zero-padding)
    LEAD_DII      = 2;      % DII lead index (1-based)
    SIGNAL_SCALE  = 1000;   % Multiply raw signal by this to get mV
    BATCH_SIZE    = 500;    % Exams per HDF5 read batch
    PROGRESS_EVERY = 1000;  % Print progress every N exams

    %% ====================================================================
    %  LOAD METADATA
    %  ====================================================================
    csv_file = fullfile(base_path, 'exams.csv');
    if ~exist(csv_file, 'file')
        error('exams.csv not found in %s', base_path);
    end

    fprintf('Loading CODE-15%% metadata...\n');
    meta = readtable(csv_file, 'TextType', 'string', 'TreatAsMissing', '', ...
                     'VariableNamingRule', 'preserve');

    % Validate required columns
    required_cols = {'exam_id', 'age', 'is_male', 'patient_id', ...
                     'normal_ecg', 'trace_file'};
    for i = 1:length(required_cols)
        if ~any(strcmp(meta.Properties.VariableNames, required_cols{i}))
            error('Required column ''%s'' not found in exams.csv', required_cols{i});
        end
    end

    n_exams_total = height(meta);
    n_patients_total = length(unique(meta.patient_id));
    fprintf('  %d exams from %d patients in metadata\n', ...
            n_exams_total, n_patients_total);

    % Standardise is_male to logical (may be string 'True'/'False' or numeric)
    if isnumeric(meta.is_male) || islogical(meta.is_male)
        meta.is_male_logical = logical(meta.is_male);
    else
        meta.is_male_logical = strcmpi(meta.is_male, "True") | ...
                               meta.is_male == "1";
    end

    %% ====================================================================
    %  SELECT FIRST EXAM PER PATIENT
    %  ====================================================================
    fprintf('\nSelecting first exam per patient...\n');

    % Sort by patient_id then exam_id to identify chronological first exam
    meta = sortrows(meta, {'patient_id', 'exam_id'}, {'ascend', 'ascend'});
    [~, first_idx] = unique(meta.patient_id, 'first');
    selected_meta = meta(first_idx, :);

    fprintf('  %d unique patients -> %d first exams selected\n', ...
            length(first_idx), height(selected_meta));

    %% ====================================================================
    %  GROUP ASSIGNMENT (from first exam's own normal_ecg flag)
    %  ====================================================================
    fprintf('Assigning groups from first-exam diagnosis...\n');

    % Parse normal_ecg flag robustly (handles numeric, logical, string)
    if isnumeric(selected_meta.normal_ecg) || islogical(selected_meta.normal_ecg)
        is_normal = (selected_meta.normal_ecg == 1);
    else
        is_normal = strcmpi(selected_meta.normal_ecg, "True") | ...
                    strcmpi(selected_meta.normal_ecg, "1");
    end

    selected_meta.group = repmat("pathological", height(selected_meta), 1);
    selected_meta.source_subset = repmat("abnormal", height(selected_meta), 1);
    selected_meta.group(is_normal) = "healthy";
    selected_meta.source_subset(is_normal) = "normal";

    n_normal = sum(is_normal);
    n_pathological = sum(~is_normal);
    fprintf('  ClinicallyNormal (first exam): %d (%.1f%%)\n', ...
            n_normal, 100 * n_normal / height(selected_meta));
    fprintf('  Pathological (first exam):     %d (%.1f%%)\n', ...
            n_pathological, 100 * n_pathological / height(selected_meta));

    % Skip exams with missing critical demographics
    valid_demo = ~isnan(selected_meta.age) & ~ismissing(selected_meta.patient_id);
    n_skipped_demo = sum(~valid_demo);
    selected_meta = selected_meta(valid_demo, :);
    if n_skipped_demo > 0
        fprintf('  Skipped (missing demographics): %d\n', n_skipped_demo);
    end

    %% ====================================================================
    %  COHORT SUMMARY
    %  ====================================================================
    fprintf('\n  --- Selected cohort ---\n');
    for g = ["healthy", "pathological"]
        mask = selected_meta.group == g;
        n_g = sum(mask);
        n_m = sum(mask & selected_meta.is_male_logical);
        n_f = sum(mask & ~selected_meta.is_male_logical);
        ages = selected_meta.age(mask);
        fprintf('  %-14s  N=%6d  (M=%d, F=%d)  Age: %.1f +/- %.1f [%d-%d]\n', ...
                g, n_g, n_m, n_f, mean(ages, 'omitnan'), std(ages, 'omitnan'), ...
                min(ages), max(ages));
    end
    fprintf('  --------------------------------\n');

    % Build lookup set of selected exam_ids for fast membership testing
    selected_exam_map = containers.Map('KeyType', 'int64', 'ValueType', 'int32');
    for j = 1:height(selected_meta)
        selected_exam_map(int64(selected_meta.exam_id(j))) = int32(j);
    end

    %% ====================================================================
    %  PROCESS SELECTED EXAMS FROM HDF5 FILES
    %  ====================================================================
    hdf5_files = unique(selected_meta.trace_file);
    hdf5_files = hdf5_files(~ismissing(hdf5_files));
    fprintf('\n  %d HDF5 part files to read\n', length(hdf5_files));

    rows = {};
    n_records_processed = 0;
    n_records_skipped = 0;
    n_beats_total = 0;
    n_failed_tend = 0;
    n_exams_processed = 0;
    n_exams_no_rpeaks = 0;
    n_missing_hdf5 = 0;
    n_exams_to_process = height(selected_meta);

    t_start = tic;

    for fi = 1:length(hdf5_files)
        hdf5_name = char(hdf5_files(fi));
        hdf5_path = fullfile(base_path, hdf5_name);

        if ~exist(hdf5_path, 'file')
            fprintf('  WARNING: %s not found — skipping\n', hdf5_name);
            n_missing_hdf5 = n_missing_hdf5 + 1;
            continue;
        end

        fprintf('\n  === %s [%d/%d] ===\n', hdf5_name, fi, length(hdf5_files));

        % Read exam IDs from HDF5
        try
            hdf5_exam_ids = h5read(hdf5_path, '/exam_id');
        catch ME
            fprintf('  ERROR reading exam_id from %s: %s\n', hdf5_name, ME.message);
            continue;
        end

        n_in_file = length(hdf5_exam_ids);

        % Pre-compute which indices in this file are selected
        is_selected = false(n_in_file, 1);
        meta_row_for = zeros(n_in_file, 1);
        for idx = 1:n_in_file
            eid = int64(hdf5_exam_ids(idx));
            if selected_exam_map.isKey(eid)
                is_selected(idx) = true;
                meta_row_for(idx) = selected_exam_map(eid);
            end
        end

        n_selected_in_file = sum(is_selected);
        fprintf('  %d/%d tracings selected\n', n_selected_in_file, n_in_file);

        if n_selected_in_file == 0
            continue;
        end

        % Process in batches
        n_batches = ceil(n_in_file / BATCH_SIZE);

        for batch = 1:n_batches
            idx_start = (batch - 1) * BATCH_SIZE + 1;
            idx_end = min(batch * BATCH_SIZE, n_in_file);
            batch_size_actual = idx_end - idx_start + 1;

            % Skip batch if no selected exams in range
            if ~any(is_selected(idx_start:idx_end))
                continue;
            end

            % Read batch of tracings
            % HDF5 stores as (N, 4096, 12) in C order; MATLAB reads
            % this as (12, 4096, N) due to column-major convention.
            tracings = h5read(hdf5_path, '/tracings', ...
                              [1, 1, idx_start], ...
                              [12, N_SAMPLES, batch_size_actual]);

            for k = 1:batch_size_actual
                linear_idx = idx_start + k - 1;

                if ~is_selected(linear_idx)
                    continue;
                end

                mi = meta_row_for(linear_idx);
                exam_id       = selected_meta.exam_id(mi);
                patient_id    = selected_meta.patient_id(mi);
                age           = selected_meta.age(mi);
                is_male       = selected_meta.is_male_logical(mi);
                group         = char(selected_meta.group(mi));
                source_subset = char(selected_meta.source_subset(mi));

                if is_male, sex_str = 'M';
                else,       sex_str = 'F';
                end

                %% Extract DII lead and convert to mV
                % tracings layout: (12 leads, 4096 samples, batch)
                ecg_raw = double(squeeze(tracings(LEAD_DII, :, k)));
                ecg_mv = ecg_raw * SIGNAL_SCALE;

                %% Strip zero-padding
                ecg_mv = strip_zero_padding(ecg_mv);

                if length(ecg_mv) < FS
                    n_records_skipped = n_records_skipped + 1;
                    continue;
                end

                %% Filter
                ecg_filt = bandpass_filter(ecg_mv, FS, 0.5, 40);
                ecg_smooth = lowpass_filter(ecg_filt, FS, 15);

                %% R-peak detection (Pan-Tompkins + refinement)
                r_peaks = detect_r_peaks(ecg_filt, FS);

                search_hw = round(0.04 * FS);
                for rp = 1:length(r_peaks)
                    s0 = max(1, r_peaks(rp) - search_hw);
                    s1 = min(length(ecg_filt), r_peaks(rp) + search_hw);
                    [~, mx_idx] = max(ecg_filt(s0:s1));
                    r_peaks(rp) = s0 + mx_idx - 1;
                end
                r_peaks = unique(r_peaks);

                if length(r_peaks) < 2
                    n_exams_no_rpeaks = n_exams_no_rpeaks + 1;
                    n_records_skipped = n_records_skipped + 1;
                    continue;
                end

                %% Detect T-ends — no filtering at export
                n_candidate = length(r_peaks) - 1;
                n_beats_record = 0;
                n_failed_record = 0;

                record_id_str = char(string(patient_id));
                recording_id_str = char(string(exam_id));

                for b = 1:n_candidate
                    [t_end, ~, ~, ~] = detect_t_end(ecg_smooth, r_peaks, b, FS);

                    if isnan(t_end)
                        n_failed_record = n_failed_record + 1;
                        continue;
                    end

                    n_beats_record = n_beats_record + 1;
                    rows{end+1} = {record_id_str, 'CODE15', source_subset, ...
                                   group, age, sex_str, FS, recording_id_str, ...
                                   n_beats_record, r_peaks(b), t_end, ...
                                   r_peaks(b+1), 'tangent_automatic'}; %#ok<AGROW>
                end

                n_beats_total = n_beats_total + n_beats_record;
                n_failed_tend = n_failed_tend + n_failed_record;
                n_records_processed = n_records_processed + 1;
                n_exams_processed = n_exams_processed + 1;

                if mod(n_exams_processed, PROGRESS_EVERY) == 0
                    elapsed = toc(t_start);
                    rate = n_exams_processed / elapsed;
                    remaining = (n_exams_to_process - n_exams_processed) / rate;
                    fprintf('    %d/%d exams (%.0f/min, ~%.0f min left, %d beats)\n', ...
                            n_exams_processed, n_exams_to_process, ...
                            rate * 60, remaining / 60, n_beats_total);
                end
            end

            clear tracings;
        end
    end

    %% ====================================================================
    %  SUMMARY AND OUTPUT
    %  ====================================================================
    elapsed_total = toc(t_start);

    fprintf('\n  ============================================================\n');
    fprintf('  CODE-15%% EXPORT SUMMARY\n');
    fprintf('  ============================================================\n');
    fprintf('  Inclusion: first exam per patient, same-exam classification\n');
    fprintf('  Patients in metadata:     %d\n', n_patients_total);
    fprintf('  Exams processed:          %d\n', n_records_processed);
    fprintf('  Exams skipped:            %d\n', n_records_skipped);
    fprintf('    (no R-peaks detected):  %d\n', n_exams_no_rpeaks);
    fprintf('  Missing HDF5 files:       %d\n', n_missing_hdf5);
    fprintf('  Total beats exported:     %d\n', n_beats_total);
    fprintf('  T-end detection failures: %d\n', n_failed_tend);
    if (n_beats_total + n_failed_tend) > 0
        fprintf('  T-end detection rate:     %.1f%%\n', ...
                100 * n_beats_total / (n_beats_total + n_failed_tend));
    end
    fprintf('  Elapsed time:             %.1f min\n', elapsed_total / 60);
    fprintf('  ============================================================\n');

    %% Build table and save
    if isempty(rows)
        warning('No beats were exported. Check that HDF5 files are present.');
        n_beats_out = 0;
        n_records_out = 0;
        return;
    end

    T = build_beats_table(rows);
    writetable(T, output_file);
    fprintf('  Saved: %s (%d rows)\n', output_file, height(T));

    %% Copy exams.csv to preprocessed directory for downstream analyses
    %  The mortality annotations (death, timey) and diagnostic flags live
    %  in exams.csv. Copying it alongside the beats CSV keeps the
    %  preprocessed directory self-contained for analysis scripts.
    [out_dir, ~, ~] = fileparts(output_file);
    exams_dst = fullfile(out_dir, 'code15_exams.csv');
    if ~strcmp(csv_file, exams_dst)
        [status, msg] = copyfile(csv_file, exams_dst);
        if status
            fprintf('  Copied: exams.csv -> %s\n', exams_dst);
        else
            warning('export_code15:copy_failed', 'Failed to copy exams.csv: %s', msg);
        end
    end    

    exams_dst = strrep(output_file, 'code15_beats.csv', 'code15_exams.csv');
    if ~strcmp(csv_file, exams_dst)
        copyfile(csv_file, exams_dst);
        fprintf('  Copied: exams.csv -> %s\n', exams_dst);
    end

    n_beats_out = n_beats_total;
    n_records_out = n_records_processed;
end

%% ========================================================================
%  CODE-15%-SPECIFIC HELPERS
%  ========================================================================

function ecg_valid = strip_zero_padding(ecg)
% STRIP_ZERO_PADDING - Remove symmetric zero-padding from CODE-15% tracings
%
% The CODE-15% dataset zero-pads short recordings (7 s / 2800 samples)
% symmetrically to reach 4096 samples. This function finds the first and
% last non-zero samples to extract the valid signal region.

    TOLERANCE = 1e-6;  % mV — well below any real ECG amplitude

    abs_ecg = abs(ecg);
    first_valid = find(abs_ecg > TOLERANCE, 1, 'first');
    last_valid  = find(abs_ecg > TOLERANCE, 1, 'last');

    if isempty(first_valid) || isempty(last_valid)
        ecg_valid = ecg;
        return;
    end

    % Small margin to avoid clipping edge values
    first_valid = max(1, first_valid - 5);
    last_valid  = min(length(ecg), last_valid + 5);

    ecg_valid = ecg(first_valid:last_valid);
end
