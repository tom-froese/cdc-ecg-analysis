function [n_beats_out, n_records_out] = export_ptb(base_path, output_file)
% EXPORT_PTB - Export PTB Diagnostic ECG beat data in unified CDC format
%
% Annotation tier: HYBRID (manual_t_auto_r)
%   T-wave endpoints are from expert clinical annotations published as
%   supplementary material to Bousseljot et al. (2006). Five independent
%   referees annotated Q-onset and T-end for each recording; the median
%   values (in milliseconds from recording start) are used here.
%   R-peaks are detected algorithmically using a Pan-Tompkins detector,
%   since the supplementary annotations provide only QRS onset/offset
%   timing in milliseconds rather than sample-level R-peak locations.
%
% Usage:
%   [n_beats, n_records] = export_ptb(base_path, output_file)
%
% Inputs:
%   base_path   - Root folder of the PTB Diagnostic ECG Database.
%                 Must contain patient subdirectories (patient001, ...) and
%                 the supplementary annotation file:
%                     12938_2006_174_MOESM1_ESM.doc
%                 This file is available from:
%                     Bousseljot et al., "Nutzung der EKG-Signaldatenbank
%                     CARDIODAT der PTB ueber das Internet",
%                     Biomedizinische Technik, 2006, Supplementary Table 1.
%                     DOI: 10.1007/s10010-006-0174-8
%
%   output_file - Path for the output CSV (unified CDC beat format).
%
% Outputs:
%   n_beats_out   - Number of beats exported.
%   n_records_out - Number of patients with at least one valid beat.
%
% Dependencies:
%   detect_r_peaks.m, build_beats_table.m (project functions)
%
% No beat-level quality filters are applied; all valid beats are exported.

    %% ====================================================================
    %  Step 1: Locate and parse the supplementary annotation document
    %  ====================================================================

    doc_file = fullfile(base_path, '12938_2006_174_MOESM1_ESM.doc');

    if ~exist(doc_file, 'file')
        error(['Annotation file not found:\n  %s\n\n' ...
               'Please download the supplementary material (.doc) from:\n' ...
               '  Bousseljot et al. (2006), DOI: 10.1007/s10010-006-0174-8\n' ...
               'and place it in the PTB database root folder.'], doc_file);
    end

    fprintf('  Parsing annotations from %s...\n', doc_file);
    T_ann = parse_ptb_annotation_doc(doc_file);
    fprintf('  Found %d annotated recordings from %d patients\n', ...
            height(T_ann), length(unique(T_ann.patient_num)));

    %% ====================================================================
    %  Step 2: Iterate over annotated recordings and export beats
    %  ====================================================================

    unique_patients = unique(T_ann.patient_num);
    n_patients = length(unique_patients);

    rows = {};
    n_patients_processed = 0;
    n_beats_total = 0;

    for p = 1:n_patients
        patient_num = unique_patients(p);
        patient_folder = sprintf('patient%03d', patient_num);
        patient_path = fullfile(base_path, patient_folder);

        if ~exist(patient_path, 'dir')
            continue;
        end

        % Get all annotated recordings for this patient
        patient_mask = (T_ann.patient_num == patient_num);
        patient_records = T_ann(patient_mask, :);
        n_recordings = height(patient_records);

        patient_had_valid_beat = false;

        for rec_idx = 1:n_recordings
            record_name_full = patient_records.record_name{rec_idx};
            record_name_base = patient_records.record_name_base{rec_idx};
            q_onset_ms = patient_records.q_onset_median(rec_idx);
            t_end_ms   = patient_records.t_end_median(rec_idx);

            % Skip recordings flagged as invalid by the referees
            if q_onset_ms == 0 || t_end_ms == 0
                continue;
            end

            % Locate the actual file on disk (handles naming variants)
            record_name = find_record_file(patient_path, ...
                                           record_name_full, record_name_base);
            if isempty(record_name)
                continue;
            end

            % Build file paths
            hea_file = fullfile(patient_path, [record_name '.hea']);
            dat_file = fullfile(patient_path, [record_name '.dat']);
            xyz_file = fullfile(patient_path, [record_name '.xyz']);

            % Read header and ECG signals
            try
                [fs, n_samples, signals] = read_ptb_header(hea_file);
                ecg_data = read_ptb_signals(dat_file, xyz_file, ...
                                            n_samples, signals);
            catch
                continue;
            end

            % Select Lead II for R-peak detection
            lead_ii_idx = find(strcmpi({signals.leads.name}, 'ii'));
            if isempty(lead_ii_idx)
                lead_ii_idx = 2;  % Fallback: second channel is typically Lead II
            end
            lead_ii = ecg_data(:, lead_ii_idx);

            % Detect R-peaks algorithmically
            r_peaks_all = detect_r_peaks(lead_ii, fs);

            % Convert manual annotations from milliseconds to samples
            q_onset_sample = round(q_onset_ms * fs / 1000);
            t_end_sample   = round(t_end_ms   * fs / 1000);

            % Find the first R-peak that falls within the annotated window
            % (between Q-onset and T-end of the annotated beat)
            valid_r_idx = find(r_peaks_all > q_onset_sample & ...
                               r_peaks_all < t_end_sample);
            if isempty(valid_r_idx)
                continue;
            end

            r_peak_idx    = valid_r_idx(1);
            r_peak_sample = r_peaks_all(r_peak_idx);

            % Require a subsequent R-peak for RR interval calculation
            if r_peak_idx >= length(r_peaks_all)
                continue;
            end
            next_r_sample = r_peaks_all(r_peak_idx + 1);

            % Extract patient demographics from the WFDB header
            age = parse_age(signals.patient_info.age);
            sex = parse_sex(signals.patient_info.sex);
            diagnosis = signals.patient_info.diagnosis;

            % Classify as healthy or pathological
            if is_healthy_diagnosis(diagnosis)
                group = 'healthy';
            else
                group = 'pathological';
            end

            % Preserve the original diagnosis as the source subset label
            if isempty(diagnosis) || isempty(strtrim(diagnosis))
                source_subset = 'Unknown';
            else
                source_subset = diagnosis;
            end

            % Append beat to output
            rows{end+1} = {patient_folder, 'PTBDB', source_subset, group, ...
                           age, sex, fs, record_name, 1, ...
                           r_peak_sample, t_end_sample, next_r_sample, ...
                           'manual_t_auto_r'}; %#ok<AGROW>

            n_beats_total = n_beats_total + 1;
            patient_had_valid_beat = true;
        end

        if patient_had_valid_beat
            n_patients_processed = n_patients_processed + 1;
        end

        if mod(p, 50) == 0
            fprintf('    Processed %d/%d patients (%d beats)\n', ...
                    p, n_patients, n_beats_total);
        end
    end

    fprintf('  Processed %d patients, %d beats\n', ...
            n_patients_processed, n_beats_total);

    %% ====================================================================
    %  Step 3: Build unified table and write to CSV
    %  ====================================================================

    T = build_beats_table(rows);
    writetable(T, output_file);
    fprintf('  Saved: %s\n', output_file);

    n_beats_out   = n_beats_total;
    n_records_out = n_patients_processed;
end


%% ========================================================================
%  ANNOTATION DOCUMENT PARSER
%  ========================================================================

function T = parse_ptb_annotation_doc(doc_file)
% PARSE_PTB_ANNOTATION_DOC - Extract fiducial point annotations from the
% supplementary .doc file published with the PTB Diagnostic ECG Database.
%
% The file is a legacy Microsoft Word binary (.doc) containing a structured
% table of expert annotations. The data section has the following layout,
% repeated for each recording:
%
%   patient<N>          (patient identifier, 1-indexed, not zero-padded)
%   s####[l|a|b]re      (recording identifier, with optional variant suffix)
%   <12 integer values>:
%       Q-onset: R1, R2, R3, R4, R5, Median   (6 values, in ms)
%       T-end:   R1, R2, R3, R4, R5, Median   (6 values, in ms)
%
% Multiple recordings per patient are listed consecutively under the same
% patient header. One recording (patient285/s0544_re) has all-zero values,
% indicating no valid ECG tracings were observed.
%
% Extraction method: The .doc binary is read as raw bytes and printable
% ASCII characters are retained. This avoids any dependency on external
% tools (antiword, LibreOffice, python-docx) while reliably recovering
% the structured numeric data embedded in the document. The approach works
% because the annotation data consists entirely of ASCII text (patient
% identifiers, record names, and integer values).
%
% Input:
%   doc_file - Path to 12938_2006_174_MOESM1_ESM.doc
%
% Output:
%   T - Table with columns:
%       patient_num       - Patient number (integer)
%       patient_folder    - Zero-padded folder name (e.g., 'patient001')
%       record_name       - Full record name with suffix (e.g., 's0014l_re')
%       record_name_base  - Base record name without suffix (e.g., 's0014_re')
%       q_onset_median    - Median Q-onset from 5 referees (ms)
%       t_end_median      - Median T-end from 5 referees (ms)

    % Read binary content and extract printable ASCII
    fid = fopen(doc_file, 'rb');
    if fid == -1
        error('Cannot open annotation file: %s', doc_file);
    end
    raw_bytes = fread(fid, inf, 'uint8');
    fclose(fid);

    ascii_mask = (raw_bytes >= 32 & raw_bytes < 127) | ...
                  raw_bytes == 10 | raw_bytes == 13 | raw_bytes == 9;
    text_chars = char(raw_bytes');
    text_chars(~ascii_mask') = ' ';
    tokens = strsplit(strtrim(text_chars));

    % Parse tokens sequentially: patient headers, record names, then values
    n_tokens = length(tokens);
    current_patient = [];
    records = {};

    i = 1;
    while i <= n_tokens
        token = tokens{i};

        % Match patient header: 'patient1', 'patient2', ..., 'patient294'
        patient_match = regexp(token, '^patient(\d+)$', 'tokens');
        if ~isempty(patient_match)
            current_patient = str2double(patient_match{1}{1});
            i = i + 1;
            continue;
        end

        % Match record name: s####re, s####lre, s####are, s####bre
        % (with or without underscores in the token)
        token_clean = strrep(token, '_', '');
        record_match = regexp(token_clean, '^(s\d{4})([lab]?)re$', 'tokens');

        if ~isempty(record_match) && ~isempty(current_patient)
            record_base   = record_match{1}{1};  % e.g., 's0014'
            record_suffix = record_match{1}{2};   % '', 'l', 'a', or 'b'

            % Construct standardised record names
            if isempty(record_suffix)
                record_name_full = [record_base '_re'];
            else
                record_name_full = [record_base record_suffix '_re'];
            end
            record_name_base = [record_base '_re'];

            % Collect the next 12 numeric tokens (6 Q-onset + 6 T-end)
            values = [];
            j = i + 1;
            while j <= n_tokens && length(values) < 12
                next_token = tokens{j};

                % Stop if we encounter a new patient or record header
                if ~isempty(regexp(next_token, '^patient\d+$', 'once'))
                    break;
                end
                next_clean = strrep(next_token, '_', '');
                if ~isempty(regexp(next_clean, '^s\d{4}[lab]?re$', 'once'))
                    break;
                end

                % Parse integer values
                val = str2double(next_token);
                if ~isnan(val) && val == round(val)
                    values(end+1) = val; %#ok<AGROW>
                end
                j = j + 1;
            end

            if length(values) == 12
                records{end+1} = struct( ...
                    'patient_num',      current_patient, ...
                    'patient_folder',   sprintf('patient%03d', current_patient), ...
                    'record_name',      record_name_full, ...
                    'record_name_base', record_name_base, ...
                    'q_onset_median',   values(6), ...
                    't_end_median',     values(12)); %#ok<AGROW>
            else
                warning('export_ptb:incomplete_annotation', ...
                        'Incomplete annotation for patient%03d/%s (found %d/12 values)', ...
                        current_patient, record_name_full, length(values));
            end

            i = j;
            continue;
        end

        i = i + 1;
    end

    % Convert struct array to table
    n = length(records);
    if n == 0
        error('export_ptb:no_annotations', ...
              'No annotations could be parsed from %s', doc_file);
    end

    patient_num      = zeros(n, 1);
    patient_folder   = cell(n, 1);
    record_name      = cell(n, 1);
    record_name_base = cell(n, 1);
    q_onset_median   = zeros(n, 1);
    t_end_median     = zeros(n, 1);

    for k = 1:n
        patient_num(k)      = records{k}.patient_num;
        patient_folder{k}   = records{k}.patient_folder;
        record_name{k}      = records{k}.record_name;
        record_name_base{k} = records{k}.record_name_base;
        q_onset_median(k)   = records{k}.q_onset_median;
        t_end_median(k)     = records{k}.t_end_median;
    end

    T = table(patient_num, patient_folder, record_name, record_name_base, ...
              q_onset_median, t_end_median);

    fprintf('  Parsed %d recordings (%d with valid annotations)\n', ...
            n, sum(q_onset_median > 0 & t_end_median > 0));
end


%% ========================================================================
%  PTB-SPECIFIC HELPERS
%  ========================================================================

function record_name = find_record_file(patient_folder, record_name_full, record_name_base)
% FIND_RECORD_FILE - Locate a PTB recording file on disk
%
% The annotation document uses naming conventions that may differ slightly
% from the actual filenames on disk. For example, the annotation may list
% 's0014lre' while the disk has 's0014l_re' or 's0014_re'. This function
% tries several naming variants to find the correct .hea file.
%
% Inputs:
%   patient_folder    - Full path to the patient directory
%   record_name_full  - Record name with suffix (e.g., 's0014l_re')
%   record_name_base  - Base record name (e.g., 's0014_re')
%
% Output:
%   record_name - Matched filename stem, or '' if not found

    record_name = '';
    variants = {
        strrep(record_name_full, '_', ''), ...  % s0014lre
        record_name_full, ...                   % s0014l_re
        strrep(record_name_base, '_', ''), ...  % s0014re
        record_name_base                        % s0014_re
    };
    for i = 1:length(variants)
        hea_file = fullfile(patient_folder, [variants{i} '.hea']);
        if exist(hea_file, 'file')
            record_name = variants{i};
            return;
        end
    end
end


function age = parse_age(age_str)
% PARSE_AGE - Extract numeric age from WFDB header string
    age = NaN;
    if isempty(age_str) || strcmp(age_str, 'Unknown'), return; end
    age_num = str2double(age_str);
    if ~isnan(age_num), age = age_num; end
end


function sex = parse_sex(sex_str)
% PARSE_SEX - Normalise sex string from WFDB header to 'M' or 'F'
    sex = '';
    if isempty(sex_str) || strcmp(sex_str, 'Unknown'), return; end
    sex_str = upper(strtrim(sex_str));
    if startsWith(sex_str, 'M'), sex = 'M';
    elseif startsWith(sex_str, 'F'), sex = 'F'; end
end


function is_healthy = is_healthy_diagnosis(diagnosis)
% IS_HEALTHY_DIAGNOSIS - Check whether a PTB diagnosis string indicates
% a healthy control subject
    is_healthy = false;
    if isempty(diagnosis), return; end
    diagnosis_lower = lower(diagnosis);
    healthy_keywords = {'healthy', 'control', 'normal'};
    for i = 1:length(healthy_keywords)
        if contains(diagnosis_lower, healthy_keywords{i})
            is_healthy = true;
            return;
        end
    end
end


%% ========================================================================
%  PTB WFDB FORMAT READERS
%  ========================================================================

function [fs, n_samples, signals] = read_ptb_header(hea_file)
% READ_PTB_HEADER - Parse a WFDB header file from the PTB database
%
% Reads the multi-line header to extract:
%   - Sampling frequency and number of samples
%   - Per-lead gain, baseline, and filename
%   - Patient metadata (age, sex, diagnosis) from comment lines
%
% Input:
%   hea_file - Path to .hea file
%
% Outputs:
%   fs        - Sampling frequency (Hz)
%   n_samples - Total number of samples in the recording
%   signals   - Struct with fields:
%       .leads        - Array of structs (name, gain, baseline, file)
%       .patient_info - Struct (age, sex, ecg_date, diagnosis)

    fid = fopen(hea_file, 'r');
    if fid == -1, error('Cannot open header file: %s', hea_file); end

    signals.leads = struct('name', {}, 'gain', {}, 'baseline', {}, 'file', {});
    signals.patient_info = struct('age', 'Unknown', 'sex', 'Unknown', ...
                                   'ecg_date', 'Unknown', 'diagnosis', '');
    sig_idx = 0;

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        if isempty(strtrim(line)), continue; end

        % Comment lines contain patient metadata
        if line(1) == '#'
            comment = strtrim(line(2:end));
            if startsWith(comment, 'age:')
                signals.patient_info.age = strtrim(comment(5:end));
            elseif startsWith(comment, 'sex:')
                signals.patient_info.sex = strtrim(comment(5:end));
            elseif startsWith(comment, 'ECG date:')
                signals.patient_info.ecg_date = strtrim(comment(10:end));
            elseif startsWith(comment, 'Reason for admission:')
                signals.patient_info.diagnosis = strtrim(comment(22:end));
            end
            continue;
        end

        parts = strsplit(line);

        % First non-comment line is the record header
        if sig_idx == 0
            fs = str2double(parts{3});
            n_samples = str2double(parts{4});
            sig_idx = sig_idx + 1;
            continue;
        end

        % Subsequent lines describe individual signal channels
        if length(parts) >= 9
            sig.file     = parts{1};
            sig.gain     = str2double(parts{3});
            sig.baseline = str2double(parts{5});
            sig.name     = lower(parts{9});
            signals.leads(end+1) = sig;
        end
    end
    fclose(fid);
end


function ecg_data = read_ptb_signals(dat_file, xyz_file, n_samples, signals)
% READ_PTB_SIGNALS - Read binary ECG data from PTB .dat and .xyz files
%
% PTB recordings split signals across two binary files:
%   .dat - Standard 12-lead ECG channels (int16, interleaved)
%   .xyz - Frank XYZ leads (int16, interleaved)
%
% Each sample is converted from raw ADC units to physical units (mV) using
% the gain and baseline specified in the header.
%
% Inputs:
%   dat_file  - Path to .dat binary file
%   xyz_file  - Path to .xyz binary file
%   n_samples - Number of samples per channel
%   signals   - Struct from read_ptb_header (contains lead metadata)
%
% Output:
%   ecg_data - [n_samples x n_channels] matrix in physical units

    dat_leads = [];
    xyz_leads = [];

    for i = 1:length(signals.leads)
        if contains(signals.leads(i).file, '.xyz')
            xyz_leads(end+1) = i; %#ok<AGROW>
        else
            dat_leads(end+1) = i; %#ok<AGROW>
        end
    end

    n_total = length(signals.leads);
    ecg_data = zeros(n_samples, n_total);

    if ~isempty(dat_leads) && exist(dat_file, 'file')
        n_dat = length(dat_leads);
        fid = fopen(dat_file, 'rb');
        raw_dat = fread(fid, [n_dat, n_samples], 'int16')';
        fclose(fid);
        for i = 1:n_dat
            idx = dat_leads(i);
            ecg_data(:, idx) = (raw_dat(:, i) - signals.leads(idx).baseline) ...
                               / signals.leads(idx).gain;
        end
    end

    if ~isempty(xyz_leads) && exist(xyz_file, 'file')
        n_xyz = length(xyz_leads);
        fid = fopen(xyz_file, 'rb');
        raw_xyz = fread(fid, [n_xyz, n_samples], 'int16')';
        fclose(fid);
        for i = 1:n_xyz
            idx = xyz_leads(i);
            ecg_data(:, idx) = (raw_xyz(:, i) - signals.leads(idx).baseline) ...
                               / signals.leads(idx).gain;
        end
    end
end