function [n_beats_out, n_records_out] = export_ptbxl(ecg_dataset, feature_dataset, output_file)
% EXPORT_PTBXL - Export PTB-XL beat data in unified CDC format
%
% Annotation tier: ECGDELI_AUTOMATIC
%   Both R-peaks and T-wave endpoints are from the ECGDeli algorithm,
%   a validated automatic fiducial point detector. Annotations are
%   provided in the PTB-XL companion feature dataset.
%
% The PTB-XL database contains ~21,800 12-lead clinical ECGs from
% ~18,800 patients with expert-validated diagnostic labels.
%
% No beat-level quality filters are applied; all valid beats are exported.

    db_file = fullfile(ecg_dataset, 'ptbxl_database.csv');
    scp_file = fullfile(ecg_dataset, 'scp_statements.csv');

    if ~exist(db_file, 'file')
        error('Database file not found: %s', db_file);
    end

    %% Load database
    fprintf('  Loading PTB-XL database...\n');
    db = load_ptbxl_database(db_file);
    fprintf('  Found %d ECG records from %d patients\n', ...
            length(db.ecg_id), length(unique(db.patient_id)));

    % Load SCP descriptions for diagnosis translation
    scp_descriptions = struct();
    if exist(scp_file, 'file')
        scp_descriptions = load_scp_descriptions(scp_file);
        fprintf('  Loaded SCP statement descriptions\n');
    end

    %% Process each ECG record
    n_records = length(db.ecg_id);
    n_records_processed = 0;
    n_records_skipped = 0;
    n_beats_total = 0;
    fs = 500;  % PTB-XL high-resolution sampling rate

    rows = {};

    fprintf('  Processing %d ECG records...\n', n_records);
    report_interval = 1000;

    for i = 1:n_records
        ecg_id = db.ecg_id(i);
        patient_id = db.patient_id(i);

        % Build annotation file path
        subfolder = sprintf('%05d', floor(ecg_id / 1000) * 1000);
        ann_file = fullfile(feature_dataset, 'fiducial_points', 'ecgdeli', subfolder, ...
                           sprintf('%05d_points_lead_II.atr', ecg_id));

        if ~exist(ann_file, 'file')
            n_records_skipped = n_records_skipped + 1;
            continue;
        end

        try
            % Read ECGDeli annotations
            annotations = read_ecgdeli_annotations(ann_file);
            [r_peaks, t_ends, next_r_peaks] = extract_ecgdeli_rt_pairs(annotations);

            n_beats = length(r_peaks);
            if n_beats < 1
                n_records_skipped = n_records_skipped + 1;
                continue;
            end

            % Patient demographics
            age = db.age(i);
            sex_code = db.sex(i);
            if sex_code == 0, sex = 'M';
            elseif sex_code == 1, sex = 'F';
            else, sex = ''; end

            % Diagnosis and group
            scp_codes = db.scp_codes{i};
            [diagnosis, group] = parse_ptbxl_diagnosis(scp_codes, scp_descriptions);

            % Identifiers
            record_id_str = sprintf('patient%05d', patient_id);
            recording_id_str = sprintf('%05d', ecg_id);

            % Store each beat
            for b = 1:n_beats
                rows{end+1} = {record_id_str, 'PTBXL', diagnosis, group, ...
                               age, sex, fs, recording_id_str, b, ...
                               r_peaks(b), t_ends(b), next_r_peaks(b), ...
                               'ecgdeli_automatic'};
                n_beats_total = n_beats_total + 1;
            end

            n_records_processed = n_records_processed + 1;

        catch ME
            n_records_skipped = n_records_skipped + 1;
            if mod(n_records_skipped, 100) == 0
                fprintf('  Warning at ECG %d: %s\n', ecg_id, ME.message);
            end
        end

        if mod(i, report_interval) == 0
            fprintf('    Processed %d/%d records (%.1f%%), %d beats so far...\n', ...
                    i, n_records, 100*i/n_records, n_beats_total);
        end
    end

    fprintf('  Processed %d records (skipped %d), %d beats\n', ...
            n_records_processed, n_records_skipped, n_beats_total);

    T = build_beats_table(rows);
    writetable(T, output_file);
    fprintf('  Saved: %s\n', output_file);

    n_beats_out = n_beats_total;
    n_records_out = n_records_processed;
end

%% ========================================================================
%  PTB-XL-SPECIFIC HELPERS
%  ========================================================================

function db = load_ptbxl_database(db_file)
    fid = fopen(db_file, 'r');
    if fid == -1, error('Cannot open database file: %s', db_file); end

    header = fgetl(fid);
    cols = strsplit(header, ',');

    ecg_id_col = find(strcmp(cols, 'ecg_id'));
    patient_id_col = find(strcmp(cols, 'patient_id'));
    age_col = find(strcmp(cols, 'age'));
    sex_col = find(strcmp(cols, 'sex'));
    scp_codes_col = find(strcmp(cols, 'scp_codes'));

    db.ecg_id = []; db.patient_id = []; db.age = []; db.sex = [];
    db.scp_codes = {};

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end

        fields = parse_csv_line(line);

        if length(fields) >= max([ecg_id_col, patient_id_col, age_col, sex_col, scp_codes_col])
            db.ecg_id(end+1) = str2double(fields{ecg_id_col});
            db.patient_id(end+1) = str2double(fields{patient_id_col});
            age_val = str2double(fields{age_col});
            db.age(end+1) = age_val;  % NaN preserved naturally
            sex_val = str2double(fields{sex_col});
            if isnan(sex_val), db.sex(end+1) = -1;
            else, db.sex(end+1) = sex_val; end
            db.scp_codes{end+1} = fields{scp_codes_col};
        end
    end
    fclose(fid);
end

function scp_desc = load_scp_descriptions(scp_file)
    scp_desc = struct();
    fid = fopen(scp_file, 'r');
    if fid == -1, return; end

    header = fgetl(fid);
    cols = strsplit(header, ',');
    desc_col_matches = find(contains(lower(cols), 'description'));
    if ~isempty(desc_col_matches), desc_col = desc_col_matches(1);
    else, desc_col = 2; end

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        fields = parse_csv_line(line);
        if length(fields) >= desc_col
            code = strtrim(fields{1});
            desc = strtrim(fields{desc_col});
            safe_code = matlab.lang.makeValidName(code);
            scp_desc.(safe_code) = desc;
        end
    end
    fclose(fid);
end

function [diagnosis, group] = parse_ptbxl_diagnosis(scp_codes_str, scp_desc)
    diagnosis = 'Unknown';
    group = 'pathological';

    if isempty(scp_codes_str), return; end

    % Check for NORM (healthy)
    if contains(scp_codes_str, 'NORM')
        diagnosis = 'Normal ECG';
        group = 'healthy';
        return;
    end

    % Extract codes
    code_pattern = '''([A-Za-z0-9_]+)''';
    matches = regexp(scp_codes_str, code_pattern, 'tokens');
    if isempty(matches), return; end

    codes_found = cellfun(@(x) x{1}, matches, 'UniformOutput', false);

    % Priority-based categorization
    mi_codes = {'IMI', 'AMI', 'LMI', 'PMI', 'ASMI', 'ILMI', 'IPLMI', 'IPMI', ...
                'INJAL', 'INJAS', 'INJIN', 'INJLA', 'INJIL'};
    if any(ismember(codes_found, mi_codes))
        diagnosis = 'Myocardial Infarction'; return;
    end

    isch_codes = {'ISC_', 'ISCA', 'ISCI', 'ISCAL', 'ISCAS', 'ISCIN', 'ISCIL', ...
                  'ISCLA', 'STTC', 'STD_', 'STE_', 'NST_', 'NDT'};
    if any(ismember(codes_found, isch_codes))
        diagnosis = 'Ischemia/ST-T Changes'; return;
    end

    block_codes = {'LBBB', 'RBBB', 'LAFB', 'LPFB', 'IRBBB', 'CLBBB', 'CRBBB', ...
                   'AVB', '1AVB', '2AVB', '3AVB', 'WPW', 'IVCD'};
    if any(ismember(codes_found, block_codes))
        diagnosis = 'Conduction Block'; return;
    end

    hyp_codes = {'LVH', 'RVH', 'LAO', 'RAO', 'LAE', 'RAE', 'SEHYP', 'HYP'};
    if any(ismember(codes_found, hyp_codes))
        diagnosis = 'Hypertrophy'; return;
    end

    af_codes = {'AFIB', 'AFLT'};
    if any(ismember(codes_found, af_codes))
        diagnosis = 'AF/Flutter'; return;
    end

    rhythm_codes = {'SBRAD', 'STACH', 'SARRH', 'SVTAC', 'PSVT', 'TRIGU', ...
                    'BIGU', 'PAC', 'PVC'};
    if any(ismember(codes_found, rhythm_codes))
        diagnosis = 'Rhythm Abnormality'; return;
    end

    axis_codes = {'LAD', 'RAD', 'APTS'};
    if any(ismember(codes_found, axis_codes))
        diagnosis = 'Axis Deviation'; return;
    end

    if any(ismember(codes_found, {'PMI', 'PACE'}))
        diagnosis = 'Paced Rhythm'; return;
    end

    % SR alone is considered healthy
    if any(ismember(codes_found, {'SR'}))
        significant_codes = setdiff(codes_found, {'SR', 'NORM'});
        if isempty(significant_codes)
            diagnosis = 'Sinus Rhythm';
            group = 'healthy';
            return;
        end
    end

    % Default
    first_code = codes_found{1};
    safe_code = matlab.lang.makeValidName(first_code);
    if isfield(scp_desc, safe_code)
        diagnosis = scp_desc.(safe_code);
    else
        diagnosis = 'Other Pathology';
    end
end

function annotations = read_ecgdeli_annotations(ann_file)
    fid = fopen(ann_file, 'rb');
    if fid == -1, error('Cannot open annotation file: %s', ann_file); end
    raw_bytes = fread(fid, inf, 'uint8');
    fclose(fid);

    annotations.samples = [];
    annotations.labels = {};
    current_sample = 0;
    pending_sample = [];

    i = 1;
    while i < length(raw_bytes)
        byte1 = raw_bytes(i);
        byte2 = raw_bytes(i + 1);
        i = i + 2;

        sample_offset = byte1 + bitand(byte2, 3) * 256;
        ann_type = bitshift(byte2, -2);

        switch ann_type
            case 59  % SKIP
                if i + 3 <= length(raw_bytes)
                    skip_bytes = raw_bytes(i:i+3);
                    skip_val = typecast(uint8(skip_bytes), 'int32');
                    current_sample = current_sample + double(skip_val);
                    i = i + 4;
                end
            case 63  % AUX
                aux_len = sample_offset;
                if aux_len > 0 && i + aux_len - 1 <= length(raw_bytes)
                    aux_bytes = raw_bytes(i:i+aux_len-1);
                    try
                        aux_text = char(aux_bytes');
                        aux_text = strtrim(aux_text(aux_text >= 32));
                    catch
                        aux_text = '';
                    end
                    i = i + aux_len;
                    if mod(aux_len, 2) == 1, i = i + 1; end

                    if ~isempty(pending_sample) && ~isempty(aux_text) && ...
                       ~startsWith(aux_text, '#')
                        annotations.samples(end+1) = pending_sample;
                        annotations.labels{end+1} = aux_text;
                        pending_sample = [];
                    end
                end
            case 22  % NOTE
                current_sample = current_sample + sample_offset;
                pending_sample = current_sample;
            case {0, 60, 61, 62}
                continue;
            otherwise
                if ann_type > 0
                    current_sample = current_sample + sample_offset;
                end
        end
    end
    annotations.samples = annotations.samples(:);
end

function [r_peaks, t_ends, next_r_peaks] = extract_ecgdeli_rt_pairs(annotations)
    r_peaks_all = [];
    t_ends_all = [];

    for i = 1:length(annotations.labels)
        label = lower(annotations.labels{i});
        samp = annotations.samples(i);
        if contains(label, 'r peak')
            r_peaks_all(end+1) = samp;
        elseif contains(label, 't-wave offset') || contains(label, 't wave offset')
            t_ends_all(end+1) = samp;
        end
    end

    r_peaks_all = r_peaks_all(:);
    t_ends_all  = t_ends_all(:);

    r_peaks = []; t_ends = []; next_r_peaks = [];

    for i = 1:length(r_peaks_all) - 1
        r = r_peaks_all(i);
        r_next = r_peaks_all(i + 1);
        t_idx = find(t_ends_all > r & t_ends_all < r_next, 1, 'first');
        if ~isempty(t_idx)
            r_peaks(end+1) = r;
            t_ends(end+1) = t_ends_all(t_idx);
            next_r_peaks(end+1) = r_next;
        end
    end

    r_peaks = r_peaks(:); t_ends = t_ends(:); next_r_peaks = next_r_peaks(:);
end

function fields = parse_csv_line(line)
    fields = {};
    current = '';
    in_quotes = false;
    for i = 1:length(line)
        c = line(i);
        if c == '"'
            in_quotes = ~in_quotes;
        elseif c == ',' && ~in_quotes
            fields{end+1} = strtrim(current);
            current = '';
        else
            current = [current c];
        end
    end
    fields{end+1} = strtrim(current);
end
