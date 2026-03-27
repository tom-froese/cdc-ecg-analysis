function [n_beats_out, n_records_out] = export_ludb(base_path, output_file)
% EXPORT_LUDB - Export LUDB beat data in unified CDC format
%
% Annotation tier: MANUAL
%   Both R-peaks and T-wave endpoints are from expert cardiologist
%   annotations provided with the Lobachevsky University database.
%
% The LUDB contains 200 subjects with 12-lead ECG recordings annotated
% by cardiologists. Fiducial points (P, QRS, T boundaries) are provided
% as WFDB annotation files per lead.
%
% No beat-level quality filters are applied; all valid beats are exported.

    n_subjects = 200;
    lead = 'ii';
    fs = 500;

    % Check for ludb.csv for diagnosis info
    csv_file = fullfile(base_path, 'ludb.csv');
    use_csv = false;
    csv_data = [];

    if exist(csv_file, 'file')
        try
            opts = detectImportOptions(csv_file, 'VariableNamingRule', 'preserve');
            opts = setvartype(opts, opts.VariableNames, 'char');
            csv_data = readtable(csv_file, opts);
            use_csv = true;
            fprintf('  Using ludb.csv for diagnosis information\n');
        catch
            fprintf('  Could not read ludb.csv, using header files\n');
        end
    end

    % Preallocate cell arrays
    rows = {};
    n_records_processed = 0;
    n_beats_total = 0;

    for subj = 1:n_subjects
        try
            hea_file = fullfile(base_path, 'data', [num2str(subj) '.hea']);
            ann_file = fullfile(base_path, 'data', [num2str(subj) '.' lead]);

            if ~exist(hea_file, 'file') || ~exist(ann_file, 'file')
                continue;
            end

            % Read demographics from header
            [age, sex, diagnoses_raw] = read_ludb_header(hea_file);

            % Determine health status
            if use_csv
                try
                    is_healthy = check_healthy_csv(csv_data, subj);
                catch
                    is_healthy = false;
                end
            else
                is_healthy = false;
            end

            % Read annotations and extract matched R-T pairs
            [ann_samples, ann_symbols] = read_wfdb_annotations(ann_file, 'ludb');
            [r_peaks, t_ends] = extract_rt_pairs(ann_samples, ann_symbols);

            n_beats = length(r_peaks);
            if n_beats < 2
                continue;
            end

            % Determine group and source_subset
            if is_healthy
                group = 'healthy';
                source_subset = 'Healthy';
            else
                group = 'pathological';
                source_subset = parse_ludb_diagnosis(diagnoses_raw);
            end

            if isempty(age) || isnan(age), age = NaN; end
            if isempty(sex), sex = ''; end

            record_id_str = num2str(subj);

            % Store each beat (except last, which has no next_r)
            for b = 1:(n_beats - 1)
                rows{end+1} = {record_id_str, 'LUDB', source_subset, group, ...
                               age, sex, fs, record_id_str, b, ...
                               r_peaks(b), t_ends(b), r_peaks(b+1), 'manual'};
                n_beats_total = n_beats_total + 1;
            end

            n_records_processed = n_records_processed + 1;

        catch ME
            fprintf('    Subject %d: Error - %s\n', subj, ME.message);
        end
    end

    fprintf('  Processed %d records, %d beats\n', n_records_processed, n_beats_total);

    % Build and save table
    T = build_beats_table(rows);
    writetable(T, output_file);
    fprintf('  Saved: %s\n', output_file);

    n_beats_out = n_beats_total;
    n_records_out = n_records_processed;
end

%% ========================================================================
%  LUDB-SPECIFIC HELPERS
%  ========================================================================

function [age, sex, diagnoses] = read_ludb_header(hea_file)
    age = NaN;
    sex = '';
    diagnoses = '';

    fid = fopen(hea_file, 'r');
    if fid == -1, return; end

    all_comments = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end

        if contains(line, '<age>:')
            parts = strsplit(line, ':');
            if length(parts) >= 2
                age = str2double(strtrim(parts{2}));
            end
        elseif contains(line, '<sex>:')
            parts = strsplit(line, ':');
            if length(parts) >= 2
                sex = upper(strtrim(parts{2}));
            end
        end

        if startsWith(line, '#')
            all_comments{end+1} = line;
        end
    end
    fclose(fid);

    diagnoses = strjoin(all_comments, ' | ');
end

function is_healthy = check_healthy_csv(csv_data, subj_id)
    id_col = csv_data{:, 1};

    if iscell(id_col)
        row_idx = find(strcmp(id_col, num2str(subj_id)), 1);
    else
        row_idx = find(id_col == subj_id, 1);
    end

    if isempty(row_idx)
        error('Subject %d not found in CSV', subj_id);
    end

    var_names = csv_data.Properties.VariableNames;
    rhythm_col = find(contains(lower(var_names), 'rhythm'), 1);
    if isempty(rhythm_col), rhythm_col = 2; end

    rhythm = csv_data{row_idx, rhythm_col};
    if iscell(rhythm), rhythm = rhythm{1}; end
    rhythm_lower = lower(char(rhythm));
    has_sinus = contains(rhythm_lower, 'sinus');

    row_data = strjoin(string(csv_data{row_idx, :}), ' ');
    row_lower = lower(char(row_data));

    pathology_keywords = {'fibrillation', 'flutter', 'block', 'hypertrophy', ...
                         'ischemia', 'infarction', 'tachycardia', 'bradycardia', ...
                         'deviation', 'abnormal', 'extrasystol', 'pacing', ...
                         'overload', 'enlargement'};

    has_pathology = false;
    for i = 1:length(pathology_keywords)
        if contains(row_lower, pathology_keywords{i})
            has_pathology = true;
            break;
        end
    end

    is_healthy = has_sinus && ~has_pathology;
end

function category = parse_ludb_diagnosis(diag_raw)
    diag = lower(diag_raw);

    if contains(diag, 'infarction') || contains(diag, ' mi ')
        category = 'Infarction';
    elseif contains(diag, 'ischemi')
        category = 'Ischemia';
    elseif contains(diag, 'hypertrophy') || contains(diag, 'enlargement') || contains(diag, 'overload')
        category = 'Hypertrophy';
    elseif contains(diag, 'fibrillation') || contains(diag, 'flutter')
        category = 'AF/Flutter';
    elseif contains(diag, 'block')
        category = 'Conduction Block';
    elseif contains(diag, 'tachycardia')
        category = 'Tachycardia';
    elseif contains(diag, 'bradycardia')
        category = 'Bradycardia';
    elseif contains(diag, 'extrasystol') || contains(diag, 'ectop')
        category = 'Ectopy';
    elseif contains(diag, 'deviation')
        category = 'Axis Deviation';
    elseif contains(diag, 'pacing') || contains(diag, 'pacemaker')
        category = 'Paced';
    elseif contains(diag, 'sinus') && ~contains(diag, 'arrhythmia')
        category = 'Sinus Rhythm Variant';
    else
        category = 'Other Pathology';
    end
end


