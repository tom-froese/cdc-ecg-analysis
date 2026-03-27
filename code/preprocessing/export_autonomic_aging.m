function [n_beats_out, n_records_out] = export_autonomic_aging(base_path, output_file)
% EXPORT_AUTONOMIC_AGING - Export Autonomic Aging beat data in unified CDC format
%
% Annotation tier: TANGENT_AUTOMATIC
%   R-peaks are detected using a Pan-Tompkins algorithm (no database
%   annotations available). T-wave endpoints are detected using the
%   tangent method (identical algorithm to the Fantasia pipeline).
%
% The Autonomic Aging database contains resting ECG recordings (~19 min
% average) from 1,121 healthy volunteers aged 18-92, recorded at Jena
% University Hospital. Age is encoded as group midpoints since the
% database provides age ranges rather than exact ages.
%
% IMPORTANT: No beat-level quality filters are applied during export.
% All beats where detect_t_end returns a valid (non-NaN) result are
% exported. Quality filtering is deferred to the analysis stage.
%
% Tom Froese, OIST Embodied Cognitive Science Unit
% Created: February 2026

    if ~exist(base_path, 'dir')
        error('Autonomic Aging database not found: %s', base_path);
    end

    records_file = fullfile(base_path, 'RECORDS');
    if ~exist(records_file, 'file')
        error('RECORDS file not found in %s', base_path);
    end

    csv_file = fullfile(base_path, 'subject-info.csv');
    if ~exist(csv_file, 'file')
        error('subject-info.csv not found in %s', base_path);
    end

    %% Read subject demographics
    subject_info = read_subject_info(csv_file);
    fprintf('  Loaded demographics for %d subjects\n', length(subject_info));

    %% Read record list
    fid = fopen(records_file, 'r');
    records = {};
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ~isempty(line) && ischar(line)
            records{end+1} = line;
        end
    end
    fclose(fid);
    fprintf('  Found %d records\n', length(records));

    %% Process each record
    rows = {};
    n_records_processed = 0;
    n_records_skipped = 0;
    n_beats_total = 0;
    n_failed_tend = 0;

    for rec_idx = 1:length(records)
        record_name = records{rec_idx};
        fprintf('  [%4d/%d] %s ... ', rec_idx, length(records), record_name);

        try
            hea_file = fullfile(base_path, [record_name '.hea']);
            dat_file = fullfile(base_path, [record_name '.dat']);

            if ~exist(hea_file, 'file') || ~exist(dat_file, 'file')
                fprintf('SKIP (missing files)\n');
                n_records_skipped = n_records_skipped + 1;
                continue;
            end

            %% Demographics from subject-info.csv
            info = lookup_subject(subject_info, record_name);
            age_midpoint = age_group_midpoint(info.age_group);
            if info.sex == 0, sex_str = 'M';
            else, sex_str = 'F'; end

            %% Read header and signal
            [fs, n_samples, n_sig, gains, baselines, sig_names, fmt, ~] = read_wfdb_header(hea_file);
            raw = read_wfdb_signal(dat_file, n_samples, n_sig, fmt);
            ecg_ch = find_ecg_channel(sig_names);
            ecg_mv = (double(raw(:, ecg_ch)) - baselines(ecg_ch)) / gains(ecg_ch);

            %% Filter
            ecg_filt = bandpass_filter(ecg_mv, fs, 0.5, 40);
            ecg_smooth = lowpass_filter(ecg_filt, fs, 15);

            %% R-peak detection (Pan-Tompkins)
            r_peaks = detect_r_peaks(ecg_filt, fs);

            % Correct to actual maxima (+/-40 ms, matching Fantasia)
            search_hw = round(0.04 * fs);
            for i = 1:length(r_peaks)
                s0 = max(1, r_peaks(i) - search_hw);
                s1 = min(length(ecg_filt), r_peaks(i) + search_hw);
                [~, idx] = max(ecg_filt(s0:s1));
                r_peaks(i) = s0 + idx - 1;
            end
            r_peaks = unique(r_peaks);

            if length(r_peaks) < 2
                fprintf('SKIP (< 2 R-peaks)\n');
                n_records_skipped = n_records_skipped + 1;
                continue;
            end

            %% Detect T-ends - NO sanity filtering
            n_total = length(r_peaks) - 1;
            n_beats_record = 0;
            n_failed_record = 0;

            for b = 1:n_total
                [t_end, ~, ~, ~] = detect_t_end(ecg_smooth, r_peaks, b, fs);

                if isnan(t_end)
                    n_failed_record = n_failed_record + 1;
                    continue;
                end

                % Export all beats where T-end was successfully detected
                n_beats_record = n_beats_record + 1;
                rows{end+1} = {record_name, 'Autonomic_Aging', 'Healthy Control', ...
                               'healthy', age_midpoint, sex_str, fs, record_name, ...
                               n_beats_record, r_peaks(b), t_end, r_peaks(b+1), ...
                               'tangent_automatic'};
            end

            n_beats_total = n_beats_total + n_beats_record;
            n_failed_tend = n_failed_tend + n_failed_record;
            n_records_processed = n_records_processed + 1;

            fprintf('OK  (beats=%d, T-end failed=%d)\n', n_beats_record, n_failed_record);

        catch ME
            fprintf('ERROR: %s\n', ME.message);
            n_records_skipped = n_records_skipped + 1;
        end
    end

    %% Summary
    fprintf('\n  Records processed: %d\n', n_records_processed);
    fprintf('  Records skipped:   %d\n', n_records_skipped);
    fprintf('  Total beats:       %d\n', n_beats_total);
    fprintf('  T-end detection failures: %d\n', n_failed_tend);
    if (n_beats_total + n_failed_tend) > 0
        fprintf('  T-end detection rate: %.1f%%\n', ...
                100 * n_beats_total / (n_beats_total + n_failed_tend));
    end

    T = build_beats_table(rows);
    writetable(T, output_file);
    fprintf('  Saved: %s\n', output_file);

    n_beats_out = n_beats_total;
    n_records_out = n_records_processed;
end

%% ========================================================================
%  AUTONOMIC AGING-SPECIFIC HELPERS
%  ========================================================================

function ch = find_ecg_channel(sig_names)
    ch = 1;
    for i = 1:length(sig_names)
        nm = upper(sig_names{i});
        if strcmp(nm, 'ECG') || strcmp(nm, 'ECG1') || contains(nm, 'ECG')
            ch = i; break;
        end
    end
end

function info = read_subject_info(csv_file)
    fid = fopen(csv_file, 'r');
    fgetl(fid);  % skip header
    info = struct('id', {}, 'age_group', {}, 'sex', {}, 'bmi', {}, ...
                  'rec_len', {}, 'device', {});
    k = 0;
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ~ischar(line) || isempty(line), continue; end
        parts = strsplit(line, ',');
        if length(parts) < 6, continue; end
        k = k + 1;
        info(k).id        = strtrim(parts{1});
        info(k).age_group = str2double(parts{2});
        info(k).sex       = str2double(parts{3});
        info(k).bmi       = str2double(parts{4});
        info(k).rec_len   = str2double(parts{5});
        info(k).device    = str2double(parts{6});
    end
    fclose(fid);
end

function info = lookup_subject(subject_info, record_name)
    id_num = str2double(record_name);
    for i = 1:length(subject_info)
        if str2double(subject_info(i).id) == id_num
            info = subject_info(i);
            return;
        end
    end
    info.id = record_name;
    info.age_group = NaN; info.sex = NaN;
    info.bmi = NaN; info.rec_len = NaN; info.device = NaN;
end

function m = age_group_midpoint(g)
    midpoints = [18.5, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 88.5];
    if g >= 1 && g <= 15, m = midpoints(g);
    else, m = NaN; end
end
