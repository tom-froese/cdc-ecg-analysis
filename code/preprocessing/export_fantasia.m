function [n_beats_out, n_records_out] = export_fantasia(base_path, output_file)
% EXPORT_FANTASIA - Export Fantasia beat data in unified CDC format
%
% Annotation tier: MANUAL_R_AUTO_T
%   R-peaks are from the database's .ecg beat annotations (automatic
%   beat detection provided with the PhysioNet record), refined to local
%   maxima within a +/-40 ms window. T-wave endpoints are detected
%   automatically using the tangent method.
%
% The Fantasia database contains 2-hour ECG recordings from 20 young
% (21-34 years) and 20 elderly (68-85 years) healthy subjects watching
% the movie Fantasia.
%
% IMPORTANT: No beat-level quality filters are applied during export.
% All beats where detect_t_end returns a valid (non-NaN) result are
% exported. Quality filtering is deferred to the analysis stage.
%
% Tom Froese, OIST Embodied Cognitive Science Unit
% Created: February 2026

    if ~exist(base_path, 'dir')
        error('Fantasia database not found: %s', base_path);
    end

    records_file = fullfile(base_path, 'RECORDS');
    if ~exist(records_file, 'file')
        error('RECORDS file not found in %s', base_path);
    end

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
        fprintf('  [%2d/%d] %s ... ', rec_idx, length(records), record_name);

        try
            hea_file = fullfile(base_path, [record_name '.hea']);
            dat_file = fullfile(base_path, [record_name '.dat']);
            ecg_file = fullfile(base_path, [record_name '.ecg']);

            if ~exist(hea_file, 'file') || ~exist(dat_file, 'file') || ~exist(ecg_file, 'file')
                fprintf('SKIP (missing files)\n');
                n_records_skipped = n_records_skipped + 1;
                continue;
            end

            %% Demographics
            age_group = classify_fantasia_group(record_name);
            if strcmp(age_group, 'young')
                source_subset = 'Young Healthy';
            else
                source_subset = 'Elderly Healthy';
            end
            [age, sex] = parse_fantasia_demographics(hea_file);

            %% Read header and signal
            [fs, n_samples, n_sig, gains, ~, sig_names, fmt, init_vals] = read_wfdb_header(hea_file);
            raw = read_wfdb_signal(dat_file, n_samples, n_sig, fmt);
            ecg_ch = find_ecg_channel(sig_names);
            ecg_mv = (double(raw(:, ecg_ch)) - init_vals(ecg_ch)) / gains(ecg_ch);

            %% Filter for T-end detection
            ecg_filt = bandpass_filter(ecg_mv, fs, 0.5, 40);
            ecg_smooth = lowpass_filter(ecg_filt, fs, 15);

            %% Read R-peaks from database annotations and refine
            [r_peaks, ~] = read_wfdb_annotations(ecg_file, 'beat');
            r_peaks = r_peaks(r_peaks >= 1 & r_peaks <= n_samples);

            % Correct to actual maxima (+/-40 ms search window)
            search_hw = round(0.04 * fs);
            for i = 1:length(r_peaks)
                s0 = max(1, r_peaks(i) - search_hw);
                s1 = min(length(ecg_filt), r_peaks(i) + search_hw);
                [~, idx] = max(ecg_filt(s0:s1));
                r_peaks(i) = s0 + idx - 1;
            end

            if length(r_peaks) < 2
                fprintf('SKIP (< 2 R-peaks)\n');
                n_records_skipped = n_records_skipped + 1;
                continue;
            end

            %% Detect T-ends for each beat - NO sanity filtering
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
                rows{end+1} = {record_name, 'Fantasia', source_subset, 'healthy', ...
                               age, sex, fs, record_name, n_beats_record, ...
                               r_peaks(b), t_end, r_peaks(b+1), 'manual_r_auto_t'};
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
%  FANTASIA-SPECIFIC HELPERS
%  ========================================================================

function g = classify_fantasia_group(name)
    if ~isempty(strfind(name, 'y')), g = 'young'; else, g = 'old'; end
end

function ch = find_ecg_channel(sig_names)
    ch = 1;
    for i = 1:length(sig_names)
        if strcmpi(sig_names{i}, 'ECG'), ch = i; break; end
    end
end

function [age, sex] = parse_fantasia_demographics(hea_file)
    fid = fopen(hea_file, 'r');
    raw = fread(fid, inf, '*char')';
    fclose(fid);
    age = NaN; sex = '';
    m = regexp(raw, 'Age:\s*(\d+)', 'tokens');
    if ~isempty(m), age = str2double(m{1}{1}); end
    m = regexp(raw, 'Sex:\s*(\w)', 'tokens');
    if ~isempty(m)
        s = upper(m{1}{1});
        if s == 'M' || s == 'F', sex = s; end
    end
end
