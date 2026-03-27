function [n_beats_out, n_records_out] = export_qtdb(base_path, output_file)
% EXPORT_QTDB - Export QT Database beat data in unified CDC format
%
% Annotation tier: MANUAL
%   Both R-peaks and T-wave endpoints are from expert manual annotations
%   (q1c annotator, with pu0/qt1/man fallbacks).
%
% The QTDB contains ECG excerpts from multiple source databases with
% expert-annotated QT intervals. Records are categorised by their
% source: arrhythmia, ST change, SVA, healthy, European ST-T,
% sudden death, and long-term databases.
%
% No beat-level quality filters are applied; all valid beats are exported.

    ann_type = 'q1c';
    fs = 250;

    % Define record groups by source database
    records_arrhythmia = {'sel100', 'sel102', 'sel103', 'sel104', 'sel114', 'sel116', ...
                          'sel117', 'sel123', 'sel213', 'sel221', 'sel223', 'sel230', ...
                          'sel231', 'sel232', 'sel233'};

    records_st_change = {'sel301', 'sel302', 'sel306', 'sel307', 'sel308', 'sel310'};

    records_sva = {'sel803', 'sel808', 'sel811', 'sel820', 'sel821', 'sel840', ...
                   'sel847', 'sel853', 'sel871', 'sel872', 'sel873', 'sel883', 'sel891'};

    records_healthy = {'sel16265', 'sel16272', 'sel16273', 'sel16420', 'sel16483', ...
                       'sel16539', 'sel16773', 'sel16786', 'sel16795', 'sel17453'};

    records_european = {'sele0104', 'sele0106', 'sele0107', 'sele0110', 'sele0111', ...
                        'sele0112', 'sele0114', 'sele0116', 'sele0121', 'sele0122', ...
                        'sele0124', 'sele0126', 'sele0129', 'sele0133', 'sele0136', ...
                        'sele0166', 'sele0170', 'sele0203', 'sele0210', 'sele0211', ...
                        'sele0303', 'sele0405', 'sele0406', 'sele0409', 'sele0411', ...
                        'sele0509', 'sele0603', 'sele0604', 'sele0606', 'sele0607', ...
                        'sele0609', 'sele0612', 'sele0704'};

    records_sudden_death = {'sel30', 'sel31', 'sel32', 'sel33', 'sel34', 'sel35', ...
                            'sel36', 'sel37', 'sel38', 'sel39', 'sel40', 'sel41', ...
                            'sel42', 'sel43', 'sel44', 'sel45', 'sel46', 'sel47', ...
                            'sel48', 'sel49', 'sel50', 'sel51', 'sel52', 'sel17152'};

    records_longterm = {'sel14046', 'sel14157', 'sel14172', 'sel15814'};

    all_records = [records_arrhythmia, records_st_change, records_sva, records_healthy, ...
                   records_european, records_sudden_death, records_longterm];

    rows = {};
    n_records_processed = 0;
    n_beats_total = 0;

    for i = 1:length(all_records)
        record = all_records{i};

        % Determine source subset and group
        if ismember(record, records_sudden_death)
            source_subset = 'BIH Sudden Death';
            group = 'sudden_death';
        elseif ismember(record, records_healthy)
            source_subset = 'MIT-BIH Normal Sinus Rhythm';
            group = 'healthy';
        elseif ismember(record, records_arrhythmia)
            source_subset = 'MIT-BIH Arrhythmia';
            group = 'pathological';
        elseif ismember(record, records_sva)
            source_subset = 'MIT-BIH SVA';
            group = 'pathological';
        elseif ismember(record, records_st_change)
            source_subset = 'MIT-BIH ST Change';
            group = 'pathological';
        elseif ismember(record, records_european)
            source_subset = 'European ST-T';
            group = 'pathological';
        elseif ismember(record, records_longterm)
            source_subset = 'MIT-BIH Long-Term';
            group = 'pathological';
        else
            source_subset = 'Unknown';
            group = 'unknown';
        end

        try
            ann_file = fullfile(base_path, [record '.' ann_type]);

            % Try alternative annotators if primary not found
            if ~exist(ann_file, 'file')
                for alt_ann = {'pu0', 'qt1', 'man'}
                    alt_file = fullfile(base_path, [record '.' alt_ann{1}]);
                    if exist(alt_file, 'file')
                        ann_file = alt_file;
                        break;
                    end
                end
            end

            if ~exist(ann_file, 'file')
                continue;
            end

            % Read annotations and extract matched R-T pairs
            [ann_samples, ann_symbols] = read_wfdb_annotations(ann_file, 'qtdb');
            [r_peaks, t_ends] = extract_rt_pairs(ann_samples, ann_symbols);

            n_beats = length(r_peaks);
            if n_beats < 2
                continue;
            end

            % Store each beat (except last)
            for b = 1:(n_beats - 1)
                rows{end+1} = {record, 'QTDB', source_subset, group, ...
                               NaN, '', fs, record, b, ...
                               r_peaks(b), t_ends(b), r_peaks(b+1), 'manual'};
                n_beats_total = n_beats_total + 1;
            end

            n_records_processed = n_records_processed + 1;

        catch ME
            fprintf('    Record %s: Error - %s\n', record, ME.message);
        end
    end

    fprintf('  Processed %d records, %d beats\n', n_records_processed, n_beats_total);

    T = build_beats_table(rows);
    writetable(T, output_file);
    fprintf('  Saved: %s\n', output_file);

    n_beats_out = n_beats_total;
    n_records_out = n_records_processed;
end
