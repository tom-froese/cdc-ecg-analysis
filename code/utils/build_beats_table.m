function T = build_beats_table(rows)
% BUILD_BEATS_TABLE - Convert cell array of beat data into a unified table
%
% Usage:
%   T = build_beats_table(rows)
%
% Input:
%   rows - Cell array where each element is a 1x13 cell:
%          {record_id, source_database, source_subset, group,
%           age, sex, fs, recording_id, beat_idx,
%           r_sample, t_end_sample, next_r_sample, annotation_method}
%
% Output:
%   T - Table with the unified CDC beat format (14 columns)

    n = length(rows);

    record_id         = cell(n,1);
    source_database   = cell(n,1);
    source_subset     = cell(n,1);
    group             = cell(n,1);
    age               = zeros(n,1);
    sex               = cell(n,1);
    fs                = zeros(n,1);
    recording_id      = cell(n,1);
    beat_idx          = zeros(n,1);
    r_sample          = zeros(n,1);
    t_end_sample      = zeros(n,1);
    next_r_sample     = zeros(n,1);
    annotation_method = cell(n,1);

    for i = 1:n
        r = rows{i};
        record_id{i}         = r{1};
        source_database{i}   = r{2};
        source_subset{i}     = r{3};
        group{i}             = r{4};
        age(i)               = r{5};
        sex{i}               = r{6};
        fs(i)                = r{7};
        recording_id{i}      = r{8};
        beat_idx(i)          = r{9};
        r_sample(i)          = r{10};
        t_end_sample(i)      = r{11};
        next_r_sample(i)     = r{12};
        annotation_method{i} = r{13};
    end

    unique_subject_id = cellfun(@(db, rid) sprintf('%s_%s', db, string(rid)), ...
                                source_database, record_id, 'UniformOutput', false);

    T = table(record_id, unique_subject_id, source_database, source_subset, group, ...
              age, sex, fs, recording_id, beat_idx, r_sample, t_end_sample, ...
              next_r_sample, annotation_method);
end