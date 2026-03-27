function [fs, n_samples, n_sig, gains, baselines, sig_names, fmt, init_vals] = read_wfdb_header(hea_file)
% READ_WFDB_HEADER - Parse a PhysioNet WFDB header (.hea) file
%
% Usage:
%   [fs, n_samples, n_sig, gains, baselines, sig_names, fmt, init_vals] = read_wfdb_header(hea_file)
%
% Inputs:
%   hea_file - Path to the .hea header file
%
% Outputs:
%   fs         - Sampling frequency (Hz)
%   n_samples  - Number of samples per signal
%   n_sig      - Number of signals
%   gains      - Gain (ADC units per mV) for each signal
%   baselines  - Baseline (ADC zero) for each signal
%   sig_names  - Cell array of signal names
%   fmt        - Storage format code for each signal
%   init_vals  - Initial value for each signal

    fid = fopen(hea_file, 'r');
    if fid == -1
        error('Cannot open header file: %s', hea_file);
    end
    raw = fread(fid, inf, '*char')';
    fclose(fid);

    lines = strsplit(raw, char(10));
    sig_names = {}; gains = []; baselines = []; fmt = []; init_vals = [];
    fs = 250; n_samples = 0; n_sig = 0;
    header_parsed = false;

    for k = 1:length(lines)
        line = strtrim(lines{k});
        if isempty(line) || line(1) == '#', continue; end

        parts = strsplit(line);

        if ~header_parsed
            if length(parts) >= 2, n_sig     = str2double(parts{2}); end
            if length(parts) >= 3
                fs_str = parts{3};
                idx = strfind(fs_str, '/');
                if ~isempty(idx), fs_str = fs_str(1:idx(1)-1); end
                fs = str2double(fs_str);
            end
            if length(parts) >= 4, n_samples = str2double(parts{4}); end
            header_parsed = true;
            continue;
        end

        % Signal line
        if length(parts) >= 2, fmt(end+1)  = str2double(parts{2}); else, fmt(end+1)  = 16; end

        if length(parts) >= 3
            gm = regexp(parts{3}, '^([\d.]+)', 'tokens');
            if ~isempty(gm), gains(end+1) = str2double(gm{1}{1}); else, gains(end+1) = 200; end
            bm = regexp(parts{3}, '\((-?\d+)\)', 'tokens');
            if ~isempty(bm), baselines(end+1) = str2double(bm{1}{1}); else, baselines(end+1) = 0; end
        else
            gains(end+1) = 200; baselines(end+1) = 0;
        end

        if length(parts) >= 5 && baselines(end) == 0
            baselines(end) = str2double(parts{5});
        end
        if length(parts) >= 6, init_vals(end+1) = str2double(parts{6}); else, init_vals(end+1) = 0; end
        if length(parts) >= 9, sig_names{end+1} = parts{9}; else, sig_names{end+1} = sprintf('Sig%d', length(sig_names)+1); end
    end

    if fs <= 0 || isnan(fs), fs = 250; end
    if n_samples <= 0 || isnan(n_samples)
        d = dir(strrep(hea_file, '.hea', '.dat'));
        if ~isempty(d), n_samples = floor(d.bytes / 2 / n_sig); end
    end
    while length(sig_names) < n_sig
        sig_names{end+1} = sprintf('Sig%d', length(sig_names)+1);
    end
end
