function raw = read_wfdb_signal(dat_file, n_samples, n_sig, fmt)
% READ_WFDB_SIGNAL - Read binary signal data from a WFDB .dat file
%
% Usage:
%   raw = read_wfdb_signal(dat_file, n_samples, n_sig, fmt)
%
% Inputs:
%   dat_file  - Path to the .dat binary signal file
%   n_samples - Number of samples per signal (from header)
%   n_sig     - Number of signals (from header)
%   fmt       - Storage format array (from header); supports 16 and 212
%
% Output:
%   raw - [n_samples x n_sig] matrix of raw ADC values

    fid = fopen(dat_file, 'rb');
    if fid == -1, error('Cannot open signal file: %s', dat_file); end

    if ~isempty(fmt) && fmt(1) == 212
        raw = read_format_212(fid, n_samples, n_sig);
    else
        raw = double(fread(fid, [n_sig, n_samples], 'int16')');
    end
    fclose(fid);
end

function data = read_format_212(fid, n_samples, n_sig)
% Decode WFDB format 212 (12-bit packed, 2 samples per 3 bytes)
    rb = fread(fid, inf, 'uint8');
    nt = floor(length(rb) / 3);

    all_samples = zeros(nt * 2, 1);
    for s = 1:nt
        bi = (s-1)*3 + 1;
        if bi+2 > length(rb), break; end
        b0 = rb(bi); b1 = rb(bi+1); b2 = rb(bi+2);
        v0 = double(b0) + double(bitand(b1,15))*256;
        if v0 >= 2048, v0 = v0 - 4096; end
        v1 = double(b2) + double(bitshift(b1,-4))*256;
        if v1 >= 2048, v1 = v1 - 4096; end
        all_samples((s-1)*2 + 1) = v0;
        all_samples((s-1)*2 + 2) = v1;
    end

    % Deinterleave into n_sig channels
    total_samples = length(all_samples);
    n_frames = min(n_samples, floor(total_samples / n_sig));
    data = zeros(n_frames, n_sig);
    for ch = 1:n_sig
        data(:, ch) = all_samples(ch:n_sig:ch + (n_frames-1)*n_sig);
    end
end
