function [samples, symbols, channels] = read_wfdb_annotations(ann_file, format)
% READ_WFDB_ANNOTATIONS - Read annotations from a WFDB annotation file
%
% Usage:
%   [samples, symbols] = read_wfdb_annotations(ann_file, 'ludb')
%   [samples, symbols, channels] = read_wfdb_annotations(ann_file, 'qtdb')
%   [samples, symbols] = read_wfdb_annotations(ann_file, 'beat')
%
% Inputs:
%   ann_file - Path to the annotation file
%   format   - Annotation format:
%              'ludb' - Lobachevsky University format (fiducial points)
%              'qtdb' - QT Database format (with channel info)
%              'beat' - Simple beat annotations (R-peaks only, e.g. Fantasia .ecg)
%
% Outputs:
%   samples  - Column vector of annotation sample numbers
%   symbols  - Cell array of annotation symbols ('N', 'p', 't', '(', ')')
%   channels - Column vector of channel indices (QTDB only; empty otherwise)

    if nargin < 2
        format = 'ludb';
    end

    fid = fopen(ann_file, 'rb');
    if fid == -1
        error('Cannot open annotation file: %s', ann_file);
    end
    raw_bytes = fread(fid, inf, 'uint8');
    fclose(fid);

    switch lower(format)
        case 'ludb'
            [samples, symbols] = parse_fiducial_annotations(raw_bytes);
            channels = [];
        case 'qtdb'
            [samples, symbols, channels] = parse_qtdb_annotations(raw_bytes);
        case 'beat'
            [samples, symbols] = parse_beat_annotations(raw_bytes);
            channels = [];
        otherwise
            error('Unknown annotation format: %s', format);
    end
end

%% ========================================================================
%  LUDB fiducial point annotations
%  ========================================================================
function [samples, symbols] = parse_fiducial_annotations(raw_bytes)
    samples = [];
    symbols = {};
    current_sample = 0;

    i = 1;
    while i < length(raw_bytes)
        byte1 = raw_bytes(i);
        byte2 = raw_bytes(i + 1);
        i = i + 2;

        sample_offset = byte1 + bitand(byte2, 3) * 256;
        ann_type = bitshift(byte2, -2);

        if ann_type == 59  % SKIP
            if i + 3 <= length(raw_bytes)
                skip_bytes = raw_bytes(i:i+3);
                skip_val = typecast(uint8(skip_bytes), 'int32');
                current_sample = current_sample + double(skip_val);
                i = i + 4;
            end
            continue;
        elseif ann_type == 63  % AUX
            aux_len = sample_offset;
            if aux_len > 0
                i = i + aux_len;
                if mod(aux_len, 2) == 1
                    i = i + 1;
                end
            end
            continue;
        elseif ann_type == 0 || ann_type == 60 || ann_type == 61 || ann_type == 62
            continue;
        end

        current_sample = current_sample + sample_offset;
        symbol = ann_type_to_symbol(ann_type);

        if ~isempty(symbol)
            samples(end+1) = current_sample;
            symbols{end+1} = symbol;
        end
    end

    samples = samples(:);
end

%% ========================================================================
%  QTDB annotations (with channel tracking)
%  ========================================================================
function [samples, symbols, channels] = parse_qtdb_annotations(raw_bytes)
    samples = [];
    symbols = {};
    channels = [];
    current_sample = 0;
    current_channel = 0;

    i = 1;
    while i + 1 <= length(raw_bytes)
        byte1 = double(raw_bytes(i));
        byte2 = double(raw_bytes(i + 1));
        i = i + 2;

        sample_offset = byte1 + bitand(byte2, 3) * 256;
        ann_type = bitshift(byte2, -2);

        if ann_type == 59  % SKIP - PDP-11 format
            if i + 3 <= length(raw_bytes)
                high_low  = double(raw_bytes(i));
                high_high = double(raw_bytes(i + 1));
                low_low   = double(raw_bytes(i + 2));
                low_high  = double(raw_bytes(i + 3));

                high_word = high_low + high_high * 256;
                low_word  = low_low  + low_high  * 256;
                skip_val  = high_word * 65536 + low_word;

                current_sample = current_sample + skip_val;
                i = i + 4;
            end
            continue;
        elseif ann_type == 62  % CHN
            current_channel = sample_offset;
            continue;
        elseif ann_type == 63  % AUX
            aux_len = sample_offset;
            if aux_len > 0
                i = i + aux_len;
                if mod(aux_len, 2) == 1
                    i = i + 1;
                end
            end
            continue;
        elseif ann_type == 0 || ann_type == 60 || ann_type == 61
            if ann_type == 0 && sample_offset == 0
                if i + 1 <= length(raw_bytes)
                    if raw_bytes(i) == 0 && raw_bytes(i+1) == 0
                        break;
                    end
                end
            end
            continue;
        end

        current_sample = current_sample + sample_offset;
        symbol = ann_type_to_symbol(ann_type);

        if ~isempty(symbol)
            samples(end+1) = current_sample;
            symbols{end+1} = symbol;
            channels(end+1) = current_channel;
        end
    end

    samples  = samples(:);
    channels = channels(:);
end

%% ========================================================================
%  Beat annotations (R-peaks only, e.g. Fantasia .ecg files)
%  ========================================================================
function [samples, symbols] = parse_beat_annotations(raw_bytes)
    samples = [];
    symbols = {};
    current_sample = 0;

    i = 1;
    while i + 1 <= length(raw_bytes)
        b1 = double(raw_bytes(i));
        b2 = double(raw_bytes(i+1));
        i = i + 2;

        so = b1 + bitand(b2, 3) * 256;
        at = bitshift(b2, -2);

        if at == 59  % SKIP
            if i + 3 <= length(raw_bytes)
                hw = double(raw_bytes(i)) + double(raw_bytes(i+1))*256;
                lw = double(raw_bytes(i+2)) + double(raw_bytes(i+3))*256;
                current_sample = current_sample + hw*65536 + lw;
                i = i + 4;
            end
            continue;
        elseif at == 62, continue;
        elseif at == 63
            if so > 0, i = i + so + mod(so, 2); end
            continue;
        elseif at == 0
            if so == 0 && i + 1 <= length(raw_bytes)
                if raw_bytes(i) == 0 && raw_bytes(i+1) == 0, break; end
            end
            continue;
        elseif at == 60 || at == 61, continue;
        end

        current_sample = current_sample + so;

        % Beat annotation types 1-11 are R-peaks
        if at >= 1 && at <= 11
            samples(end+1) = current_sample;
            symbols{end+1} = 'N';
        end
    end

    samples = samples(:);
end

%% ========================================================================
%  Annotation type to symbol mapping
%  ========================================================================
function symbol = ann_type_to_symbol(ann_type)
    switch ann_type
        case 1,  symbol = 'N';   % Normal beat (R-peak)
        case 24, symbol = 'p';   % P-wave peak
        case 27, symbol = 't';   % T-wave peak
        case 28, symbol = 'u';   % U-wave
        case 39, symbol = '(';   % Waveform onset
        case 40, symbol = ')';   % Waveform offset
        otherwise, symbol = '';
    end
end
