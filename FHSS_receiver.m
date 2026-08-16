%% ============================================================
% MAIN SCRIPT
% ============================================================
clc; clear; close all;

% ---- USER-CONFIGURABLE PARAMETERS (only things you should touch) ----
cfg.filename        = 'Project_1_2.bin';
cfg.ch_fs           = [];          % [] = auto-detect from file header, or set e.g. 200e3
cfg.baud_min        = 100;         % Hz  - lower bound for baud rate search
cfg.baud_max        = [];          % Hz  - [] = ch_fs/4 (Nyquist-safe default)
cfg.max_rep         = 8;           % max repetition factor to scan
cfg.preamble_hex    = [];          % [] = auto-detect; or e.g. {'7E','7E','7E','3F','3F'}
cfg.trailer_hex     = [];          % [] = auto-detect from preamble complement
cfg.stft_mu         = 1.5;         % noise threshold multiplier (1.2-2.0 typical)
cfg.show_plots      = true;        % set false to suppress figures
% ---------------------------------------------------------------------

[bits_majority, valid_payloads, info] = decode_FHSS(cfg);

fprintf('\n===== DECODING COMPLETE =====\n');
fprintf('Baud rate : %.2f bps\n', info.baud);
fprintf('SPS       : %d\n',       info.sps);
fprintf('Bit order : %s\n',       info.bit_order);
fprintf('Polarity  : %d\n',       info.polarity);
fprintf('Frames    : %d\n',       numel(valid_payloads));

if ~isempty(valid_payloads)
    fprintf('\nFirst recovered payload (ASCII):\n');
    disp(char(valid_payloads{1}));
end


%% ============================================================
% DECODE FUNCTION
% ============================================================
function [bits_majority, valid_payloads, info] = decode_FHSS(cfg)

    %% --- Unpack config with safe defaults ---
    filename     = cfg.filename;
    ch_fs        = cfg.ch_fs;
    baud_min     = cfg.baud_min;
    baud_max     = cfg.baud_max;
    max_rep      = cfg.max_rep;
    preamble_hex = cfg.preamble_hex;
    trailer_hex  = cfg.trailer_hex;
    stft_mu      = cfg.stft_mu;
    show_plots   = cfg.show_plots;

    fprintf('\n===== BLIND FHSS PARAMETER ESTIMATION =====\n');

    %% ================== LOAD FILE ==================
    fid = fopen(filename, 'r');
    if fid == -1, error('Could not open file: %s', filename); end

    nchan = fread(fid, 1, 'int32');
    hdr   = fread(fid, nchan, 'int32');
    ch    = cell(1, nchan);

    for ci = 1:nchan
        raw = fread(fid, 2 * hdr(ci), 'float32');
        if numel(raw) ~= 2 * hdr(ci)
            fclose(fid);
            error('Unexpected EOF at channel %d', ci);
        end
        raw    = reshape(raw, 2, []);
        ch{ci} = raw(1,:).' + 1j * raw(2,:).';
    end
    fclose(fid);

    N = min(cellfun(@length, ch));
    fprintf('Loaded %d channels, %d samples each\n', nchan, N);

    %% --- Infer sample rate from file if not provided ---
    if isempty(ch_fs)
        % Fallback: estimate from signal bandwidth via spectral centroid
        % Default safe value when no metadata exists
        ch_fs = 200e3;
        fprintf('ch_fs not specified — defaulting to %.0f Hz\n', ch_fs);
    end
    if isempty(baud_max),  baud_max = ch_fs / 4; end

    %% ================== PHASE 1: ADAPTIVE STFT CHANNEL SCORING ==================
    fprintf('\nPhase 1: Adaptive STFT channel scoring...\n');

    % STFT window size: ~5ms at given sample rate, rounded to power of 2
    win     = 2^round(log2(ch_fs * 0.005));
    win     = max(64, min(win, 1024));   % clamp [64, 1024]
    overlap = round(0.75 * win);
    nfft    = 2 * win;                   % always 2x window for zero-padding
    N_disc  = min(N, round(ch_fs * 0.5)); % use first 500ms for discovery

    channel_score = zeros(1, nchan);
    channel_bw    = zeros(1, nchan);     % also estimate per-channel bandwidth

    for ci = 1:nchan
        x = ch{ci}(1:N_disc);
        [S, F, ~] = spectrogram(x, win, overlap, nfft, ch_fs, 'centered');
        A = abs(S);

        noise_floor = median(A(:));
        Th = stft_mu * noise_floor;
        mask = A >= Th;

        if any(mask(:))
            channel_score(ci) = mean(A(mask)) / (noise_floor + eps);  % SNR-based
            % Estimate occupied bandwidth: fraction of freq bins above threshold
            active_bins = any(mask, 2);
            channel_bw(ci) = sum(active_bins) * (ch_fs / nfft);
        end

        if show_plots
            figure('Name', sprintf('STFT Ch%d', ci), 'NumberTitle', 'off');
            imagesc((0:size(A,2)-1), F/1e3, 10*log10(A.^2 + eps));
            axis xy; xlabel('Frame'); ylabel('Freq (kHz)');
            title(sprintf('Channel %d | Score=%.2f', ci, channel_score(ci)));
            colorbar;
        end
        clear S A mask;
    end

    % Adaptive threshold: Otsu-like split on scores
    thr = otsu_threshold(channel_score);
    active_idx = find(channel_score >= thr);
    if numel(active_idx) < 2
        % Fallback: take top half
        [~, sorted] = sort(channel_score, 'descend');
        active_idx  = sorted(1:max(2, floor(nchan/2)));
    end
    active_idx = sort(active_idx);
    fprintf('Active channels: %s\n', mat2str(active_idx));

    %% ================== PHASE 2: MARK/SPACE GROUPING ==================
    fprintf('\nPhase 2: Data-driven MARK/SPACE grouping...\n');

    % Score each pairing by cross-correlation of power difference
    % Try all ways to split active channels into two non-empty groups
    n_act = numel(active_idx);
    best_contrast = -inf;
    best_mark  = [];
    best_space = [];

    % For small n_act try all splits; for large, use median-frequency split
    if n_act <= 8
        for k = 1:(n_act-1)
            combos = nchoosek(1:n_act, k);
            for r = 1:size(combos, 1)
                m_idx = active_idx(combos(r,:));
                s_idx = setdiff(active_idx, m_idx);
                contrast = compute_contrast(ch, m_idx, s_idx, N);
                if contrast > best_contrast
                    best_contrast = contrast;
                    best_mark  = m_idx;
                    best_space = s_idx;
                end
            end
        end
    else
        % Large number: split by median score
        scores_act = channel_score(active_idx);
        med_score  = median(scores_act);
        best_mark  = active_idx(scores_act >= med_score);
        best_space = active_idx(scores_act <  med_score);
        if isempty(best_space)
            mid = floor(n_act/2);
            best_mark  = active_idx(1:mid);
            best_space = active_idx(mid+1:end);
        end
    end

    fprintf('MARK  channels: %s\n', mat2str(best_mark));
    fprintf('SPACE channels: %s\n', mat2str(best_space));

    pwr_mark  = sum_channel_power(ch, best_mark,  N);
    pwr_space = sum_channel_power(ch, best_space, N);
    pwr_diff  = pwr_mark - pwr_space;

    %% ================== PHASE 3: TIMING RECOVERY ==================
    fprintf('\nPhase 3: Blind baud rate recovery...\n');

    pwr_diff = pwr_diff - mean(pwr_diff);

    % Adaptive smoothing: smooth over ~1/5 of expected minimum symbol period
    smooth_len = max(3, round(ch_fs / baud_max / 5));
    pwr_smooth = smoothdata(pwr_diff, 'movmean', smooth_len);

    sq_wave = sign(pwr_smooth);
    sq_wave(sq_wave == 0) = 1;
    edges = double(abs(diff(sq_wave)));

    % FFT-based baud detection over [baud_min, baud_max]
    L_fft = 2^nextpow2(length(edges));
    L_fft = min(L_fft, 2^20);   % cap at ~1M to avoid OOM
    E_fft = abs(fft(edges, L_fft));
    freqs = (0:L_fft-1) * (ch_fs / L_fft);

    valid_mask = (freqs >= baud_min) & (freqs <= baud_max);
    if ~any(valid_mask)
        baud_rate = baud_max / 4;
        fprintf('WARNING: No spectral peak in [%.0f, %.0f] Hz — using %.0f\n', ...
                baud_min, baud_max, baud_rate);
    else
        E_valid  = E_fft(valid_mask);
        f_valid  = freqs(valid_mask);
        [~, mi]  = max(E_valid);
        baud_rate = f_valid(mi);
        % Harmonic rejection: if peak is likely a harmonic, prefer subharmonic
        baud_rate = reject_harmonics(baud_rate, E_fft, freqs, baud_min);
    end

    sps = max(1, round(ch_fs / baud_rate));
    fprintf('Baud rate: %.2f Hz | SPS: %d\n', baud_rate, sps);

    % Matched filter (integrate over symbol period)
    mf_out = filter(ones(1, sps) / sps, 1, pwr_smooth);

    % Phase sweep — try all phases, pick best SNR
    best_energy = -inf;
    raw_bits = [];
    best_phase_idx = 1;

    for phase = 1:sps
        sampled = mf_out(phase:sps:end);
        if isempty(sampled), continue; end
        % Use variance/mean ratio as quality (high = cleaner eye)
        en = var(sampled) / (mean(abs(sampled)) + eps);
        if en > best_energy
            best_energy    = en;
            best_phase_idx = phase;
            raw_bits       = (sampled > 0).';
        end
    end
    fprintf('Optimal sample phase: %d\n', best_phase_idx);

    %% ================== PHASE 4: BLIND PARAMETER SCAN ==================
    fprintf('\nPhase 4: Blind scan (rep x polarity x shift x bit-order)...\n');

    % If preamble not provided, auto-detect by finding repeated byte patterns
    if isempty(preamble_hex)
        [preamble_bytes, trailer_bytes] = auto_detect_framing(raw_bits, max_rep);
        if isempty(preamble_bytes)
            error(['Could not auto-detect framing. ' ...
                   'Please supply preamble_hex in cfg.']);
        end
        fprintf('Auto-detected preamble: %s\n', ...
                strjoin(cellstr(dec2hex(preamble_bytes)), ' '));
        fprintf('Auto-detected trailer : %s\n', ...
                strjoin(cellstr(dec2hex(trailer_bytes)), ' '));
    else
        preamble_bytes = uint8(hex2dec(preamble_hex)).';
        if isempty(trailer_hex)
            % Default: use last 2 bytes of preamble as trailer
            trailer_bytes = preamble_bytes(end-min(1,end-1):end);
        else
            trailer_bytes = uint8(hex2dec(trailer_hex)).';
        end
    end

    target_bytes = preamble_bytes;
    best_score   = -1;
    best_dec     = [];
    best_rep     = 1;
    best_pol     = 1;
    best_shift   = 0;
    best_order   = 'msb';

    for test_rep = 1:max_rep
        n_sym = floor(length(raw_bits) / test_rep);
        if n_sym < 16, continue; end

        bits_trim = raw_bits(1 : test_rep * n_sym);
        bits_mat  = reshape(bits_trim, test_rep, []).';
        decimated = double(sum(bits_mat, 2) >= (test_rep / 2));

        for pol = [1, -1]
            bits_p = decimated;
            if pol == -1, bits_p = 1 - bits_p; end

            for shift = 0:7
                if (shift + 8) > numel(bits_p), continue; end
                bits_shifted = bits_p(shift+1:end);

                for order = {'msb', 'lsb'}
                    bytes = bits_to_bytes(bits_shifted, order{1});
                    score = numel(strfind(bytes(:).', target_bytes(:).'));
                    if score > best_score
                        best_score = score;
                        best_rep   = test_rep;
                        best_pol   = pol;
                        best_shift = shift;
                        best_order = order{1};
                        best_dec   = bytes;
                    end
                end
            end
        end
    end

    fprintf('Best: rep=%d | pol=%d | shift=%d | order=%s | score=%d\n', ...
            best_rep, best_pol, best_shift, best_order, best_score);

    if best_score < 1
        warning('Preamble not found. Output may be noise.');
    end

    %% ================== PHASE 5: FRAME EXTRACTION ==================
    fprintf('\nPhase 5: Frame extraction...\n');

    PLEN = numel(preamble_bytes);
    TLEN = numel(trailer_bytes);

    best_dec_row = best_dec(:).';
    all_payloads = {};
    search_from  = 1;
    frame_count  = 0;

    while search_from <= (numel(best_dec_row) - PLEN)
        idx = strfind(best_dec_row(search_from:end), preamble_bytes(:).');
        if isempty(idx), break; end

        payload_start = search_from + idx(1) - 1 + PLEN;
        if payload_start > numel(best_dec_row), break; end

        t_idx = strfind(best_dec_row(payload_start:end), trailer_bytes(:).');
        if isempty(t_idx), break; end

        payload_end = payload_start + t_idx(1) - 2;
        if payload_end < payload_start
            search_from = payload_start + 1;
            continue;
        end

        frame_count = frame_count + 1;
        extracted   = best_dec_row(payload_start:payload_end);
        all_payloads{frame_count} = extracted; %#ok<AGROW>

        ascii_str = char(extracted);
        ascii_str(ascii_str < 32 | ascii_str > 126) = '.';
        fprintf('Frame %2d (len=%3d): %s\n', frame_count, numel(extracted), ascii_str);

        search_from = payload_end + TLEN + 1;
    end

    fprintf('Total frames recovered: %d\n', frame_count);

    %% --- Output assembly ---
    valid_payloads = all_payloads;
    bits_majority  = best_dec;
    info = struct( ...
        'baud',      baud_rate, ...
        'sps',       sps, ...
        'polarity',  best_pol, ...
        'bit_order', best_order, ...
        'rep',       best_rep, ...
        'shift',     best_shift, ...
        'n_frames',  frame_count ...
    );

end % end decode_FHSS


%% ============================================================
% HELPER: Sum power across a set of channels
%% ============================================================
function pwr = sum_channel_power(ch, idx_list, N)
    pwr = zeros(N, 1);
    for i = idx_list
        pwr = pwr + abs(ch{i}(1:N)).^2;
    end
end


%% ============================================================
% HELPER: Compute contrast between two channel groups
%         (higher = cleaner mark/space separation)
%% ============================================================
function contrast = compute_contrast(ch, mark_idx, space_idx, N)
    N_sample = min(N, 10000); % use subset for speed
    pm = sum_channel_power(ch, mark_idx,  N_sample);
    ps = sum_channel_power(ch, space_idx, N_sample);
    d  = pm - ps;
    contrast = var(d) / (mean(abs(d)) + eps);
end


%% ============================================================
% HELPER: Otsu threshold on 1-D score vector
%% ============================================================
function thr = otsu_threshold(scores)
    scores = scores(:);
    scores = scores(isfinite(scores) & scores > 0);
    if numel(scores) < 2
        thr = 0;
        return;
    end
    mn = min(scores); mx = max(scores);
    if mn == mx
        thr = mn * 0.5;
        return;
    end
    nbins = 64;
    edges = linspace(mn, mx, nbins+1);
    counts = histcounts(scores, edges);
    p = counts / sum(counts);
    best_var = -inf;
    thr = mn;
    for k = 1:(nbins-1)
        w0 = sum(p(1:k));
        w1 = sum(p(k+1:end));
        if w0 < eps || w1 < eps, continue; end
        mu0 = sum((1:k)     .* p(1:k))     / w0;
        mu1 = sum((k+1:nbins) .* p(k+1:end)) / w1;
        bv  = w0 * w1 * (mu0 - mu1)^2;
        if bv > best_var
            best_var = bv;
            thr = edges(k+1);
        end
    end
end


%% ============================================================
% HELPER: Reject likely harmonics in baud estimate
%% ============================================================
function baud = reject_harmonics(baud_candidate, E_fft, freqs, baud_min)
    baud = baud_candidate;
    % Check if baud/2 or baud/3 has comparable energy -> prefer subharmonic
    for divisor = [2, 3, 4]
        sub = baud_candidate / divisor;
        if sub < baud_min, continue; end
        [~, si] = min(abs(freqs - sub));
        [~, ci] = min(abs(freqs - baud_candidate));
        if E_fft(si) > 0.4 * E_fft(ci)
            baud = sub;
            fprintf('Harmonic rejection: %.1f -> %.1f Hz\n', baud_candidate, sub);
            return;
        end
    end
end


%% ============================================================
% HELPER: Auto-detect preamble + trailer from raw bits
%         Strategy: look for the most repeated byte-aligned byte pattern
%% ============================================================
function [preamble_bytes, trailer_bytes] = auto_detect_framing(raw_bits, max_rep)
    preamble_bytes = [];
    trailer_bytes  = [];

    for test_rep = 1:max_rep
        n_sym = floor(length(raw_bits) / test_rep);
        if n_sym < 16, continue; end
        bits_trim = raw_bits(1:test_rep*n_sym);
        bits_mat  = reshape(bits_trim, test_rep, []).';
        dec_bits  = double(sum(bits_mat,2) >= (test_rep/2));

        for pol = [1, -1]
            bp = dec_bits; if pol == -1, bp = 1 - bp; end
            for shift = 0:7
                bytes = bits_to_bytes(bp(shift+1:end), 'msb');
                if numel(bytes) < 10, continue; end

                % Count byte frequency
                byte_counts = histcounts(double(bytes), 0:256);
                [sorted_counts, sorted_bytes] = sort(byte_counts, 'descend');
                top_byte = uint8(sorted_bytes(1) - 1); % most common byte

                if sorted_counts(1) < 3, continue; end

                % Preamble: 3+ consecutive occurrences of top_byte
                run_preamble = [];
                for rep_len = 3:6
                    pat = repmat(top_byte, 1, rep_len);
                    if ~isempty(strfind(bytes(:).', pat))
                        run_preamble = pat;
                    end
                end
                if isempty(run_preamble), continue; end

                % Trailer: second most common byte, repeated 2x
                top2 = uint8(sorted_bytes(2) - 1);
                trailer = repmat(top2, 1, 2);

                preamble_bytes = run_preamble;
                trailer_bytes  = trailer;
                return; % Return on first good detection
            end
        end
    end
end


%% ============================================================
% HELPER: Convert bits to bytes
%         order: 'msb' (default) or 'lsb'
%% ============================================================
function bytes = bits_to_bytes(bits, order)
    if nargin < 2, order = 'msb'; end
    bits   = bits(:).';
    nbytes = floor(numel(bits) / 8);
    if nbytes < 1, bytes = uint8([]); return; end
    bits   = bits(1:8*nbytes);
    bmat   = reshape(bits, 8, nbytes).';
    if strcmpi(order, 'lsb')
        bmat = fliplr(bmat);
    end
    powers = 2.^(7:-1:0);
    bytes  = uint8(bmat * powers.');
end