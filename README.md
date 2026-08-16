# Blind FHSS-FSK Receiver

A MATLAB receiver for **Frequency Hopping Spread Spectrum** signals carrying binary FSK data, built for ECE 303P. The receiver recovers the transmitted text from a channelized I/Q capture **without being told the baud rate, framing, bit order, polarity, or which channels carry MARK and SPACE** — every one of those parameters is estimated blindly from the data.

The repository also contains the guided-study research document covering BER/Eb-N₀, PRBS-based BER measurement, FHSS operation and motivation, processing gain, fast vs. slow hopping, and hop-rate hardware limits.

---

## Contents

- [What the receiver does](#what-the-receiver-does)
- [Input file format](#input-file-format)
- [Processing pipeline](#processing-pipeline)
- [Configuration](#configuration)
- [Running it](#running-it)
- [Output](#output)
- [Function reference](#function-reference)
- [Guided study document](#guided-study-document)
- [Requirements](#requirements)
- [Status and known gaps](#status-and-known-gaps)

---

## What the receiver does

The assignment asks for a receiver that can determine the hopping frequency range, determine and classify the hopping rate, demodulate the FSK signal, recover the text, and identify errors. This implementation approaches all of that **blindly** — rather than being handed the hopping sequence and bit duration, it works them out from the capture:

| Unknown | How it is recovered |
|---|---|
| Which channels carry signal | Adaptive STFT energy scoring + Otsu threshold |
| Which are MARK vs SPACE | Exhaustive group split maximizing power-difference contrast |
| Baud rate / bit duration | FFT of the symbol-transition edge train, with harmonic rejection |
| Symbol timing phase | Matched-filter phase sweep over all `sps` offsets |
| Repetition (spreading) factor | Scanned 1…8, majority-vote decimation |
| Bit polarity | Both tried, scored against preamble |
| Byte alignment | All 8 shifts tried |
| Bit order (MSB/LSB) | Both tried |
| Preamble / trailer bytes | Auto-detected as the most repeated byte patterns |

The winning parameter combination is chosen by a single objective: **how many times the preamble appears in the decoded byte stream.**

---

## Input file format

The receiver reads a custom binary container (`Project_1_2.bin` by default) holding the signal already split into per-hop-channel baseband streams:

```
int32                    nchan            % number of channels
int32   × nchan          hdr(1..nchan)    % samples per channel
float32 × 2·hdr(ci)      per channel      % interleaved I,Q,I,Q,...
```

Each channel is reassembled as `I + jQ` into a complex column vector. Because the data arrives pre-channelized, the receiver operates on **per-channel envelope power** rather than on a single wideband stream — which is exactly how a non-coherent FSK detector behaves after the dehopper.

---

## Processing pipeline

### Phase 1 — Adaptive STFT channel scoring

Each channel's first 500 ms is transformed with a `spectrogram` whose window is sized to ≈5 ms of data (rounded to a power of two, clamped to [64, 1024]), 75 % overlap, and 2× zero-padded FFT. The median magnitude is taken as the noise floor; bins above `stft_mu × noise_floor` form an activity mask. Each channel gets:

- **score** = mean masked magnitude ÷ noise floor (an SNR proxy)
- **bandwidth** = number of active frequency bins × bin spacing

An **Otsu threshold** on the score vector splits active from idle channels, with a top-half fallback if fewer than two channels survive. This is the step that answers *"over what frequency range is the signal hopping?"*

### Phase 2 — MARK/SPACE grouping

For ≤ 8 active channels, every non-empty two-way split is enumerated with `nchoosek`, and each is scored by the contrast of the group power difference:

```
contrast = var(P_mark − P_space) / mean(|P_mark − P_space|)
```

A clean binary FSK split produces a strongly bimodal difference signal and therefore the highest contrast. For larger channel counts the code falls back to a median-score split. The winning split yields `pwr_diff = P_mark − P_space`, a soft-decision waveform where sign carries the bit.

### Phase 3 — Blind timing recovery

`pwr_diff` is mean-removed, smoothed with a moving average sized to ≈1/5 of the shortest plausible symbol period, then hard-limited to a square wave. The absolute first difference of that square wave is an **edge train** — an impulse at every symbol transition — whose FFT peaks at the baud rate.

The peak is searched only within `[baud_min, baud_max]`, then passed through `reject_harmonics`: if the candidate divided by 2, 3 or 4 still carries >40 % of the peak energy, the subharmonic is preferred (guarding against locking onto 2× the true rate).

Samples per symbol follow as `sps = round(ch_fs / baud)`. A boxcar matched filter integrates over one symbol, and all `sps` sampling phases are tried — the phase maximising `var/mean` of the sampled values (a cheap eye-opening metric) wins. Bits are sliced as `sampled > 0`.

### Phase 4 — Blind parameter scan

A four-way nested search over

- repetition factor `1 … max_rep` (majority-vote decimation of repeated bits),
- polarity (normal / inverted),
- byte-alignment shift `0 … 7`,
- bit order (`msb` / `lsb`),

packs bits into bytes and counts preamble occurrences. Highest count wins. If no preamble was supplied, `auto_detect_framing` first infers one by finding the most frequently repeated byte (3–6 consecutive occurrences) and takes the second most common byte, doubled, as the trailer.

### Phase 5 — Frame extraction

The decoded byte stream is walked with `strfind`, extracting every payload bounded by a preamble/trailer pair. Each frame is printed with non-printable bytes replaced by `.`, which makes residual bit errors immediately visible against otherwise clean ASCII.

---

## Configuration

All user-facing settings live in the `cfg` struct at the top of the script:

| Field | Default | Meaning |
|---|---|---|
| `filename` | `'Project_1_2.bin'` | Input capture |
| `ch_fs` | `[]` | Per-channel sample rate; `[]` falls back to 200 kHz |
| `baud_min` | `100` | Lower bound of the baud search (Hz) |
| `baud_max` | `[]` | Upper bound; `[]` → `ch_fs/4` |
| `max_rep` | `8` | Largest repetition factor scanned |
| `preamble_hex` | `[]` | e.g. `{'7E','7E','7E','3F','3F'}`; `[]` → auto-detect |
| `trailer_hex` | `[]` | `[]` → derived from the preamble |
| `stft_mu` | `1.5` | Noise-threshold multiplier (1.2–2.0 typical) |
| `show_plots` | `true` | Per-channel STFT figures |

Nothing below that block needs editing for a normal run.

---

## Running it

```matlab
% 1. Place the capture next to the script
% 2. Adjust cfg.filename (and cfg.ch_fs if you know the true rate)
% 3. Run
>> FHSS_receiver
```

If the true per-channel sample rate is known, **set `cfg.ch_fs` explicitly** — every downstream estimate (baud rate, `sps`, channel bandwidth) scales directly with it, so the 200 kHz fallback will silently skew all reported numbers if it is wrong.

If auto-detection fails with *"Could not auto-detect framing"*, supply `cfg.preamble_hex` manually and rerun.

---

## Output

```
===== BLIND FHSS PARAMETER ESTIMATION =====
Loaded N channels, M samples each

Phase 1: Adaptive STFT channel scoring...
Active channels: [...]

Phase 2: Data-driven MARK/SPACE grouping...
MARK  channels: [...]
SPACE channels: [...]

Phase 3: Blind baud rate recovery...
Baud rate: ... Hz | SPS: ...
Optimal sample phase: ...

Phase 4: Blind scan (rep x polarity x shift x bit-order)...
Auto-detected preamble: ...
Best: rep=... | pol=... | shift=... | order=... | score=...

Phase 5: Frame extraction...
Frame  1 (len= ...): <recovered ASCII>
Total frames recovered: ...

===== DECODING COMPLETE =====
```

Returned by `decode_FHSS`:

- `bits_majority` — the decoded byte stream under the winning parameter set
- `valid_payloads` — cell array of extracted frame payloads
- `info` — struct with `baud`, `sps`, `polarity`, `bit_order`, `rep`, `shift`, `n_frames`

---

## Function reference

| Function | Role |
|---|---|
| `decode_FHSS` | Full five-phase pipeline |
| `sum_channel_power` | Sums \|x\|² across a channel group |
| `compute_contrast` | Scores a candidate MARK/SPACE split |
| `otsu_threshold` | 1-D Otsu split on the channel score vector |
| `reject_harmonics` | Prefers a subharmonic when the baud peak looks like a harmonic |
| `auto_detect_framing` | Infers preamble/trailer from repeated byte patterns |
| `bits_to_bytes` | MSB- or LSB-first bit packing |

---

## Guided study document

`FHSS_Guided_Study_Answers.docx` responds to the six research questions:

1. **BER vs SER, and Eb/N₀** — why BER = SER for binary FSK, the Gray-coded M-ary approximation BER ≈ SER/log₂(M), and why Eb/N₀ is preferred over SNR as a bandwidth-independent metric.
2. **Live BER measurement** — PRBS transmission, receiver-side sequence synchronization, XOR comparison, and why an m-sequence beats a fixed pattern (spectral flatness, ISI stress, two-valued autocorrelation, reproducibility).
3. **FHSS operation and motivation** — hop-by-hop physical behaviour, plus four problems solved with a deployed system for each: anti-jam (HAVE QUICK II), coexistence (Bluetooth), low probability of intercept (Link-16), and multipath diversity (GSM slow frequency hopping).
4. **The two bandwidths** — instantaneous B\_ch ≈ 2Δf + R\_b vs total hopping bandwidth W\_ss, processing gain G\_p = W\_ss/B\_ch = N, and the resulting jammer power dilution.
5. **Fast vs slow hopping** — defined by hop rate against bit rate, burst-error behaviour under slow hopping vs per-bit diversity combining under fast hopping.
6. **Hop rate and bandwidth** — why neither bandwidth depends directly on hop rate, spectral splatter at each transition, blanking overhead scaling as T\_blank × R\_h, and synthesizer settling time as the binding hardware constraint (R\_h ≤ 1/(T\_s + T\_d)).

---

## Requirements

- MATLAB R2019b or later
- **Signal Processing Toolbox** — `spectrogram`, `smoothdata`
- No hardware required; the receiver runs entirely offline from the capture file

---

## Status and known gaps

Honest accounting of what is and is not implemented against the assignment's five stated abilities:

| Requirement | Status |
|---|---|
| (a) Determine the hopping frequency range | **Partial.** Active channels are identified and per-channel bandwidth is computed, but `channel_bw` is never printed and no aggregate hop-band figure is reported. |
| (b) Determine the hopping rate; classify slow vs fast | **Not implemented.** Phase 3 recovers the *baud rate*, not the hop rate. Classification needs the hop rate compared against the bit rate (slow: R\_h < R\_b; fast: R\_h ≥ R\_b). |
| (c) Demodulate given the sequence and bit duration | **Exceeded.** Both are estimated blindly rather than supplied. |
| (d) Determine the received text | **Done.** Frames are extracted and printed as ASCII. |
| (e) Identify errors in the text | **Partial.** Repetition majority voting corrects some errors and non-printables are flagged with `.`, but there is no explicit error count, BER estimate, or per-frame cross-check between repeated frames. |

Other items worth cleaning up:

- `cfg.ch_fs` is documented as "auto-detect from file header", but the container carries no sample-rate field — the code simply falls back to 200 kHz.
- `channel_bw` is computed in Phase 1 and never used.
- `bits_majority` is named for bits but returns the decoded **byte** stream.
- Timing recovery is one-shot; there is no tracking loop, so a small baud estimation error accumulates as symbol slip over a long capture.
