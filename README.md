# HK 20-inch PMT Pre-calibration Analysis (Korean group)

Analysis code for the 20-inch PMT pre-calibration setup at Korea. The DAQ
writes RAW waveforms; the pipeline below turns them into per-run results and
the plots used in the campaign reports.

## Pipeline

The three stages `script_v7.sh` launches, in order:

| Stage | File | What it does |
|---|---|---|
| 1. Acquire | `ADC_test8_20260808.cpp` | CAEN digitizer readout (built by `Makefile` into `execute_DAQ_v2`) |
| 2. Produce | `prod_ntp_v7.C` | RAW waveforms -> per-event histograms (charge, pulse height, timing), per-run ROOT file |
| 3. Read | `read_ntp_v7.C` | Fits those histograms -> QE, SPE gain, TTS/FWHM, dark rate, per-run result tree |
| 4. Contour | `Draw_Contour_v3.C` | Per-run angular contour plot |

Shared headers: `path_builder.h`, `path_builder2.h`, `angle_convert.h`,
`repro_common.h`. ROOT session defaults live in `rootlogon.C`.

### Current analysis settings

- **Threshold**: 1.5 mV (`Analysis_Threshold_mV` in `prod_ntp_v7.C`).
  Lowered from 3 mV on 2026-09-01; it feeds the amplitude cut, the timing
  crossing, and the dark count together, so the campaign was reprocessed
  from RAW after the change.
- **QE definition**: from the pulse-height distribution, with the timing cut
  applied as a logical AND -- not from the charge fit.
- **Timing cut**: peak +- 20 ns, re-centred per run on the channel's own
  exGaus fit mean. Applied uniformly to all three channels (before
  2026-08-31 only the Monitor channel was peak-centred while Rot1/Rot2 used
  a fixed absolute window, whose margin shrank to ~3.5 samples at high HV).

## Layout

```
.                  pipeline stages, shared headers, build
reports/           the standard report plots (uniformity, stability,
                   reproducibility, ACF, per-PMT summaries)
studies/
  bfield/          B-field campaign: pulse shape / distributions, centre and Y-40deg
  hv_scan/         HV scan: pulse shape, time-walk, gain curves
  threshold/       3 mV -> 1.5 mV threshold study (prod/read variants + comparisons)
  int_window/      charge integration-window study
tools/             run-info dumps, waveform/TQ viewers, laser and monitoring scripts
  bfield_sim/      coil field simulation, compared against the measured field
reprocess/         one-off batch reprocessing runs, kept as an audit trail of
                   which reprocessing was applied when
logs/              run/taking logs referenced by the analysis
legacy/            superseded sources kept for reference
```

## Running it

```bash
make                                  # build the DAQ binary
root -l -b -q 'prod_ntp_v7.C(<run>)'  # stage 2
root -l -b -q 'read_ntp_v7.C(<run>)'  # stage 3
```

Most macros under `reports/` and `studies/` are run the same way and read the
per-run outputs produced above; each carries a `// Usage:` line at the top.
