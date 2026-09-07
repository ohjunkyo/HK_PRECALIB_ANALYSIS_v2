# PMT Analysis Code — Guide

ROOT analysis code for PMT photocathode uniformity & QE measurement,
run on CAEN DT5730 waveform data.

---

## 1. Pipeline overview

```
  RAW (.root)          PRODUCED (prd .root)         RESULT (.root + .png)
  waveforms      →     per-event ntuple +     →     SPE fit + QE + plots
  (ADC samples)        histograms (per ch)          (per ch)

  precal_raw_*    prod_ntp_v7.C                read_ntp_v7.C
```

Three stages: (1) RAW holds raw ADC waveforms from the digitizer;
(2) `prod_ntp_v7.C` reduces each event to scalar features (charge, timing, …) and
fills per-channel histograms → "produced" (prd) file; (3) `read_ntp_v7.C` fits the
single-photoelectron (SPE) charge spectrum and computes Quantum Efficiency (QE).

---

## 2. `prod_ntp_v7.C` — Produce step

Input: RAW file. Output: prd file (per-channel TTree + histograms + analysis parameters).

### Per-event features
| Branch | Meaning |
|--------|---------|
| `pedestal` | pre-signal baseline level (ADC) |
| `pico` | integrated charge in the signal window |
| `max` | peak pulse height |
| `falltime` | sample where the pulse crosses threshold |
| `diff` | `falltime − triggerTime` (timing) |
| `LiveTime` | accumulated live time (from TriggerTimeTag deltas) |

### Key logic
- **Pedestal** — `GetPedestal()` over a window before the signal. Laser mode uses a
  wide window (signal fixed near sample ~610); Dark mode keeps it narrow to avoid
  overlapping the PostTrigger-based signal window.
- **Threshold** — `threshold = pedestal − Analysis_Threshold_mV / ADC_to_mV`
  (pulses are negative-going, so threshold sits a fixed mV below the pedestal).
- **Charge (`pico`)** — `GetCharge()` integrates the ADC below the pedestal over the
  signal window. A *dynamic* window (`falltime−10 … +40`) is also filled into the 2D
  `QT` (charge-vs-time) histogram.
- **Dark / noise hits** — for non-trigger channels, `GetTimesBelowThreshold()` counts
  threshold crossings *outside* the signal window → `NoiseCount`, `NoiseCountRate`
  (per channel). Filled into `NoiseTime` (arrival time) and `NoiseMT`
  (amplitude-vs-time).
  > Past bug: a single shared counter made every channel report an identical dark
  > rate; fixed with per-channel `noiseHits[i]`.
- **Angle tagging** — `angle_convert.h::GetHamamatsuAngle()` converts the stage
  rotate/tilt + cable direction into the PMT incidence angle, used in every
  histogram title.

Run:
```bash
root -l 'prod_ntp_v7.C(800, "Data/RAW/Laser/precal_raw_kor_run_20260622_800.root")'
```

---

## 3. `read_ntp_v7.C` — Analysis step (SPE fit + QE)

Input: prd file. Output: result file + charge/timing plot PNG.

### (A) SPE charge fit
`FitGaussian()` fits the charge spectrum with a **Poisson–Gaussian** model
(7 parameters): a pedestal Gaussian + a sum over n-p.e. peaks weighted by
`Poisson(n, μ)`, plus an Erf-based backscattering tail.

| Parameter | Meaning |
|-----------|---------|
| `Mu (μ)` | mean number of photoelectrons |
| `SPE_Gain` | single-p.e. charge centroid |
| `Ped/SPE_Sigma` | pedestal / SPE widths |
| `Back_Weight` | backscattering fraction |

An empty channel (e.g. disconnected) skips the fit and is flagged **"Fit Empty"**
(ROOT errors suppressed).

### (B) Two QE methods

QE here is a **relative QE** (the laser intensity is fixed, so the measured signal
level is proportional to the photocathode QE).

**1. Counting QE (PHC = Pulse-Height + Timing Cut)**
```
raw   QE = N_phc  / N_total                       (PHC cut only)
k_fact   = (N_total · timing_window) / TotalNoiseTime
real_sig = N_all − k_fact · NoiseCount            (dark subtraction)
final QE = real_sig / N_total                     (PHC + timing cut, dark-subtracted)
```
Count events whose pulse passes the pulse-height cut (`N_phc`) and, additionally, the
timing cut (`N_all`). The dark (noise) contribution expected inside the timing window
is scaled by `k_fact` from the out-of-window noise rate and subtracted.

**2. Poisson QE**
```
raw   QE = μ_fit                                  (mean p.e. from the SPE fit)
μ_dark   = NoiseRate · SignalWindow
final QE = μ_fit − μ_dark                         (dark-subtracted)
```
Uses Poisson statistics: the probability of zero photoelectrons is `P(0) = e^(−μ)`, so
the fitted `μ` *is* the mean detected p.e. per trigger, i.e. the relative QE. The
expected dark p.e. in the signal window is subtracted.

> On the charge plot: `QE_PHC: Final (raw)` (green), `QE_Poisson: Final (raw)`
> (magenta). The console prints a summary table with Timing (mean·TTS), SPE Mean,
> χ²/ndf, and both QE methods (raw vs dark-subtracted).

Run:
```bash
root -l 'read_ntp_v7.C(800, "Data/production/precal_prd_kor_run_20260622_800.root")'
```

---

## 4. Supporting macros

| File | Role |
|------|------|
| `RateScan_v7.C` | dark-rate vs threshold scan (0.5–5.0 mV) |
| `Draw_Contour_v3.C` | waveform persistence/contour (Voltage vs Sample) |
| `Draw_Uniformity_Norm_v7.C` | QE uniformity vs angle |
| `Draw_Overlay_Uniformity_v7.C` | overlay uniformity of multiple PMTs |
| `run_info.C` | per-date run summary table (SN/HV/angle/laser) |
| `check_laser_runs.py` | verify the laser LD was emitting during runs (+plot) |
| `dump_runinfo.C` | dump RunInfo branches |
| `angle_convert.h` | cable direction → PMT incidence angle (shared header) |

---

## 5. Key concepts

- **Pedestal** — baseline; charge is measured relative to it.
- **SPE (single photoelectron)** — the unit-gain charge peak; calibrates gain.
- **μ (Poisson mean)** — mean detected p.e. per trigger; relative QE ∝ μ.
- **Dark subtraction** — subtract thermionic/noise hits (`NoiseCount`/`NoiseRate`)
  measured outside the signal window.
- **ADC_to_mV** = 0.1220703125 mV/ADC (= 2000 mV / 16384; 14-bit, 2 Vpp).

Config: `/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h`
(`TriggerCh=3`, paths). Helpers: `Base/analysisCode/Analysis.cpp`.
