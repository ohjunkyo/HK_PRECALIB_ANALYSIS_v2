# Analysis pipeline -- README

This covers the *analysis* side of `ADC_test` (raw file -> plots), separate
from `README`/`README_latest` which cover the DAQ hardware/PC setup.

All commands below assume you're in `/home/precalkor/ADC/ADC_test`.

## Quick start: `analyze.sh`

Every analysis macro can be run directly from ROOT (`root -l -b -q
'Macro.C(args)'`), but the argument syntax (TString quoting, `{}`-vector
literals for multi-tag macros) is easy to get wrong by hand. `./analyze.sh`
wraps all of the macros below in a normal CLI:

```bash
./analyze.sh help
```

lists every subcommand with its arguments. This is the intended way for a
second person to run analysis without going through the GUI at all -- the
GUI's own pipeline buttons (Produce/Analysis/Contour/Waveform) call the exact
same macros via `run_cpp_script_v2.sh`, just with the tilt/rot/wavelength
already filled in from that scan point; nothing about running these from
`analyze.sh` conflicts with also running scans from the GUI at the same time,
since they're independent read/write operations on already-recorded files.

## GUI vs. CLI -- what needs the GUI, what doesn't

The GUI (`DAQ_Control_SW/main.py`) is required for anything that talks to
**hardware**: starting a DAQ acquisition, moving the rotation stage, firing
the laser, reading HV/B-field/temperature live. None of that has a CLI
equivalent, on purpose -- hardware control needs the GUI's interlocks and
state tracking (is_running, is_moving, the Console job-slot bookkeeping) to
stay safe.

Everything downstream of a **finished** raw/produced file is pure data
analysis and has no hardware dependency, so all of it is now reachable from
`analyze.sh`:

| Task | GUI button | CLI |
|---|---|---|
| raw -> produced | "3. Produce" | `./analyze.sh prod <run>` |
| produced -> analyzed | (part of Produce chain) | `./analyze.sh read <run>` |
| 2D contour | "6. Waveform (2D Contour)" | `./analyze.sh contour <run>` |
| waveform inspection | "5. Waveform inspection" | `./analyze.sh waveform <run>` |
| Uniformity graphs + report | "7. Uniformity (tag, runs)" | `./analyze.sh uniformity <tag> <start> <end>` |
| Overlay multiple tags | (menu item) | `./analyze.sh overlay` (reads `Data/UNIFORMITY/overlay_tags.txt`) |
| Stability report | (menu item) | `./analyze.sh stability <tag> <start> <end>` |
| Reproducibility (2 repeats) | -- (analysis-only) | `./analyze.sh repro <tag1> <tag2>` |
| Reproducibility (N sets vs a reference) | -- (analysis-only) | `./analyze.sh repron <refTag> <tag2> ...` |
| Reproducibility (3+ repeats, covariance) | -- (analysis-only) | `./analyze.sh repro3 <tag1> <tag2> <tag3> ...` |
| Per-PMT cathode uniformity (Hamamatsu style) | -- (analysis-only) | `./analyze.sh hama <tag> [--separate-axes]` |
| Per-PMT representative values | "PMT Info" tab | `./analyze.sh pmtrep [metric]` |

So: a second person analyzing alongside the DAQ operator can do everything
in this table over SSH/CLI without needing the GUI open at all, as long as
the raw files they need already exist on disk.

## Stability

**What it answers:** "does this PMT's QE/Gain/TTS/ChargeResolution drift
over hours at a fixed angle, and is any drift beyond counting statistics?"

```bash
./analyze.sh stability 20260813 800 899
```

`tag` + `run_start`/`run_end` select the FinalResult files
(`precal_result_kor_run_<tag>_<run>.root`) that make up one Stability block.
Output: `Data/image/Stability/Stability_Report_<tag>_<start>_<end>.pdf`, with
per-channel QE/Gain/TTS/ChargeResolution vs. real acquisition time, plus
B-field and Dark Box temperature/humidity pages pulled from the same window
(reuses `Draw_Overlay_Uniformity_v7.C`'s `DataCondition*` helpers -- don't
reimplement B-field parsing elsewhere, see that macro if you need to extend
it).

Known finding (see `Draw_ACF_Diagnostic.C` for the underlying autocorrelation
analysis): Stability data shows a real, reproducible ~40.2 min QE
periodicity, phase-locked (with a stable ~24 min lag) to the Dark Box ambient
temperature -- not a DAQ/software artifact (checked: CPU/IO were idle
throughout via `sar`), and not the laser diode's own TEC temperature either
(checked directly against `~/ADC/ADC_test/LOG/LASER/laser_data_<wl>_<date>.csv`'s
`temp_c` column -- flat, no periodicity). Most likely explanation: ambient
temperature cycling (HVAC duty cycle) reaching the optical path through some
channel other than the diode's own regulated junction temperature.

## Reproducibility

**What it answers:** "if I scan the same angles twice (or more), how much do
the two measurements disagree, and is that disagreement bigger than
statistics alone would predict?"

Needs `Data/UNIFORMITY/Graphs_Uniformity_<tag>.root` for each tag involved --
build those first with `./analyze.sh uniformity <tag> <start> <end>` if they
don't already exist.

**Two repeats** (`repro`, `Draw_Reproducibility.C`): point-by-point pull test,
`RMS(pull) ~ 1` means the two scans agree to within their own stated
statistical error; `> 1` means there's a real systematic on top of that. This
is the trustworthy metric for comparing any two General Scan sets.

**How the systematic is quoted (changed 2026-08-25).** It used to print
`sqrt(max(0, sigma_obs^2 - sigma_exp^2))`, which reports exactly `0.00 %`
whenever the observed scatter lands at or below the statistical expectation --
routine with N=23 points, and not at all the same statement as "there is no
systematic". RMS(pull) itself carries an uncertainty of `R/sqrt(2N)` (~15% at
N=23), so a real systematic of a couple of percent is perfectly compatible
with measuring `RMS(pull) < 1`. The macro now reports a one-sided 95% CL
upper limit instead:

    R_UL    = R * (1 + 1.645/sqrt(2N))
    sys_UL  = sqrt(max(0, R_UL^2 - 1)) * sigma_exp

and only quotes a *value* when `R - 1 > 1.645 * R/sqrt(2N)` (i.e. the
systematic is actually resolved). If even `R_UL < 1`, the per-point errors are
larger than the scatter they describe -- they are overestimated, a limit built
on them would be meaningless, and the macro falls back to the always-valid
bound `sys <= sigma_obs` and says so.

**Multi-set against one reference** (`repron`, `Draw_Reproducibility_Multi.C`):
with N sets there are N(N-1)/2 pairs and no order to read them in, so drift
over time is invisible. Pinning the first tag as reference gives N-1
comparisons on a shared baseline, plus a summary table, a combined pull page,
and a page listing the coil current / |B| each set was taken under (mean +-
stddev over that set's own window) -- because "were these even taken under the
same conditions?" is the first question any reproducibility difference raises.

```bash
./analyze.sh repro 20260815_0_45 20260816_0_45
```

**Three-or-more repeats** (`repro3`, `Draw_Reproducibility_Covariance3.C`):
estimates an empirical point-to-point covariance matrix across all repeats
and tests each one against it. **Caveat, tested and confirmed 2026-08:** with
N=23 scan points and only R=3-4 repeats, this test is close to
self-referential (the covariance is built from the same handful of vectors
it's testing) and gives a uniform, uninformative "consistent" verdict
regardless of the actual data -- confirmed by running it on real R=3 and R=4
General Scan sets and getting an identical chi2/N for every single
repeat/axis/PMT both times. Don't trust `repro3`'s verdict until R is much
larger (dozens of repeats) relative to N; use pairwise `repro` instead in the
meantime.

## Uniformity

`Draw_Uniformity_Norm_v7.C` reads the FinalResult files for one scan block,
Monitor-normalizes, and writes both a PDF report and the
`Graphs_Uniformity_<tag>.root` file that `repro`/`repro3`/`overlay` all read.
Run this first for any new scan block before doing reproducibility on it.

```bash
./analyze.sh uniformity 20260817 0 45          # channels 0,1,2 (default)
./analyze.sh uniformity 20260817 0 45 0,1      # only channels 0,1
```

## Overlay

`Draw_Overlay_Uniformity_v7.C` overlays multiple `Graphs_Uniformity_<tag>.root`
tags on the same plot (e.g. comparing wavelengths, or comparing dates). It
takes no CLI args -- edit `Data/UNIFORMITY/overlay_tags.txt` (one tag per
line, `tag,label` for a custom legend label) and optionally
`Data/UNIFORMITY/overlay_channels.txt` (comma-separated channel indices),
then:

```bash
./analyze.sh overlay
```

## Single-run pipeline (prod / read / contour / waveform)

These operate on one run number at a time and normally run automatically as
part of a scan (`script_v7.sh`'s background analysis chain). To re-run any
of them manually (e.g. to regenerate a plot, or to process a run that failed
partway through):

```bash
./analyze.sh prod 145                 # raw -> produced (.prd)
./analyze.sh read 145                 # produced -> analyzed result
./analyze.sh contour 145              # 2D contour PNG/PDF
./analyze.sh waveform 145             # waveform inspection
```

Pass an explicit file path as the 2nd arg if the run number alone doesn't
resolve to the right file (e.g. a run that was reprocessed into a different
directory).

## Per-PMT results

`hama` (`Draw_HamamatsuStyle_Uniformity.C`) draws the Hamamatsu-slide-style
cathode uniformity page per PMT (`--separate-axes` puts X and Y on their own
pages) and, as a side effect, appends that PMT's Center/Side/C-S-ratio/sigma
to two accumulating stores:

- `Data/summary/pmt_uniformity_summary.csv`
- `Data/summary/pmt_uniformity_summary.root`  (TTree `pmt_summary`)

Both are written in LONG format -- one row per (PMT, tag, metric, **axis**)
rather than one wide row with `cs_X`/`cs_Y`/`cs_pooled` columns. That is what
makes "store together now, split X vs Y later" a filter rather than a schema
change:

```bash
root -l Data/summary/pmt_uniformity_summary.root
root [1] pmt_summary->Scan("serial:tag:axis:cs_ratio", "axis==\"X\"")
```

`pmtrep` (`PMT_Representative.C`) turns those rows into one representative
number per PMT/axis. A single scan gives one value; several give a spread, and
that spread ACROSS repeats -- not the per-point statistical error -- is the
dominant uncertainty here:

    value      = (1/R) * sum_r x_r
    sig_repro  = sqrt( 1/(R-1) * sum_r (x_r - value)^2 )     [R >= 2]
    err(value) = sig_repro / sqrt(R)

With R = 1 there is no reproducibility estimate at all, and that is reported
as "single meas." rather than "+- 0" -- same reasoning as the systematic
limits above. It also prints, per measurement window, the mean +- stddev of
Dark Box T/H, laser TEC temperature and drive, coil current and |B|, so two
PMTs can be checked for having been measured under comparable conditions
(previously only a campaign-wide average existed, which cannot answer that).

Outputs: `Data/summary/pmt_representative.csv` and
`Data/image/Summary/PMT_Representative.pdf`. The GUI's **PMT Info** tab reads
the same CSVs and plots each PMT's repeat curves with their mean +-
reproducibility band; it is still marked WIP.

## Angle conventions

Graphs are plotted against the Hamamatsu incidence angle by default. Pass
`xaxis="rawstage"` to `Draw_Uniformity_Norm_v7` to plot against the stage
angle the scan actually commanded instead (what the DAQ GUI's Live Scan view
is keyed on):

```bash
root -l -b -q 'Draw_Uniformity_Norm_v7.C("20260817",0,45,"0,1,2","rawstage")'
```

Each output file is stamped with an `AngleAxis` TNamed recording which
convention it used, and a `gr_<serial>_<axis>_RawStageAngle` companion graph
carrying the raw stage tilt of every point. Overlay reads that stamp and
refuses to combine tags that don't share a convention, rather than silently
drawing two different quantities on one axis.

The GUI's live view uses `DAQ_Control_SW/angle_convert.py`, a line-for-line
port of `angle_convert.h` (verified identical over all 8 cable directions x
5 tilts x both axes). If the polynomial or sign rules change in the header,
change them in the port too.

## Housekeeping -- Underway 

- **ACLiC build artifacts** (`*_C.so`, `*_C.d`, `*_C_ACLiC_dict_rdict.pcm`)
  used to get written next to the source `.C` files whenever a macro was run
  with `.C+` (compiled) instead of plain interpreted. `rootlogon.C` in this
  directory now redirects ACLiC's build dir to `.aclic_cache/`, so this
  doesn't happen anymore as long as you start ROOT from this directory (or
  `analyze.sh`, which always does). Safe to `rm -rf .aclic_cache` any time.
- If you see stray `*_C.so`/`*_C.d` files at the top level anyway, they were
  built from a different working directory (rootlogon.C only applies to the
  directory ROOT was started in) -- just delete them.
