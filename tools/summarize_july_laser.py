#!/usr/bin/env python3
"""
Summarize all July 2026 runs: laser wavelength, Pulse/Bias current (mA),
Charge (SPE, pC) and QE per PMT channel. Built for the laser-intensity
optimization study (which wavelength / how much drive current gives a good
delivered intensity).

Key columns
-----------
  run, datetime            : run id + acquisition time (from RAW file mtime)
  wavelength_nm            : laser wavelength (RunInfo.Wavelength)
  laser_mA                 : drive current recorded at acquisition (RunInfo.Laser_mA == Pulse)
  pulse_mA, bias_mA        : Pulse / Bias from the matching laser CSV row (nearest in time)
  temp_c                   : laser TEC temperature at that moment
  ch, label, SN, HV        : channel (0/1/2), cable label (Mon./Rot#1/Rot#2), PMT serial, HV
  rot, tilt                : rotation-stage angles for that channel
  charge_pC                : SPE mean charge  (tree_chN.spe_mean)
  poisson_qe, relative_qe  : QE metrics       (tree_chN.poisson_qe / relativeQE)

Monitor channel (ch0 / "Mon.") sits at a fixed position, so its charge_pC is the
cleanest angle-independent proxy for the delivered laser intensity.
"""
import os, glob, csv, bisect
from datetime import datetime
import uproot

FINAL_DIR = "/home/precalkor/ADC/ADC_test/Data/FinalResult"
RAW_DIRS  = ["/home/precalkor/ADC/ADC_test/Data/RAW",
             "/home/precalkor/external_HDD_1_4T/Data_Backup/RAW"]
LASER_DIR = "/home/precalkor/ADC/ADC_test/LOG/LASER"
OUT_DIR   = "/home/precalkor/ADC/ADC_test/Data/summary"
os.makedirs(OUT_DIR, exist_ok=True)

CH_LABEL = {0: "Mon.", 1: "Rot#1", 2: "Rot#2"}

# ---- 1. RAW acquisition timestamps: {(date,run): datetime} ------------------
raw_time = {}
for d in RAW_DIRS:
    for f in glob.glob(os.path.join(d, "**", "precal_raw_kor_run_202607*.root"), recursive=True):
        base = os.path.basename(f)
        try:
            _, _, _, _, date, run = base.replace(".root", "").split("_")
        except ValueError:
            continue
        key = (date, run)
        mt = datetime.fromtimestamp(os.path.getmtime(f))
        # prefer the earliest mtime seen (original acquisition, not a backup copy)
        if key not in raw_time or mt < raw_time[key]:
            raw_time[key] = mt

# ---- 2. Laser CSV rows per wavelength, time-sorted for nearest lookup -------
#   files: laser_data_<wl>nm_YYYYMMDD.csv  with timestamp,temp_c,pulse_ma,bias_ma,...
laser = {}  # wl(int) -> (sorted_times[list], rows[list of (pulse,bias,temp)])
for f in glob.glob(os.path.join(LASER_DIR, "laser_data_*nm_202607*.csv")):
    wl = int(os.path.basename(f).split("_")[2].replace("nm", ""))
    laser.setdefault(wl, [])
    with open(f) as fh:
        for row in csv.DictReader(fh):
            try:
                ts = datetime.fromisoformat(row["timestamp"])
                laser[wl].append((ts, float(row["pulse_ma"]),
                                  float(row.get("bias_ma", 0) or 0),
                                  float(row.get("temp_c", 0) or 0)))
            except (ValueError, KeyError):
                continue
for wl in laser:
    laser[wl].sort(key=lambda r: r[0])
laser_times = {wl: [r[0] for r in rows] for wl, rows in laser.items()}

def match_laser(wl, when):
    """Nearest laser CSV row (pulse,bias,temp) to `when` for wavelength wl."""
    rows = laser.get(wl)
    if not rows or when is None:
        return (None, None, None)
    times = laser_times[wl]
    i = bisect.bisect_left(times, when)
    cands = []
    if i < len(rows):   cands.append(rows[i])
    if i > 0:           cands.append(rows[i-1])
    best = min(cands, key=lambda r: abs((r[0]-when).total_seconds()))
    # ignore matches more than 1 hour away (wrong day / gap)
    if abs((best[0]-when).total_seconds()) > 3600:
        return (None, None, None)
    return (best[1], best[2], best[3])

# ---- 3. Walk every July result file ----------------------------------------
files = sorted(glob.glob(os.path.join(FINAL_DIR, "precal_result_kor_run_202607*.root")))
rows_out = []
skipped = 0
for n, path in enumerate(files, 1):
    base = os.path.basename(path)
    try:
        _, _, _, _, date, run = base.replace(".root", "").split("_")
    except ValueError:
        continue
    try:
        f = uproot.open(path)
        ri = f["RunInfo"]
        g = lambda k: ri[k].array()[0]
        wl = int(g("Wavelength")); laser_mA = int(g("Laser_mA"))
        run_mode = str(g("RunMode"))
        sn = {0: str(g("SN1")), 1: str(g("SN2")), 2: str(g("SN3"))}
        hv = {0: int(g("HV1")), 1: int(g("HV2")), 2: int(g("HV3"))}
        rot = {0: 0, 1: int(g("RawRotateAngle2")), 2: int(g("RawRotateAngle3"))}
        tilt = {0: 0, 1: int(g("RawTiltAngle2")), 2: int(g("RawTiltAngle3"))}
    except Exception:
        skipped += 1
        continue

    when = raw_time.get((date, run))
    pulse, bias, temp = match_laser(wl, when)
    # Pulse falls back to RunInfo Laser_mA when no CSV match
    if pulse is None:
        pulse = laser_mA

    for ch in (0, 1, 2):
        try:
            t = f[f"tree_ch{ch}"]
            charge = float(t["spe_mean"].array()[0])
            pqe = float(t["poisson_qe"].array()[0])
            rqe = float(t["relativeQE"].array()[0])
            mu = float(t["poisson_mu"].array()[0])
        except Exception:
            continue
        rows_out.append(dict(
            date=date, run=int(run),
            datetime=when.isoformat(sep=" ") if when else "",
            run_mode=run_mode, wavelength_nm=wl, laser_mA=laser_mA,
            pulse_mA=round(pulse, 2) if pulse is not None else "",
            bias_mA=round(bias, 2) if bias is not None else "",
            temp_c=round(temp, 2) if temp is not None else "",
            ch=ch, label=CH_LABEL[ch], SN=sn[ch], HV=hv[ch],
            rot=rot[ch], tilt=tilt[ch],
            charge_pC=round(charge, 4), poisson_mu=round(mu, 5),
            poisson_qe=round(pqe, 4), relative_qe=round(rqe, 4),
        ))
    if n % 100 == 0:
        print(f"  ...{n}/{len(files)} files")

# ---- 4. Write master CSV ----------------------------------------------------
cols = ["date","run","datetime","run_mode","wavelength_nm","laser_mA",
        "pulse_mA","bias_mA","temp_c","ch","label","SN","HV","rot","tilt",
        "charge_pC","poisson_mu","poisson_qe","relative_qe"]
master = os.path.join(OUT_DIR, "july_laser_summary_master.csv")
with open(master, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
    for r in rows_out:
        w.writerow(r)
print(f"\n[OK] Master rows: {len(rows_out)}  (skipped {skipped} files) -> {master}")

# ---- 5. Intensity-optimization pivot ---------------------------------------
# Monitor channel = angle-independent delivered-intensity proxy.
# Group by (wavelength, laser_mA) -> mean Mon charge, and mean test-PMT QE.
from collections import defaultdict
import statistics as st
QE_MAX = 10.0   # poisson_qe above this is a fit blow-up, not a real value

agg = defaultdict(lambda: defaultdict(list))
for r in rows_out:
    key = (r["wavelength_nm"], r["laser_mA"])
    agg[key][f"charge_ch{r['ch']}"].append(r["charge_pC"])
    if r["ch"] in (1, 2) and 0 < r["poisson_qe"] < QE_MAX:
        agg[key][f"pqe_ch{r['ch']}"].append(r["poisson_qe"])

def med(v): return round(st.median(v), 4) if v else ""

pivot = os.path.join(OUT_DIR, "july_intensity_pivot.csv")
with open(pivot, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["wavelength_nm","laser_mA","n_runs",
                "Mon_charge_pC(med)","Rot1_charge_pC(med)","Rot2_charge_pC(med)",
                "Rot1_poissonQE(med)","Rot2_poissonQE(med)"])
    for key in sorted(agg):
        wl, mA = key
        a = agg[key]
        nruns = max(len(a.get(f"charge_ch{c}",[])) for c in (0,1,2)) if a else 0
        w.writerow([wl, mA, nruns, med(a.get("charge_ch0",[])),
                    med(a.get("charge_ch1",[])), med(a.get("charge_ch2",[])),
                    med(a.get("pqe_ch1",[])), med(a.get("pqe_ch2",[]))])
print(f"[OK] Intensity pivot -> {pivot}")

# ---- 6. Console preview: wavelength x laser_mA (median, outlier-robust) -----
print("\n===== Wavelength x Drive current -> delivered intensity (medians) =====")
print("(Mon = fixed monitor PMT = angle-independent intensity proxy; QE fit blow-ups >10 excluded)")
print(f"{'WL(nm)':>7} {'Laser_mA':>9} {'n':>4} {'Mon Q[pC]':>11} {'Rot1 Q':>9} {'Rot2 Q':>9} {'Rot1 pQE':>9} {'Rot2 pQE':>9}")
for key in sorted(agg):
    wl, mA = key; a = agg[key]
    nruns = max(len(a.get(f"charge_ch{c}",[])) for c in (0,1,2))
    def mm(v): return f"{st.median(v):.3f}" if v else "   -"
    print(f"{wl:>7} {mA:>9} {nruns:>4} {mm(a.get('charge_ch0',[])):>11} "
          f"{mm(a.get('charge_ch1',[])):>9} {mm(a.get('charge_ch2',[])):>9} "
          f"{mm(a.get('pqe_ch1',[])):>9} {mm(a.get('pqe_ch2',[])):>9}")

# ---- 7. Plot: Monitor charge (delivered intensity) vs wavelength ------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Per-wavelength: all Mon-channel charges, to show delivered-intensity spread
    by_wl = defaultdict(list)
    for r in rows_out:
        if r["ch"] == 0:
            by_wl[r["wavelength_nm"]].append(r["charge_pC"])
    wls = sorted(by_wl)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = {375: "#7b2ff7", 405: "#2e6fdb", 450: "#1c9e4a", 473: "#e8a400"}
    data = [by_wl[wl] for wl in wls]
    bp = ax.boxplot(data, positions=range(len(wls)), widths=0.5,
                    patch_artist=True, showfliers=False)
    for patch, wl in zip(bp["boxes"], wls):
        patch.set_facecolor(colors.get(wl, "#888")); patch.set_alpha(0.55)
    for i, wl in enumerate(wls):
        vals = by_wl[wl]
        ax.text(i, st.median(vals), f" med={st.median(vals):.2f}\n n={len(vals)}",
                va="center", fontsize=8)
    ax.set_xticks(range(len(wls)))
    ax.set_xticklabels([f"{wl} nm" for wl in wls])
    ax.set_ylabel("Monitor PMT Charge [pC]  (delivered intensity)")
    ax.set_title("July delivered-intensity uniformity across laser wavelengths\n"
                 "(Monitor channel = fixed position, angle-independent)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    outpng = "/home/precalkor/ADC/ADC_test/Data/image/july_intensity_by_wavelength.png"
    os.makedirs(os.path.dirname(outpng), exist_ok=True)
    fig.savefig(outpng, dpi=130)
    print(f"[OK] Plot -> {outpng}")
except Exception as e:
    print(f"[WARN] plot skipped: {e}")

# =============================================================================
# 8. PER-PMT optimal-intensity view
# -----------------------------------------------------------------------------
# For SPE calibration the decision variable is occupancy mu (mean p.e./pulse):
# too high -> pile-up spoils the single-p.e. peak; too low -> poor statistics.
# Each PMT has different QE/gain, so the same laser current yields different mu.
# Within a uniformity scan the ON-AXIS (peak) point sees the most light, so the
# PEAK mu is the constraint that must stay inside the SPE window.
#
# Target on-axis window (adjust to your calibration preference):
MU_LO, MU_HI = 0.05, 0.15
# =============================================================================
# group -> {(label, SN, wavelength, laser_mA): [(mu, charge, qe), ...]}
# A point is a "valid fit" only when both mu and qe are physical; the Poisson
# fit occasionally rails to a ceiling (mu ~ 2.5, qe > 10) at low-stat off-axis
# angles -- those must not be mistaken for a high on-axis occupancy.
MU_CEIL = 0.5   # mu above this is a railed fit, not real occupancy
pmt = defaultdict(list)
for r in rows_out:
    if r["run_mode"] == "Dark":
        continue
    pmt[(r["label"], r["SN"], r["wavelength_nm"], r["laser_mA"])].append(
        (r["poisson_mu"], r["charge_pC"], r["poisson_qe"]))

pmt_csv = os.path.join(OUT_DIR, "july_per_pmt_optimal_intensity.csv")
rows_pmt = []
for (label, sn, wl, mA), vals in pmt.items():
    # valid-fit occupancies only (drop railed fits)
    mus = [m for m, _, q in vals if 0 < m < MU_CEIL and 0 < q < QE_MAX]
    if not mus:
        # fall back to any physical mu if no jointly-valid point exists
        mus = [m for m, _, _ in vals if 0 < m < MU_CEIL]
    if not mus:
        continue
    peak_mu = max(mus)          # on-axis (max illumination) among valid fits
    med_mu = st.median(mus)
    peak_q = max(c for _, c, _ in vals)
    qes = [q for _, _, q in vals if 0 < q < QE_MAX]
    if peak_mu < MU_LO:
        verdict = "UNDER (raise mA)"
    elif peak_mu > MU_HI:
        verdict = "OVER (lower mA)"
    else:
        verdict = "OK"
    rows_pmt.append(dict(label=label, SN=sn, wavelength_nm=wl, laser_mA=mA,
                         n=len(vals), peak_mu=round(peak_mu, 4),
                         median_mu=round(med_mu, 4), peak_charge_pC=round(peak_q, 3),
                         median_qe=round(st.median(qes), 3) if qes else "",
                         verdict=verdict))

# sort: PMT (SN) then wavelength then current
rows_pmt.sort(key=lambda r: (r["SN"], r["wavelength_nm"], r["laser_mA"]))
with open(pmt_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["label","SN","wavelength_nm","laser_mA",
        "n","peak_mu","median_mu","peak_charge_pC","median_qe","verdict"])
    w.writeheader()
    for r in rows_pmt:
        w.writerow(r)
print(f"[OK] Per-PMT optimal-intensity table -> {pmt_csv}")

# Console: one block per PMT
print(f"\n===== Per-PMT optimal intensity  (SPE target on-axis mu {MU_LO}-{MU_HI}) =====")
cur_sn = None
for r in rows_pmt:
    if r["SN"] != cur_sn:
        cur_sn = r["SN"]
        print(f"\n  ## {r['label']}  ({r['SN']})")
        print(f"     {'WL':>5} {'mA':>5} {'n':>4} {'peak_mu':>8} {'med_mu':>8} {'peakQ':>7} {'medQE':>7}  verdict")
    print(f"     {r['wavelength_nm']:>5} {r['laser_mA']:>5} {r['n']:>4} "
          f"{r['peak_mu']:>8} {r['median_mu']:>8} {r['peak_charge_pC']:>7} "
          f"{str(r['median_qe']):>7}  {r['verdict']}")

# Per-PMT plot: peak mu vs wavelength, one marker per drive current, target band
try:
    sns_order = []
    for r in rows_pmt:
        if (r["label"], r["SN"]) not in sns_order:
            sns_order.append((r["label"], r["SN"]))
    fig, axes = plt.subplots(1, len(sns_order), figsize=(5.2*len(sns_order), 5.2),
                             sharey=True)
    if len(sns_order) == 1:
        axes = [axes]
    for ax, (label, sn) in zip(axes, sns_order):
        pts = [r for r in rows_pmt if r["SN"] == sn]
        xs = [r["wavelength_nm"] for r in pts]
        ys = [r["peak_mu"] for r in pts]
        sizes = [max(20, min(200, r["n"]*3)) for r in pts]
        ax.axhspan(MU_LO, MU_HI, color="#2e9e4a", alpha=0.15, label=f"SPE window {MU_LO}-{MU_HI}")
        ax.scatter(xs, ys, s=sizes, c=[colors.get(x, "#555") for x in xs],
                   edgecolor="k", zorder=3)
        for r in pts:
            ax.annotate(f"{r['laser_mA']}mA", (r["wavelength_nm"], r["peak_mu"]),
                        fontsize=7, xytext=(4, 4), textcoords="offset points")
        ax.set_yscale("log")
        ax.set_title(f"{label}\n{sn}")
        ax.set_xlabel("Wavelength [nm]")
        ax.set_xticks([375, 405, 450, 473])
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("Peak (on-axis) occupancy  $\\mu$  [p.e./pulse]")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.suptitle("Per-PMT optimal laser intensity across wavelengths\n"
                 "(marker = one drive-current setting; inside green band = well-tuned)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    outpng2 = "/home/precalkor/ADC/ADC_test/Data/image/july_per_pmt_optimal_mu.png"
    fig.savefig(outpng2, dpi=130)
    print(f"[OK] Per-PMT plot -> {outpng2}")
except Exception as e:
    print(f"[WARN] per-PMT plot skipped: {e}")
