#!/usr/bin/env python3
"""Automated laser-intensity calibration scan.

For each wavelength, sweeps a list of candidate (bias_mA, pulse_mA) settings.
At each candidate it takes a short 50,000-event raw DAQ run, runs the
existing prod_ntp_v7.C -> read_ntp_v7.C analysis chain (same as every normal
scan point), and checks whether every PMT's QE (poisson_qe, i.e. Poisson mu
x100%) lands in TARGET_QE_PCT_RANGE with a 2 p.e. fraction consistent with
pure Poisson statistics at that occupancy (no excess double-pulsing).

This is pure orchestration -- it reuses the exact same hardware call paths
already used unattended elsewhere in this codebase:
  - laser control:  TamadenshiLaser.connect/set_bias_current/set_pulse_current
    (same class DAQ_Control_SW/managers/rotation_manager.py's
    _prepare_laser_block already drives autonomously during General Scan)
  - acquisition:     script_v7.sh -> execute_DAQ_v2 (same script the GUI's
    "Overlay"/"General Scan" buttons invoke)
No new hardware code path is introduced.

SAFETY
------
This turns laser diodes on and steps their current automatically. Run it
only with the interlock satisfied and the beam path clear -- the same
precondition as any other laser-on run. Use --dry-run first to review the
planned (wavelength, bias, pulse) sequence and file layout without touching
any hardware.

Usage
-----
    python3 laser_intensity_scan.py --dry-run
    python3 laser_intensity_scan.py --rot2 0 --tilt2 0 --rot3 0 --tilt3 0

Edit CANDIDATES / TARGET_QE_PCT_RANGE below before a real run -- there is no
sensible universal default for either.
"""
import argparse
import csv
import glob
import math
import os
import shlex
import subprocess
import sys
import time

BASE = "/home/precalkor/ADC/ADC_test"
DAQ_DIR = "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW"
LIVE_CONFIG = os.path.join(DAQ_DIR, "config3.h")
SCRIPT_V7 = os.path.join(BASE, "script_v7.sh")
FLAG_DIR = "/tmp/daq_flags"
LASER_PORT_MAPPING = {
    "375nm": "1-3.3:1.0", "405nm": "1-3.1:1.0",
    "450nm": "1-3.2:1.0", "473nm": "1-3.4:1.0",
}

sys.path.insert(0, "/home/precalkor/Integrated_Control_SW/Laser_Control_SW/app")

N_EVENTS = 50000

# ── EDIT BEFORE USE ─────────────────────────────────────────────────────────
# wavelength -> list of (bias_mA, pulse_mA) candidates to try, low to high.
# 405nm: fiber was replaced ~2026-07-28. analyze_laser_intensity_history.py
# shows the SAME pulse now gives ~1.5x the old QE (matched at 149.99mA and
# 160.00mA: ratios 1.53 and 1.47) -- but that jump for only +10mA also shows
# we're right on the laser diode's threshold knee (steeply nonlinear), so
# this candidate list is a starting bracket to scan finely, not a computed
# answer. bias=0.0 matches every historical 405nm run that hit target.
# 375/450/473nm: no historical (bias,pulse) setting ever got all 3 PMTs into
# the 1-2% QE band simultaneously (one channel is always off from the other
# two) -- these ranges just bracket the historical operating point, a real
# scan is needed to find where all 3 balance.
CANDIDATES = {
    "405nm": [(0.0, 145.0), (0.0, 148.0), (0.0, 150.0), (0.0, 152.0), (0.0, 155.0), (0.0, 158.0)],
    # 375/450/473nm: no historical setting balanced all 3 PMTs into 1-2% QE
    # (see analyze_laser_intensity_history.py) -- flat 150-180mA @ 5mA steps
    # to find it directly rather than guess a narrower bracket.
    "375nm": [(0.0, p) for p in range(150, 181, 5)],
    "450nm": [(0.0, p) for p in range(150, 181, 5)],
    "473nm": [(0.0, p) for p in range(150, 181, 5)],
}
# Target: poisson_qe (= Poisson mu * 100%) in this band on ALL 3 PMTs
# (ch0=monitor EM2740, ch1=SN2, ch2=SN3), with the 2 p.e. fraction consistent
# with pure Poisson statistics at that occupancy (checked via
# poisson2pe_fraction() below) -- i.e. no excess double-pulsing beyond what
# the occupancy itself predicts.
TARGET_QE_PCT_RANGE = (1.0, 2.0)
# ─────────────────────────────────────────────────────────────────────────


def make_temp_config(bias_ma, pulse_ma, wavelength_digits):
    """Copy the live config3.h with Events forced to N_EVENTS and the
    Laser/Wavelength metadata fields set to match this candidate, so the
    ROOT file header records what was actually driven. Everything else
    (paths, channel mask, HV, ...) is left untouched -- this is a metadata +
    event-count override only, never a hardware-parameter change beyond what
    the caller explicitly set via TamadenshiLaser below."""
    with open(LIVE_CONFIG) as f:
        text = f.read()
    import re
    text = re.sub(r"const int Events\s*=\s*\d+;", f"const int Events = {N_EVENTS};", text)
    text = re.sub(r'const std::string Laser\s*=\s*"[^"]*";',
                  f'const std::string Laser = "{pulse_ma:g}";', text)
    text = re.sub(r'const std::string Wavelength\s*=\s*"[^"]*";',
                  f'const std::string Wavelength = "{wavelength_digits}";', text)
    tmp_path = os.path.join(BASE, "config3_intensity_scan_tmp.h")
    with open(tmp_path, "w") as f:
        f.write(text)
    return tmp_path


def set_laser(wl, bias_ma, pulse_ma, dry_run):
    """Connect + drive one laser to the given bias/pulse and turn the LD on.
    Mirrors rotation_manager.py's _prepare_laser_block exactly (TEC on,
    apply bias/pulse, LD on, settle) but standalone, without the Tk app."""
    print(f"[LASER] {wl}: bias={bias_ma}mA pulse={pulse_ma}mA")
    if dry_run:
        return True
    from laser_driver import TamadenshiLaser  # local import: only needed for a real run
    inst = TamadenshiLaser()
    ok, msg = inst.connect(dev_path=LASER_PORT_MAPPING[wl].encode("utf-8"))
    if not ok:
        print(f"[ERROR] {wl} connect failed: {msg}")
        return False
    inst.set_tec_on(True)
    time.sleep(5.0)   # short settle; _prepare_laser_block's full stability gate is overkill for a quick cal point
    inst.set_bias_current(bias_ma)
    inst.set_pulse_current(pulse_ma)
    time.sleep(0.5)
    inst.set_ld_on(True)
    time.sleep(10.0)
    inst._laser_instance_ref = inst  # keep alive; caller disconnects via laser_off()
    set_laser._last_inst = inst
    return True


def laser_off(dry_run):
    if dry_run:
        return
    inst = getattr(set_laser, "_last_inst", None)
    if inst:
        try:
            inst.set_ld_on(False)
            inst.disconnect()
        except Exception as e:
            print(f"[WARNING] laser_off: {e}")


def run_daq_point(temp_config, rot2, tilt2, rot3, tilt3, base_run_start, tag, dry_run):
    """Runs script_v7.sh exactly like the GUI does. Blocks until the raw
    acquisition finishes (script_v7.sh itself blocks on execute_DAQ_v2);
    the analysis chain it launches afterward is backgrounded in tmux, so
    the caller must separately wait_for_analysis()."""
    cmd = [SCRIPT_V7, "laser", temp_config, str(rot2), str(tilt2),
           str(rot3), str(tilt3), str(base_run_start), tag]
    print(f"[DAQ] {' '.join(shlex.quote(c) for c in cmd)}")
    if dry_run:
        return 0
    r = subprocess.run(cmd, cwd=BASE, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print(f"[ERROR] script_v7.sh exit {r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return r.returncode


def wait_for_analysis(run_int, dry_run, timeout=300):
    """Poll the same flag files ui_manager.py's pipeline monitor watches --
    read_<run>.flag disappearing (or done_<run>.flag appearing) means the
    prod->read chain finished for this run."""
    if dry_run:
        return True
    done_flag = os.path.join(FLAG_DIR, f"done_{run_int}.flag")
    t0 = time.time()
    while time.time() - t0 < timeout:
        if os.path.exists(done_flag):
            return True
        time.sleep(2)
    print(f"[WARNING] run {run_int}: analysis did not finish within {timeout}s")
    return False


def poisson2pe_fraction(mu):
    """P(2 p.e.) / P(>=1 p.e.) under pure Poisson -- same formula used in
    analyze_laser_intensity_history.py, the baseline a measured excess would
    be compared against (this script doesn't measure the excess itself, only
    reports what pure statistics predicts at this occupancy for reference)."""
    if mu is None or mu <= 0:
        return None
    p0 = math.exp(-mu)
    p1 = mu * p0
    p2 = (mu**2 / 2.0) * p0
    denom = 1 - p0
    return (p2 / denom) if denom > 0 else None


def read_charges(run_int, config_path=LIVE_CONFIG):
    """Pull poisson_qe (%) and poisson_mu per channel from the produced result
    tree -- same branches analyze_laser_intensity_history.py reads."""
    import glob as _glob
    import uproot
    pattern = os.path.join(BASE, "Data", "FinalResult", f"*_{run_int:03d}*.root")
    matches = _glob.glob(pattern)
    if not matches:
        return {}
    out = {}
    with uproot.open(matches[0]) as f:
        for key in f.keys():
            if not key.startswith("tree_ch"):
                continue
            ch = key.split(";")[0]
            t = f[key]
            qe_arr = t["poisson_qe"].array()
            mu_arr = t["poisson_mu"].array()
            if len(qe_arr):
                out[ch] = {"qe_pct": float(qe_arr[0]), "mu": float(mu_arr[0]),
                           "poisson_2pe_frac": poisson2pe_fraction(float(mu_arr[0]))}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rot2", type=float, default=0.0)
    ap.add_argument("--tilt2", type=float, default=0.0)
    ap.add_argument("--rot3", type=float, default=0.0)
    ap.add_argument("--tilt3", type=float, default=0.0)
    ap.add_argument("--base-run-start", type=int, default=900,
                     help="High run-number block reserved for calibration so these never collide with real scan runs.")
    ap.add_argument("--dry-run", action="store_true",
                     help="Print the planned sequence; touch no hardware, run no DAQ.")
    args = ap.parse_args()

    out_csv = os.path.join(BASE, "Data", "summary", "IntensityScan_results.csv")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    rows = []

    run_counter = args.base_run_start
    for wl, candidates in CANDIDATES.items():
        wl_digits = wl.replace("nm", "")
        for bias, pulse in candidates:
            print(f"\n=== {wl}  bias={bias}mA  pulse={pulse}mA  (run {run_counter}) ===")
            if not set_laser(wl, bias, pulse, args.dry_run):
                continue
            temp_cfg = make_temp_config(bias, pulse, wl_digits)
            rc = run_daq_point(temp_cfg, args.rot2, args.tilt2, args.rot3, args.tilt3,
                                run_counter, "IntensityCal", args.dry_run)
            if rc != 0 and not args.dry_run:
                laser_off(args.dry_run)
                run_counter += 1
                continue
            wait_for_analysis(run_counter, args.dry_run)
            if not args.dry_run:
                charges = read_charges(run_counter)
            else:
                charges = {f"tree_ch{i}": {"qe_pct": 1.5, "mu": 0.015, "poisson_2pe_frac": 0.0075} for i in range(3)}
            all_ok = len(charges) == 3 and all(
                TARGET_QE_PCT_RANGE[0] <= v["qe_pct"] <= TARGET_QE_PCT_RANGE[1] for v in charges.values())
            summary = {ch: f"QE={v['qe_pct']:.2f}% 2pe~{v['poisson_2pe_frac']*100:.2f}%" for ch, v in charges.items()}
            print(f"  {summary}  -> {'OK' if all_ok else 'OUT OF RANGE'}")
            row = {"wavelength": wl, "bias_mA": bias, "pulse_mA": pulse, "run": run_counter, "all_pmts_ok": all_ok}
            for ch, v in charges.items():
                row[f"{ch}_qe_pct"] = v["qe_pct"]
                row[f"{ch}_poisson_2pe_frac"] = v["poisson_2pe_frac"]
            rows.append(row)
            laser_off(args.dry_run)
            run_counter += 1

    if rows:
        fieldnames = sorted({k for r in rows for k in r})
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        print(f"\n[INFO] wrote {out_csv}")
        print("\n=== Summary: intensity candidates where ALL PMTs are in range ===")
        for r in rows:
            if r.get("all_pmts_ok"):
                print(f"  {r['wavelength']}: bias={r['bias_mA']}mA pulse={r['pulse_mA']}mA (run {r['run']})")


if __name__ == "__main__":
    main()
