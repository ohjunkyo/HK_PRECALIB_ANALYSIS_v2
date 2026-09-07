#!/usr/bin/env python3
"""Compares the calc_bfield simulation output (run with today's real synced
coil currents) against the 4 real magnetometer readings from log_bfield_*.csv.

NOTE on the current mapping: the real hardware log only has 4 channels
(I1..I4) while the simulation wants 7 (X1,X2,X3,Y1,Y2,Z1,Z2). This script
assumes X1=X2=X3=I1, Y1=Y2=I2, Z1=I3, Z2=I4 (i.e. each axis's coils are
driven in series from one supply) -- this has NOT been confirmed against the
real wiring, so treat the absolute B-field values as indicative, not exact,
until that mapping is verified.
"""
import os
import uproot
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIM_ROOT = "/home/precalkor/ADC/ADC_test/bfield_sim/results/real_test.root"
REAL_CSV = "/home/precalkor/ADC/ADC_test/log_bfield_20260627.csv"
OUT_PATH = "/home/precalkor/ADC/ADC_test/bfield_sim/sim_vs_real.png"

with uproot.open(SIM_ROOT) as f:
    t = f["tree"]
    pm = t["pm"].array(library="np")
    Bmag = t["Bmag"].array(library="np")
    Bx = t["Bx"].array(library="np")
    By = t["By"].array(library="np")
    Bz = t["Bz"].array(library="np")

real = pd.read_csv(REAL_CSV).iloc[-1]
real_currents = {f"I{n}": float(real[f"I{n}"]) for n in (1, 2, 3, 4)}
real_mG = {}
for n in (1, 2, 3, 4):
    xyz = [float(str(real[f"M{n}_{ax}"])) if str(real[f"M{n}_{ax}"]) != "N/A" else None for ax in "XYZ"]
    real_mG[f"M{n}"] = None if any(v is None for v in xyz) else (sum(v * v for v in xyz) ** 0.5) * 1000.0

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

ax = axes[0]
labels = [f"PMT {i+1}" for i in pm]
ax.bar(labels, Bmag, color="#4c78a8")
ax.set_ylabel("Simulated |B| [mG]", fontsize=12)
ax.set_title(f"calc_bfield: |B| at 6 PMT positions\n(currents: I1={real_currents['I1']:.2f}A I2={real_currents['I2']:.2f}A "
             f"I3={real_currents['I3']:.2f}A I4={real_currents['I4']:.2f}A)", fontsize=10)
ax.tick_params(labelsize=10)
ax.grid(alpha=0.3, axis="y")
for i, v in enumerate(Bmag):
    ax.text(i, v + 0.3, f"{v:.1f}", ha="center", fontsize=9)

ax = axes[1]
mnames = [k for k, v in real_mG.items() if v is not None]
mvals = [real_mG[k] for k in mnames]
ax.bar(mnames, mvals, color="#e45756")
ax.set_ylabel("Measured |B| [mG] (magnetometers)", fontsize=12)
ax.set_title(f"Real log_bfield_*.csv, last sample\n{real['timestamp']}", fontsize=10)
ax.tick_params(labelsize=10)
ax.grid(alpha=0.3, axis="y")
for i, v in enumerate(mvals):
    ax.text(i, v + max(mvals) * 0.02, f"{v:.1f}", ha="center", fontsize=9)

fig.suptitle("Simulated field at PMT positions (left) vs. real magnetometer readings (right)\n"
             "-- NOT directly comparable: magnetometers M1-4 are not confirmed to sit at the PMT positions", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT_PATH, dpi=150)
print(f"[INFO] wrote {OUT_PATH}")
print("\nSimulated per-PMT field:")
for i, bm, bx, by, bz in zip(pm, Bmag, Bx, By, Bz):
    print(f"  PMT {i+1}: Bx={bx:+.2f} By={by:+.2f} Bz={bz:+.2f} |B|={bm:.2f} mG")
print("\nReal magnetometers:", real_mG)
