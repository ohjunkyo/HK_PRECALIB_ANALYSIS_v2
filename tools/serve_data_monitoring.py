#!/usr/bin/env python3
"""Live DAQ data-monitoring dashboard.

Serves an interactive table of every production run's metadata (RunInfo +
dark-rate QC + which Uniformity tag covers it), rescanning in the background so
it stays current with no manual step.

Run:      python3 serve_data_monitoring.py [port]
Managed by the systemd --user unit daq-data-monitoring.service.
Open http://<machine-IP>:<port>/  from any device that can reach this host.

Search / sort / filter all run client-side against a JSON snapshot (/api/runs),
so interacting with the table never waits on the server.
"""
import bisect
import csv
import datetime
import glob
import json
import os
import re
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import uproot

BASE = "/home/precalkor/ADC/ADC_test"
PROD_GLOB = os.path.join(BASE, "Data", "production", "precal_prd_kor_run_*.root")
RAW_GLOB = os.path.join(BASE, "Data", "RAW", "Laser", "precal_raw_kor_run_*.root")
FINAL_GLOB = os.path.join(BASE, "Data", "FinalResult", "precal_result_kor_run_*.root")
UNI_GLOB = os.path.join(BASE, "Data", "UNIFORMITY", "Graphs_Uniformity_*.root")
TAKINGLOG_GLOB = os.path.join(BASE, "LOG", "DAQ", "TakingLog_*.txt")
BFIELD_GLOB = os.path.join(BASE, "log_bfield_*.csv")
SCANMAP_GLOB = "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/LOG/ScanHistory/scanmap_*.json"

# A contiguous run block touching at least this many distinct tilt angles is a
# General Scan sweep; manual work sits at one or two fixed angles.
SCAN_MIN_DISTINCT_TILTS = 6

FNAME_RE = re.compile(r"precal_prd_kor_run_(\d{8})_(\d+)\.root$")
RUN_RE = re.compile(r"_(\d{8})_(\d+)\.root$")
# Uniformity tags: "<date>_<runStart>_<runEnd>" (range) or bare "<date>".
UNI_RANGE_RE = re.compile(r"^(\d{8})_(\d+)_(\d+)$")
# Acquisition timestamps come from either line the DAQ writes per run:
#   "... | [INFO] Generated File: precal_raw_kor_run_20260731_321.root"
#   "... | Launching background analysis for precal_raw_kor_run_20260731_825.root..."
# Some runs only ever get the second one (e.g. 20260731_825/827), so matching
# just "Generated File" silently lost their timestamps -- and with them the
# B-field lookup. Both lines are written as the run finishes, so either is a
# valid acquisition time for a 10-min-sampled magnet log.
TAKING_RE = re.compile(
    r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}).*?"
    r"(?:Generated File:|Launching background analysis for)\s*"
    r"precal_raw_kor_run_(\d{8})_(\d+)\.root")

# The B-field logger samples every ~10 min, so a run (~5 min) is matched to the
# nearest sample. Beyond this the reading is too far away to attribute to the
# run, and the B-field columns are left blank rather than guessed.
BFIELD_MATCH_TOLERANCE_S = 15 * 60
# Total coil current above which the magnet counts as energized (same 0.4 A
# threshold DataConditionBfieldStatus() uses in Draw_Overlay_Uniformity_v7.C).
BFIELD_ON_AMPS = 0.4

RESCAN_INTERVAL_S = 20
# A production file still being written has a partial/empty RunInfo. Skip
# anything touched this recently and retry it on the next sweep, so a run that
# is mid-write never lands in the table as a blank row.
MIN_FILE_AGE_S = 15
DEFAULT_PORT = 8090

FIELDS = ["SN1", "HV1", "Direction1", "SN2", "HV2", "Direction2",
          "RawRotateAngle2", "RawTiltAngle2", "SN3", "HV3", "Direction3",
          "RawRotateAngle3", "RawTiltAngle3", "Wavelength", "Laser_mA",
          "NOTE", "RunMode", "Expert", "Shifter"]

_lock = threading.Lock()
_rows = {}              # (date, run) -> row dict
_scanned = set()        # basenames successfully read
_aux = {"raw": set(), "final": set(), "uni": [], "taken": {}, "bf_t": [], "bf_v": [],
        "scanmap": set()}
_scan_info = {"time": None, "new": 0, "total": 0, "skipped_young": 0, "error": None}


def _runinfo(f):
    out = {}
    try:
        ri = f["RunInfo"]
        names = set(ri.keys())
    except Exception:
        return None
    for field in FIELDS:
        if field not in names:
            out[field] = None
            continue
        try:
            arr = ri[field].array(library="np")
            out[field] = arr[0] if len(arr) else None
        except Exception:
            out[field] = None
    return out


def _param(f, name):
    try:
        return f[name].member("fVal")
    except Exception:
        return None


def _jsonable(v):
    if v is None:
        return None
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if hasattr(v, "item"):          # numpy scalar
        try:
            return v.item()
        except Exception:
            return str(v)
    return v


def scan_aux():
    """Cheap filesystem-only side scans: which runs have RAW / FinalResult, and
    the list of Uniformity tag ranges (for the per-run 'covered by tag' column)."""
    raw, final = set(), set()
    for pattern, dest in ((RAW_GLOB, raw), (FINAL_GLOB, final)):
        for p in glob.glob(pattern):
            m = RUN_RE.search(os.path.basename(p))
            if m:
                dest.add((m.group(1), int(m.group(2))))

    uni = []
    for p in glob.glob(UNI_GLOB):
        tag = os.path.basename(p)[len("Graphs_Uniformity_"):-len(".root")]
        m = UNI_RANGE_RE.match(tag)
        if m:
            uni.append((m.group(1), int(m.group(2)), int(m.group(3)), tag))
        else:
            uni.append((tag[:8] if tag[:8].isdigit() else None, None, None, tag))

    taken = load_taking_log()
    bf_t, bf_v = load_bfield()
    scanmap = load_scanmap_runs()
    with _lock:
        _aux["raw"], _aux["final"], _aux["uni"] = raw, final, uni
        _aux["taken"], _aux["bf_t"], _aux["bf_v"] = taken, bf_t, bf_v
        _aux["scanmap"] = scanmap


def load_scanmap_runs():
    """(date, run) pairs the scan manager recorded as General-Scan points.

    NOT a complete history: the map is keyed by (axis, tilt, wavelength), so a
    second scan on the same date overwrites the first one's entries. Presence
    therefore proves 'General Scan', but absence proves nothing -- which is why
    classify_scan_type() falls back to the tilt-sweep shape below."""
    runs = set()
    for path in glob.glob(SCANMAP_GLOB):
        try:
            with open(path, errors="replace") as f:
                data = json.load(f)
        except Exception:
            continue
        for entry in (data or {}).values():
            m = RUN_RE.search(os.path.basename(entry.get("file", "") or ""))
            if m:
                runs.add((m.group(1), int(m.group(2))))
    return runs


def classify_scan_type(rows_by_key):
    """-> {(date, run): 'General Scan' | 'Manual'}

    A General Scan steps TILT systematically across its range, so a contiguous
    block of runs from one belongs to many distinct tilts. Manual runs sit at
    one or two fixed angles (e.g. the 800-829 motor test: 30 runs, 2 tilts).
    Counting distinct tilts per contiguous run block separates the two cleanly,
    and unlike scanmap it works for every date in the archive."""
    out = {}
    by_date = {}
    for (date, run) in rows_by_key:
        by_date.setdefault(date, []).append(run)

    for date, runs in by_date.items():
        runs.sort()
        block = [runs[0]]
        blocks = []
        for prev, cur in zip(runs, runs[1:]):
            # <=2 tolerates a skipped/failed run inside one scan block
            (block.append(cur) if cur - prev <= 2
             else (blocks.append(block), block := [cur]))
        blocks.append(block)

        for block in blocks:
            tilts = {rows_by_key[(date, r)].get("Tilt2") for r in block}
            tilts.discard(None)
            is_scan = len(tilts) >= SCAN_MIN_DISTINCT_TILTS
            for r in block:
                # scanmap presence is authoritative and overrides the shape test
                definite = (date, r) in _aux.get("scanmap", set())
                out[(date, r)] = "General Scan" if (definite or is_scan) else "Manual"
    return out


def load_taking_log():
    """(date, run) -> datetime the RAW file was produced, from LOG/DAQ/TakingLog_*.

    This is the ACQUISITION time. A production file's mtime is only when it was
    processed (reprocessing rewrites it), so it can't be used to look up what
    the magnet was doing while the run was actually being taken."""
    taken = {}
    for path in glob.glob(TAKINGLOG_GLOB):
        try:
            with open(path, errors="replace") as f:
                for line in f:
                    m = TAKING_RE.match(line.strip())
                    if not m:
                        continue
                    try:
                        ts = datetime.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M")
                    except ValueError:
                        continue
                    taken[(m.group(2), int(m.group(3)))] = ts   # last entry wins (re-takes)
        except OSError:
            continue
    return taken


def load_bfield():
    """Whole B-field log as two parallel, time-sorted lists (timestamps, values)
    so each run can be matched with a bisect instead of a linear scan."""
    rows = []
    for path in glob.glob(BFIELD_GLOB):
        try:
            with open(path, newline="", errors="replace") as f:
                for rec in csv.DictReader(f):
                    raw_ts = (rec.get("timestamp") or "").strip()
                    if not raw_ts:
                        continue
                    try:
                        ts = datetime.datetime.fromisoformat(raw_ts)
                    except ValueError:
                        continue
                    rows.append((ts, rec))
        except OSError:
            continue
    rows.sort(key=lambda r: r[0])
    return [r[0] for r in rows], [r[1] for r in rows]


def _fnum(rec, key):
    """CSV cell -> float, treating 'N/A'/blank as missing (a disconnected probe
    writes N/A; coercing that to 0.0 would fake a real zero-field reading)."""
    v = (rec.get(key) or "").strip()
    if not v or v.upper() == "N/A":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def bfield_for(run_time):
    """B-field/current readings nearest to run_time, or {} when the log has no
    sample close enough (logger stopped, or a date range it never covered)."""
    times, vals = _aux["bf_t"], _aux["bf_v"]
    if not times or run_time is None:
        return {}
    i = bisect.bisect_left(times, run_time)
    best, best_dt = None, None
    for j in (i - 1, i):
        if 0 <= j < len(times):
            dt = abs((times[j] - run_time).total_seconds())
            if best_dt is None or dt < best_dt:
                best, best_dt = j, dt
    if best is None or best_dt > BFIELD_MATCH_TOLERANCE_S:
        return {}

    rec = vals[best]
    out = {"B_MatchMin": round(best_dt / 60.0, 1)}
    currents = []
    for n in (1, 2, 3, 4):
        amp = _fnum(rec, f"I{n}")
        out[f"Coil_I{n}_A"] = None if amp is None else round(amp, 3)
        if amp is not None:
            currents.append(abs(amp))
    out["Coil_Isum_A"] = round(sum(currents), 3) if currents else None
    out["Bfield"] = ("" if not currents else
                     ("ON" if sum(currents) > BFIELD_ON_AMPS else "OFF"))
    for n in (1, 2, 3, 4):
        xyz = [_fnum(rec, f"M{n}_{ax}") for ax in ("X", "Y", "Z")]
        out[f"B_M{n}_mG"] = (None if any(c is None for c in xyz) else
                             round((sum(c * c for c in xyz) ** 0.5) * 1000.0, 1))
    return out


def uni_tags_for(date, run):
    tags = [t for (d, a, b, t) in _aux["uni"]
            if d == date and a is not None and a <= run <= b]
    return ", ".join(sorted(tags))


def scan_one(path):
    base = os.path.basename(path)
    m = FNAME_RE.match(base)
    if not m:
        return None
    date, run = m.group(1), int(m.group(2))
    try:
        f = uproot.open(path)
    except Exception:
        return None                     # unreadable (mid-write) -> retry later
    info = _runinfo(f)
    if info is None or info.get("SN1") is None:
        return None                     # RunInfo not written yet -> retry later
    g = lambda k: _jsonable(info.get(k))
    return {
        "Date": date, "Run": run,
        "RunMode": g("RunMode"),
        "Wavelength": g("Wavelength"), "Laser_mA": g("Laser_mA"),
        "SN1": g("SN1"), "HV1": g("HV1"), "Dir1": g("Direction1"),
        "SN2": g("SN2"), "HV2": g("HV2"), "Dir2": g("Direction2"),
        "Rot2": g("RawRotateAngle2"), "Tilt2": g("RawTiltAngle2"),
        "SN3": g("SN3"), "HV3": g("HV3"), "Dir3": g("Direction3"),
        "Rot3": g("RawRotateAngle3"), "Tilt3": g("RawTiltAngle3"),
        "Shifter": g("Shifter"), "Expert": g("Expert"), "NOTE": g("NOTE"),
        "Dark0": _param(f, "NoiseCountRate_ch0"),
        "Dark1": _param(f, "NoiseCountRate_ch1"),
        "Dark2": _param(f, "NoiseCountRate_ch2"),
        "Processed": time.strftime("%Y-%m-%d %H:%M:%S",
                                    time.localtime(os.path.getmtime(path))),
    }


def rescan():
    try:
        scan_aux()
        now = time.time()
        new = skipped = 0
        for path in sorted(glob.glob(PROD_GLOB)):
            base = os.path.basename(path)
            if base in _scanned:
                continue
            try:
                if now - os.path.getmtime(path) < MIN_FILE_AGE_S:
                    skipped += 1        # still being written; pick it up next sweep
                    continue
            except OSError:
                continue
            row = scan_one(path)
            if row is None:
                continue                # left out of _scanned on purpose -> retried
            with _lock:
                _rows[(row["Date"], row["Run"])] = row
                _scanned.add(base)
            new += 1

        with _lock:
            scan_types = classify_scan_type(_rows)
            for key, row in _rows.items():
                row["ScanType"] = scan_types.get(key, "")
                row["HasRAW"] = "Y" if key in _aux["raw"] else ""
                row["HasResult"] = "Y" if key in _aux["final"] else ""
                taken = _aux["taken"].get(key)
                row["TakenAt"] = taken.strftime("%Y-%m-%d %H:%M") if taken else ""
                row.update({"Bfield": "", "Coil_Isum_A": None, "B_MatchMin": None,
                            **{f"Coil_I{n}_A": None for n in (1, 2, 3, 4)},
                            **{f"B_M{n}_mG": None for n in (1, 2, 3, 4)}})
                row.update(bfield_for(taken))
                row["UniTag"] = uni_tags_for(*key)
            _scan_info.update(time=time.strftime("%Y-%m-%d %H:%M:%S"), new=new,
                              total=len(_rows), skipped_young=skipped, error=None)
    except Exception as e:
        with _lock:
            _scan_info["error"] = str(e)


def loop():
    while True:
        rescan()
        time.sleep(RESCAN_INTERVAL_S)


PAGE = r"""<!doctype html><html><head><meta charset="utf-8">
<title>DAQ Data Monitoring</title>
<style>
:root{--bg:#0e1117;--panel:#161b26;--line:#2a3142;--fg:#e6e6e6;--dim:#9aa0a6;--acc:#4d9fff;}
*{box-sizing:border-box}
body{font-family:-apple-system,Segoe UI,Arial,sans-serif;margin:0;background:var(--bg);color:var(--fg)}
header{padding:14px 18px;border-bottom:1px solid var(--line);background:var(--panel);position:sticky;top:0;z-index:10}
h1{margin:0 0 8px;font-size:17px}
.meta{color:var(--dim);font-size:12px;margin-bottom:10px}
.controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center}
input,select{background:#0e1117;color:var(--fg);border:1px solid var(--line);border-radius:6px;
  padding:7px 10px;font-size:13px;font-family:inherit}
input:focus,select:focus{outline:none;border-color:var(--acc)}
#q{min-width:280px}
.count{color:var(--dim);font-size:12px;margin-left:auto}
button{background:#222a3a;color:var(--fg);border:1px solid var(--line);border-radius:6px;
  padding:7px 12px;font-size:13px;cursor:pointer}
button:hover{border-color:var(--acc)}
.wrap{overflow:auto;height:calc(100vh - 150px)}
table{border-collapse:collapse;width:100%;font-size:12px}
th,td{border-bottom:1px solid var(--line);padding:5px 9px;text-align:left;white-space:nowrap}
th{background:#1c2333;position:sticky;top:0;cursor:pointer;user-select:none;z-index:5}
th:hover{color:var(--acc)}
th .ar{color:var(--acc);font-size:10px}
tr:hover td{background:#222a3a}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.hot{color:#ff6b6b;font-weight:600}
mark{background:#4d9fff44;color:#cfe4ff;border-radius:2px}
.tag{color:#8bd17c}
</style></head><body>
<header>
  <h1>DAQ Production Run Monitoring</h1>
  <div class="meta" id="meta">loading…</div>
  <div class="controls">
    <input id="q" placeholder="Search anything (e.g. 135, EM5370, motor on, 405)  — space = AND">
    <select id="fScan"><option value="">Scan + Manual</option><option>General Scan</option><option>Manual</option></select>
    <select id="fMode"><option value="">All modes</option></select>
    <select id="fDate"><option value="">All dates</option></select>
    <select id="fWl"><option value="">All wavelengths</option></select>
    <select id="fPmt"><option value="">All PMTs</option></select>
    <select id="fBf"><option value="">B-field: any</option><option>ON</option><option>OFF</option></select>
    <select id="fUni"><option value="">All Uniformity tags</option></select>
    <button id="reset">Reset</button>
    <span class="count" id="count"></span>
  </div>
</header>
<div class="wrap"><table>
  <thead><tr id="hrow"></tr></thead>
  <tbody id="body"></tbody>
</table></div>
<script>
const COLS=[
 {k:"Date",t:"Date"},{k:"Run",t:"Run",n:1},{k:"ScanType",t:"Scan type"},{k:"RunMode",t:"Mode"},
 {k:"Wavelength",t:"λ(nm)",n:1},{k:"Laser_mA",t:"Laser(mA)",n:1},
 {k:"SN1",t:"Monitor SN"},{k:"HV1",t:"HV1",n:1},
 {k:"SN2",t:"SN2"},{k:"HV2",t:"HV2",n:1},{k:"Rot2",t:"Rot2",n:1},{k:"Tilt2",t:"Tilt2",n:1},
 {k:"SN3",t:"SN3"},{k:"HV3",t:"HV3",n:1},{k:"Rot3",t:"Rot3",n:1},{k:"Tilt3",t:"Tilt3",n:1},
 {k:"Shifter",t:"Shifter"},{k:"Expert",t:"Expert"},{k:"NOTE",t:"Note"},
 {k:"Dark0",t:"Dark0(Hz)",n:1},{k:"Dark1",t:"Dark1(Hz)",n:1},{k:"Dark2",t:"Dark2(Hz)",n:1},
 {k:"Bfield",t:"B-field"},{k:"Coil_Isum_A",t:"ΣI(A)",n:1},
 {k:"Coil_I1_A",t:"I1(A)",n:1},{k:"Coil_I2_A",t:"I2(A)",n:1},
 {k:"Coil_I3_A",t:"I3(A)",n:1},{k:"Coil_I4_A",t:"I4(A)",n:1},
 {k:"B_M1_mG",t:"|B|M1(mG)",n:1},{k:"B_M2_mG",t:"|B|M2(mG)",n:1},{k:"B_M3_mG",t:"|B|M3(mG)",n:1},
 {k:"B_MatchMin",t:"Blog Δmin",n:1},
 {k:"UniTag",t:"Uniformity tag"},{k:"HasRAW",t:"RAW"},{k:"HasResult",t:"Result"},
 {k:"TakenAt",t:"Taken at"},{k:"Processed",t:"Processed at"}];

let DATA=[], sortKey="__chrono", sortDir=-1;   // default: newest run first
const $=id=>document.getElementById(id);

function chrono(r){ return (r.Date||"")+String(r.Run??0).padStart(6,"0"); }

function opts(sel,vals,label){
  const cur=sel.value;
  sel.innerHTML=`<option value="">${label}</option>`+
    vals.map(v=>`<option${v===cur?" selected":""}>${v}</option>`).join("");
}

function buildFilters(){
  const uniq=k=>[...new Set(DATA.map(r=>r[k]).filter(v=>v!==null&&v!==""))].sort();
  opts($("fMode"),uniq("RunMode"),"All modes");
  opts($("fDate"),uniq("Date").reverse(),"All dates");
  opts($("fWl"),uniq("Wavelength"),"All wavelengths");
  opts($("fPmt"),[...new Set(DATA.flatMap(r=>[r.SN2,r.SN3]).filter(Boolean))].sort(),"All PMTs");
  const tags=[...new Set(DATA.flatMap(r=>(r.UniTag||"").split(", ")).filter(Boolean))].sort().reverse();
  opts($("fUni"),tags,"All Uniformity tags");
}

function rowText(r){ return COLS.map(c=>r[c.k]??"").join(" ").toLowerCase(); }

function apply(){
  const terms=$("q").value.toLowerCase().split(/\s+/).filter(Boolean);
  const fm=$("fMode").value, fd=$("fDate").value, fw=$("fWl").value,
        fp=$("fPmt").value, fu=$("fUni").value, fb=$("fBf").value, fs=$("fScan").value;
  let out=DATA.filter(r=>{
    if(fs&&r.ScanType!==fs) return false;
    if(fm&&r.RunMode!==fm) return false;
    if(fd&&r.Date!==fd) return false;
    if(fw&&String(r.Wavelength)!==fw) return false;
    if(fp&&r.SN2!==fp&&r.SN3!==fp) return false;
    if(fb&&r.Bfield!==fb) return false;
    if(fu&&!(r.UniTag||"").split(", ").includes(fu)) return false;
    if(terms.length){ const t=rowText(r); return terms.every(x=>t.includes(x)); }
    return true;
  });
  const col=COLS.find(c=>c.k===sortKey);
  out.sort((a,b)=>{
    let x,y;
    if(sortKey==="__chrono"){ x=chrono(a); y=chrono(b); }
    else if(col&&col.n){ x=a[sortKey]??-Infinity; y=b[sortKey]??-Infinity; }
    else { x=String(a[sortKey]??"").toLowerCase(); y=String(b[sortKey]??"").toLowerCase(); }
    return (x<y?-1:x>y?1:0)*sortDir;
  });
  render(out,terms);
}

function hl(v,terms){
  let s=String(v??"");
  if(!s||!terms.length) return s.replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));
  s=s.replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));
  for(const t of terms){
    if(!t) continue;
    s=s.replace(new RegExp("("+t.replace(/[.*+?^${}()|[\]\\]/g,"\\$&")+")","ig"),"<mark>$1</mark>");
  }
  return s;
}

function render(rows,terms){
  $("hrow").innerHTML=COLS.map(c=>{
    const ar=sortKey===c.k?(sortDir>0?" ▲":" ▼"):"";
    return `<th data-k="${c.k}">${c.t}<span class="ar">${ar}</span></th>`;}).join("");
  $("hrow").querySelectorAll("th").forEach(th=>th.onclick=()=>{
    const k=th.dataset.k;
    if(sortKey===k) sortDir*=-1; else { sortKey=k; sortDir=1; }
    apply();
  });
  const MAX=1500;
  $("body").innerHTML=rows.slice(0,MAX).map(r=>"<tr>"+COLS.map(c=>{
    let v=r[c.k];
    if(c.n&&typeof v==="number"){
      if(c.k.startsWith("Dark")||c.k.startsWith("B_M")) v=v.toFixed(0);
      else if(c.k.startsWith("Coil_")) v=v.toFixed(2);
    }
    const cls=[c.n?"num":"", (c.k.startsWith("Dark")&&typeof r[c.k]==="number"&&r[c.k]>30000)?"hot":"",
               c.k==="UniTag"?"tag":""].filter(Boolean).join(" ");
    return `<td class="${cls}">${hl(v,terms)}</td>`;}).join("")+"</tr>").join("");
  $("count").textContent=`${rows.length} rows${rows.length>MAX?` (showing ${MAX})`:""} of ${DATA.length}`;
}

async function refresh(){
  const r=await fetch("/api/runs");const j=await r.json();
  DATA=j.rows; buildFilters();
  $("meta").textContent=`Last scan: ${j.info.time} | new: ${j.info.new} | `+
    `mid-write skipped: ${j.info.skipped_young} | total runs: ${j.info.total} | `+
    `auto-refresh 20s | error: ${j.info.error}`;
  apply();
}
["q","fScan","fMode","fDate","fWl","fPmt","fBf","fUni"].forEach(id=>{
  $(id).addEventListener(id==="q"?"input":"change",apply);});
$("reset").onclick=()=>{["q","fScan","fMode","fDate","fWl","fPmt","fBf","fUni"].forEach(id=>$(id).value="");
  sortKey="__chrono";sortDir=-1;apply();};
refresh(); setInterval(refresh,20000);
</script></body></html>"""


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _send(self, body, ctype):
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path.startswith("/api/runs"):
            with _lock:
                payload = json.dumps({"rows": list(_rows.values()),
                                       "info": dict(_scan_info)}).encode()
            self._send(payload, "application/json")
        else:
            self._send(PAGE.encode("utf-8"), "text/html; charset=utf-8")


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PORT
    print(f"[INFO] Initial scan (~30 s for ~3000 files)…", flush=True)
    rescan()
    print(f"[INFO] Loaded {len(_rows)} runs.", flush=True)
    threading.Thread(target=loop, daemon=True).start()
    print(f"[INFO] Serving on http://0.0.0.0:{port}/", flush=True)
    ThreadingHTTPServer(("0.0.0.0", port), Handler).serve_forever()


if __name__ == "__main__":
    main()
