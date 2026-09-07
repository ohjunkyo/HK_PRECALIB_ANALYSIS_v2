#!/bin/bash
# Re-run read_ntp_v7.C for the 2026-08-31 Timing Cut re-centering fix
# (ch1/ch2 now re-center tLow/tHigh on their own meanFit, matching ch0's
# pre-existing behavior -- see project memory "pending-adaptive-threshold-fix").
# prod_ntp_v7.C is untouched, so all existing Data/production/*.root files
# are reused -- read-only re-run, no RAW needed for any block here.
#
# Scope (user-approved 2026-08-31, "8.13~HV Scan까지 그리고 지금 받고 있는것"):
#   - Same scope as the earlier SPE-error-fix reprocess (HV campaign 8/27-8/30,
#     Stability 8/13-14, recent General Scans 8/26/8/17/8/16x2)
#   - PLUS today's live General Scan, 20260831 block 000-039 (run 040 still
#     being written as of launch -- excluded; 800-range runs are laser-calib,
#     not scan angle points -- excluded)
#
# Run with low CPU/IO priority (nice/ionice) since today's General Scan is
# still live (currently on run 040) -- user confirmed OK to run concurrently.
set -u

RUN_CPP="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/run_cpp_script_v2.sh"
READ_C="/home/precalkor/ADC/ADC_test/read_ntp_v7.C"
CONFIG="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"
PROD_DIR="/home/precalkor/ADC/ADC_test/Data/production"

LOG="/home/precalkor/ADC/ADC_test/logs/reprocess_read_20260831_timingfix_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# date lo hi
BLOCKS="
20260827 0 45
20260827 100 145
20260828 0 45
20260828 100 145
20260829 0 45
20260829 100 145
20260829 200 245
20260830 0 45
20260815 0 45
20260813 800 899
20260814 800 899
20260826 0 4
20260817 0 45
20260816 0 45
20260816 100 145
20260831 0 39
"

echo "$BLOCKS" | while read -r date lo hi; do
    [ -z "$date" ] && continue
    log "=== Block ${date}_${lo}~${hi} ==="
    for i in $(seq "$lo" "$hi"); do
        run=$(printf '%03d' "$i")
        prod="${PROD_DIR}/precal_prd_kor_run_${date}_${run}.root"
        if [ ! -f "$prod" ]; then
            log "  MISSING production for ${date}_${run} -- skipping."
            continue
        fi
        bash "$RUN_CPP" "$READ_C" "$CONFIG" "$i" "\"$prod\"" >> "$LOG" 2>&1
        log "  read: ${date}_${run} done"
    done
done

log "===== reprocess_read (timing fix) finished ====="

cd /home/precalkor/ADC/ADC_test
log "Running Draw_Overlay_Uniformity_v7.C ..."
root -l -b -q Draw_Overlay_Uniformity_v7.C >> "$LOG" 2>&1
log "Running Draw_HVScan_Report.C ..."
root -l -b -q Draw_HVScan_Report.C >> "$LOG" 2>&1

log "===== ALL DONE ====="
log "Log: $LOG"
