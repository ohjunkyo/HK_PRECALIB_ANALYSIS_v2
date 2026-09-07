#!/bin/bash
# FULL campaign reprocess from RAW after Analysis_Threshold_mV 3->1.5mV
# change in prod_ntp_v7.C (2026-09-01). Runs prod_ntp_v7.C (RAW->production)
# THEN read_ntp_v7.C (production->FinalResult) for every run below.
# 20260815_0-45 already done via the diagnostic test run and copied into
# Data/production//Data/FinalResult directly -- skipped here.
set -u

RUN_CPP="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/run_cpp_script_v2.sh"
PROD_C="/home/precalkor/ADC/ADC_test/prod_ntp_v7.C"
READ_C="/home/precalkor/ADC/ADC_test/read_ntp_v7.C"
CONFIG="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"

LOCAL_RAW_DIR="/home/precalkor/ADC/ADC_test/Data/RAW/Laser"
EXT_RAW_DIR="/media/precalkor/HD-EDS-E/Data_Backup/RAW/Laser"
PROD_DIR="/home/precalkor/ADC/ADC_test/Data/production"

LOG="/home/precalkor/ADC/ADC_test/logs/reprocess_full_20260901_thr1p5_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# date lo hi
BLOCKS="
20260813 800 899
20260814 800 899
20260816 0 45
20260816 100 145
20260817 0 45
20260826 0 4
20260827 0 45
20260827 100 145
20260828 0 45
20260828 100 145
20260829 0 45
20260829 100 145
20260829 200 245
20260830 0 45
20260831 0 46
20260831 100 145
"

total=0; done_count=0; missing=0
FINAL_DIR="/home/precalkor/ADC/ADC_test/Data/FinalResult"
# Resume must NOT trust "file exists" -- almost every run in this campaign
# already had a FinalResult from months of normal (3mV) operation, long
# before Analysis_Threshold_mV was changed to 1.5 (found 2026-09-01: a run
# resumed this way silently skipped ~400 stale-3mV files and finished in 5
# minutes claiming "done"). Only skip a run if its FinalResult is NEWER than
# the threshold-change commit -- that's the only reliable signal it was
# actually reprocessed under the new threshold.
THRESHOLD_CHANGE_EPOCH=$(date -d "2026-09-01 14:50:35" +%s)

# Self-throttle: never touch the External-HDD RAW file while execute_DAQ_v2
# is actually mid-acquisition (that's the USB-bus contention that caused the
# 2026-09-01 motor-comm-timeout incident). Polls every 3s and waits out the
# whole acquisition -- General Scan points run for several minutes each, so
# a fixed sleep would either be too short (still racing) or too long
# (wastes idle time in the motor-move gaps between points).
wait_for_daq_idle() {
    while pgrep -x execute_DAQ_v2 > /dev/null 2>&1; do
        sleep 3
    done
    sleep 2   # small buffer after DAQ exits before we start our own read
}

echo "$BLOCKS" | while read -r date lo hi; do
    [ -z "$date" ] && continue
    log "=== Block ${date}_${lo}~${hi} ==="
    for i in $(seq "$lo" "$hi"); do
        run=$(printf '%03d' "$i")
        final="${FINAL_DIR}/precal_result_kor_run_${date}_${run}.root"
        if [ -f "$final" ]; then
            final_mtime=$(stat -c %Y "$final" 2>/dev/null || echo 0)
            if [ "$final_mtime" -ge "$THRESHOLD_CHANGE_EPOCH" ]; then
                continue   # genuinely reprocessed under the new (1.5mV) threshold -- skip
            fi
        fi
        raw="${LOCAL_RAW_DIR}/precal_raw_kor_run_${date}_${run}.root"
        if [ ! -f "$raw" ]; then
            ext="${EXT_RAW_DIR}/precal_raw_kor_run_${date}_${run}.root"
            if [ -f "$ext" ]; then
                raw="$ext"
            else
                log "  MISSING RAW for ${date}_${run} (checked local + external) -- skipping."
                continue
            fi
        fi
        wait_for_daq_idle
        log "  prod: ${date}_${run}"
        bash "$RUN_CPP" "$PROD_C" "$CONFIG" "$i" "\"$raw\"" >> "$LOG" 2>&1
        prod="${PROD_DIR}/precal_prd_kor_run_${date}_${run}.root"
        if [ ! -f "$prod" ]; then
            log "  MISSING production output for ${date}_${run} -- skipping read."
            continue
        fi
        log "  read: ${date}_${run}"
        bash "$RUN_CPP" "$READ_C" "$CONFIG" "$i" "\"$prod\"" >> "$LOG" 2>&1
    done
done

log "===== reprocess_full (1.5mV) finished ====="

cd /home/precalkor/ADC/ADC_test
log "Running Draw_Overlay_Uniformity_v7.C ..."
root -l -b -q Draw_Overlay_Uniformity_v7.C >> "$LOG" 2>&1
log "Running Draw_HVScan_Report.C ..."
root -l -b -q Draw_HVScan_Report.C >> "$LOG" 2>&1

log "===== ALL DONE ====="
log "Log: $LOG"
