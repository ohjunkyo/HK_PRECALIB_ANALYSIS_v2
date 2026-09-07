#!/bin/bash
# Diagnostic-only reprocess of the 8/15 Recom. block (46 runs) with the
# shared Analysis_Threshold_mV constant set to 1.5mV instead of 3mV --
# affects BOTH the amplitude cut AND the diff/timing-crossing calculation
# (see project memory "pending-adaptive-threshold-fix", 2026-09-01 finding).
# Writes to separate Data/production_thr1p5/ and Data/FinalResult_thr1p5/
# directories -- does NOT touch real production/FinalResult.
set -u

RUN_CPP="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/run_cpp_script_v2.sh"
PROD_C="/home/precalkor/ADC/ADC_test/prod_ntp_v7_thr1p5.C"
READ_C="/home/precalkor/ADC/ADC_test/read_ntp_v7_thr1p5.C"
CONFIG="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"
EXT_RAW_DIR="/media/precalkor/HD-EDS-E/Data_Backup/RAW/Laser"

LOG="/home/precalkor/ADC/ADC_test/logs/thr1p5_test_20260815_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

for i in $(seq 0 45); do
    run=$(printf '%03d' "$i")
    raw="${EXT_RAW_DIR}/precal_raw_kor_run_20260815_${run}.root"
    if [ ! -f "$raw" ]; then log "  MISSING RAW for $run -- skip"; continue; fi
    log "prod(1.5mV): 20260815_${run}"
    bash "$RUN_CPP" "$PROD_C" "$CONFIG" "$i" "\"$raw\"" >> "$LOG" 2>&1
    prod="./Data/production_thr1p5/precal_prd_kor_run_20260815_${run}.root"
    if [ ! -f "$prod" ]; then log "  MISSING prod output for $run -- skip read"; continue; fi
    log "read(1.5mV): 20260815_${run}"
    bash "$RUN_CPP" "$READ_C" "$CONFIG" "$i" "\"$prod\"" >> "$LOG" 2>&1
done

log "===== DONE ====="
log "Log: $LOG"
