#!/bin/bash
# Reprocess the 2026-08-29 MORNING General Scan (runs 000-045, 07:22~11:27,
# HV2=1579 / HV3=1650, "HV scan, recom. value - 100V").
#
# Why this exists: on 2026-08-29 three General Scans ran, and the run-number
# allocator handed out 000-045 TWICE -- once to the morning scan and again to
# the evening scan (17:44~21:50, HV2=1529). The evening scan therefore
# overwrote the morning scan's FinalResult files. The morning RAW files
# survived only on the external HDD, because a manual backup had already
# reclaimed them from local before the collision happened.
#
# Recovery (2026-08-30): the evening scan was renamed to 200-245 (RAW and
# FinalResult), the morning RAW was copied back from external to local as
# 000-045, and this script regenerates the morning FinalResult files that
# were lost.
#
# Runs the same three-stage chain script_v7.sh uses per point:
#   prod_ntp_v7.C  ->  read_ntp_v7.C  ->  Draw_Contour_v3.C
#
# Safe to re-run: any run that already has a FinalResult is skipped, so an
# interrupted pass can simply be started again.
set -u

RUN_CPP="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/run_cpp_script_v2.sh"
PROD_C="/home/precalkor/ADC/ADC_test/prod_ntp_v7.C"
READ_C="/home/precalkor/ADC/ADC_test/read_ntp_v7.C"
CONTOUR_C="/home/precalkor/ADC/ADC_test/Draw_Contour_v3.C"
CONFIG="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"

RAW_DIR="/home/precalkor/ADC/ADC_test/Data/RAW/Laser"
PROD_DIR="/home/precalkor/ADC/ADC_test/Data/production"
RESULT_DIR="/home/precalkor/ADC/ADC_test/Data/FinalResult"

DATE_TAG="20260829"
RUN_LO=0
RUN_HI=45

# The morning scan's own HV, used to verify we are reprocessing the right
# files and not the evening scan by mistake.
EXPECT_HV2=1579

LOG="/home/precalkor/ADC/ADC_test/logs/reprocess_${DATE_TAG}_morning_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "===== reprocess ${DATE_TAG} morning (runs ${RUN_LO}-${RUN_HI}) started ====="

processed=0
skipped_done=0
failed=0
wrong_hv=0

for i in $(seq "$RUN_LO" "$RUN_HI"); do
    run=$(printf '%03d' "$i")
    raw="${RAW_DIR}/precal_raw_kor_run_${DATE_TAG}_${run}.root"
    prd="${PROD_DIR}/precal_prd_kor_run_${DATE_TAG}_${run}.root"
    result="${RESULT_DIR}/precal_result_kor_run_${DATE_TAG}_${run}.root"

    if [ ! -f "$raw" ]; then
        log "  MISSING RAW for run ${run} -- skipping."
        failed=$((failed + 1))
        continue
    fi

    if [ -f "$result" ]; then
        log "  run ${run}: FinalResult already exists -- skipping."
        skipped_done=$((skipped_done + 1))
        continue
    fi

    # Guard: confirm this RAW really is the morning scan before spending
    # ~minutes analysing it. A mismatch means the wrong file is in place and
    # the operator needs to look before anything else runs.
    hv2=$(python3 -c "
import uproot,sys
try:
    print(int(uproot.open('$raw')['RunInfo']['HV2'].array()[0]))
except Exception:
    print(-1)
" 2>/dev/null)
    if [ "$hv2" != "$EXPECT_HV2" ]; then
        log "  run ${run}: HV2=${hv2} but expected ${EXPECT_HV2} -- NOT the morning scan. Skipping."
        wrong_hv=$((wrong_hv + 1))
        continue
    fi

    log "Processing run ${DATE_TAG}_${run} (HV2=${hv2})"
    bash "$RUN_CPP" "$PROD_C"    "$CONFIG" "$i" "\"$raw\"" >> "$LOG" 2>&1
    bash "$RUN_CPP" "$READ_C"    "$CONFIG" "$i" "\"$prd\"" >> "$LOG" 2>&1
    bash "$RUN_CPP" "$CONTOUR_C" "$CONFIG" "$i" "\"$raw\"" >> "$LOG" 2>&1

    if [ -f "$result" ]; then
        log "  -> DONE: FinalResult produced for run ${run}."
        processed=$((processed + 1))
    else
        log "  -> WARN: FinalResult still missing for run ${run} -- check log above."
        failed=$((failed + 1))
    fi
done

log "===== finished: processed=${processed} skipped(already done)=${skipped_done} wrong_hv=${wrong_hv} failed=${failed} ====="
log "Log written to ${LOG}"
