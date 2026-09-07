#!/bin/bash
# Waits for the currently-running (or an about-to-start scheduled) General
# Scan to genuinely finish, then runs the offline analysis chain
# (prod_ntp_v7 -> read_ntp_v7 -> Draw_Contour_v3) for every RAW run that
# doesn't already have a FinalResult -- catching up whatever the normal
# per-point tmux analysis chain in script_v7.sh missed.
#
# Why this exists: on 2026-08-12 the DAQ watchdog killed run 100 (fixed-size
# 600s timeout, too short for a 700000-event point) before the scan ever
# reached its per-point analysis launch step, and every run after it was
# never auto-analyzed either. Requested by the operator to run once, chained
# after the 2026-08-13 03:00 JST scheduled scan (see queued_schedules.json),
# via cron: `0 3 * * * /home/precalkor/ADC/ADC_test/auto_catchup_analysis.sh`
#
# General on purpose -- does NOT hardcode a run-number range, since the
# 03:00 schedule starts a fresh General Scan whose run numbers aren't known
# ahead of time. Instead it waits for genuine scan-idle, then processes
# every RAW file that has no matching FinalResult yet, regardless of date.
set -u

RUN_CPP="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/run_cpp_script_v2.sh"
PROD_C="/home/precalkor/ADC/ADC_test/prod_ntp_v7.C"
READ_C="/home/precalkor/ADC/ADC_test/read_ntp_v7.C"
CONTOUR_C="/home/precalkor/ADC/ADC_test/Draw_Contour_v3.C"
CONFIG="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"

RAW_DIR="/home/precalkor/ADC/ADC_test/Data/RAW"
PROD_DIR="/home/precalkor/ADC/ADC_test/Data/production"
RESULT_DIR="/home/precalkor/ADC/ADC_test/Data/FinalResult"

LOG="/home/precalkor/ADC/ADC_test/logs/auto_catchup_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "===== auto_catchup_analysis.sh started ====="

# ---- 1) Wait for the scan to actually be idle, not just between two points ----
# execute_DAQ_v2 is absent for a few seconds between every scan point (motor
# move + settle), so a single absence check would fire mid-scan. Require
# IDLE_STREAK_NEEDED consecutive absence checks, CHECK_INTERVAL_S apart, before
# concluding the scan is truly done -- comfortably longer than any normal
# inter-point gap.
CHECK_INTERVAL_S=60
IDLE_STREAK_NEEDED=5          # 5 * 60s = 5 minutes of sustained idle
MAX_WAIT_HOURS=14             # safety cap so this can't wait forever
max_checks=$(( MAX_WAIT_HOURS * 3600 / CHECK_INTERVAL_S ))

idle_streak=0
checks=0
while [ "$checks" -lt "$max_checks" ]; do
    if pgrep -x execute_DAQ_v2 > /dev/null 2>&1; then
        if [ "$idle_streak" -gt 0 ]; then
            log "DAQ active again (idle streak reset)."
        fi
        idle_streak=0
    else
        idle_streak=$((idle_streak + 1))
        log "No execute_DAQ_v2 running (idle streak ${idle_streak}/${IDLE_STREAK_NEEDED})."
        if [ "$idle_streak" -ge "$IDLE_STREAK_NEEDED" ]; then
            log "Scan appears finished (idle for $((IDLE_STREAK_NEEDED * CHECK_INTERVAL_S / 60)) min). Proceeding."
            break
        fi
    fi
    checks=$((checks + 1))
    sleep "$CHECK_INTERVAL_S"
done

if [ "$checks" -ge "$max_checks" ]; then
    log "WARNING: gave up waiting after ${MAX_WAIT_HOURS}h without a clean idle window. Proceeding anyway."
fi

# ---- 2) Process every RAW run missing a FinalResult ----
processed=0
skipped_done=0
failed=0

for raw in "$RAW_DIR"/Laser/precal_raw_kor_run_*.root "$RAW_DIR"/Dark/precal_raw_kor_run_*.root; do
    [ -f "$raw" ] || continue
    base=$(basename "$raw" .root)                    # precal_raw_kor_run_YYYYMMDD_NNN
    core=${base#precal_raw_kor_run_}                  # YYYYMMDD_NNN
    date_tag=${core%_*}
    run=${core##*_}
    run_int=$((10#$run))                              # strip leading zeros for the macro's int arg

    result="${RESULT_DIR}/precal_result_kor_run_${date_tag}_${run}.root"
    prd="${PROD_DIR}/precal_prd_kor_run_${date_tag}_${run}.root"

    if [ -f "$result" ]; then
        skipped_done=$((skipped_done + 1))
        continue
    fi

    log "Processing run ${date_tag}_${run} (raw: $raw)"

    bash "$RUN_CPP" "$PROD_C" "$CONFIG" "$run_int" "\"$raw\"" >> "$LOG" 2>&1
    bash "$RUN_CPP" "$READ_C" "$CONFIG" "$run_int" "\"$prd\"" >> "$LOG" 2>&1
    bash "$RUN_CPP" "$CONTOUR_C" "$CONFIG" "$run_int" "\"$raw\"" >> "$LOG" 2>&1

    if [ -f "$result" ]; then
        log "  -> DONE: FinalResult produced for run ${date_tag}_${run}."
        processed=$((processed + 1))
    else
        log "  -> WARN: FinalResult still missing for run ${date_tag}_${run} -- check log above."
        failed=$((failed + 1))
    fi
done

log "===== auto_catchup_analysis.sh finished: processed=$processed skipped(already done)=$skipped_done failed=$failed ====="

# ---- 3) Self-remove from crontab -- this was a one-shot request (2026-08-13 only) ----
crontab -l 2>/dev/null | grep -v "auto_catchup_analysis.sh" | grep -v "^# One-shot (2026-08-13 only)" | grep -v "^# catch up any un-analyzed runs\. Self-removes" | crontab -
log "Removed self from crontab (one-shot done)."
