#!/bin/bash
# Pipelined FULL-campaign reprocess, producing BOTH thresholds in one pass:
#   1.5mV (production)  -> Data/production        -> Data/FinalResult
#   3.0mV (comparison)  -> Data/production_thr3_compare -> Data/FinalResult_thr3compare
#
# Replaces the strictly-serial reprocess_full_20260901_thr1p5.sh, which spent
# ~85% of its wall time idle: it ran ONE run at a time and, before each one,
# waited out the entire General Scan acquisition (322s cycle for 45s of real
# work). Measured bottlenecks (2026-09-02):
#   external HDD : 26 MB/s  (HD-EDS-E is on a 480Mbps USB2 root hub, SHARED
#                            with the CAEN DT5xxx digitizer -- that shared bus
#                            is why the serial script had to throttle against
#                            execute_DAQ_v2 at all)
#   local NVMe   : 505 MB/s
#   prod CPU     : ~19s/run, read_ntp ~5s/run
# So the only part that must be serialised against the DAQ is pulling bytes off
# the external disk. Everything downstream is local and CPU-bound.
#
# Design: two decoupled stages.
#   STAGER  - one process, copies RAW external->local NVMe, DAQ-throttled,
#             keeps at most STAGE_BUFFER files parked locally.
#   WORKERS - N parallel chains on locally-staged files. Never touch the
#             external disk, so they never contend with the DAQ and never wait.
# Because each RAW is staged ONCE and both thresholds are derived from that one
# local copy, the 3mV comparison costs only CPU -- which hides entirely behind
# the USB2 staging rate. Expected ~3h for 506 runs x 2 thresholds, against ~38h
# for the serial script doing 1.5mV alone.
set -u

NWORKERS="${NWORKERS:-8}"
STAGE_BUFFER="${STAGE_BUFFER:-24}"     # max RAW files parked locally (~15GB)

RUN_CPP="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/run_cpp_script_v2.sh"
CONFIG="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"
PROD_C="/home/precalkor/ADC/ADC_test/prod_ntp_v7.C"
READ_C="/home/precalkor/ADC/ADC_test/read_ntp_v7.C"
PROD3_C="/home/precalkor/ADC/ADC_test/prod_ntp_v7_thr3.C"
READ3_C="/home/precalkor/ADC/ADC_test/read_ntp_v7_thr3compare.C"

LOCAL_RAW_DIR="/home/precalkor/ADC/ADC_test/Data/RAW/Laser"
EXT_RAW_DIR="/media/precalkor/HD-EDS-E/Data_Backup/RAW/Laser"
PROD_DIR="/home/precalkor/ADC/ADC_test/Data/production"
FINAL_DIR="/home/precalkor/ADC/ADC_test/Data/FinalResult"
PROD3_DIR="/home/precalkor/ADC/ADC_test/Data/production_thr3_compare"
FINAL3_DIR="/home/precalkor/ADC/ADC_test/Data/FinalResult_thr3compare"
STAGE_DIR="/home/precalkor/ADC/ADC_test/Data/RAW/_stage"

LOGDIR="/home/precalkor/ADC/ADC_test/logs"
LOG="$LOGDIR/reprocess_pipelined_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$LOGDIR" "$STAGE_DIR" "$PROD3_DIR" "$FINAL3_DIR"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# Runs needing the 1.5mV chain: the whole campaign.
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

# Runs needing the 3mV chain: only the 9 HV points x 46 angles that feed the
# threshold-comparison figure. 20260815 (Recom.) is here but NOT in BLOCKS --
# its 1.5mV results were produced by the earlier diagnostic pass, so it needs
# the 3mV chain only.
HV_BLOCKS="
20260829 1 -200V
20260829 2 -150V
20260829 0 -100V
20260827 0 -50V
20260815 0 Recom
20260827 1 +50V
20260830 0 +100V
20260828 0 +150V
20260828 1 +200V
"

# A 1.5mV FinalResult only counts as done if it is NEWER than the 3mV->1.5mV
# threshold change. Nearly every run already had a stale 3mV FinalResult from
# months of normal operation, and a plain "file exists" check silently skipped
# ~400 of them once already.
THRESHOLD_CHANGE_EPOCH=$(date -d "2026-09-01 14:50:35" +%s)

# ------------------------------------------------------------- TASK LIST ----
# One line per run: "<date> <run> <need15> <need3>". Built once up front so the
# stager pulls each RAW exactly once even when both chains want it.
TASKS="$STAGE_DIR/.tasks"
build_tasks() {
    local tmp="$STAGE_DIR/.tasks.tmp"; : > "$tmp"
    echo "$BLOCKS" | while read -r date lo hi; do
        [ -z "$date" ] && continue
        for i in $(seq "$lo" "$hi"); do printf '%s %03d\n' "$date" "$i"; done
    done >> "$tmp"
    echo "$HV_BLOCKS" | while read -r date blk label; do
        [ -z "$date" ] && continue
        for a in $(seq 0 45); do printf '%s %03d\n' "$date" "$((blk*100+a))"; done
    done >> "$tmp"

    sort -u "$tmp" | while read -r date run; do
        [ -z "$date" ] && continue
        local n15=0 n3=0
        local f15="${FINAL_DIR}/precal_result_kor_run_${date}_${run}.root"
        local f3="${FINAL3_DIR}/precal_result_kor_run_${date}_${run}.root"
        if [ -f "$f15" ]; then
            local m; m=$(stat -c %Y "$f15" 2>/dev/null || echo 0)
            [ "$m" -ge "$THRESHOLD_CHANGE_EPOCH" ] || n15=1
        else
            n15=1
        fi
        # The 3mV chain is wanted only for the HV-scan runs; everything else
        # stays 1.5mV-only, so an absent 3mV result there is not "missing".
        if echo "$HV_WANTED" | grep -qx "${date}_${run}"; then
            [ -f "$f3" ] || n3=1
        fi
        [ "$n15" = 0 ] && [ "$n3" = 0 ] && continue
        echo "$date $run $n15 $n3"
    done > "$TASKS"
    rm -f "$tmp"
}

HV_WANTED=$(echo "$HV_BLOCKS" | while read -r date blk label; do
    [ -z "$date" ] && continue
    for a in $(seq 0 45); do printf '%s_%03d\n' "$date" "$((blk*100+a))"; done
done)
export HV_WANTED

# ---------------------------------------------------------------- STAGER ----
# Only this stage reads the external disk, so only this stage throttles.
stage_all() {
    while read -r date run n15 n3; do
        [ -z "$date" ] && continue
        local_raw="${LOCAL_RAW_DIR}/precal_raw_kor_run_${date}_${run}.root"
        ext_raw="${EXT_RAW_DIR}/precal_raw_kor_run_${date}_${run}.root"
        staged="${STAGE_DIR}/precal_raw_kor_run_${date}_${run}.root"

        # Back-pressure: never let the staging dir outgrow the buffer.
        while [ "$(ls -1 "$STAGE_DIR"/precal_raw_kor_run_*.root 2>/dev/null | wc -l)" -ge "$STAGE_BUFFER" ]; do
            sleep 5
        done

        echo "$n15 $n3" > "$STAGE_DIR/.want_${date}_${run}"

        if [ -f "$local_raw" ]; then
            ln -sf "$local_raw" "$staged"          # already local, no copy
            continue
        fi
        if [ ! -f "$ext_raw" ]; then
            log "  MISSING RAW ${date}_${run}"
            rm -f "$STAGE_DIR/.want_${date}_${run}"
            continue
        fi

        # Wait out any live acquisition before touching the shared USB2 bus.
        while pgrep -x execute_DAQ_v2 > /dev/null 2>&1; do sleep 3; done

        cp "$ext_raw" "$staged.part" 2>>"$LOG" && mv "$staged.part" "$staged" \
            || { log "  STAGE FAIL ${date}_${run}"; rm -f "$staged.part" "$STAGE_DIR/.want_${date}_${run}"; }
    done < "$TASKS"
    touch "$STAGE_DIR/.staging_done"
}

# --------------------------------------------------------------- WORKERS ----
process_one() {
    local date="$1" run="$2" idx="$3"
    local staged="${STAGE_DIR}/precal_raw_kor_run_${date}_${run}.root"
    local prod="${PROD_DIR}/precal_prd_kor_run_${date}_${run}.root"
    local prod3="${PROD3_DIR}/precal_prd_kor_run_${date}_${run}.root"
    local want="${STAGE_DIR}/.want_${date}_${run}"
    local n15=1 n3=0
    [ -f "$want" ] && read -r n15 n3 < "$want"

    if [ "$n15" = 1 ]; then
        bash "$RUN_CPP" "$PROD_C" "$CONFIG" "$idx" "\"$staged\"" >> "$LOG.w" 2>&1
        if [ -f "$prod" ]; then
            bash "$RUN_CPP" "$READ_C" "$CONFIG" "$idx" "\"$prod\"" >> "$LOG.w" 2>&1
        else
            log "  MISSING production ${date}_${run} -- skip read"
        fi
    fi

    if [ "$n3" = 1 ]; then
        bash "$RUN_CPP" "$PROD3_C" "$CONFIG" "$idx" "\"$staged\"" >> "$LOG.w" 2>&1
        if [ -f "$prod3" ]; then
            bash "$RUN_CPP" "$READ3_C" "$CONFIG" "$idx" "\"$prod3\"" >> "$LOG.w" 2>&1
        else
            log "  MISSING production_thr3 ${date}_${run} -- skip read3"
        fi
    fi

    # Free the buffer slot (symlinks to the live local RAW are only unlinked).
    rm -f "$staged" "$want"
    rmdir "$STAGE_DIR/.claim_precal_raw_kor_run_${date}_${run}" 2>/dev/null
    log "  done ${date}_${run} (15mV=$n15 3mV=$n3)"
}
export -f process_one log
export STAGE_DIR PROD_DIR PROD3_DIR RUN_CPP PROD_C READ_C PROD3_C READ3_C CONFIG LOG

consume() {
    while true; do
        for staged in "$STAGE_DIR"/precal_raw_kor_run_*.root; do
            [ -e "$staged" ] || continue
            base=$(basename "$staged" .root)
            # Atomic claim: mkdir succeeds for exactly one consumer, so a file
            # can never be dispatched twice, and nothing is lost if the stager
            # adds files while we iterate.
            mkdir "$STAGE_DIR/.claim_$base" 2>/dev/null || continue
            d=$(echo "$base" | sed -E 's/.*_run_([0-9]{8})_([0-9]{3})$/\1/')
            r=$(echo "$base" | sed -E 's/.*_run_([0-9]{8})_([0-9]{3})$/\2/')
            echo "$d $r $((10#$r))"
        done | xargs -r -P "$NWORKERS" -n 3 bash -c 'process_one "$0" "$1" "$2"'

        if [ -f "$STAGE_DIR/.staging_done" ]; then
            ls -1 "$STAGE_DIR"/precal_raw_kor_run_*.root >/dev/null 2>&1 || break
        fi
        sleep 3
    done
}

# ------------------------------------------------------------------ MAIN ----
rm -rf "$STAGE_DIR"/.claim_* "$STAGE_DIR"/.want_* "$STAGE_DIR/.staging_done"

log "===== pipelined reprocess start (workers=$NWORKERS buffer=$STAGE_BUFFER) ====="
log "Building task list ..."
build_tasks
NTASK=$(wc -l < "$TASKS")
N15=$(awk '$3==1' "$TASKS" | wc -l)
N3=$(awk '$4==1' "$TASKS" | wc -l)
log "  $NTASK runs to stage  (1.5mV chain: $N15, 3mV chain: $N3)"

stage_all &
STAGER_PID=$!
consume
wait $STAGER_PID 2>/dev/null

log "===== reprocess finished -- drawing reports ====="
cd /home/precalkor/ADC/ADC_test
log "Draw_Overlay_Uniformity_v7.C ..."
root -l -b -q Draw_Overlay_Uniformity_v7.C >> "$LOG" 2>&1
log "Draw_HVScan_Report.C ..."
root -l -b -q Draw_HVScan_Report.C >> "$LOG" 2>&1
log "Draw_Conclusion_ThresholdStages.C (all HV x all angles) ..."
root -l -b -q Draw_Conclusion_ThresholdStages.C >> "$LOG" 2>&1
log "===== ALL DONE ====="
log "Log: $LOG"
