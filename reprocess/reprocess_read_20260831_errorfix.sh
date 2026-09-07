#!/bin/bash
# Re-run read_ntp_v7.C (SPE fit -> spe_mean/spe_mean_error, QE, etc.) for the
# blocks affected by the 2026-08-31 signalError fix (covariance-corrected
# instead of GetParError(3) alone -- see project memory
# "gain-curve-spe-fit-decision"). prod_ntp_v7.C itself was NOT changed, so
# existing Data/production/*.root files are reused -- only 20260829_200~245
# needs the full RAW->prod->read chain, because its production/ files were
# lost in the 8/29 run-number-collision incident (see project memory
# "run-number-collision-20260829") and never regenerated.
#
# Scope (user-approved 2026-08-31):
#   - HV scan campaign: 8/27-8/30, 9 blocks
#   - Most recent 2 Stability runs: 8/13, 8/14 (block 800)
#   - Most recent 4 General Scans outside the HV campaign: 8/26, 8/17, 8/16x2
#
# After this: Draw_Overlay_Uniformity_v7.C reruns to pick up corrected errors
# in the Overlay report.
set -u

RUN_CPP="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/run_cpp_script_v2.sh"
PROD_C="/home/precalkor/ADC/ADC_test/prod_ntp_v7.C"
READ_C="/home/precalkor/ADC/ADC_test/read_ntp_v7.C"
CONFIG="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"

RAW_DIR="/home/precalkor/Data/RAW/Laser"
EXT_RAW_DIR="/media/precalkor/HD-EDS-E/Data_Backup/RAW/Laser"
PROD_DIR="/home/precalkor/ADC/ADC_test/Data/production"

LOG="/home/precalkor/ADC/ADC_test/logs/reprocess_read_20260831_errorfix_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# date lo hi needs_prod(1/0)
BLOCKS="
20260827 0 45 0
20260827 100 145 0
20260828 0 45 0
20260828 100 145 0
20260829 0 45 0
20260829 100 145 0
20260829 200 245 1
20260830 0 45 0
20260815 0 45 0
20260813 800 899 0
20260814 800 899 0
20260826 0 4 0
20260817 0 45 0
20260816 0 45 0
20260816 100 145 0
"

processed=0; failed=0; skipped=0

echo "$BLOCKS" | while read -r date lo hi needs_prod; do
    [ -z "$date" ] && continue
    log "=== Block ${date}_${lo}~${hi} (needs_prod=${needs_prod}) ==="
    for i in $(seq "$lo" "$hi"); do
        run=$(printf '%03d' "$i")
        prod="${PROD_DIR}/precal_prd_kor_run_${date}_${run}.root"

        if [ "$needs_prod" = "1" ] && [ ! -f "$prod" ]; then
            raw="${RAW_DIR}/precal_raw_kor_run_${date}_${run}.root"
            if [ ! -f "$raw" ]; then
                ext="${EXT_RAW_DIR}/precal_raw_kor_run_${date}_${run}.root"
                if [ -f "$ext" ]; then
                    cp "$ext" "$raw"
                else
                    log "  MISSING RAW for ${date}_${run} (checked local + external) -- skipping."
                    continue
                fi
            fi
            log "  prod: ${date}_${run}"
            bash "$RUN_CPP" "$PROD_C" "$CONFIG" "$i" "\"$raw\"" >> "$LOG" 2>&1
        fi

        if [ ! -f "$prod" ]; then
            log "  MISSING production for ${date}_${run} -- skipping read."
            continue
        fi

        bash "$RUN_CPP" "$READ_C" "$CONFIG" "$i" "\"$prod\"" >> "$LOG" 2>&1
        log "  read: ${date}_${run} done"
    done
done

log "===== reprocess_read finished ====="

# Regenerate the Overlay Uniformity report so it reflects the corrected errors.
cd /home/precalkor/ADC/ADC_test
log "Running Draw_Overlay_Uniformity_v7.C ..."
root -l -b -q Draw_Overlay_Uniformity_v7.C >> "$LOG" 2>&1
log "Running Draw_GainCurve_v1.C ..."
root -l -b -q Draw_GainCurve_v1.C >> "$LOG" 2>&1

log "===== ALL DONE ====="
log "Log: $LOG"
