#!/bin/bash
# Re-run read_ntp_v7.C over the 9-block HV-scan campaign only, to add the
# new fwhm_exG/sigma_exG branches (2026-09-03, 장지승 박사님 feedback via
# chat: "FWHM으로 같은 방식으로 봐보세" -- computed alongside TTS in the same
# exGaus fit, just never saved before). prod_ntp_v7.C itself is unchanged,
# so Data/production/*.root is reused as-is -- read-only pass, no RAW I/O,
# so no DAQ self-throttle needed (unlike prod reprocessing).
set -u

READ_C="/home/precalkor/ADC/ADC_test/read_ntp_v7.C"
PROD_DIR="/home/precalkor/ADC/ADC_test/Data/production"
LOG="/home/precalkor/ADC/ADC_test/logs/reprocess_read_20260903_fwhmsigma_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# date lo hi
BLOCKS="
20260829 100 145
20260829 200 245
20260829 0 45
20260827 0 45
20260815 0 45
20260827 100 145
20260830 0 45
20260828 0 45
20260828 100 145
"

processed=0; failed=0; skipped=0
cd /home/precalkor/ADC/ADC_test || exit 1

echo "$BLOCKS" | while read -r date lo hi; do
    [ -z "$date" ] && continue
    log "=== Block ${date}_${lo}~${hi} ==="
    for i in $(seq "$lo" "$hi"); do
        run=$(printf '%03d' "$i")
        prod="${PROD_DIR}/precal_prd_kor_run_${date}_${run}.root"
        if [ ! -f "$prod" ]; then
            log "  SKIP ${date}_${run}: no production file"
            continue
        fi
        if root -l -b -q "read_ntp_v7.C(${i},\"${prod}\")" >> "$LOG" 2>&1; then
            log "  OK   ${date}_${run}"
        else
            log "  FAIL ${date}_${run}"
        fi
    done
done

log "=== DONE ==="
