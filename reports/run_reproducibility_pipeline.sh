#!/usr/bin/env bash
# Runs the full reproducibility-analysis chain built up during the QE
# reproducibility study (2026-08): observed-vs-expected reproducibility,
# angle-trend check + bootstrap uncertainty, and the Pull distribution.
# Each stage is a separate, independently-runnable macro (Draw_Reproducibility.C
# is the pre-existing production macro; PlotTrendCheck.C and PlotPullDist.C
# were written during this analysis) -- this script just calls them in order
# with the same two tags so a single command reproduces the whole set of
# Q&A slide figures instead of remembering four separate invocations.
#
# Usage:
#   ./run_reproducibility_pipeline.sh <tag1> <tag2> ["EM5370","EL9590"] [metric]
#
# Example (the Rep.1 vs Rep.2 study used throughout this session):
#   ./run_reproducibility_pipeline.sh 20260803_100_145 20260803_200_245

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

TAG1="${1:?Usage: $0 <tag1> <tag2> [pmt_list] [metric]}"
TAG2="${2:?Usage: $0 <tag1> <tag2> [pmt_list] [metric]}"
PMTS="${3:-{\"EM5370\",\"EL9590\"}}"
METRIC="${4:-MonNorm_CorrRawQE}"

SUMMARY="./Data/UNIFORMITY/Pipeline_Summary_${TAG1}_vs_${TAG2}.txt"
mkdir -p ./Data/UNIFORMITY
> "$SUMMARY"

echo "=================================================="
echo " Reproducibility pipeline: $TAG1 vs $TAG2"
echo "=================================================="

echo ""
echo "[1/3] Draw_Reproducibility.C -- observed/expected/systematic per axis+PMT"
root -l -b -q "Draw_Reproducibility.C(\"$TAG1\", \"$TAG2\", $PMTS, \"$METRIC\", 1.0)" \
    2>&1 | tee -a "$SUMMARY" | grep -v "^Warning\|^Info"

echo ""
echo "[2/3] PlotTrendCheck.C -- angle-trend fit + bootstrap uncertainty on sigma_obs"
root -l -b -q "PlotTrendCheck.C(\"$TAG1\", \"$TAG2\", \"$METRIC\", 1.0)" \
    2>&1 | tee -a "$SUMMARY" | grep -v "^Warning\|^Info"

echo ""
echo "[3/3] PlotPullDist.C -- pull distribution vs N(0,1)"
root -l -b -q "PlotPullDist.C(\"$TAG1\", \"$TAG2\", \"$METRIC\", 1.0)" \
    2>&1 | tee -a "$SUMMARY" | grep -v "^Warning\|^Info"

echo ""
echo "=================================================="
echo " Pipeline complete."
echo "=================================================="
echo " Console log       : $SUMMARY"
echo " Reproducibility   : ./Data/image/Uniformity/Reproducibility_${TAG1}_vs_${TAG2}.pdf"
echo " Trend + bootstrap : ./Data/UNIFORMITY/TrendCheck_RatioVsAngle.png"
echo "                     ./Data/UNIFORMITY/TrendCheck_Bootstrap.png"
echo " Pull distribution : ./Data/UNIFORMITY/PullDistribution_byGroup.png"
echo "                     ./Data/UNIFORMITY/PullDistribution_combined.png"
