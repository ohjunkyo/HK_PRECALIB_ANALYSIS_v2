#!/bin/bash
# analyze.sh -- standalone CLI entry point for the analysis macros, so
# someone can run the same analyses the GUI triggers without touching the
# GUI at all (needed once a second person works on analysis alongside the
# DAQ operator). Run with no args (or 'help') for usage.
#
# All macros are invoked plain-interpreted (no ACLiC '+'), matching how
# run_cpp_script_v2.sh already calls them from the GUI -- same behavior,
# just reachable directly from a shell.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

usage() {
    cat <<'EOF'
Usage: ./analyze.sh <command> [args...]

Pipeline (single run):
  prod      <run> [rawFilePath]                 -- prod_ntp_v7: raw -> produced (.prd)
  read      <run> [processedFilePath]            -- read_ntp_v7: produced -> analyzed
  contour   <run> [file] [y_lo] [y_hi] [x_start] [x_end]
                                                  -- Draw_Contour_v3: 2D contour plot
  waveform  <run> [file]                          -- Draw_waveform_v2: waveform inspection

Uniformity / scan-level:
  uniformity <tag> <run_start> <run_end> [chsel]  -- Draw_Uniformity_Norm_v7:
                                                     builds Graphs_Uniformity_<tag>.root
                                                     + a PDF report. chsel default "0,1,2".
  overlay                                          -- Draw_Overlay_Uniformity_v7: overlays
                                                     the tag list from
                                                     Data/UNIFORMITY/overlay_tags.txt
                                                     (edit that file to choose tags/labels)

Stability (fixed-angle repeat runs):
  stability <tag> <run_start> <run_end> [resultDir]
                                                  -- Draw_Stability_v1: QE/Gain/TTS/etc
                                                     stability report + B-field/env pages

Reproducibility (needs Uniformity graphs already built via 'uniformity'):
  repro     <tag1> <tag2> [serials] [metric]      -- Draw_Reproducibility: pairwise
                                                     Rep.1 vs Rep.2 pull test
  repro3    <tag1> <tag2> <tag3> [tag4...] [--serials=A,B] [--metric=M]
                                                  -- Draw_Reproducibility_Covariance3:
                                                     needs >=3 tags. See NOTE below.

Per-PMT cathode uniformity (Hamamatsu slide style):
  hama      <tag> [tag2...] [--serials=A,B] [--metric=M]
            [--separate-axes] [--center-normalize]
                                                  -- Draw_HamamatsuStyle_Uniformity:
                                                     one page per PMT (or per PMT+axis
                                                     with --separate-axes), Center/Side/
                                                     C-S-Ratio/sigma box, and appends a
                                                     summary row per PMT to
                                                     Data/summary/pmt_uniformity_summary.csv
  pmtrep    [metric]                              -- PMT_Representative: mean +
                                                     reproducibility per PMT across
                                                     repeats, plus the monitoring
                                                     conditions of each window

  serials default: EM6400,EL5150   metric default: MonNorm_CorrRawQE

Examples:
  ./analyze.sh prod 145
  ./analyze.sh uniformity 20260817 0 45
  ./analyze.sh repro 20260815_0_45 20260816_0_45
  ./analyze.sh repro3 20260815_0_45 20260816_0_45 20260816_100_145 20260817_0_45
  ./analyze.sh stability 20260813 800 899

NOTE on repro3's small-R limitation: with R repeats and N scan points, the
covariance estimate is only meaningfully conditioned once R is large relative
to N (N=23 typically needs dozens of repeats) -- at R=3 or 4 the per-repeat
chi2/N test is close to self-referential and its "consistent" verdict isn't
strong evidence either way. The 'repro' pairwise pull test (RMS(pull)~1 means
statistics-only, >1 means a real systematic) is the trustworthy metric until
enough repeats accumulate for repro3 to mean something.
EOF
}

run_root() {
    echo "[analyze.sh] root -l -b -q '$1'"
    root -l -b -q "$1"
}

# Comma-join args 2.. for a ROOT call, e.g. join_args a b c -> a,b,c
join_csv() { local IFS=,; echo "$*"; }

cmd="${1:-help}"; shift || true

case "$cmd" in
    prod)
        run="$1"; file="${2:-}"
        [ -z "$run" ] && { echo "usage: ./analyze.sh prod <run> [rawFilePath]"; exit 1; }
        run_root "prod_ntp_v7.C($run,\"$file\")"
        ;;
    read)
        run="$1"; file="${2:-}"
        [ -z "$run" ] && { echo "usage: ./analyze.sh read <run> [processedFilePath]"; exit 1; }
        run_root "read_ntp_v7.C($run,\"$file\")"
        ;;
    contour)
        run="$1"; file="${2:-}"; ylo="${3:-10.0}"; yhi="${4:-3.0}"; xs="${5:--1}"; xe="${6:--1}"
        [ -z "$run" ] && { echo "usage: ./analyze.sh contour <run> [file] [y_lo] [y_hi] [x_start] [x_end]"; exit 1; }
        run_root "Draw_Contour_v3.C($run,\"$file\",$ylo,$yhi,$xs,$xe)"
        ;;
    waveform)
        run="$1"; file="${2:-}"
        [ -z "$run" ] && { echo "usage: ./analyze.sh waveform <run> [file]"; exit 1; }
        run_root "Draw_waveform_v2.C($run,\"$file\")"
        ;;
    uniformity)
        tag="$1"; rs="$2"; re="$3"; chsel="${4:-0,1,2}"
        [ -z "$re" ] && { echo "usage: ./analyze.sh uniformity <tag> <run_start> <run_end> [chsel]"; exit 1; }
        run_root "Draw_Uniformity_Norm_v7.C(\"$tag\",$rs,$re,\"$chsel\")"
        ;;
    overlay)
        run_root "Draw_Overlay_Uniformity_v7.C()"
        ;;
    stability)
        tag="$1"; rs="$2"; re="$3"; resultdir="${4:-./Data/FinalResult}"
        [ -z "$re" ] && { echo "usage: ./analyze.sh stability <tag> <run_start> <run_end> [resultDir]"; exit 1; }
        run_root "Draw_Stability_v1.C(\"$tag\",$rs,$re,\"$resultdir\")"
        ;;
    repro)
        t1="$1"; t2="$2"; serials="${3:-EM6400,EL5150}"; metric="${4:-MonNorm_CorrRawQE}"
        [ -z "$t2" ] && { echo "usage: ./analyze.sh repro <tag1> <tag2> [serials] [metric]"; exit 1; }
        IFS=',' read -ra sarr <<< "$serials"
        serialsCpp="{"; for s in "${sarr[@]}"; do serialsCpp+="\"$s\","; done; serialsCpp="${serialsCpp%,}}"
        run_root "Draw_Reproducibility.C(\"$t1\",\"$t2\",$serialsCpp,\"$metric\",1.0)"
        ;;
    repro3)
        tags=(); serials="EM6400,EL5150"; metric="MonNorm_CorrRawQE"
        for a in "$@"; do
            case "$a" in
                --serials=*) serials="${a#--serials=}" ;;
                --metric=*)  metric="${a#--metric=}" ;;
                *) tags+=("$a") ;;
            esac
        done
        [ "${#tags[@]}" -lt 3 ] && { echo "usage: ./analyze.sh repro3 <tag1> <tag2> <tag3> [tag4...] [--serials=A,B] [--metric=M]  (need >=3 tags)"; exit 1; }
        tagsCpp="{"; for t in "${tags[@]}"; do tagsCpp+="\"$t\","; done; tagsCpp="${tagsCpp%,}}"
        IFS=',' read -ra sarr <<< "$serials"
        serialsCpp="{"; for s in "${sarr[@]}"; do serialsCpp+="\"$s\","; done; serialsCpp="${serialsCpp%,}}"
        run_root "Draw_Reproducibility_Covariance3.C($tagsCpp,$serialsCpp,\"$metric\",1.0)"
        ;;
    repron)
        tags=(); serials="EM6400,EL5150"; metric="MonNorm_CorrRawQE"; labels=""
        for a in "$@"; do
            case "$a" in
                --serials=*) serials="${a#--serials=}" ;;
                --metric=*)  metric="${a#--metric=}" ;;
                --labels=*)  labels="${a#--labels=}" ;;
                *) tags+=("$a") ;;
            esac
        done
        [ "${#tags[@]}" -lt 2 ] && { echo "usage: ./analyze.sh repron <refTag> <tag2> [tag3...] [--serials=A,B] [--metric=M] [--labels=Rep.1,Rep.2,...]"; exit 1; }
        tagsCpp="{"; for t in "${tags[@]}"; do tagsCpp+="\"$t\","; done; tagsCpp="${tagsCpp%,}}"
        IFS=',' read -ra sarr <<< "$serials"
        serialsCpp="{"; for s in "${sarr[@]}"; do serialsCpp+="\"$s\","; done; serialsCpp="${serialsCpp%,}}"
        labelsCpp="{}"
        if [ -n "$labels" ]; then
            IFS=',' read -ra larr <<< "$labels"
            labelsCpp="{"; for l in "${larr[@]}"; do labelsCpp+="\"$l\","; done; labelsCpp="${labelsCpp%,}}"
        fi
        run_root "Draw_Reproducibility_Multi.C($tagsCpp,$serialsCpp,\"$metric\",1.0,$labelsCpp)"
        ;;
    pmtrep)
        metric="${1:-MonNorm_CorrRawQE}"
        run_root "PMT_Representative.C(\"$metric\")"
        ;;
    hama)
        tags=(); serials="EM6400,EL5150"; metric="MonNorm_CorrRawQE"
        centernorm=0; separate=0
        for a in "$@"; do
            case "$a" in
                --serials=*) serials="${a#--serials=}" ;;
                --metric=*)  metric="${a#--metric=}" ;;
                --center-normalize) centernorm=1 ;;
                --separate-axes)    separate=1 ;;
                *) tags+=("$a") ;;
            esac
        done
        [ "${#tags[@]}" -lt 1 ] && { echo "usage: ./analyze.sh hama <tag> [tag2...] [--serials=A,B] [--metric=M] [--separate-axes] [--center-normalize]"; exit 1; }
        tagsCpp="{"; for t in "${tags[@]}"; do tagsCpp+="\"$t\","; done; tagsCpp="${tagsCpp%,}}"
        IFS=',' read -ra sarr <<< "$serials"
        serialsCpp="{"; for s in "${sarr[@]}"; do serialsCpp+="\"$s\","; done; serialsCpp="${serialsCpp%,}}"
        run_root "Draw_HamamatsuStyle_Uniformity.C($tagsCpp,$serialsCpp,\"$metric\",$centernorm,$separate)"
        ;;
    help|-h|--help|"")
        usage
        ;;
    *)
        echo "Unknown command: $cmd"; echo
        usage
        exit 1
        ;;
esac
