#!/bin/bash
# ==============================================================================
# Script Name: script_v7.sh
# Author: Junkyo Oh
# Date: 2026-8-xx
# Description: DAQ execution script with auto-launching analysis (v7).
# ==============================================================================

set -e

cleanup_flags() {
    log_and_print "[WARNING] DAQ script interrupted or closed. Purging runtime flags..."
    if [ -n "$RUN_INT" ]; then
        rm -f "/tmp/daq_flags/daq_${RUN_INT}.flag"
        rm -f "/tmp/daq_flags/prod_${RUN_INT}.flag"
        rm -f "/tmp/daq_flags/read_${RUN_INT}.flag"
        rm -f "/tmp/daq_flags/contour_${RUN_INT}.flag"
    fi
    exit 1
}
trap cleanup_flags SIGINT SIGTERM SIGHUP

ANALYSIS_CODE_PATH="/home/precalkor/ADC/ADC_test/prod_ntp_v7.C"
READ_CODE_PATH="/home/precalkor/ADC/ADC_test/read_ntp_v7.C"
CONTOUR_CODE_PATH="/home/precalkor/ADC/ADC_test/Draw_Contour_v3.C"
RUN_CPP_SCRIPT_PATH="/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/run_cpp_script_v2.sh"
# =============================================

run_mode=$1
CONFIG_FILE=$2
ARG_ROT2=$3
ARG_TILT2=$4
ARG_ROT3=$5
ARG_TILT3=$6
BASE_RUN_START=${7:-0}
FILE_TAG=${8:-""}

if [ "$#" -lt 2 ]; then
    echo "Error: This script requires <run_mode> and <config_file_path>"
    exit 1
fi

parse_config() {
    local var_name=$1
    grep "const .* $var_name" "$CONFIG_FILE" | sed -n 's/.*= *"\([^"]*\)".*/\1/p; s/.*= *\([0-9]\+\);.*/\1/p' | tr -d '\r' | xargs
}

log_and_print() {
    local message="$1"
    local LOG_DIR=$(parse_config "LogDir")
    if [ -z "$LOG_DIR" ]; then LOG_DIR="./LOG/DAQ"; fi
    local LOG_FILE="${LOG_DIR}/TakingLog_$(date +"%Y%m%d").txt"
    mkdir -p "$(dirname "$LOG_FILE")"
    echo "$(date +"%Y-%m-%d %H:%M") | $message" | tee -a "$LOG_FILE"
}

get_next_run_filepath() {
    local base_path=$1
    local run_number=1
    while true; do
        local formatted_run_number=$(printf "%04d" $run_number)
        local next_filepath="${base_path}.${formatted_run_number}"
        if [ ! -f "${next_filepath}.root" ]; then echo "$next_filepath"; return; fi
        run_number=$((run_number + 1))
    done
}

DAQ_PROGRAM_PATH=$(parse_config "DaqProgramPath")
DATA_DIR=$(parse_config "RawDataPath")
PROCESSED_DIR=$(parse_config "ProcessedDataPath")
NUM_SEQUENCES=$(parse_config "NumSequences")
INTERVAL_TIME=$(parse_config "IntervalTime")
EVENTS=$(parse_config "Events")
TIME_WINDOW=$(parse_config "TimeWindow")
CHANNEL=$(parse_config "ChannelMask")
POST_TRIGGER=$(parse_config "PostTrigger")

SN2=$(parse_config "SN2"); SN3=$(parse_config "SN3")
ROT2=${ARG_ROT2:-$(parse_config "RotateAngle2")}
TILT2=${ARG_TILT2:-$(parse_config "TiltAngle2")}
ROT3=${ARG_ROT3:-$(parse_config "RotateAngle3")}
TILT3=${ARG_TILT3:-$(parse_config "TiltAngle3")}

LASER_MA=$(parse_config "Laser")
WAVELENGTH=$(parse_config "Wavelength")
HV1=$(parse_config "HV1")
HV2=$(parse_config "HV2")
HV3=$(parse_config "HV3")
SHIFTER=$(parse_config "Shift_worker")
EXPERT=$(parse_config "Expert")
NOTE_VAL=$(parse_config "NOTE")
# boost::program_options aborts on an empty adjacent value ("--note=" ->
# "argument should follow immediately after the equal sign"), so an empty
# NOTE/Shifter/Expert in config3.h must be replaced with a placeholder.
[ -z "$SHIFTER" ]  && SHIFTER="-"
[ -z "$EXPERT" ]   && EXPERT="-"
[ -z "$NOTE_VAL" ] && NOTE_VAL="-"

#--------------------------
if [ "${FILE_FORMAT,,}" == "csv" ]; then
    ext="csv"
    FORMAT_FLAG="" 
else
    ext="root"
    FORMAT_FLAG="-r"
fi

if [ "${run_mode,,}" == "dark" ]; then
    final_data_dir="${DATA_DIR}/Dark"
else
    final_data_dir="${DATA_DIR}/Laser"
fi
mkdir -p "$final_data_dir"

#START_DATE=${SCAN_START_DATE:-$(date +"%Y%m%d")}
if [ -n "$SCAN_START_DATE" ]; then
    START_DATE="$SCAN_START_DATE"
else
    START_DATE=$(date +"%Y%m%d")
fi
LOCATION="kor"
BASE_PREFIX="precal_raw_${LOCATION}_run_${START_DATE}"

#UPPER_BOUND=$((BASE_RUN_START + 99)){{{
#LAST_NUM=$(ls ${final_data_dir}/${BASE_PREFIX}_*.root 2>/dev/null | awk -F'_' '{print $NF}' | sed 's/.root//' | awk -v min="$BASE_RUN_START" -v max="$UPPER_BOUND" '{num=$1+0; if(num>=min && num<=max) print num}' | sort -n | tail -1)}}}
UPPER_BOUND=$((BASE_RUN_START + 99))
LAST_NUM=$(ls ${final_data_dir}/${BASE_PREFIX}_*.root ${final_data_dir}/${BASE_PREFIX}_*.csv 2>/dev/null | awk -F'_' '{print $NF}' | sed -E 's/\.(root|csv)//' | awk -v min="$BASE_RUN_START" -v max="$UPPER_BOUND" '{num=$1+0; if(num>=min && num<=max) print num}' | sort -n | tail -1)

#if [ "${run_mode,,}" == "dark" ]; then{{{
#    LAST_NUM=$(ls ${final_data_dir}/${BASE_PREFIX}_*_dark*.root ${final_data_dir}/${BASE_PREFIX}_*_dark*.csv 2>/dev/null | awk -F'_' '{print $NF}' | sed -E 's/\.(root|csv)//' | awk -v min="$BASE_RUN_START" -v max="$UPPER_BOUND" '{num=$1+0; if(num>=min && num<=max) print num}' | sort -n | tail -1)
#else
#    LAST_NUM=$(ls ${final_data_dir}/${BASE_PREFIX}_*_laser*.root ${final_data_dir}/${BASE_PREFIX}_*_laser*.csv 2>/dev/null | awk -F'_' '{print $NF}' | sed -E 's/\.(root|csv)//' | awk -v min="$BASE_RUN_START" -v max="$UPPER_BOUND" '{num=$1+0; if(num>=min && num<=max) print num}' | sort -n | tail -1)
#fi
# }}}
if [ -z "$LAST_NUM" ]; then
    NEXT_NUM=$(printf "%03d" $BASE_RUN_START)
else
    NEXT_NUM=$(printf "%03d" $((10#$LAST_NUM + 1)))
fi

filename_core="${BASE_PREFIX}_${NEXT_NUM}"
log_and_print "=========================================================="
log_and_print "[INFO] Run Type      : ${run_mode^}"
log_and_print "[INFO] Generated File: ${filename_core}.${ext}" 
log_and_print "[INFO] Target Dir    : ${final_data_dir}"

FIRST_FILE=""; LAST_FILE=""

for (( seq=1; seq<=NUM_SEQUENCES; seq++ )); do
    if [ "$seq" -gt 1 ]; then
        NEXT_NUM=$(printf "%03d" $((10#$NEXT_NUM + 1)))
        filename_core="${BASE_PREFIX}_${NEXT_NUM}"
    fi

    output_filepath="${final_data_dir}/${filename_core}"
    actual_filename="${filename_core}.${ext}"

    [ $seq -eq 1 ] && FIRST_FILE="$actual_filename"
    LAST_FILE="$actual_filename"

    DAQ_EXEC="${DAQ_PROGRAM_PATH}execute_DAQ_v2"
    
    if [ "$run_mode" == "dark" ] || [ "$run_mode" == "Dark" ]; then
        MODE_FLAG="-y 1 -z 1"
    else
        MODE_FLAG="-y 1 -e"
    fi

    DAQ_COMMAND="$DAQ_EXEC -c $CHANNEL -o $output_filepath $FORMAT_FLAG -n $EVENTS -d $TIME_WINDOW $MODE_FLAG -p $POST_TRIGGER -M ${run_mode^} --rot2=$ROT2 --tilt2=$TILT2 --rot3=$ROT3 --tilt3=$TILT3 --laser=$LASER_MA --wavelength=$WAVELENGTH --hv1=$HV1 --hv2=$HV2 --hv3=$HV3 --shifter=\"$SHIFTER\" --expert=\"$EXPERT\" --note=\"$NOTE_VAL\""
    
    RUN_INT=$((10#$NEXT_NUM))

    FLAG_DIR="/tmp/daq_flags"
    mkdir -p "$FLAG_DIR"
    sync 
    sleep 0.1

    DAQ_FLAG="${FLAG_DIR}/daq_${RUN_INT}.flag"
    PROD_FLAG="${FLAG_DIR}/prod_${RUN_INT}.flag"
    READ_FLAG="${FLAG_DIR}/read_${RUN_INT}.flag"
    CONT_FLAG="${FLAG_DIR}/contour_${RUN_INT}.flag"
    DONE_FLAG="${FLAG_DIR}/done_${RUN_INT}.flag"

    touch "$DAQ_FLAG"

    log_and_print "Sequence $seq/$NUM_SEQUENCES: Executing DAQ..."
    eval "$DAQ_COMMAND"

    rm -f "$DAQ_FLAG"
    
    if [ "$ext" == "root" ]; then
        # prod_ntp_v7 reads the RAW file and writes the PRODUCED (prd) file into
        # ProcessedDataPath together with the analysis params (Config_Threshold_mV,
        # NoiseCountRate, ...). read_ntp_v7 NEEDS those params, so it must be given the
        # prd file, NOT the raw file (otherwise it reads an empty/wrong result).
        # prod/contour read the raw file directly, which is why only read was failing.
        prd_core="${filename_core/raw/prd}"
        prd_filepath="${PROCESSED_DIR%/}/${prd_core}.root"

        PROD_CMD="${RUN_CPP_SCRIPT_PATH} ${ANALYSIS_CODE_PATH} ${CONFIG_FILE} ${RUN_INT} \\\"${output_filepath}.root\\\""
        READ_CMD="${RUN_CPP_SCRIPT_PATH} ${READ_CODE_PATH} ${CONFIG_FILE} ${RUN_INT} \\\"${prd_filepath}\\\""
        CONTOUR_CMD="${RUN_CPP_SCRIPT_PATH} ${CONTOUR_CODE_PATH} ${CONFIG_FILE} ${RUN_INT} \\\"${output_filepath}.root\\\""

        FULL_ANA_CMD="touch ${PROD_FLAG} && bash ${PROD_CMD} && rm -f ${PROD_FLAG} && touch ${READ_FLAG} && bash ${READ_CMD} && rm -f ${READ_FLAG} && touch ${CONT_FLAG} && bash ${CONTOUR_CMD} && rm -f ${CONT_FLAG} && touch ${DONE_FLAG} ; echo 'Analysis sequence completed.' ; exec bash"

        log_and_print "Launching background analysis for $actual_filename..."
        if command -v tmux &> /dev/null && [ -n "$TMUX" ]; then
            tmux new-window -d -n "Ana_Seq$seq" "bash -c '$FULL_ANA_CMD'"
        else
            tmux new-session -d -s "Ana_Seq$seq" "bash -c '$FULL_ANA_CMD'"
        fi
    else
        log_and_print "[INFO] CSV format detected. Skipping background analysis chain."
    fi

    if [ "$seq" -lt "$NUM_SEQUENCES" ]; then
        log_and_print "Waiting $INTERVAL_TIME sec..."
        sleep "$INTERVAL_TIME"
    fi
done

log_and_print "=========================================================="
