// path_builder.h

#ifndef PATH_BUILDER1_H
#define PATH_BUILDER1_H

#include <string>
#include <vector>
#include <TString.h>
#include <stdexcept>
#include <TSystem.h>
#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config2.h"

std::string format_angle(const std::string& prefix, const std::string& angle_str) {
	if (angle_str.empty()) return "";
	try {
		int angle_val = std::stoi(angle_str);
		char sign = (angle_val < 0) ? 'M' : 'P';
		if (angle_val == 0) {
			return prefix + "00";
		}
		return prefix + sign + std::to_string(std::abs(angle_val));
	} catch (const std::invalid_argument& e) {
		return prefix + "ERR"; 
	}
}

namespace PathUtils {
	inline std::string get_filename_core() {
        std::vector<std::string> serials      = {SN1, SN2, SN3};
        std::vector<std::string> directions   = {direction1, direction2, direction3};
        std::vector<std::string> hvs          = {HV1, HV2, HV3};
        std::vector<std::string> rotateAngles = {RotateAngle1, RotateAngle2, RotateAngle3}; 
        std::vector<std::string> tiltAngles   = {TiltAngle1, TiltAngle2, TiltAngle3};   
        //  {RotateAngle1, RotateAngle2, RotateAngle3}
        //  {TiltAngle1, TiltAngle2, TiltAngle3}

        std::string core = "";
        for (size_t i = 0; i < serials.size(); ++i) {
            if (!serials[i].empty()) {
                if (!core.empty()) core += "_";

                std::string rot_part = format_angle("R", rotateAngles[i]);
                std::string tilt_part = format_angle("T", tiltAngles[i]);

                core += serials[i] + directions[i] + "_hv" + hvs[i] + "_" + rot_part + "_" + tilt_part;
            }
        }
        return core;
    }
}

/*{{{
  namespace PathUtils {
  inline std::string get_filename_core() {
  std::vector<std::string> serials    = {SN1, SN2, SN3};
  std::vector<std::string> directions = {direction1, direction2, direction3};
  std::vector<std::string> hvs        = {HV1, HV2, HV3};

  std::string core = "";
  for (size_t i = 0; i < serials.size(); ++i) {
  if (!serials[i].empty()) {
  if (!core.empty()) core += "_";
  core += serials[i] + directions[i] + "_hv" + hvs[i];
  }
  }
  return core;
  }
  }
  *//*}}}*/
enum FileType {
    RAW_DATA,
    PROCESSED_DATA,
    FINAL_RESULT
};

std::string build_filepath(int run, FileType type, bool is_dark_mode) {
    std::string filename_core = PathUtils::get_filename_core();
    TString run_str_formatted;
    run_str_formatted.Form("%04d", run);

    std::string note_suffix = NOTE.empty() ? "" : "_" + NOTE;

    std::string mode_part;
    if (is_dark_mode) {
        mode_part = "_dark" + note_suffix;
    } else {
        mode_part = "_laser" + Laser + note_suffix;
    }

    TString path;
    switch (type) {
        case RAW_DATA:
            path = TString::Format("%s%s/%s%s.%s.root",
                    RawDataPath.c_str(), (is_dark_mode ? "Dark" : "Laser"),
                    filename_core.c_str(), mode_part.c_str(), run_str_formatted.Data());
            break;
        case PROCESSED_DATA:
            path = TString::Format("%sprd_%s%s.%s.root",
                    ProcessedDataPath.c_str(), filename_core.c_str(), mode_part.c_str(), run_str_formatted.Data());
            break;
        case FINAL_RESULT:
            path = TString::Format("%sresult_%s%s.%s.root",
                    FinalResultPath.c_str(), filename_core.c_str(), mode_part.c_str(), run_str_formatted.Data());
            break;
    }
    return std::string(path.Data());
}

/*{{{*/
/*
   std::string build_filepath(int run, FileType type, bool is_dark_mode) {
   std::string filename_core = PathUtils::get_filename_core();
   TString run_str_formatted;
   run_str_formatted.Form("%04d", run);
   std::string note_suffix = NOTE.empty() ? "" : "_" + NOTE;

   std::string mode_part;
   if (is_dark_mode) {
   mode_part = "_dark" + note_suffix;
   } else {
   mode_part = "_laser" + Laser + "_rot" + RotateAngle + "_tilt" + TiltAngle + note_suffix;
   }

   TString path;
   switch (type) {
   case RAW_DATA:
   path = TString::Format("%s%s/%s%s.%s.root",
   RawDataPath.c_str(), (is_dark_mode ? "Dark" : "Laser"),
   filename_core.c_str(), mode_part.c_str(), run_str_formatted.Data());
   break;
   case PROCESSED_DATA:
   path = TString::Format("%sprd_%s%s.%s.root",
   ProcessedDataPath.c_str(), filename_core.c_str(), mode_part.c_str(), run_str_formatted.Data());
   break;
   case FINAL_RESULT:
   path = TString::Format("%sresult_%s%s.%s.root",
   FinalResultPath.c_str(), filename_core.c_str(), mode_part.c_str(), run_str_formatted.Data());
   break;
   }
   return std::string(path.Data());
   }
   */
/*}}}*/
std::string find_filepath(int run, FileType type) {
    std::string laser_path_str = build_filepath(run, type, false);
    std::string dark_path_str = build_filepath(run, type, true);

    if (!gSystem->AccessPathName(laser_path_str.c_str())) {
        return laser_path_str;
    }
    if (!gSystem->AccessPathName(dark_path_str.c_str())) {
        return dark_path_str;
    }
    return laser_path_str;
}

#endif // PATH_BUILDER_H
