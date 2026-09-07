// path_builder2.h
#ifndef PATH_BUILDER_H
#define PATH_BUILDER_H

#include <string>
#include <TString.h>
#include <TSystem.h>
#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"

/*
namespace PathUtils {
    inline std::string find_raw_file(int run, bool is_dark_mode) {
        TString dir = RawDataPath.c_str();
        dir += (is_dark_mode ? "Dark/" : "Laser/");
        
        TString run_str;
        run_str.Form("%03d", run); 

        TString cmd = TString::Format("ls %sprecal_raw_kor_run_*_%s.root 2>/dev/null | tail -1", dir.Data(), run_str.Data());
        TString result = gSystem->GetFromPipe(cmd);
        result.ReplaceAll("\n", ""); 
        
        return std::string(result.Data());
    }
*/

namespace PathUtils {
    inline std::string find_raw_file(int run, bool is_dark_mode) {
        TString run_str;
        run_str.Form("%03d", run);

        TString dir_local = RawDataPath.c_str();
        dir_local += (is_dark_mode ? "Dark/" : "Laser/");
        TString cmd = TString::Format("ls %s*%s*.root 2>/dev/null | tail -1", dir_local.Data(), run_str.Data());
        TString result = gSystem->GetFromPipe(cmd);
        result.ReplaceAll("\n", "");

        if (result.Length() == 0) {
            TString dir_ext = ExternalPath.c_str(); 
            dir_ext += (is_dark_mode ? "Dark/" : "Laser/");
            cmd = TString::Format("ls %s*%s*.root 2>/dev/null | tail -1", dir_ext.Data(), run_str.Data());
            result = gSystem->GetFromPipe(cmd);
            result.ReplaceAll("\n", "");
        }

        return std::string(result.Data());
    }
    inline std::string get_base_name(const std::string& raw_path, int run) {
        if (raw_path.empty()) {
            TString fallback;
            fallback.Form("RunNotFound_%03d", run);
            return std::string(fallback.Data());
        }
        TString p(raw_path.c_str());
        int slash_idx = p.Last('/');
        if (slash_idx != kNPOS) p = p(slash_idx + 1, p.Length());
        p.ReplaceAll(".root", "");
        return std::string(p.Data());
    }
}

enum FileType {
    RAW_DATA,
    PROCESSED_DATA,
    FINAL_RESULT
};


std::string build_filepath(int run, FileType type, bool is_dark_mode) {
    std::string raw_path = PathUtils::find_raw_file(run, is_dark_mode);

    if (type == RAW_DATA) {
        if (raw_path.empty()) {
            TString fallback;
            fallback.Form("%s%s/NotFound_%03d.root", RawDataPath.c_str(), (is_dark_mode ? "Dark" : "Laser"), run);
            return std::string(fallback.Data());
        }
        return raw_path;
    }

    std::string base_name = PathUtils::get_base_name(raw_path, run);
    TString bName(base_name.c_str());
    TString path = "";

    switch (type) {
        case PROCESSED_DATA:
            bName.ReplaceAll("raw", "prd");
            if (!bName.Contains("prd")) bName.Prepend("prd_");
            path = TString::Format("%s%s.root", ProcessedDataPath.c_str(), bName.Data());
            break;
        case FINAL_RESULT:
            bName.ReplaceAll("raw", "result");
            bName.ReplaceAll("prd", "result");
            if (!bName.Contains("result")) bName.Prepend("result_");
            path = TString::Format("%s%s.root", FinalResultPath.c_str(), bName.Data());
            break;
        default:
            break;
    }

    return std::string(path.Data());
}

std::string find_filepath(int run, FileType type) {
    std::string laser_path_str = build_filepath(run, type, false);
    std::string dark_path_str = build_filepath(run, type, true);

    if (!gSystem->AccessPathName(laser_path_str.c_str())) return laser_path_str;
    if (!gSystem->AccessPathName(dark_path_str.c_str())) return dark_path_str;

    return laser_path_str; 
}

#endif
