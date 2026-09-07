// Compare_DarkWindow_SignalWindow.C
//
// Standalone check (separate from prod_ntp_v7.C / the Uniformity pipeline):
// for a run range taken with the laser physically off (RunMode=="Laser",
// i.e. external-trigger acquisition, just with no light), count how many
// dark hits land in the standard dark-search window (noiseStart..noiseEnd)
// vs. the standard signal integration window (sigStart..sigEnd), per PMT,
// per angle. Uses the exact same window/threshold formulas as prod_ntp_v7.C's
// laser-mode branch so the comparison is apples-to-apples with production.
//
// Usage:
//   root -l -b -q 'Compare_DarkWindow_SignalWindow.C("20260728", {{0,22},{100,122}}, "black_sheet_check")'

#include <vector>
#include <utility>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TMath.h>
#include <iostream>
#include <cmath>

#include "./Base/analysisCode/Analysis.cpp"
#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"
#include "angle_convert.h"

static const char* kRawDir = "/home/precalkor/Data/RAW/Laser";

struct ChPoint {
    double angle = 0;
    double darkRate = 0, darkErr = 0;
    double sigRate = 0, sigErr = 0;
};

// Process one run: replicate prod_ntp_v7.C's LASER-mode window/threshold
// formulas (these are "Laser" RunMode files, laser physically off) and count
// hits per channel in noiseStart..noiseEnd vs sigStart..sigEnd.
static bool ProcessOneRun(const TString& date, int run, std::vector<TString>& serials,
                          std::vector<double>& angleOut, std::vector<long long>& darkHits,
                          std::vector<long long>& sigHits, double& darkWinLiveTime, double& sigWinLiveTime) {
    TString path = Form("%s/precal_raw_kor_run_%s_%03d.root", kRawDir, date.Data(), run);
    if (gSystem->AccessPathName(path)) { std::cout << "[WARN] missing " << path << std::endl; return false; }
    TFile* file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) { std::cout << "[WARN] could not open " << path << std::endl; return false; }

    char sn1[100] = "Mon", sn2[100] = "PMT1", sn3[100] = "PMT2";
    char dir1[8] = "", dir2[8] = "", dir3[8] = "";
    int tilt2 = 0, rot2 = 0, tilt3 = 0, rot3 = 0;
    TTree* infoTree = (TTree*)file->Get("RunInfo");
    if (infoTree) {
        if (infoTree->GetBranch("SN1")) infoTree->SetBranchAddress("SN1", sn1);
        if (infoTree->GetBranch("SN2")) infoTree->SetBranchAddress("SN2", sn2);
        if (infoTree->GetBranch("SN3")) infoTree->SetBranchAddress("SN3", sn3);
        if (infoTree->GetBranch("Direction1")) infoTree->SetBranchAddress("Direction1", dir1);
        if (infoTree->GetBranch("Direction2")) infoTree->SetBranchAddress("Direction2", dir2);
        if (infoTree->GetBranch("Direction3")) infoTree->SetBranchAddress("Direction3", dir3);
        if (infoTree->GetBranch("RawTiltAngle2")) infoTree->SetBranchAddress("RawTiltAngle2", &tilt2);
        if (infoTree->GetBranch("RawRotateAngle2")) infoTree->SetBranchAddress("RawRotateAngle2", &rot2);
        if (infoTree->GetBranch("RawTiltAngle3")) infoTree->SetBranchAddress("RawTiltAngle3", &tilt3);
        if (infoTree->GetBranch("RawRotateAngle3")) infoTree->SetBranchAddress("RawRotateAngle3", &rot3);
        infoTree->GetEntry(0);
    }
    auto dirOf = [&](int ch) -> char {
        const char* d = (ch == 0) ? dir1 : (ch == 1) ? dir2 : dir3;
        return (d && d[0] != '\0') ? d[0] : DirForCh(ch);
    };

    TTree* T = (TTree*)file->Get("T");
    if (!T) { file->Close(); return false; }
    long long NEntry = T->GetEntries();
    unsigned int RecordLength = 1024, PostTrigger = 60;
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("RecordLength", 1); T->SetBranchStatus("PostTrigger", 1);
    T->SetBranchAddress("RecordLength", &RecordLength);
    T->SetBranchAddress("PostTrigger", &PostTrigger);
    T->GetEntry(0);
    int nSamples = RecordLength;

    UShort_t* ADC_buffer = new UShort_t[8 * nSamples];
    if (!T->GetBranch("ADC")) { std::cout << "[WARN] no ADC branch in " << path << std::endl; delete[] ADC_buffer; file->Close(); return false; }
    T->SetBranchStatus("ADC", 1);
    T->SetBranchAddress("ADC", ADC_buffer);

    std::vector<int> active_channels;
    for (int c = 0; c < 8; ++c) if (T->GetBranch(Form("OffsetValue%d", c))) active_channels.push_back(c);
    int channelCount = active_channels.size();
    if (channelCount == 0) { delete[] ADC_buffer; file->Close(); return false; }

    // Cable-shortening date split, same as prod_ntp_v7.C.
    bool useShortCable = (date.Atoi() >= 20260720);
    int sigStart = useShortCable ? 595 : 610;
    int sigEnd = (sigStart + 50 > nSamples) ? nSamples : sigStart + 50;
    int trgPoint = (int)(nSamples * (1.0 - (PostTrigger / 100.0)));
    int noiseStart = 50;
    int noiseEnd = trgPoint - 50;
    const double thresholdMv = 3.0;   // matches prod_ntp_v7.C's Analysis_Threshold_mV

    serials.assign(channelCount, "");
    angleOut.assign(channelCount, 0.0);
    darkHits.assign(channelCount, 0LL);
    sigHits.assign(channelCount, 0LL);

    // Rate denominator: the TOTAL TIME ACTUALLY SEARCHED, not the run's real
    // wall-clock elapsed time. Each event only has a narrow slice
    // (noiseStart..noiseEnd or sigStart..sigEnd, a few hundred ns) inspected,
    // and this "Laser"/external-trigger acquisition runs at ~1kHz -- e.g. a
    // 100s run has NEntry~100000 events but each dark-window look is only
    // ~618ns, so the observed duty cycle is ~0.06%. Dividing hit counts by
    // the full 100s (as an earlier version of this script did) undercounts
    // the true rate by ~1600x. Matches prod_ntp_v7.C's own totalNoiseLiveTime
    // convention (NEntry * window_width, not TriggerTimeTag-based elapsed
    // time) for exactly this reason.
    darkWinLiveTime = (double)NEntry * (noiseEnd - noiseStart) * 2.0 * 1e-9;
    sigWinLiveTime  = (double)NEntry * (sigEnd - sigStart) * 2.0 * 1e-9;

    std::vector<unsigned int> ADC(nSamples, 0);
    for (long long e = 0; e < NEntry; ++e) {
        T->GetEntry(e);

        for (int i = 0; i < channelCount; ++i) {
            int ch = active_channels[i];
            if (ch == TriggerCh) continue;
            for (int s = 0; s < nSamples; ++s) ADC[s] = ADC_buffer[ch * nSamples + s];
            double pedestal = GetPedestal(ADC.data(), 50, 550, nSamples);
            double threshold = pedestal - (thresholdMv / Config::ADC_to_mV);
            darkHits[i] += (long long)GetTimesBelowThreshold(ADC.data(), threshold, nSamples, noiseStart, noiseEnd).size();
            sigHits[i]  += (long long)GetTimesBelowThreshold(ADC.data(), threshold, nSamples, sigStart, sigEnd).size();
        }
    }

    for (int i = 0; i < channelCount; ++i) {
        int ch = active_channels[i];
        serials[i] = (ch == 0) ? sn1 : (ch == 1) ? sn2 : (ch == 2) ? sn3 : Form("CH%d", ch);
        double tilt_val = (ch == 1) ? (double)tilt2 : (ch == 2) ? (double)tilt3 : 0.0;
        double rot_val  = (ch == 1) ? (double)rot2  : (ch == 2) ? (double)rot3  : 0.0;
        const char* axisLabel = "?-axis";
        angleOut[i] = GetHamamatsuAngle(dirOf(ch), tilt_val, rot_val, axisLabel);
    }

    delete[] ADC_buffer;
    file->Close();
    return true;
}

void Compare_DarkWindow_SignalWindow(TString date, std::vector<std::pair<int,int>> ranges, TString outTag = "") {
    if (outTag.IsNull()) outTag = date;

    std::vector<int> runs;
    for (auto& r : ranges) for (int run = r.first; run <= r.second; ++run) runs.push_back(run);
    if (runs.empty()) { std::cout << "[ERROR] no runs in range" << std::endl; return; }

    // Per-channel-index (0=Monitor,1=PMT1,2=PMT2) accumulated points.
    std::vector<TString> serialByCh(3, "");
    std::vector<std::vector<ChPoint>> pts(3);

    for (int run : runs) {
        std::vector<TString> serials; std::vector<double> angle;
        std::vector<long long> darkHits, sigHits;
        double darkWinLiveTime = 0, sigWinLiveTime = 0;
        if (!ProcessOneRun(date, run, serials, angle, darkHits, sigHits, darkWinLiveTime, sigWinLiveTime)) continue;
        if (darkWinLiveTime <= 0 || sigWinLiveTime <= 0) { std::cout << "[WARN] run " << run << " zero window livetime, skipping" << std::endl; continue; }

        for (size_t i = 0; i < serials.size() && i < 3; ++i) {
            if ((int)i == TriggerCh) continue;
            if (serialByCh[i].IsNull()) serialByCh[i] = serials[i];
            ChPoint p;
            p.angle = angle[i];
            p.darkRate = darkHits[i] / darkWinLiveTime;
            p.darkErr  = (darkHits[i] > 0) ? std::sqrt((double)darkHits[i]) / darkWinLiveTime : 0.0;
            p.sigRate  = sigHits[i] / sigWinLiveTime;
            p.sigErr   = (sigHits[i] > 0) ? std::sqrt((double)sigHits[i]) / sigWinLiveTime : 0.0;
            pts[i].push_back(p);
        }
        std::cout << Form("[INFO] run %03d: darkWinLiveTime=%.6fs sigWinLiveTime=%.6fs", run, darkWinLiveTime, sigWinLiveTime) << std::endl;
    }

    std::cout << "\n===== Summary (mean +/- std over all angles) =====" << std::endl;
    for (int i = 0; i < 3; ++i) {
        if (pts[i].empty()) continue;
        double dSum = 0, sSum = 0; int n = pts[i].size();
        for (auto& p : pts[i]) { dSum += p.darkRate; sSum += p.sigRate; }
        double dMean = dSum / n, sMean = sSum / n;
        double dVar = 0, sVar = 0;
        for (auto& p : pts[i]) { dVar += (p.darkRate - dMean) * (p.darkRate - dMean); sVar += (p.sigRate - sMean) * (p.sigRate - sMean); }
        double dStd = std::sqrt(dVar / n), sStd = std::sqrt(sVar / n);
        std::cout << Form("  %-8s : Dark search window = %6.2f +/- %5.2f Hz  |  Signal range window = %6.2f +/- %5.2f Hz  (N=%d angles)",
                          serialByCh[i].Data(), dMean, dStd, sMean, sStd, n) << std::endl;
    }
    std::cout << "====================================================\n" << std::endl;

    TString pdfPath = "./Data/image/Uniformity/DarkVsSigWindow_" + outTag + ".pdf";
    gSystem->mkdir("./Data/image/Uniformity", kTRUE);
    // Channels with data, known up front -- printing "(" on the first and ")"
    // on the last avoids a separate trailing "close" canvas leaving a blank
    // extra page at the end.
    std::vector<int> chsWithData;
    for (int i = 0; i < 3; ++i) if (!pts[i].empty()) chsWithData.push_back(i);
    if (chsWithData.empty()) { std::cout << "[WARN] No data found for any channel/run in range -- nothing saved." << std::endl; return; }

    for (size_t idx = 0; idx < chsWithData.size(); ++idx) {
        int i = chsWithData[idx];
        int n = pts[i].size();
        std::vector<double> ang(n), darkR(n), darkE(n), sigR(n), sigE(n), zero(n, 0.0);
        for (int k = 0; k < n; ++k) {
            ang[k] = pts[i][k].angle; darkR[k] = pts[i][k].darkRate; darkE[k] = pts[i][k].darkErr;
            sigR[k] = pts[i][k].sigRate; sigE[k] = pts[i][k].sigErr;
        }
        TGraphErrors* gDark = new TGraphErrors(n, ang.data(), darkR.data(), zero.data(), darkE.data());
        TGraphErrors* gSig  = new TGraphErrors(n, ang.data(), sigR.data(), zero.data(), sigE.data());
        gDark->SetMarkerStyle(21); gDark->SetMarkerColor(kBlack); gDark->SetLineColor(kBlack);
        gSig->SetMarkerStyle(20); gSig->SetMarkerColor(kRed + 1); gSig->SetLineColor(kRed + 1);

        TCanvas* c = new TCanvas(Form("c_%d", i), serialByCh[i], 1200, 800);
        c->SetGrid();
        TMultiGraph* mg = new TMultiGraph();
        mg->Add(gDark, "P"); mg->Add(gSig, "P");
        mg->SetTitle(Form("%s: Dark-Search-Window vs Signal-Window Rate;Position angle [degree];Rate [Hz]", serialByCh[i].Data()));
        mg->Draw("A");
        // Middle band, not the top-right corner -- with the dark-search-window
        // cluster riding high and the signal-window cluster sitting low, the
        // corner overlapped the dark points while the gap between the two
        // clusters was empty.
        TLegend* leg = new TLegend(0.20, 0.68, 0.42, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(gDark, "Dark search window", "lp");
        leg->AddEntry(gSig, "Signal range window", "lp");
        leg->Draw();
        TString suffix = (idx == 0) ? "(" : (idx == chsWithData.size() - 1) ? ")" : "";
        c->Print(pdfPath + suffix, Form("pdf Title:%s", serialByCh[i].Data()));
        c->SaveAs(Form("./Data/image/Uniformity/DarkVsSigWindow_%s_%s.png", outTag.Data(), serialByCh[i].Data()));
    }
    std::cout << "[INFO] Saved " << pdfPath << std::endl;
}
