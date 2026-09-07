// RateScan_v7.C
// ----------------------------------------------------------------------------
// Dark-mode threshold scan: from a SINGLE self-triggered dark run, sweep the
// software discrimination threshold (mV) and count the dark-pulse rate per
// channel at each threshold. No re-acquisition is needed — every threshold is
// re-applied in software to the already-recorded waveforms.
//
// Produces a "Rate [kHz] vs Threshold [mV]" plot (log-y) for all active
// channels, including the trigger channel, and saves it as a PNG that the
// GUI Image Viewer can display.
//
// Usage (same arg convention as prod_ntp_v7 so it plugs into run_cpp_script_v2.sh):
//   root -l -b -q 'RateScan_v7.C(RUN, "/path/to/raw_file.root")'
// ----------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TLatex.h>

#include "./Base/analysisCode/Analysis.cpp"
#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"
#include "path_builder2.h"

// ---- Sweep configuration --------------------------------------------------
static const double kThrMin_mV  = 2.0;   // lowest threshold to scan (hw trigger ~3 mV; below is undersampled)
static const double kThrMax_mV  = 5.0;   // highest threshold to scan
static const double kThrStep_mV = 0.2;   // step

void RateScan_v7(int run, const char* raw_file_path = "",
                 double thr_min_mV = kThrMin_mV, double thr_max_mV = kThrMax_mV,
                 double thr_step_mV = kThrStep_mV) {
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);

    // ---- Resolve input file -------------------------------------------------
    std::string OpenPath;
    if (strlen(raw_file_path) > 0) {
        OpenPath = raw_file_path;
    } else {
        OpenPath = find_filepath(run, RAW_DATA);
    }

    TFile *file = TFile::Open(OpenPath.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Could not open raw data file: " << OpenPath << std::endl;
        return;
    }

    // ---- Run mode (warn if not a dark run) ----------------------------------
    char runMode[100] = "Unknown";
    // [DEVELOP] double elapsedTime_s = 0.0;   // for N/T trigger-rate anchor validation
    TTree* infoTree = (TTree*)file->Get("RunInfo");
    if (infoTree && infoTree->GetBranch("RunMode")) {
        infoTree->SetBranchAddress("RunMode", runMode);
    }
    // [DEVELOP] if (infoTree && infoTree->GetBranch("ElapsedTime")) {
    // [DEVELOP]     infoTree->SetBranchAddress("ElapsedTime", &elapsedTime_s);
    // [DEVELOP] }
    if (infoTree) infoTree->GetEntry(0);
    if (TString(runMode) != "Dark" && TString(runMode) != "DARK") {
        std::cout << "[WARN] Run mode is '" << runMode << "', not Dark. "
                     "Rate scan is intended for self-triggered dark runs." << std::endl;
    }

    TTree *T = (TTree*)file->Get("T");
    if (!T) { std::cerr << "[ERROR] Tree 'T' not found." << std::endl; return; }
    long long NEntry = T->GetEntries();

    unsigned int RecordLength = 1024, PostTrigger = 60, TriggerTimeTag = 0;
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("RecordLength", 1);
    T->SetBranchStatus("PostTrigger", 1);
    T->SetBranchStatus("TriggerTimeTag", 1);
    T->SetBranchAddress("RecordLength", &RecordLength);
    T->SetBranchAddress("PostTrigger", &PostTrigger);
    T->SetBranchAddress("TriggerTimeTag", &TriggerTimeTag);
    T->GetEntry(0);

    int nSamples = RecordLength;

    UShort_t* ADC_buffer = new UShort_t[8 * nSamples];
    if (T->GetBranch("ADC")) {
        T->SetBranchStatus("ADC", 1);
        T->SetBranchAddress("ADC", ADC_buffer);
    } else {
        std::cerr << "[ERROR] 'ADC' branch not found." << std::endl;
        delete[] ADC_buffer; return;
    }

    // ---- Active channels (exclude trigger/reference channel from config) ----
    std::vector<int> active_channels;
    for (int i = 0; i < 8; ++i) {
        if (i == TriggerCh) continue;   // skip FG/trigger reference channel
        if (T->GetBranch(Form("OffsetValue%d", i))) active_channels.push_back(i);
    }
    int channelCount = active_channels.size();
    if (channelCount == 0) { std::cerr << "[ERROR] No active channels.\n"; delete[] ADC_buffer; return; }
    std::cout << "[INFO] TriggerCh=" << TriggerCh << " excluded from rate scan." << std::endl;

    // ---- Windows (same convention as prod_ntp_v7) ---------------------------
    int trgPoint   = (int)(nSamples * (1.0 - (PostTrigger / 100.0)));
    int noiseStart = 50;
    int noiseEnd   = trgPoint - 50;
    if (noiseEnd <= noiseStart) noiseEnd = nSamples;   // fallback
    double winDuration_s = (double)(noiseEnd - noiseStart) * 2.0 * 1e-9;  // 2 ns/sample

    int pedStart = 100, pedEnd = 300;

    // ---- Threshold grid -----------------------------------------------------
    std::vector<double> thr_mV;
    for (double t = thr_min_mV; t <= thr_max_mV + 1e-9; t += thr_step_mV)
        thr_mV.push_back(t);
    int nThr = (int)thr_mV.size();

    // hits[ch_index][thr_index]
    std::vector<std::vector<double>> hits(channelCount, std::vector<double>(nThr, 0.0));
    double totalLiveTime_s = 0.0;

    std::cout << "===== Rate Scan: " << NEntry << " events, "
              << channelCount << " ch, " << nThr << " thresholds ("
              << thr_min_mV << "-" << thr_max_mV << " mV) =====" << std::endl;

    std::vector<unsigned int> ADC(nSamples, 0);

    // ---- Event loop ---------------------------------------------------------
    for (long long iEntry = 0; iEntry < NEntry; ++iEntry) {
        T->GetEntry(iEntry);
        if (iEntry % 5000 == 0)
            std::cout << "Processing... " << (iEntry * 100 / NEntry) << "%\r" << std::flush;

        // Each recorded event contributes one noise-window observation per channel.
        totalLiveTime_s += winDuration_s;

        for (int i = 0; i < channelCount; ++i) {
            int ch = active_channels[i];
            for (int s = 0; s < nSamples; ++s) ADC[s] = ADC_buffer[ch * nSamples + s];

            double pedestal = GetPedestal(ADC.data(), pedStart, pedEnd, nSamples);

            for (int it = 0; it < nThr; ++it) {
                double threshold_adc = pedestal - (thr_mV[it] / Config::ADC_to_mV);
                std::vector<int> pulses =
                    GetTimesBelowThreshold(ADC.data(), threshold_adc, nSamples, noiseStart, noiseEnd);
                hits[i][it] += (double)pulses.size();
            }
        }
    }
    std::cout << "Processing... 100%   " << std::endl;

    if (totalLiveTime_s <= 0) {
        std::cerr << "[ERROR] Zero live time — cannot compute rate." << std::endl;
        delete[] ADC_buffer; return;
    }

    // ---- Build graphs (Rate in kHz) ----------------------------------------
    // Distinct colors per channel (matplotlib-like tab palette).
    int palette[8] = {kAzure+1, kOrange+1, kGreen+2, kRed+1, kViolet+1, kSpring-6, kPink+7, kGray+2};

    TCanvas *c = new TCanvas("cRateScan", "Dark Rate vs Threshold", 1000, 750);
    c->SetGrid();
    c->SetLogy();

    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.72, 0.15, 0.88, 0.15 + 0.06 * channelCount);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);

    for (int i = 0; i < channelCount; ++i) {
        int ch = active_channels[i];
        TGraph *g = new TGraph();
        int np = 0;
        for (int it = 0; it < nThr; ++it) {
            double rate_kHz = (hits[i][it] / totalLiveTime_s) / 1000.0;
            if (rate_kHz > 0) {            // log-y: only plot positive rates
                g->SetPoint(np++, thr_mV[it], rate_kHz);
            }
        }
        int col = palette[ch % 8];
        g->SetLineColor(col);
        g->SetMarkerColor(col);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.8);
        g->SetLineWidth(2);
        mg->Add(g, "LP");
        leg->AddEntry(g, Form("ch%d", ch), "lp");
    }

    mg->SetTitle(Form("Dark Rate vs Threshold (Run %d);Threshold [mV];Rate [kHz]", run));
    mg->Draw("A");
    mg->GetXaxis()->SetLimits(thr_min_mV - 0.2, thr_max_mV + 0.2);

    // [DEVELOP] True trigger-rate anchor: draw N/T as horizontal dashed line.
    // Requires ElapsedTime in RunInfo (ADC_test7.cpp) and TriggerCh excluded
    // from ChannelMask so FG does not dominate the trigger rate.
    // Uncomment the [DEVELOP] blocks here and in ADC_test7.cpp to enable.

    leg->Draw();
    c->Modified();
    c->Update();

    // ---- Save PNG into image dir --------------------------------------------
    std::string imageDir = std::string(ImagePath) + "Rate/";
    gSystem->mkdir(imageDir.c_str(), kTRUE);

    TString baseName = gSystem->BaseName(OpenPath.c_str());
    baseName.ReplaceAll(".root", "");
    std::string outPng = imageDir + std::string(baseName.Data()) + "_RateScan.png";
    c->SaveAs(outPng.c_str());

    std::cout << "==============================================================" << std::endl;
    std::cout << "  Rate scan complete." << std::endl;
    std::cout << "  - Input  : " << OpenPath << std::endl;
    std::cout << "  - Output : " << outPng << std::endl;
    std::cout << "  - Live time: " << totalLiveTime_s << " s per channel" << std::endl;
    std::cout << "==============================================================" << std::endl;

    delete[] ADC_buffer;
}
