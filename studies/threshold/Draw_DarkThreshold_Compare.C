// Usage:  root -l -b -q Draw_DarkThreshold_Compare.C
//
// Dark/noise rate compared between two threshold definitions, analogous to
// Draw_QEThreshold_Compare.C but for the DARK side of the same threshold
// question (user: "Dark도 비교해봐야되네", 2026-08-31):
//   1. Fixed 3mV  (what prod_ntp_v7.C actually uses today, Analysis_Threshold_mV)
//   2. Pedestal mu+5sigma -- HV-independent noise-scaled threshold
//
// Unlike the QE comparison, NoiseCountRate is NOT a post-hoc histogram
// re-cut: prod_ntp_v7.C counts below-threshold SAMPLE crossings live during
// its event loop (GetTimesBelowThreshold over noiseStart..noiseEnd), so this
// macro re-loops the same RAW waveforms and counts hits under BOTH
// thresholds in a single pass (no new DAQ needed -- RAW files already exist).
//
// DIAGNOSTIC ONLY -- does not touch prod_ntp_v7.C / production files. See
// project memory "pending-adaptive-threshold-fix".

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include <vector>
#include <map>
#include <iostream>
#include <cmath>

#include "./Base/analysisCode/Analysis.cpp"
#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"

const double Analysis_Threshold_mV = 3.0;

struct HVPoint { const char* date; const char* run; const char* label; };
static const std::vector<HVPoint> kPoints = {
    {"20260829", "111", "-200V"},
    {"20260815", "011", "Recom."},
    {"20260828", "111", "+200V"},
};

// Re-loops one RAW file and returns, per non-trigger channel: dark hit
// counts under the fixed 3mV threshold and under a pedestal-mu+5sigma
// threshold (sigma estimated from the pedestal-window sample RMS, averaged
// over the first 2000 events), plus the live time used for both rates.
struct DarkResult { double rate3mv = 0, rate5sig = 0, thr5sigMv = 0; };

std::map<int, DarkResult> ComputeDarkRates(const char* rawPath) {
    std::map<int, DarkResult> out;
    TFile* file = TFile::Open(rawPath, "READ");
    if (!file || file->IsZombie()) { if (file) file->Close(); return out; }

    bool darkmode = false;
    char runMode[100] = "Unknown";
    TTree* infoTree = (TTree*)file->Get("RunInfo");
    if (infoTree) {
        if (infoTree->GetBranch("RunMode")) infoTree->SetBranchAddress("RunMode", runMode);
        infoTree->GetEntry(0);
        if (TString(runMode) == "Dark" || TString(runMode) == "DARK") darkmode = true;
    }

    TTree* T = (TTree*)file->Get("T");
    if (!T) { file->Close(); return out; }
    long long NEntry = T->GetEntries();
    unsigned int RecordLength = 1024;
    bool constsOnT = (T->GetBranch("RecordLength") != nullptr);
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("TriggerTimeTag", 1);
    unsigned int TriggerTimeTag = 0;
    T->SetBranchAddress("TriggerTimeTag", &TriggerTimeTag);
    if (constsOnT) {
        T->SetBranchStatus("RecordLength", 1);
        T->SetBranchAddress("RecordLength", &RecordLength);
    }
    T->GetEntry(0);
    int nSamples = RecordLength;

    UShort_t* ADC_buffer = new UShort_t[8 * nSamples];
    T->SetBranchStatus("ADC", 1);
    T->SetBranchAddress("ADC", ADC_buffer);

    std::vector<int> activeChannels;
    for (int i = 0; i < 8; ++i) if (T->GetBranch(Form("OffsetValue%d", i))) activeChannels.push_back(i);
    if (activeChannels.empty() && infoTree) {
        if (infoTree->GetBranch("ChannelMask")) {
            char chMask[16] = {0};
            infoTree->SetBranchAddress("ChannelMask", chMask);
            infoTree->GetEntry(0);
            infoTree->ResetBranchAddresses();
            int len = (int)strlen(chMask);
            for (int i = 0; i < len && i < 8; ++i)
                if (chMask[len - 1 - i] == '1') activeChannels.push_back(i);
        }
    }
    if (activeChannels.empty()) { delete[] ADC_buffer; file->Close(); return out; }

    int trgPoint = (int)(nSamples * (1.0 - (60 / 100.0)));  // PostTrigger fallback ~60
    int pedStart = darkmode ? 100 : 50;
    int pedEnd   = darkmode ? 300 : 550;
    int sigStart = 580;  // approximation, second-shift-era cable; only noise window matters here
    int noiseStart = 50;
    int noiseEnd = darkmode ? std::max(noiseStart + 10, sigStart - 20) : (trgPoint - 50);
    double liveTimePerEvt_s = (double)(noiseEnd - noiseStart) * 2.0 * 1e-9;

    // Pass 1: estimate pedestal sigma per channel from the first 2000 events'
    // sample-to-sample spread within the pedestal window (baseline noise RMS).
    long long sigmaScanN = std::min((long long)2000, NEntry);
    std::map<int, double> sumSq, sumV; std::map<int, long long> nSamp;
    std::vector<unsigned int> ADC(nSamples);
    for (long long e = 0; e < sigmaScanN; ++e) {
        T->GetEntry(e);
        for (int ch : activeChannels) {
            for (int s = 0; s < nSamples; ++s) ADC[s] = ADC_buffer[ch * nSamples + s];
            double ped = GetPedestal(ADC.data(), pedStart, pedEnd, nSamples);
            for (int s = pedStart; s < pedEnd; ++s) {
                double v = (double)ADC[s] - ped;
                sumSq[ch] += v * v; nSamp[ch]++;
            }
        }
    }
    std::map<int, double> sigmaAdc;
    for (int ch : activeChannels) sigmaAdc[ch] = (nSamp[ch] > 0) ? std::sqrt(sumSq[ch] / nSamp[ch]) : 0.0;

    // Pass 2: count below-threshold sample-crossings under both thresholds.
    std::map<int, long long> hits3mv, hits5sig;
    double LiveTime = 0;
    unsigned int prevTag = 0;
    for (long long e = 0; e < NEntry; ++e) {
        T->GetEntry(e);
        if (e > 0) {
            double diffTag = (double)TriggerTimeTag - (double)prevTag;
            if (diffTag < 0) diffTag += TMath::Power(2, 31);
            LiveTime += diffTag * 8.0 * 1e-9;
        }
        prevTag = TriggerTimeTag;
        for (int ch : activeChannels) {
            if (ch == TriggerCh) continue;
            for (int s = 0; s < nSamples; ++s) ADC[s] = ADC_buffer[ch * nSamples + s];
            double ped = GetPedestal(ADC.data(), pedStart, pedEnd, nSamples);
            double thr3 = ped - (Analysis_Threshold_mV / Config::ADC_to_mV);
            double thr5 = ped - (5.0 * sigmaAdc[ch]);
            hits3mv[ch]  += (long long)GetTimesBelowThreshold(ADC.data(), thr3, nSamples, noiseStart, noiseEnd).size();
            hits5sig[ch] += (long long)GetTimesBelowThreshold(ADC.data(), thr5, nSamples, noiseStart, noiseEnd).size();
        }
    }

    for (int ch : activeChannels) {
        if (ch == TriggerCh) continue;
        DarkResult r;
        r.rate3mv  = (LiveTime > 0) ? hits3mv[ch]  / LiveTime : 0.0;
        r.rate5sig = (LiveTime > 0) ? hits5sig[ch] / LiveTime : 0.0;
        r.thr5sigMv = 5.0 * sigmaAdc[ch] * Config::ADC_to_mV;
        out[ch] = r;
    }

    delete[] ADC_buffer;
    file->Close();
    return out;
}

void DrawChannelPanel(TVirtualPad* pad, int ch, const char* chLabel) {
    pad->cd();
    pad->SetGridy();
    pad->SetLeftMargin(0.16);

    const int n = (int)kPoints.size();
    TH1D* h3mv  = new TH1D(Form("hDark3mv_ch%d", ch), "", n, 0, n);
    TH1D* h5sig = new TH1D(Form("hDark5sig_ch%d", ch), "", n, 0, n);

    for (int i = 0; i < n; ++i) {
        const auto& p = kPoints[i];
        TString path = Form("/home/precalkor/ADC/ADC_test/Data/RAW/Laser/precal_raw_kor_run_%s_%s.root", p.date, p.run);
        if (gSystem->AccessPathName(path)) {
            path = Form("/media/precalkor/HD-EDS-E/Data_Backup/RAW/Laser/precal_raw_kor_run_%s_%s.root", p.date, p.run);
        }
        auto rates = ComputeDarkRates(path);
        double r3 = 0, r5 = 0, thr5mv = 0;
        if (rates.count(ch)) { r3 = rates[ch].rate3mv; r5 = rates[ch].rate5sig; thr5mv = rates[ch].thr5sigMv; }

        h3mv->SetBinContent(i + 1, r3);
        h5sig->SetBinContent(i + 1, r5);
        h3mv->GetXaxis()->SetBinLabel(i + 1, p.label);
        std::cout << Form("%-18s ch%d  3mV_rate=%.1fHz  mu+5sig_rate=%.1fHz (thr=%.2fmV)",
                          p.label, ch, r3, r5, thr5mv) << std::endl;
    }

    h3mv->SetBarWidth(0.4);  h3mv->SetBarOffset(0.05);
    h3mv->SetFillColor(kGray + 1); h3mv->SetLineColor(kGray + 2);
    h5sig->SetBarWidth(0.4); h5sig->SetBarOffset(0.5);
    h5sig->SetFillColor(kRed); h5sig->SetLineColor(kRed + 2);

    double ymax = std::max(h3mv->GetMaximum(), h5sig->GetMaximum());
    h3mv->SetMaximum(ymax * 1.25); h3mv->SetMinimum(0);
    h3mv->SetTitle(Form("%s : Dark/Noise Rate, 3mV vs #mu+5#sigma;;Rate [Hz]", chLabel));
    h3mv->GetXaxis()->SetLabelFont(132); h3mv->GetXaxis()->SetLabelSize(0.05);
    h3mv->GetYaxis()->SetLabelFont(132); h3mv->GetYaxis()->SetTitleFont(132);
    h3mv->GetYaxis()->SetTitleSize(0.045); h3mv->GetYaxis()->SetLabelSize(0.04);

    h3mv->Draw("bar0");
    h5sig->Draw("bar0 same");

    for (int i = 0; i < n; ++i) {
        double x3 = h3mv->GetXaxis()->GetBinCenter(i + 1) - 0.2;
        double x5 = h5sig->GetXaxis()->GetBinCenter(i + 1) + 0.2;
        TLatex lat; lat.SetTextFont(132); lat.SetTextSize(0.03); lat.SetTextAlign(21);
        lat.DrawLatex(x3, h3mv->GetBinContent(i + 1) + ymax * 0.02, Form("%.0f", h3mv->GetBinContent(i + 1)));
        lat.DrawLatex(x5, h5sig->GetBinContent(i + 1) + ymax * 0.02, Form("%.0f", h5sig->GetBinContent(i + 1)));
    }

    TLegend* leg = new TLegend(0.18, 0.78, 0.6, 0.88);
    leg->SetTextFont(132); leg->SetTextSize(0.03); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(h3mv, "Fixed 3mV threshold", "f");
    leg->AddEntry(h5sig, "Pedestal #mu+5#sigma threshold", "f");
    leg->Draw();
}

void Draw_DarkThreshold_Compare() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");
    gStyle->SetTitleFontSize(0.05);

    TCanvas* c = new TCanvas("cDarkThreshold", "Dark rate: 3mV vs mu+5sigma", 1400, 650);
    c->Divide(2, 1);
    DrawChannelPanel(c->cd(1), 1, "Rot1 (EM6400)");
    DrawChannelPanel(c->cd(2), 2, "Rot2 (EL5150)");

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + "DarkThreshold_Compare_3mV_vs_5sigma.pdf";
    c->SaveAs(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
