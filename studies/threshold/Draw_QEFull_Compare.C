// Usage:  root -l -b -q Draw_QEFull_Compare.C
//
// The REAL comparison (user, 2026-09-01, correcting an earlier version that
// only re-cut a histogram): counting QE computed exactly the way production
// does -- amplitude threshold + Timing cut + Dark subtraction -- once with
// the fixed 3mV threshold (= what read_ntp_v7.C actually uses today) and
// once with 1.5mV (proxy for the per-run pedestal mu+5sigma threshold,
// already shown to be numerically equivalent -- see project memory
// "pending-adaptive-threshold-fix"). Everything else (Timing cut window,
// dark subtraction) is IDENTICAL between the two -- only the amplitude
// threshold differs, so this isolates that one variable properly (unlike
// relativeQE_raw, which also drops the Timing cut).
//
// Re-loops each run's own per-event tree (production file's tree_ch{ch}:
// max/diff/pedestal/LiveTime), replicating read_ntp_v7.C's N_all + dark-sub
// logic verbatim. DIAGNOSTIC ONLY -- does not touch read_ntp_v7.C.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include <cmath>

const double ADC_TO_MV = 0.1220703125;

// Whole block (all ~46 angles), not just the center run -- stacking all
// angles' N_all/N_total/dark counts gives one aggregate QE per HV point with
// far higher statistics than a single center run (user, 2026-09-01: "모든
// 각도에 대해서 히스토그램을 쌓아서 Overlay해줘"). Block/label mapping
// verified against real RunInfo HV2/HV3 (see the -200V/-100V swap bug found
// earlier the same day, project memory "pending-adaptive-threshold-fix").
struct HVBlock { const char* date; int lo, hi; const char* label; };
static const std::vector<HVBlock> kBlocks = {
    {"20260829", 100, 145, "-200V"},
    {"20260815", 0,   45,  "Recom."},
    {"20260828", 100, 145, "+200V"},
};

// Replicates read_ntp_v7.C's N_all + dark-subtraction block exactly, for a
// given amplitude threshold (mV). Timing cut window uses the same fixed
// 20-sample width, recentered on the Diff histogram's own peak bin (a close
// stand-in for the ex-Gaussian meanFit -- the window is wide enough that a
// few-sample difference in centering barely moves N_all).
// One run's contribution to the aggregate: N_all (at this threshold),
// N_total, NoiseCount, TotalNoiseTime, winSeqS (same Timing-cut width for
// every run so winSeqS is really a per-run constant, but keep it per-call
// for clarity) -- summed across the whole block by the caller.
bool CountingQE_OneRun(const char* prodPath, int ch, double thrMv,
                        double& N_all, double& N_total, double& NoiseCount,
                        double& TotalNoiseTime, double& winSeqS) {
    TFile* f = TFile::Open(prodPath, "READ");
    if (!f || f->IsZombie()) { if (f) f->Close(); return false; }

    TTree* tr = (TTree*)f->Get(Form("tree_ch%d", ch));
    TH1D* hDiff = (TH1D*)f->Get(Form("Diff_ch%d", ch));
    TParameter<Long64_t>* pNC = (TParameter<Long64_t>*)f->Get(Form("NoiseCount_ch%d", ch));
    TParameter<double>* pTNT = (TParameter<double>*)f->Get(Form("TotalNoiseLiveTime_ch%d", ch));
    if (!tr || !hDiff || !pNC || !pTNT) { f->Close(); return false; }

    double meanMax = hDiff->GetXaxis()->GetBinCenter(hDiff->GetMaximumBin());
    double tLow = meanMax - 10.0, tHigh = meanMax + 10.0;   // fixed 20-sample width, recentered

    double tMax = 0, tDiff = 0;
    tr->SetBranchAddress("max", &tMax);
    tr->SetBranchAddress("diff", &tDiff);
    long long nEntries = tr->GetEntries();

    double ADCthreshold = thrMv / ADC_TO_MV;
    N_total = nEntries; N_all = 0;
    for (long long ev = 0; ev < nEntries; ++ev) {
        tr->GetEntry(ev);
        if (tMax > ADCthreshold && (tDiff > tLow && tDiff < tHigh)) N_all++;
    }

    NoiseCount = (double)pNC->GetVal();
    TotalNoiseTime = pTNT->GetVal();
    winSeqS = (tHigh - tLow) * 2.0 * 1e-9;

    f->Close();
    return true;
}

// Aggregate over every run (angle) in the block: sum N_all/N_total across
// all angles, and sum the dark-subtraction terms the SAME way production
// combines them for one run (expectedNoise = N_total*winSeqS*NoiseRate,
// NoiseRate = NoiseCount/TotalNoiseTime) -- summing the per-run expected
// noise counts (not just NoiseCount) keeps unequal per-run live times
// correctly weighted.
double CountingQE_Block(const char* date, int lo, int hi, int ch, double thrMv) {
    double sumNall = 0, sumNtotal = 0, sumExpectedNoise = 0;
    for (int run = lo; run <= hi; ++run) {
        TString prodPath = Form("./Data/production/precal_prd_kor_run_%s_%03d.root", date, run);
        double N_all, N_total, NoiseCount, TotalNoiseTime, winSeqS;
        if (!CountingQE_OneRun(prodPath, ch, thrMv, N_all, N_total, NoiseCount, TotalNoiseTime, winSeqS)) continue;
        double noiseRate = (TotalNoiseTime > 0) ? NoiseCount / TotalNoiseTime : 0.0;
        double expectedNoise = N_total * winSeqS * noiseRate;
        sumNall += N_all; sumNtotal += N_total; sumExpectedNoise += expectedNoise;
    }
    double realSig = sumNall - sumExpectedNoise; if (realSig < 0) realSig = 0;
    return (sumNtotal > 0) ? (realSig / sumNtotal) * 100.0 : 0.0;
}

// Per-angle QE (single run, own dark-subtraction/timing-cut) -- same formula
// as CountingQE_Block but for exactly one run, not summed over the block.
double CountingQE_OneRunFull(const char* date, int run, int ch, double thrMv) {
    TString prodPath = Form("./Data/production/precal_prd_kor_run_%s_%03d.root", date, run);
    double N_all, N_total, NoiseCount, TotalNoiseTime, winSeqS;
    if (!CountingQE_OneRun(prodPath, ch, thrMv, N_all, N_total, NoiseCount, TotalNoiseTime, winSeqS)) return -1;
    double noiseRate = (TotalNoiseTime > 0) ? NoiseCount / TotalNoiseTime : 0.0;
    double expectedNoise = N_total * winSeqS * noiseRate;
    double realSig = N_all - expectedNoise; if (realSig < 0) realSig = 0;
    return (N_total > 0) ? (realSig / N_total) * 100.0 : 0.0;
}

// One page per HV block: histogram of the 46 per-angle QE values (X=QE[%],
// Y=Entries), 3mV vs 1.5mV overlaid (user, 2026-09-01: "그냥 X축을 QE, Y축을
// Entry로 하고... 3mV와 1.5mV 방법 두개에 대한 히스토그램").
void DrawQEDistPage(TCanvas* c, const HVBlock& b, int ch, const char* chLabel) {
    std::vector<double> qe3vals, qe15vals;
    for (int run = b.lo; run <= b.hi; ++run) {
        double qe3 = CountingQE_OneRunFull(b.date, run, ch, 3.0);
        double qe15 = CountingQE_OneRunFull(b.date, run, ch, 1.5);
        if (qe3 >= 0) qe3vals.push_back(qe3);
        if (qe15 >= 0) qe15vals.push_back(qe15);
    }
    if (qe3vals.empty()) return;

    double lo = 1e9, hi = -1e9;
    for (double v : qe3vals) { lo = std::min(lo, v); hi = std::max(hi, v); }
    for (double v : qe15vals) { lo = std::min(lo, v); hi = std::max(hi, v); }
    double pad = (hi - lo) * 0.1 + 0.01;
    lo -= pad; hi += pad;

    TH1D* h3mv = new TH1D(Form("hQEdist3mv_%s_%d_%d", b.date, ch, b.lo), "", 40, lo, hi);
    TH1D* h15  = new TH1D(Form("hQEdist15_%s_%d_%d", b.date, ch, b.lo), "", 40, lo, hi);
    for (double v : qe3vals) h3mv->Fill(v);
    for (double v : qe15vals) h15->Fill(v);

    c->cd();
    gPad->SetGridy();
    h3mv->SetLineColor(kRed + 1); h3mv->SetLineWidth(2); h3mv->SetFillColorAlpha(kRed, 0.4);
    h15->SetLineColor(kBlue + 1); h15->SetLineWidth(2); h15->SetFillColorAlpha(kBlue, 0.4);
    h3mv->SetTitle(Form("%s : %s -- per-angle QE distribution;QE [%%];Entries (angles)", chLabel, b.label));
    h3mv->GetXaxis()->SetLabelFont(132); h3mv->GetXaxis()->SetTitleFont(132);
    h3mv->GetYaxis()->SetLabelFont(132); h3mv->GetYaxis()->SetTitleFont(132);
    double ymax = std::max(h3mv->GetMaximum(), h15->GetMaximum());
    h3mv->SetMaximum(ymax * 1.3);
    h3mv->Draw("HIST");
    h15->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.62, 0.75, 0.92, 0.90);
    leg->SetTextFont(132); leg->SetTextSize(0.03); leg->SetBorderSize(1); leg->SetFillColor(kWhite);
    leg->AddEntry(h3mv, Form("3mV (mean %.3f%%)", h3mv->GetMean()), "f");
    leg->AddEntry(h15, Form("1.5mV (mean %.3f%%)", h15->GetMean()), "f");
    leg->Draw();

    std::cout << Form("%-8s ch%d  n=%zu  3mV mean=%.4f std=%.4f  |  1.5mV mean=%.4f std=%.4f",
                      b.label, ch, qe3vals.size(), h3mv->GetMean(), h3mv->GetStdDev(),
                      h15->GetMean(), h15->GetStdDev()) << std::endl;
}

void Draw_QEFull_Compare() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");
    gStyle->SetTitleFontSize(0.05);

    TCanvas* c = new TCanvas("cQEDist", "Per-angle QE distribution", 1100, 650);
    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + "QEDist_byAngle_3mV_vs_1p5mV.pdf";

    bool firstPage = true;
    for (const auto& b : kBlocks) {
        for (int ch : {1, 2}) {
            const char* chLabel = (ch == 1) ? "Rot1 (EM6400)" : "Rot2 (EL5150)";
            c->Clear();
            DrawQEDistPage(c, b, ch, chLabel);
            if (firstPage) { c->Print(outPath + "("); firstPage = false; }
            else c->Print(outPath);
        }
    }
    c->Print(outPath + ")");
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
