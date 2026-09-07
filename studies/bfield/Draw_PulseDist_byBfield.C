// Usage:  root -l -b -q Draw_PulseDist_byBfield.C
//
// Same "Pulse Height / Time / Charge, overlaid across campaign points"
// distribution page as Draw_HVScan_Report.C's BuildDistributionPage(), just
// pointed at the B-field campaign (2026-08-31~09-01) instead of the 9-point
// HV scan -- 2026-09-05, user: "그 Center에 대한 Distribution Overlay한
// 그것도 동일하게 진행해줘". CenterRun(b)=b.lo+11 still applies unchanged:
// these are the same 46-run angle-sweep block layout as the HV campaign
// (unlike the later single-run Gain-curve test blocks, which needed their
// own CenterRun-free script).

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include <vector>
#include <iostream>

struct BfieldBlock {
    const char* date;
    int lo, hi;
    const char* label;
};

// Rot1's "B=0" label is per calibration intent (known current for B=0 at
// Rot1's position), not the mis-oriented sensor's own readout -- see
// project memory / 2026-09-04 chat with user.
static const std::vector<BfieldBlock> kBlocks = {
    {"20260831", 0,   45,  "Nominal"},
    {"20260831", 100, 145, "B=0 (Rot1 target)"},
    {"20260901", 100, 145, "Bx=-70mG"},
    {"20260901", 200, 245, "Bx=-150mG"},
};

int kColors[] = {kBlack, kRed + 1, kAzure + 2, kViolet + 1};

TString BlockPath(const BfieldBlock& b, int run) {
    return Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", b.date, run);
}
int CenterRun(const BfieldBlock& b) { return b.lo + 11; }

struct DistDef { const char* branch; const char* label; double xlo, xhi; bool showThreshold; bool normalizePeak; };
static const std::vector<DistDef> kDists = {
    {"Max_ch",  "Pulse Height [ADC Counts]", -20, 400, true,  false},
    {"Diff_ch", "Time [ADC Samples] (peak-norm.)", 170, 220, false, true},
    {"Pico_ch", "Charge [pC]",               -2,  10,  false, false},
};

void DrawDistPad(TVirtualPad* pad, const DistDef& d, int ch, const char* chLabel) {
    pad->cd();
    if (!d.normalizePeak) pad->SetLogy();
    pad->SetGridx(); pad->SetGridy();
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(1);

    TLegend* leg = new TLegend(0.58, 0.55, 0.94, 0.90);
    leg->SetTextFont(132);
    leg->SetTextSize(0.028);
    leg->SetBorderSize(1);

    struct Loaded { TH1D* h; int blockIdx; };
    std::vector<Loaded> loaded;
    double autoLo = 1e18, autoHi = -1e18;

    for (size_t i = 0; i < kBlocks.size(); ++i) {
        const auto& b = kBlocks[i];
        int centerRun = CenterRun(b);
        if (centerRun > b.hi) continue;
        TString path = BlockPath(b, centerRun);
        if (gSystem->AccessPathName(path)) continue;
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        TH1D* h = (TH1D*)f->Get(Form("%s%d", d.branch, ch));
        if (!h) { f->Close(); continue; }
        h->SetDirectory(0);
        f->Close();

        h->GetYaxis()->UnZoom();

        h->SetLineColor(kColors[i % 4]);
        h->SetLineWidth(1);
        h->SetTitle(Form("%s : %s;%s;Entries", chLabel, d.label, d.label));
        h->GetXaxis()->SetLabelFont(132); h->GetXaxis()->SetTitleFont(132);
        h->GetYaxis()->SetLabelFont(132); h->GetYaxis()->SetTitleFont(132);
        h->GetXaxis()->SetLabelSize(0.038); h->GetXaxis()->SetTitleSize(0.042);
        h->GetYaxis()->SetLabelSize(0.038); h->GetYaxis()->SetTitleSize(0.042);

        if (d.xlo < d.xhi) {
            h->GetXaxis()->SetRangeUser(d.xlo, d.xhi);
        } else {
            int peakBin = h->GetMaximumBin();
            double peakX = h->GetXaxis()->GetBinCenter(peakBin);
            autoLo = std::min(autoLo, peakX - 30);
            autoHi = std::max(autoHi, peakX + 30);
            h->GetYaxis()->SetTitle("Normalized [a.u.]");
        }

        if (d.normalizePeak && h->GetMaximum() > 0) h->Scale(1.0 / h->GetMaximum());

        loaded.push_back({h, (int)i});
    }

    if (d.xlo >= d.xhi && autoHi > autoLo) {
        for (auto& L : loaded) L.h->GetXaxis()->SetRangeUser(autoLo, autoHi);
    }

    double drawMax = 0;
    for (auto& L : loaded) if (L.h->GetMaximum() > drawMax) drawMax = L.h->GetMaximum();

    std::vector<Loaded> drawOrder(loaded.rbegin(), loaded.rend());
    for (size_t k = 0; k < drawOrder.size(); ++k) {
        if (k == 0 && drawMax > 0) {
            drawOrder[k].h->SetMaximum(drawMax * (d.normalizePeak ? 1.25 : 3.0));
            if (d.normalizePeak) drawOrder[k].h->SetMinimum(0);
        }
        drawOrder[k].h->DrawCopy(k == 0 ? "HIST" : "HIST SAME");
    }

    for (auto& L : loaded) leg->AddEntry(L.h, kBlocks[L.blockIdx].label, "l");

    leg->Draw();
}

TCanvas* BuildDistributionPage() {
    TCanvas* c = new TCanvas("cDistBfield", "Pulse/Time/Charge, B-field campaign", 2000, 1900);
    c->Divide(3, 3);
    const char* chLabels[3] = {"Mon (EM2740)", "Rot1 (EM6400)", "Rot2 (EL5150)"};
    int padIdx = 1;
    for (const auto& d : kDists) {
        for (int ch = 0; ch < 3; ++ch) {
            DrawDistPad(c->cd(padIdx++), d, ch, chLabels[ch]);
        }
    }
    return c;
}

void Draw_PulseDist_byBfield() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFont(132, "");
    gStyle->SetTitleFontSize(0.040);

    TCanvas* c = BuildDistributionPage();
    TString outPath = "./Data/image/ScanReport/PulseDist_byBfield.pdf";
    c->Print(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
