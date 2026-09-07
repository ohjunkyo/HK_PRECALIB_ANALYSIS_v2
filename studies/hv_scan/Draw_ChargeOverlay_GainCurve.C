// Usage:  root -l -b -q Draw_ChargeOverlay_GainCurve.C
//
// Charge (Pico_ch) distribution of the exact center runs used by the Gain
// Curve page in Draw_HVScan_Report.C, one page, Rot1/Rot2 only (Monitor
// excluded per user request, 2026-09-01).

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include <vector>

struct HVPoint { const char* date; const char* run; const char* label; int color; };
static const std::vector<HVPoint> kPoints = {
    // 2026-09-01: -200V/-100V were swapped here (same bug found earlier in
    // Draw_GainCurve_byAngle.C -- verified against real RunInfo HV2/HV3, see
    // project memory "pending-adaptive-threshold-fix"). run_111 is the true
    // -200V (HV2=1479, offset -201V from Recom.=1680), run_011 is -100V
    // (HV2=1579, offset -101V).
    {"20260829", "111", "-200V", kRed},
    {"20260829", "211", "-150V", kOrange+1},
    {"20260829", "011", "-100V", kPink+7},
    {"20260827", "011", "-50V",  kGreen+2},
    {"20260815", "011", "Recom.", kBlack},
    {"20260827", "111", "+50V",  kCyan+2},
    {"20260830", "011", "+100V", kAzure+2},
    {"20260828", "011", "+150V", kViolet+1},
    {"20260828", "111", "+200V", kBlue+1},
};

// Grid version: one pad per HV point (3x3), each showing that run's own
// Charge distribution alone -- not overlaid (user, 2026-09-01: "패드가
// 여러개인거야 저 상태 그대로", i.e. one pad per run like the individual
// per-run PNGs, not a single overlaid plot).
void DrawChargeGrid(int ch, const char* chLabel) {
    TCanvas* c = new TCanvas(Form("cChargeGrid_ch%d", ch), Form("Charge Distribution ch%d", ch), 1500, 1400);
    c->Divide(3, 3);

    for (size_t i = 0; i < kPoints.size(); ++i) {
        const auto& p = kPoints[i];
        TVirtualPad* pad = c->cd(i + 1);
        pad->SetLogy();
        pad->SetGridx(); pad->SetGridy();
        pad->SetBottomMargin(0.13); pad->SetLeftMargin(0.13);

        TString path = Form("./Data/FinalResult/precal_result_kor_run_%s_%s.root", p.date, p.run);
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        TH1D* h = (TH1D*)f->Get(Form("Pico_ch%d", ch));
        if (!h) { f->Close(); continue; }
        h->SetDirectory(0);
        h->GetYaxis()->UnZoom();
        f->Close();

        h->SetLineColor(p.color); h->SetLineWidth(2);
        h->SetTitle(Form("%s (%s_%s);Charge [pC];Entries", p.label, p.date, p.run));
        h->GetXaxis()->SetRangeUser(-2, 10);
        h->GetXaxis()->SetLabelFont(132); h->GetXaxis()->SetTitleFont(132);
        h->GetYaxis()->SetLabelFont(132); h->GetYaxis()->SetTitleFont(132);
        h->GetXaxis()->SetLabelSize(0.05); h->GetYaxis()->SetLabelSize(0.05);
        h->DrawCopy("HIST");
    }

    c->cd(0);
    TLatex title;
    title.SetTextFont(132); title.SetTextSize(0.025); title.SetTextAlign(21);
    title.DrawLatexNDC(0.5, 0.985, Form("%s: Charge Distribution per HV point (Gain Curve center runs)", chLabel));

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + Form("ChargeGrid_GainCurve_ch%d.pdf", ch);
    c->SaveAs(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}

void Draw_ChargeOverlay_GainCurve() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");

    DrawChargeGrid(1, "Rot1 (EM6400)");
    DrawChargeGrid(2, "Rot2 (EL5150)");
}
