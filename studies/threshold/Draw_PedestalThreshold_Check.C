// Usage:  root -l -b -q Draw_PedestalThreshold_Check.C

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include <vector>
#include <cmath>

const double ADC_TO_MV = 0.1220703125;
const double THRESHOLD_MV = 3.0;

struct HVPoint { const char* date; const char* run; const char* label; int color; };

static const std::vector<HVPoint> kPoints = {
    {"20260829", "111", "-200V", kRed},
    {"20260829", "211", "-150V", kOrange + 1},
    {"20260829", "011", "-100V", kMagenta + 1},
    {"20260827", "011", "-50V",  kGreen + 2},
    {"20260815", "011", "Recom.", kBlack},
    {"20260827", "111", "+50V",  kCyan + 2},
    {"20260830", "011", "+100V", kAzure + 2},
    {"20260828", "011", "+150V", kViolet + 1},
    {"20260828", "111", "+200V", kBlue + 1},
};

void DrawPedestalPanel(TVirtualPad* pad, int ch, const char* chLabel) {
    pad->cd();
    pad->SetLogy();
    pad->SetGridx(); pad->SetGridy();
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(0.15);
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(1);

    TLegend* leg = new TLegend(0.55, 0.52, 0.94, 0.92);
    leg->SetTextFont(132); leg->SetTextSize(0.024); leg->SetBorderSize(1);

    bool first = true;
    double drawMax = 0;
    double sum5sig = 0;
    int n5sig = 0;
    std::vector<TH1D*> keep;

    for (const auto& p : kPoints) {
        TString path = Form("./Data/FinalResult/precal_result_kor_run_%s_%s.root", p.date, p.run);
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        TH1D* h = (TH1D*)f->Get(Form("Max_ch%d", ch));
        h->SetDirectory(0);
        f->Close();

        h->SetLineColor(p.color);
        h->SetLineWidth(1);
        h->SetTitle(Form("%s: Pedestal Threshold Check;ADC Counts;Entries", chLabel));
        h->GetXaxis()->SetRangeUser(-20, 200);
        h->GetXaxis()->SetLabelFont(132); h->GetXaxis()->SetTitleFont(132);
        h->GetYaxis()->SetLabelFont(132); h->GetYaxis()->SetTitleFont(132);
        h->DrawCopy(first ? "HIST" : "HIST SAME");
        first = false;
        if (h->GetMaximum() > drawMax) drawMax = h->GetMaximum();
        leg->AddEntry(h, p.label, "l");
        keep.push_back(h);

        TF1* fit = new TF1(Form("pedfit_%s_ch%d", p.label, ch), "gaus", -10, 8);
        h->Fit(fit, "RQN0", "", -10, 8);
        sum5sig += fit->GetParameter(1) + 5 * std::abs(fit->GetParameter(2));
        n5sig++;
        delete fit;
    }
    if (keep.empty()) return;
    keep[0]->SetMaximum(drawMax * 3.0);

    double avg5sigADC = sum5sig / n5sig;

    TLine* l5sig = new TLine(avg5sigADC, 0.5, avg5sigADC, drawMax * 3.0);
    l5sig->SetLineColor(kRed + 2); l5sig->SetLineStyle(1); l5sig->SetLineWidth(2);
    l5sig->Draw();

    double thrADC = THRESHOLD_MV / ADC_TO_MV;
    TLine* thr3 = new TLine(thrADC, 0.5, thrADC, drawMax * 3.0);
    thr3->SetLineColor(kGray + 2); thr3->SetLineStyle(2); thr3->SetLineWidth(2);
    thr3->Draw();

    leg->AddEntry(l5sig, Form("Pedestal #mu+5#sigma (%.2fmV)", avg5sigADC * ADC_TO_MV), "l");
    leg->AddEntry(thr3, Form("3mV threshold (%.1fmV)", THRESHOLD_MV), "l");
    leg->Draw();
}

void Draw_PedestalThreshold_Check() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");

    // 2 panels side by side: Rot1 (EM6400) and Rot2 (EL5150) -- was Rot1-only
    // (2026-08-31, user: "Pedstal Check에서 Rot2도 추가해줘").
    TCanvas* c = new TCanvas("cPedestalThreshold", "Pedestal+5sigma vs 3mV Threshold", 2000, 850);
    c->Divide(2, 1);
    DrawPedestalPanel(c->cd(1), 1, "Rot1 (EM6400)");
    DrawPedestalPanel(c->cd(2), 2, "Rot2 (EL5150)");

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + "Pedestal5Sigma_vs_Valley.pdf";
    c->SaveAs(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
