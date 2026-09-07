// Usage:  root -l -b -q Draw_QEThreshold_Compare.C
//
// Counting QE compared between two threshold definitions:
//   1. Fixed 3mV (what prod_ntp_v7.C/read_ntp_v7.C actually use today)
//   2. Pedestal mu+5sigma (HV-independent, ~1.5mV -- guards against noise
//      only, doesn't chase the SPE peak the way an adaptive valley-finder
//      would, but that fragility -- valley detection latching onto a
//      statistical dip in a low-count bin -- is exactly why it was dropped
//      in favor of this simpler, more robust definition; see project memory
//      "gain-curve-spe-fit-decision". A 3rd "adaptive valley" bar was tried
//      and compared here too, 2026-08-31, but the user decided mu+5sigma
//      is the better of the two and to drop the valley variant.)
//
// This is a DIAGNOSTIC comparison, not a production change -- see project
// memory "pending-adaptive-threshold-fix" for the (still deferred) plan to
// actually change read_ntp_v7.C/prod_ntp_v7.C.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include <cmath>

const double ADC_TO_MV = 0.1220703125;

struct HVPoint {
    const char* date;
    const char* run;
    const char* label;
};

static const std::vector<HVPoint> kPoints = {
    {"20260829", "111", "-200V"},
    {"20260815", "011", "Recom."},
    {"20260828", "111", "+200V"},
};

// Same simple counting as the existing 3mV logic, just at a different fixed
// ADC threshold -- used for the 1.5mV comparison bar (2026-09-01, user:
// "QE 구한 1.5mV로 해야하지않을까?" -- the per-run mu+5sigma averages to
// ~1.5mV and barely varies across HV, see project memory
// "pending-adaptive-threshold-fix", so this checks whether just fixing the
// threshold at that average value gives basically the same answer as the
// per-run fit).
double CountingQE_FixedMv(TH1D* h, double thrMv) {
    double thrADC = thrMv / ADC_TO_MV;
    int thrBin = h->GetXaxis()->FindBin(thrADC);
    double nAbove = h->Integral(thrBin, h->GetNbinsX());
    double nTotal = h->Integral(1, h->GetNbinsX());
    return (nTotal > 0) ? nAbove / nTotal * 100.0 : 0.0;
}

double CountingQE_Mu5Sigma(TH1D* h, double& thrMvOut) {
    TF1* fit = new TF1("pedfit", "gaus", -10, 8);
    h->Fit(fit, "RQN0", "", -10, 8);
    double mu = fit->GetParameter(1), sigma = std::abs(fit->GetParameter(2));
    double thrADC = mu + 5 * sigma;
    thrMvOut = thrADC * ADC_TO_MV;

    int thrBin = h->GetXaxis()->FindBin(thrADC);
    double nAbove = h->Integral(thrBin, h->GetNbinsX());
    double nTotal = h->Integral(1, h->GetNbinsX());
    delete fit;
    return (nTotal > 0) ? nAbove / nTotal * 100.0 : 0.0;
}

void DrawChannelPanel(TVirtualPad* pad, int ch, const char* chLabel) {
    pad->cd();
    pad->SetGridy();
    pad->SetLeftMargin(0.14);

    const int n = (int)kPoints.size();
    TH1D* h3mv = new TH1D(Form("h3mv_ch%d", ch), "", n, 0, n);
    TH1D* h15  = new TH1D(Form("h15_ch%d", ch), "", n, 0, n);
    TH1D* h5sig = new TH1D(Form("h5sig_ch%d", ch), "", n, 0, n);

    for (int i = 0; i < n; ++i) {
        const auto& p = kPoints[i];
        TString path = Form("./Data/FinalResult/precal_result_kor_run_%s_%s.root", p.date, p.run);
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }

        TTree* tr = (TTree*)f->Get(Form("tree_ch%d", ch));
        double qe3mv = 0;
        if (tr) { tr->SetBranchAddress("relativeQE_raw", &qe3mv); tr->GetEntry(0); }

        TH1D* hMax = (TH1D*)f->Get(Form("Max_ch%d", ch));
        double thrMv = 0;
        double qe5sig = hMax ? CountingQE_Mu5Sigma(hMax, thrMv) : 0.0;
        double qe15 = hMax ? CountingQE_FixedMv(hMax, 1.5) : 0.0;
        f->Close();

        h3mv->SetBinContent(i + 1, qe3mv);
        h15->SetBinContent(i + 1, qe15);
        h5sig->SetBinContent(i + 1, qe5sig);
        h3mv->GetXaxis()->SetBinLabel(i + 1, p.label);

        std::cout << Form("%-18s ch%d  3mV_QE=%.3f%%  1.5mV(fixed)_QE=%.3f%%  mu+5sig_QE=%.3f%% (thr=%.2fmV)",
                          p.label, ch, qe3mv, qe15, qe5sig, thrMv) << std::endl;
    }

    h3mv->SetBarWidth(0.27);
    h3mv->SetBarOffset(0.02);
    h3mv->SetFillColor(kGray + 1);
    h3mv->SetLineColor(kGray + 2);

    h15->SetBarWidth(0.27);
    h15->SetBarOffset(0.36);
    h15->SetFillColor(kAzure + 1);
    h15->SetLineColor(kAzure + 2);

    h5sig->SetBarWidth(0.27);
    h5sig->SetBarOffset(0.70);
    h5sig->SetFillColor(kRed);
    h5sig->SetLineColor(kRed + 2);

    double ymax = std::max({h3mv->GetMaximum(), h15->GetMaximum(), h5sig->GetMaximum()});
    h3mv->SetMaximum(ymax * 1.25);
    h3mv->SetMinimum(0);
    h3mv->SetTitle(Form("%s : Counting QE, 3mV vs 1.5mV(fixed) vs #mu+5#sigma;;QE [%%]", chLabel));
    h3mv->GetXaxis()->SetLabelFont(132); h3mv->GetXaxis()->SetLabelSize(0.05);
    h3mv->GetYaxis()->SetLabelFont(132); h3mv->GetYaxis()->SetTitleFont(132);
    h3mv->GetYaxis()->SetTitleSize(0.045); h3mv->GetYaxis()->SetLabelSize(0.04);

    h3mv->Draw("bar0");
    h15->Draw("bar0 same");
    h5sig->Draw("bar0 same");

    for (int i = 0; i < n; ++i) {
        double x3 = h3mv->GetXaxis()->GetBinCenter(i + 1) - 0.34;
        double x15 = h15->GetXaxis()->GetBinCenter(i + 1) + 0.0;
        double x5 = h5sig->GetXaxis()->GetBinCenter(i + 1) + 0.34;
        TLatex lat;
        lat.SetTextFont(132); lat.SetTextSize(0.026); lat.SetTextAlign(21);
        lat.DrawLatex(x3, h3mv->GetBinContent(i + 1) + ymax * 0.02,
                     Form("%.2f%%", h3mv->GetBinContent(i + 1)));
        lat.DrawLatex(x15, h15->GetBinContent(i + 1) + ymax * 0.02,
                     Form("%.2f%%", h15->GetBinContent(i + 1)));
        lat.DrawLatex(x5, h5sig->GetBinContent(i + 1) + ymax * 0.02,
                     Form("%.2f%%", h5sig->GetBinContent(i + 1)));
    }

    TLegend* leg = new TLegend(0.15, 0.74, 0.62, 0.90);
    leg->SetTextFont(132); leg->SetTextSize(0.028); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(h3mv, "Fixed 3mV threshold", "f");
    leg->AddEntry(h15, "Fixed 1.5mV threshold", "f");
    leg->AddEntry(h5sig, "Pedestal #mu+5#sigma threshold", "f");
    leg->Draw();
}

void Draw_QEThreshold_Compare() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");
    gStyle->SetTitleFontSize(0.05);

    TCanvas* c = new TCanvas("cQEThreshold", "Counting QE: 3mV vs mu+5sigma", 1400, 650);
    c->Divide(2, 1);
    DrawChannelPanel(c->cd(1), 1, "Rot1 (EM6400)");
    DrawChannelPanel(c->cd(2), 2, "Rot2 (EL5150)");

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + "QEThreshold_Compare_3mV_vs_5sigma.pdf";
    c->SaveAs(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
