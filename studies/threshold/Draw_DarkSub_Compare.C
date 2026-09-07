// Usage:  root -l -b -q Draw_DarkSub_Compare.C
//
// One dataset (Recom. block, 20260815_000~045), comparing, with the SAME
// Timing cut applied in both cases (isolates dark-subtraction only --
// relativeQE_raw is NOT usable for this since it also drops the Timing cut
// entirely, confounding two different effects; see project memory
// "pending-adaptive-threshold-fix", 2026-09-01 correction):
//   1. "Existing": relativeQE = (N_all - dark)/N_total*100 (Timing cut ON, dark-subtracted)
//   2. "No dark sub": N_all/N_total*100, recovered as
//      relativeQE / (1 - dark_frac_phc/100) since dark_frac_phc is already
//      stored as (subtracted dark count)/N_all*100 -- Timing cut ON, no
//      dark subtraction.
// vs angle, for Mon/Rot1/Rot2, plus the Mon-Norm ratio (Test/Mon) for both.
//
// DIAGNOSTIC ONLY.

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TString.h"
#include <vector>
#include <iostream>

struct Pt { double angle; double qeCorr, qeCorrErr, qeRaw, qeRawErr; };

std::vector<Pt> LoadChannel(const char* date, int lo, int hi, int ch, bool useRot2Tilt) {
    std::vector<Pt> pts;
    for (int i = lo; i <= hi; ++i) {
        TString path = Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", date, i);
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }

        TTree* info = (TTree*)f->Get("RunInfo");
        int tilt2 = 0; double angle = 0;
        if (info) { info->SetBranchAddress("RawTiltAngle2", &tilt2); info->GetEntry(0); angle = tilt2; }

        TTree* tr = (TTree*)f->Get(Form("tree_ch%d", ch));
        double qeCorr = 0, qeCorrErr = 0, darkFrac = 0;
        if (tr) {
            tr->SetBranchAddress("relativeQE", &qeCorr);
            tr->SetBranchAddress("relativeQE_err", &qeCorrErr);
            tr->SetBranchAddress("dark_frac_phc", &darkFrac);
            tr->GetEntry(0);
        }
        f->Close();
        // "No dark sub" but SAME Timing cut: N_all/N_total*100, recovered
        // from relativeQE (= (N_all-dark)/N_total*100) and dark_frac_phc
        // (= dark/N_all*100) as relativeQE / (1 - darkFrac/100). Error
        // scales by the same factor (dark_frac_phc's own error not stored,
        // treated as exact here -- fine for a diagnostic comparison).
        double keepFrac = 1.0 - darkFrac / 100.0;
        double qeRaw = (keepFrac > 0) ? qeCorr / keepFrac : 0.0;
        double qeRawErr = (keepFrac > 0) ? qeCorrErr / keepFrac : 0.0;
        pts.push_back({angle, qeCorr, qeCorrErr, qeRaw, qeRawErr});
    }
    return pts;
}

void Draw_DarkSub_Compare() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");

    const char* date = "20260815";
    int lo = 0, hi = 45;

    auto mon  = LoadChannel(date, lo, hi, 0, true);
    auto rot1 = LoadChannel(date, lo, hi, 1, true);
    auto rot2 = LoadChannel(date, lo, hi, 2, true);

    int n = mon.size();
    std::vector<double> ang(n), angErr(n, 0.0);
    std::vector<double> mon1c(n), mon1r(n), rot1c(n), rot1r(n), rot2c(n), rot2r(n);
    std::vector<double> rot1cErr(n), rot1rErr(n), rot2cErr(n), rot2rErr(n);
    std::vector<double> norm1c(n), norm1r(n), norm2c(n), norm2r(n);
    std::vector<double> norm1cErr(n), norm1rErr(n), norm2cErr(n), norm2rErr(n);
    for (int i = 0; i < n; ++i) {
        ang[i] = mon[i].angle;
        mon1c[i]  = mon[i].qeCorr;  mon1r[i]  = mon[i].qeRaw;
        rot1c[i]  = rot1[i].qeCorr; rot1r[i]  = rot1[i].qeRaw;
        rot2c[i]  = rot2[i].qeCorr; rot2r[i]  = rot2[i].qeRaw;
        rot1cErr[i] = rot1[i].qeCorrErr; rot1rErr[i] = rot1[i].qeRawErr;
        rot2cErr[i] = rot2[i].qeCorrErr; rot2rErr[i] = rot2[i].qeRawErr;

        norm1c[i] = (mon[i].qeCorr != 0) ? rot1[i].qeCorr / mon[i].qeCorr : 0;
        norm1r[i] = (mon[i].qeRaw  != 0) ? rot1[i].qeRaw  / mon[i].qeRaw  : 0;
        norm2c[i] = (mon[i].qeCorr != 0) ? rot2[i].qeCorr / mon[i].qeCorr : 0;
        norm2r[i] = (mon[i].qeRaw  != 0) ? rot2[i].qeRaw  / mon[i].qeRaw  : 0;

        auto ratioErr = [](double num, double numErr, double den, double denErr) {
            if (num == 0 || den == 0) return 0.0;
            double relN = numErr / num, relD = denErr / den;
            return std::abs(num / den) * std::sqrt(relN*relN + relD*relD);
        };
        norm1cErr[i] = ratioErr(rot1[i].qeCorr, rot1[i].qeCorrErr, mon[i].qeCorr, mon[i].qeCorrErr);
        norm1rErr[i] = ratioErr(rot1[i].qeRaw,  rot1[i].qeRawErr,  mon[i].qeRaw,  mon[i].qeRawErr);
        norm2cErr[i] = ratioErr(rot2[i].qeCorr, rot2[i].qeCorrErr, mon[i].qeCorr, mon[i].qeCorrErr);
        norm2rErr[i] = ratioErr(rot2[i].qeRaw,  rot2[i].qeRawErr,  mon[i].qeRaw,  mon[i].qeRawErr);

        std::cout << Form("angle=%5.1f  Rot1 corr=%.3f%% raw=%.3f%% (diff %.2f%%)   Rot2 corr=%.3f%% raw=%.3f%% (diff %.2f%%)",
                          ang[i], rot1c[i], rot1r[i], 100.0*(rot1r[i]-rot1c[i])/rot1c[i],
                          rot2c[i], rot2r[i], 100.0*(rot2r[i]-rot2c[i])/rot2c[i]) << std::endl;
    }

    TCanvas* c = new TCanvas("cDarkSub", "Dark-subtracted vs raw QE", 1400, 950);
    TPad* padGrid = new TPad("padGrid", "", 0, 0.08, 1, 1);
    TPad* padLeg  = new TPad("padLeg", "", 0, 0, 1, 0.08);
    padGrid->Draw(); padLeg->Draw();
    padGrid->Divide(2, 2);

    TGraphErrors *gcSaved = nullptr, *grSaved = nullptr;

    auto mkPad = [&](TVirtualPad* pad, const char* title, const char* yTitle,
                      std::vector<double>& yc, std::vector<double>& ycErr,
                      std::vector<double>& yr, std::vector<double>& yrErr) {
        pad->cd();
        TGraphErrors* gc = new TGraphErrors(n, ang.data(), yc.data(), angErr.data(), ycErr.data());
        TGraphErrors* gr = new TGraphErrors(n, ang.data(), yr.data(), angErr.data(), yrErr.data());
        gc->SetLineColor(kBlue+1); gc->SetMarkerColor(kBlue+1); gc->SetMarkerStyle(20); gc->SetLineWidth(1);
        gr->SetLineColor(kRed+1);  gr->SetMarkerColor(kRed+1);  gr->SetMarkerStyle(21); gr->SetLineWidth(1);
        TMultiGraph* mg = new TMultiGraph();
        mg->Add(gc, "P"); mg->Add(gr, "P");
        mg->SetTitle(Form("%s;Tilt Angle [deg];%s", title, yTitle));
        mg->Draw("A");
        gcSaved = gc; grSaved = gr;
    };

    // Top row: counting-QE itself, in percent (N_phc[-dark]/N_total*100).
    // Bottom row: that same QE normalized to the Monitor(ch0)'s QE at the
    // same angle, i.e. QE_{Rot}/QE_{Mon} -- removes common-mode drift
    // (laser intensity, etc.), same quantity as MonNorm_CorrRawQE in
    // Draw_Overlay_Uniformity_v7.C, just Existing vs No-dark-sub here.
    mkPad(padGrid->cd(1), "Rot1 QE", "Counting QE [%]", rot1c, rot1cErr, rot1r, rot1rErr);
    mkPad(padGrid->cd(2), "Rot2 QE", "Counting QE [%]", rot2c, rot2cErr, rot2r, rot2rErr);
    mkPad(padGrid->cd(3), "Rot1 Mon-Norm QE", "QE_{Rot1}/QE_{Mon} [a.u.]", norm1c, norm1cErr, norm1r, norm1rErr);
    mkPad(padGrid->cd(4), "Rot2 Mon-Norm QE", "QE_{Rot2}/QE_{Mon} [a.u.]", norm2c, norm2cErr, norm2r, norm2rErr);

    // Single shared legend, bottom-center of the whole canvas, outside any
    // data pad so it never overlaps points.
    padLeg->cd();
    TLegend* legAll = new TLegend(0.3, 0.05, 0.7, 0.95);
    legAll->SetNColumns(2);
    legAll->SetBorderSize(0); legAll->SetFillStyle(0); legAll->SetTextFont(132); legAll->SetTextSize(0.5);
    legAll->AddEntry(gcSaved, "Existing (Timing cut + dark-subtracted)", "p");
    legAll->AddEntry(grSaved, "Timing cut only, no dark sub", "p");
    legAll->Draw();

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + Form("DarkSub_Compare_%s_%03d_%03d.pdf", date, lo, hi);
    c->SaveAs(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
