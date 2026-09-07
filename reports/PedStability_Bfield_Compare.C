// PedStability_Bfield_Compare.C -- standalone, non-invasive. Compares
// Pedestal Mean+-Sigma (mV) vs raw equipment tilt between two run blocks
// (B-field coil ON vs OFF), reading the PedMean_mV_ch%d / PedSigma_mV_ch%d
// TParameters. Drawn as TGraphErrors (Y=mean, error=sigma) styled like the
// Overlay's PedStability page: one pad PER PMT, each with its own tightly
// autoscaled Y-axis. Earlier versions shared one axis across Rot#1 (~1838 mV)
// and Rot#2 (~1852 mV) -- a 14 mV gap that swallowed the ~0.03-0.25 mV sigma
// differences we actually care about, exactly the "range too different"
// problem. Splitting by PMT (not by mounting plane) fixes that.
//
// Run:  root -l -b -q 'PedStability_Bfield_Compare.C("20260715",0,45,"ON","20260716",0,45,"OFF")'
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <vector>
#include <iostream>

struct PedPt { int tilt; double mean0, mean1, mean2, sig0, sig1, sig2; };

static std::vector<PedPt> LoadPedBlock(const TString& tag, int runA, int runB) {
    std::vector<PedPt> out;
    for (int r = runA; r <= runB; ++r) {
        TString rp = Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", tag.Data(), r);
        TFile* fr = TFile::Open(rp); if (!fr || fr->IsZombie()) { if (fr) fr->Close(); continue; }
        TTree* ri = (TTree*)fr->Get("RunInfo");
        int tilt = 0;
        if (ri) { ri->SetBranchAddress("RawTiltAngle2", &tilt); ri->GetEntry(0); }
        fr->Close();

        TString pp = Form("./Data/production/precal_prd_kor_run_%s_%03d.root", tag.Data(), r);
        TFile* fp = TFile::Open(pp); if (!fp || fp->IsZombie()) { if (fp) fp->Close(); continue; }
        auto getP = [&](const char* name) -> double {
            TParameter<double>* p = (TParameter<double>*)fp->Get(name);
            return p ? p->GetVal() : 0.0;
        };
        PedPt p; p.tilt = tilt;
        p.mean0 = getP("PedMean_mV_ch0"); p.mean1 = getP("PedMean_mV_ch1"); p.mean2 = getP("PedMean_mV_ch2");
        p.sig0  = getP("PedSigma_mV_ch0"); p.sig1  = getP("PedSigma_mV_ch1"); p.sig2  = getP("PedSigma_mV_ch2");
        fp->Close();
        out.push_back(p);
    }
    return out;
}

void PedStability_Bfield_Compare(const char* tagA = "20260715", int runA0 = 0, int runA1 = 45, const char* labelA = "B-field ON",
                                  const char* tagB = "20260716", int runB0 = 0, int runB1 = 45, const char* labelB = "B-field OFF") {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132); gStyle->SetLabelFont(132, "xyz"); gStyle->SetTitleFont(132, "xyz");

    auto A = LoadPedBlock(tagA, runA0, runA1);
    auto B = LoadPedBlock(tagB, runB0, runB1);
    std::cout << "[INFO] " << labelA << " (" << tagA << "): " << A.size() << " runs\n";
    std::cout << "[INFO] " << labelB << " (" << tagB << "): " << B.size() << " runs\n";

    // TGraphErrors: Y = Pedestal Mean, error = Pedestal Sigma -- same
    // mean+-sigma-in-one-point convention as the Overlay PedStability page.
    TGraphErrors *gA0 = new TGraphErrors(), *gB0 = new TGraphErrors();
    TGraphErrors *gA1 = new TGraphErrors(), *gB1 = new TGraphErrors();
    TGraphErrors *gA2 = new TGraphErrors(), *gB2 = new TGraphErrors();

    auto addPt = [](TGraphErrors* g, double x, double y, double ey) { int n = g->GetN(); g->SetPoint(n, x, y); g->SetPointError(n, 0, ey); };

    for (auto& p : A) { addPt(gA0, p.tilt, p.mean0, p.sig0); addPt(gA1, p.tilt, p.mean1, p.sig1); addPt(gA2, p.tilt, p.mean2, p.sig2); }
    for (auto& p : B) { addPt(gB0, p.tilt, p.mean0, p.sig0); addPt(gB1, p.tilt, p.mean1, p.sig1); addPt(gB2, p.tilt, p.mean2, p.sig2); }

    auto style = [](TGraphErrors* g, int col, int mk) { g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(1.3); g->SetLineWidth(2); };
    // ON = filled marker, OFF = open marker, same color per PMT.
    style(gA0, kGreen + 2, 20); style(gB0, kGreen + 2, 24);
    style(gA1, kRed + 1, 20);   style(gB1, kRed + 1, 24);
    style(gA2, kBlue + 1, 20);  style(gB2, kBlue + 1, 24);

    auto avgSig = [](TGraphErrors* g) -> double {
        if (!g || g->GetN() == 0) return 0;
        double s = 0; for (int i = 0; i < g->GetN(); ++i) s += g->GetEY()[i];
        return s / g->GetN();
    };
    auto pctDiff = [&](TGraphErrors* on, TGraphErrors* off) -> double {
        double mOn = avgSig(on), mOff = avgSig(off);
        return (mOff != 0) ? (mOn - mOff) / mOff * 100.0 : 0;
    };
    double sigPct0 = pctDiff(gA0, gB0), sigPct1 = pctDiff(gA1, gB1), sigPct2 = pctDiff(gA2, gB2);
    std::cout << Form("[INFO] Pedestal Sigma ON vs OFF:  Mon. %+.1f%%  |  Rot#1 %+.1f%%  |  Rot#2 %+.1f%%",
                       sigPct0, sigPct1, sigPct2) << std::endl;

    struct PadDef { const char* name; TGraphErrors* on; TGraphErrors* off; double pct; };
    PadDef pads[3] = {
        {"Monitor, EM2740", gA0, gB0, sigPct0},
        {"Rot#1, EM5370",   gA1, gB1, sigPct1},
        {"Rot#2, EL9590",   gA2, gB2, sigPct2},
    };

    // 3 columns (Monitor / Rot#1 / Rot#2) x 2 rows (ON on top, OFF on bottom).
    // ON and OFF get their OWN pad -- overlaying filled vs open markers on one
    // pad made them hard to tell apart. Each COLUMN shares one Y-range (taken
    // from both ON and OFF including error bars) so a straight up-down glance
    // compares the same PMT's stability with the coil on vs off on an identical
    // axis.
    auto rangeOf = [](TGraphErrors* g, double& lo, double& hi) {
        for (int k = 0; k < g->GetN(); ++k) {
            double y = g->GetY()[k], ey = g->GetEY()[k];
            if (y + ey > hi) hi = y + ey;
            if (y - ey < lo) lo = y - ey;
        }
    };
    double colLo[3], colHi[3];
    for (int i = 0; i < 3; ++i) {
        colLo[i] = 1e9; colHi[i] = -1e9;
        rangeOf(pads[i].on, colLo[i], colHi[i]);
        rangeOf(pads[i].off, colLo[i], colHi[i]);
        double diff = colHi[i] - colLo[i];
        colLo[i] -= diff * 0.15; colHi[i] += diff * 0.15;
    }

    TCanvas* c = new TCanvas("c_pbc", "Pedestal Stability Bfield Compare", 2100, 1300);
    c->Divide(3, 2, 0.010, 0.012);

    // row 0 = ON, row 1 = OFF
    for (int row = 0; row < 2; ++row) {
        bool isOn = (row == 0);
        const char* rowLabel = isOn ? labelA : labelB;
        for (int i = 0; i < 3; ++i) {
            c->cd(row * 3 + i + 1); gPad->SetGrid(); gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.11);
            TGraphErrors* g = isOn ? pads[i].on : pads[i].off;
            TMultiGraph* mg = new TMultiGraph();
            if (g->GetN()) mg->Add(g, "PZ");
            mg->SetTitle(Form("%s  |  %s;Raw stage tilt #theta [degree];Pedestal Mean #pm #sigma [mV]", pads[i].name, rowLabel));
            mg->Draw("A");
            if (colHi[i] > colLo[i]) { mg->SetMinimum(colLo[i]); mg->SetMaximum(colHi[i]); }
            mg->GetXaxis()->SetTitleSize(0.05); mg->GetYaxis()->SetTitleSize(0.05);
            mg->GetXaxis()->SetLabelSize(0.045); mg->GetYaxis()->SetLabelSize(0.045);
            gPad->Modified(); gPad->Update();
            // Print the sigma ON-vs-OFF % diff once per column, on the OFF (bottom)
            // pad so it sits with the reference series.
            if (!isOn) {
                TLatex* t = new TLatex(0.18, 0.83, Form("#sigma ON vs OFF (mean): %+.1f%%", pads[i].pct));
                t->SetNDC(); t->SetTextSize(0.040); t->SetTextFont(132); t->SetTextColor(kGray + 3); t->Draw();
            }
        }
    }

    gSystem->mkdir("./Data/image/Noise", kTRUE);
    TString out = Form("./Data/image/Noise/PedStability_Bfield_Compare_%s_vs_%s.png", tagA, tagB);
    c->SaveAs(out);
    std::cout << "[OK] " << out << std::endl;
}
