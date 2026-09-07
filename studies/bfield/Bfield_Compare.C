// Bfield_Compare.C -- standalone, non-invasive. Compares dark rate vs RAW
// equipment tilt between two run blocks (e.g. B-field coil ON vs OFF), reusing
// the already-computed production dark counts (fast: no RAW re-scan needed).
// Split by mounting plane (rot) exactly like Noise_vs_Angle.C, since raw tilt
// (not theta_Ham) is the real physical scan axis.
//
// Run:  root -l -b -q 'Bfield_Compare.C("20260715",0,45,"ON","20260716",0,45,"OFF")'
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <map>
#include <vector>
#include <iostream>

struct Pt { int tilt, rot; double rate0, rate1, rate2, cnt0, cnt1, cnt2; };

static std::vector<Pt> LoadBlock(const TString& tag, int runA, int runB) {
    std::vector<Pt> out;
    for (int r = runA; r <= runB; ++r) {
        TString rp = Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", tag.Data(), r);
        TFile* fr = TFile::Open(rp); if (!fr || fr->IsZombie()) { if (fr) fr->Close(); continue; }
        TTree* ri = (TTree*)fr->Get("RunInfo");
        int tilt = 0, rot = 0;
        ri->SetBranchAddress("RawTiltAngle2", &tilt); ri->SetBranchAddress("RawRotateAngle2", &rot);
        ri->GetEntry(0);
        fr->Close();

        TString pp = Form("./Data/production/precal_prd_kor_run_%s_%03d.root", tag.Data(), r);
        TFile* fp = TFile::Open(pp); if (!fp || fp->IsZombie()) { if (fp) fp->Close(); continue; }
        TParameter<double>* r0 = (TParameter<double>*)fp->Get("NoiseCountRate_ch0");
        TParameter<double>* r1 = (TParameter<double>*)fp->Get("NoiseCountRate_ch1");
        TParameter<double>* r2 = (TParameter<double>*)fp->Get("NoiseCountRate_ch2");
        TParameter<Long64_t>* c0 = (TParameter<Long64_t>*)fp->Get("NoiseCount_ch0");
        TParameter<Long64_t>* c1 = (TParameter<Long64_t>*)fp->Get("NoiseCount_ch1");
        TParameter<Long64_t>* c2 = (TParameter<Long64_t>*)fp->Get("NoiseCount_ch2");
        Pt p; p.tilt = tilt; p.rot = rot;
        p.rate0 = r0 ? r0->GetVal() : 0; p.rate1 = r1 ? r1->GetVal() : 0; p.rate2 = r2 ? r2->GetVal() : 0;
        p.cnt0 = c0 ? (double)c0->GetVal() : 0; p.cnt1 = c1 ? (double)c1->GetVal() : 0; p.cnt2 = c2 ? (double)c2->GetVal() : 0;
        fp->Close();
        out.push_back(p);
    }
    return out;
}

void Bfield_Compare(const char* tagA = "20260715", int runA0 = 0, int runA1 = 45, const char* labelA = "B-field ON",
                    const char* tagB = "20260716", int runB0 = 0, int runB1 = 45, const char* labelB = "B-field OFF") {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132); gStyle->SetLabelFont(132, "xyz"); gStyle->SetTitleFont(132, "xyz");

    auto A = LoadBlock(tagA, runA0, runA1);
    auto B = LoadBlock(tagB, runB0, runB1);
    std::cout << "[INFO] " << labelA << " (" << tagA << "): " << A.size() << " runs\n";
    std::cout << "[INFO] " << labelB << " (" << tagB << "): " << B.size() << " runs\n";

    // Split by rot plane (first two distinct rot values encountered in A).
    std::map<int, int> rotPlane; int rotLabel[2] = {0, 0};
    auto planeOf = [&](int rot) -> int {
        if (!rotPlane.count(rot)) {
            if ((int)rotPlane.size() >= 2) return -1;
            int p = rotPlane.size(); rotPlane[rot] = p; rotLabel[p] = rot; return p;
        }
        return rotPlane[rot];
    };
    for (auto& p : A) planeOf(p.rot);   // seed plane map from the reference (ON) block

    TGraph *gA1[2], *gA2[2], *gB1[2], *gB2[2];
    for (int p = 0; p < 2; ++p) { gA1[p] = new TGraph(); gA2[p] = new TGraph(); gB1[p] = new TGraph(); gB2[p] = new TGraph(); }

    for (auto& p : A) { int pl = planeOf(p.rot); if (pl < 0) continue; gA1[pl]->SetPoint(gA1[pl]->GetN(), p.tilt, p.rate1); gA2[pl]->SetPoint(gA2[pl]->GetN(), p.tilt, p.rate2); }
    for (auto& p : B) { int pl = planeOf(p.rot); if (pl < 0) continue; gB1[pl]->SetPoint(gB1[pl]->GetN(), p.tilt, p.rate1); gB2[pl]->SetPoint(gB2[pl]->GetN(), p.tilt, p.rate2); }

    // Monitor (ch0) doesn't rotate, so it has no "plane" -- plot it once, all
    // runs together, vs the rotating stage's tilt (a control channel: does the
    // coil / stage angle also perturb the fixed monitor PMT's dark rate?).
    TGraph *gA0 = new TGraph(), *gB0 = new TGraph();
    for (auto& p : A) gA0->SetPoint(gA0->GetN(), p.tilt, p.rate0);
    for (auto& p : B) gB0->SetPoint(gB0->GetN(), p.tilt, p.rate0);

    auto style = [](TGraph* g, int col, int mk) { g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(1.3); g->SetLineWidth(2); };
    // ON = filled + solid line, OFF = open + dashed line, same color per PMT.
    for (int p = 0; p < 2; ++p) {
        style(gA1[p], kRed + 1, 20);  style(gA2[p], kBlue + 1, 20);
        style(gB1[p], kRed + 1, 24);  style(gB2[p], kBlue + 1, 24);
        gB1[p]->SetLineStyle(2); gB2[p]->SetLineStyle(2);
    }
    style(gA0, kGreen + 2, 20); style(gB0, kGreen + 2, 24); gB0->SetLineStyle(2);

    // Mean-rate %% difference (ON vs OFF), averaged over all points in each
    // series -- a single number to go with the per-angle scatter above.
    auto avgY = [](TGraph* g) -> double {
        if (!g || g->GetN() == 0) return 0;
        double s = 0; for (int i = 0; i < g->GetN(); ++i) s += g->GetY()[i];
        return s / g->GetN();
    };
    auto pctDiff = [&](TGraph* on, TGraph* off) -> double {
        double mOn = avgY(on), mOff = avgY(off);
        return (mOff != 0) ? (mOn - mOff) / mOff * 100.0 : 0;
    };
    double pct0 = pctDiff(gA0, gB0);
    double pct1[2], pct2[2];
    for (int p = 0; p < 2; ++p) { pct1[p] = pctDiff(gA1[p], gB1[p]); pct2[p] = pctDiff(gA2[p], gB2[p]); }
    std::cout << Form("[INFO] Mean-rate ON vs OFF:  Mon. %+.1f%%   |  rot=%d plane: Rot#1 %+.1f%%, Rot#2 %+.1f%%"
                       "  |  rot=%d plane: Rot#1 %+.1f%%, Rot#2 %+.1f%%",
                       pct0, rotLabel[0], pct1[0], pct2[0], rotLabel[1], pct1[1], pct2[1]) << std::endl;

    TCanvas* c = new TCanvas("c_bfc", "Bfield Compare", 2100, 700);
    c->Divide(3, 1, 0.010, 0.001);

    c->cd(1); gPad->SetGrid(); gPad->SetLeftMargin(0.12); gPad->SetTopMargin(0.10);
    TMultiGraph* mg0 = new TMultiGraph();
    if (gA0->GetN()) mg0->Add(gA0, "P"); if (gB0->GetN()) mg0->Add(gB0, "P");
    mg0->SetTitle("Dark Rate vs raw tilt  (Monitor, EM2740);Raw stage tilt #theta [degree];Dark Rate [Hz]");
    mg0->Draw("A");
    mg0->GetXaxis()->SetTitleSize(0.05); mg0->GetYaxis()->SetTitleSize(0.05);
    mg0->GetXaxis()->SetLabelSize(0.045); mg0->GetYaxis()->SetLabelSize(0.045);
    // Sits in the gap between the ON and OFF traces (not the top, which the ON
    // series occupies), so it doesn't overlap either line.
    TLegend* l0 = new TLegend(0.14, 0.46, 0.60, 0.60); l0->SetTextSize(0.036); l0->SetBorderSize(0); l0->SetFillStyle(0);
    if (gA0->GetN()) l0->AddEntry(gA0, Form("Mon. %s", labelA), "p");
    if (gB0->GetN()) l0->AddEntry(gB0, Form("Mon. %s", labelB), "p");
    l0->Draw();
    { TLatex* t = new TLatex(0.14, 0.40, Form("ON vs OFF (mean): %+.1f%%", pct0)); t->SetNDC(); t->SetTextSize(0.038); t->SetTextFont(132); t->Draw(); }

    for (int p = 0; p < 2; ++p) {
        c->cd(p + 2); gPad->SetGrid(); gPad->SetLeftMargin(0.12); gPad->SetTopMargin(0.10);
        TMultiGraph* mg = new TMultiGraph();
        if (gA1[p]->GetN()) mg->Add(gA1[p], "P"); if (gA2[p]->GetN()) mg->Add(gA2[p], "P");
        if (gB1[p]->GetN()) mg->Add(gB1[p], "P"); if (gB2[p]->GetN()) mg->Add(gB2[p], "P");
        mg->SetTitle(Form("Dark Rate vs raw tilt  (rot=%d plane);Raw stage tilt #theta [degree];Dark Rate [Hz]", rotLabel[p]));
        mg->Draw("A");
        mg->GetXaxis()->SetTitleSize(0.05); mg->GetYaxis()->SetTitleSize(0.05);
        mg->GetXaxis()->SetLabelSize(0.045); mg->GetYaxis()->SetLabelSize(0.045);
        TLegend* l = new TLegend(0.14, 0.68, 0.75, 0.89); l->SetNColumns(2); l->SetTextSize(0.032); l->SetBorderSize(0); l->SetFillStyle(0);
        if (gA1[p]->GetN()) l->AddEntry(gA1[p], Form("Rot#1 %s", labelA), "p");
        if (gA2[p]->GetN()) l->AddEntry(gA2[p], Form("Rot#2 %s", labelA), "p");
        if (gB1[p]->GetN()) l->AddEntry(gB1[p], Form("Rot#1 %s", labelB), "p");
        if (gB2[p]->GetN()) l->AddEntry(gB2[p], Form("Rot#2 %s", labelB), "p");
        l->Draw();
        TLatex* t1 = new TLatex(0.14, 0.60, Form("Rot#1 ON vs OFF (mean): %+.1f%%", pct1[p])); t1->SetNDC();
        t1->SetTextSize(0.030); t1->SetTextFont(132); t1->SetTextColor(kRed + 1); t1->Draw();
        TLatex* t2 = new TLatex(0.14, 0.55, Form("Rot#2 ON vs OFF (mean): %+.1f%%", pct2[p])); t2->SetNDC();
        t2->SetTextSize(0.030); t2->SetTextFont(132); t2->SetTextColor(kBlue + 1); t2->Draw();
    }
    gSystem->mkdir("./Data/image/Noise", kTRUE);
    TString out = Form("./Data/image/Noise/Bfield_Compare_%s_vs_%s.png", tagA, tagB);
    c->SaveAs(out);
    std::cout << "[OK] " << out << std::endl;
}
