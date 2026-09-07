// Point-by-point reproducibility between two repeat scans of the same PMT,
// split by scan axis (X / Y) the same way Draw_Overlay_Uniformity_v7.C is --
// one page per axis, one column per PMT, Rep.1/Rep.2 overlaid on top and their
// ratio underneath.
//
// Question this answers
// ---------------------
//   "If I measure the same PMT at the same angle twice, how much do the two
//    numbers differ -- and is that difference just counting statistics, or is
//    there a systematic floor I cannot go below?"
//
// Method
// ------
// For every angle present in BOTH scans (matched within angleTol), take the two
// Monitor-normalized values v1 +- e1 and v2 +- e2 and form:
//
//   ratio = v2 / v1                        <- how reproducible, in relative terms
//   pull  = (v2 - v1) / sqrt(e1^2 + e2^2)  <- difference in units of its own expected error
//
// The pull distribution is the key object:
//   RMS(pull) ~ 1  -> the two scans agree to within counting statistics; the
//                     per-point statistical error IS the total uncertainty.
//   RMS(pull) > 1  -> the repeats disagree by MORE than statistics allows, i.e.
//                     there is a systematic (setup/stage/laser-state) component.
//
// The systematic component is then extracted the same way as in
// Draw_MonitorStability.C -- by subtracting the statistical variance from the
// observed variance of the ratio:
//
//   sigma_sys = sqrt(max(0, RMS(ratio-1)^2 - <expected relative error>^2))
//
// which is the number to quote as "reproducibility of this measurement".
//
// X and Y are kept separate on purpose: pooling them puts two independent
// measurements at (nominally) the same angle on top of each other, which makes
// the scatter look larger than it is for either axis alone.

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "repro_common.h"

void Draw_Reproducibility(TString tag1 = "20260803_100_145", TString tag2 = "20260803_200_245",
                           std::vector<TString> serialList = {"EM5370", "EL9590"},
                           TString metric = "MonNorm_CorrRawQE", double angleTol = 1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(22, "");
    gStyle->SetTitleFont(132, "XY"); gStyle->SetLabelFont(132, "XY");
    // Grey grid, matching Draw_Overlay_Uniformity_v7 / Draw_Stability_v1 --
    // ROOT's default grid is as dark as the data and competes with the
    // markers (2026-08-25).
    gStyle->SetGridColor(kGray + 1);

    const char* axes[2] = {"X", "Y"};
    int nSer = serialList.size();

    TString outPath = Form("./Data/image/Uniformity/Reproducibility_%s_vs_%s.pdf", tag1.Data(), tag2.Data());
    TCanvas* c = new TCanvas("c_repro", "Reproducibility", 1600, 950);
    c->Print(outPath + "[");

    TFile* f1 = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag1.Data()));
    TFile* f2 = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag2.Data()));
    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
        std::cout << "[ERROR] cannot open one of the tag files." << std::endl; return;
    }

    // Collected for the final pull-distribution page (all axes/PMTs together).
    std::vector<double> allPulls;

    for (int a = 0; a < 2; ++a) {
        c->Clear(); c->cd();

        TLatex* pageTitle = new TLatex(0.5, 0.965, Form("%s-axis   Reproducibility   ( %s  vs  %s )",
                                                          axes[a], tag1.Data(), tag2.Data()));
        pageTitle->SetNDC(); pageTitle->SetTextAlign(22); pageTitle->SetTextFont(22);
        pageTitle->SetTextSize(0.036); pageTitle->Draw();

        bool drewAny = false;

        for (int s = 0; s < nSer; ++s) {
            TString serial = serialList[s];

            TGraphErrors* g1 = (TGraphErrors*)f1->Get(Form("gr_%s_%s_%s", serial.Data(), axes[a], metric.Data()));
            TGraphErrors* g2 = (TGraphErrors*)f2->Get(Form("gr_%s_%s_%s", serial.Data(), axes[a], metric.Data()));
            if (!g1 || !g2) { std::cout << "[WARN] missing graph: " << serial << " " << axes[a] << std::endl; continue; }

            std::vector<MatchedPoint> pts;
            for (int i = 0; i < g1->GetN(); ++i) {
                double an1 = g1->GetX()[i], v1 = g1->GetY()[i], e1 = g1->GetEY()[i];
                if (v1 <= 0) continue;
                for (int j = 0; j < g2->GetN(); ++j) {
                    double an2 = g2->GetX()[j], v2 = g2->GetY()[j], e2 = g2->GetEY()[j];
                    if (v2 <= 0) continue;
                    if (std::abs(an1 - an2) > angleTol) continue;

                    MatchedPoint mp;
                    mp.angle = an1; mp.v1 = v1; mp.e1 = e1; mp.v2 = v2; mp.e2 = e2;
                    mp.ratio = v2 / v1;
                    double r1 = e1 / v1, r2 = e2 / v2;
                    mp.ratioErr = mp.ratio * std::sqrt(r1 * r1 + r2 * r2);
                    double combErr = std::sqrt(e1 * e1 + e2 * e2);
                    mp.pull = (combErr > 0) ? (v2 - v1) / combErr : 0;
                    pts.push_back(mp);
                    break;
                }
            }
            if (pts.size() < 2) { std::cout << "[WARN] too few matched points: " << serial << " " << axes[a] << std::endl; continue; }
            std::sort(pts.begin(), pts.end(), [](const MatchedPoint& x, const MatchedPoint& y){ return x.angle < y.angle; });
            for (auto& p : pts) allPulls.push_back(p.pull);

            ReproStats st = ComputeRepro(pts);
            drewAny = true;

            // ---- console ----
            std::cout << "\n\033[1;33m[ " << serial << " , " << axes[a] << "-axis ]\033[0m" << std::endl;
            std::cout << "  ----------------------------------------------------" << std::endl;
            std::cout << Form("  Matched points            : %d", st.n) << std::endl;
            std::cout << Form("  Stat. err/pt  Rep.1       : %.2f %%", st.meanRelErr1 * 100) << std::endl;
            std::cout << Form("  Stat. err/pt  Rep.2       : %.2f %%", st.meanRelErr2 * 100) << std::endl;
            std::cout << Form("  Mean ratio (Rep2/Rep1)    : %.4f  (bias %+.2f %%)", st.meanRatio, (st.meanRatio - 1) * 100) << std::endl;
            std::cout << Form("  Observed scatter of ratio : %.2f %%", st.obsSpread * 100) << std::endl;
            std::cout << Form("  Expected from statistics  : %.2f %%", st.expSpread * 100) << std::endl;
            std::cout << Form("  RMS(pull)                 : %.2f +- %.2f",
                              st.pullRMS, st.pullRMS / std::sqrt(2.0 * st.n)) << std::endl;
            if (st.sysResolved)
                std::cout << Form("  \033[1;36m=> Systematic            : %.2f %%\033[0m", st.sysComp * 100) << std::endl;
            else if (st.statOverestimated)
                std::cout << Form("  \033[1;36m=> Systematic            : < %.2f %% (bounded by observed scatter;"
                                  " stat. err appear overestimated)\033[0m", st.sysUL * 100) << std::endl;
            else
                std::cout << Form("  \033[1;36m=> Systematic            : < %.2f %% (95%% CL, not resolved)\033[0m",
                                  st.sysUL * 100) << std::endl;

            // ---- column layout: value pad on top, ratio pad below ----
            double xLoPad = (double)s / nSer, xHiPad = (double)(s + 1) / nSer;
            double headTop = 0.945, ratioH = 0.24;

            c->cd();
            TPad* colPad = new TPad(Form("col_%d_%d", a, s), "", xLoPad, 0.0, xHiPad, headTop);
            colPad->Draw(); colPad->cd();

            TPad* pV = new TPad(Form("pv_%d_%d", a, s), "", 0.0, ratioH, 1.0, 1.0);
            TPad* pR = new TPad(Form("pr_%d_%d", a, s), "", 0.0, 0.0, 1.0, ratioH);
            pV->SetLeftMargin(0.155); pV->SetRightMargin(0.04); pV->SetTopMargin(0.07); pV->SetBottomMargin(0.02);
            pR->SetLeftMargin(0.155); pR->SetRightMargin(0.04); pR->SetTopMargin(0.03); pR->SetBottomMargin(0.34);
            pV->SetGridx(); pV->SetGridy(); pR->SetGridx(); pR->SetGridy();
            pV->Draw(); pR->Draw();

            int n = st.n;
            std::vector<double> ang, r1v, r1e, r2v, r2e, ratioV, ratioE;
            for (auto& p : pts) {
                ang.push_back(p.angle);
                r1v.push_back(p.v1); r1e.push_back(p.e1);
                r2v.push_back(p.v2); r2e.push_back(p.e2);
                ratioV.push_back(p.ratio); ratioE.push_back(p.ratioErr);
            }

            // ---- value pad ----
            pV->cd();
            TGraphErrors* gr1 = new TGraphErrors(n, &ang[0], &r1v[0], 0, &r1e[0]);
            TGraphErrors* gr2 = new TGraphErrors(n, &ang[0], &r2v[0], 0, &r2e[0]);
            gr1->SetTitle(Form("%s, %s-axis;;Monitor-Normalized QE [a.u.]", serial.Data(), axes[a]));
            // Small markers so the +-1sigma bars stay legible (a large glyph hides a ~2% bar).
            // Open square vs filled circle: shape distinguishes the two repeats even in
            // greyscale, where blue-vs-red alone collapses (2026-08-25).
            gr1->SetMarkerStyle(25); gr1->SetMarkerColor(kBlue + 2); gr1->SetLineColor(kBlue + 2); gr1->SetMarkerSize(0.8);
            gr2->SetMarkerStyle(20); gr2->SetMarkerColor(kRed + 1);  gr2->SetLineColor(kRed + 1);  gr2->SetMarkerSize(0.7);
            gr1->Draw("APZ"); gr2->Draw("PZ same");
            gr1->GetYaxis()->SetTitleSize(0.052); gr1->GetYaxis()->SetLabelSize(0.045); gr1->GetYaxis()->SetTitleOffset(1.35);
            gr1->GetXaxis()->SetLabelSize(0);
            gr1->GetYaxis()->SetNdivisions(508);

            // Reserve empty band at the top for the legend and at the bottom
            // for the stats box, instead of letting either land on the curve.
            {
                double lo = 1e30, hi = -1e30;
                for (int i = 0; i < n; ++i) {
                    lo = std::min(lo, std::min(r1v[i] - r1e[i], r2v[i] - r2e[i]));
                    hi = std::max(hi, std::max(r1v[i] + r1e[i], r2v[i] + r2e[i]));
                }
                double span = (hi > lo) ? (hi - lo) : std::max(1e-6, std::abs(hi) * 0.1);
                gr1->SetMinimum(lo - span * 0.45);   // room for the stats box
                gr1->SetMaximum(hi + span * 0.35);   // room for the legend
            }

            TLegend* leg = new TLegend(0.66, 0.79, 0.95, 0.93);
            leg->SetBorderSize(1); leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
            leg->SetTextFont(132); leg->SetTextSize(0.042);
            leg->AddEntry(gr1, "Rep.1", "lp"); leg->AddEntry(gr2, "Rep.2", "lp");
            leg->Draw();

            TPaveText* box = new TPaveText(0.185, 0.045, 0.62, 0.33, "NDC");
            box->SetBorderSize(2); box->SetFillColor(kWhite); box->SetFillStyle(1001);
            box->SetTextFont(132); box->SetTextSize(0.039); box->SetTextAlign(12);
            box->AddText(Form("N = %d,  stat/pt = %.2f %%", st.n, (st.meanRelErr1 + st.meanRelErr2) / 2 * 100));
            box->AddText(Form("#sigma_{obs} = %.2f %%,  #sigma_{exp} = %.2f %%", st.obsSpread * 100, st.expSpread * 100));
            box->AddText(Form("RMS(pull) = %.2f #pm %.2f", st.pullRMS, st.pullRMS / std::sqrt(2.0 * st.n)));
            if (st.sysResolved)
                box->AddText(Form("#bf{syst. = %.2f %%}", st.sysComp * 100));
            else if (st.statOverestimated)
                box->AddText(Form("#bf{syst. < %.2f %%}  (stat. err overest.)", st.sysUL * 100));
            else
                box->AddText(Form("#bf{syst. < %.2f %% (95%% CL)}", st.sysUL * 100));
            box->Draw();

            // ---- ratio pad ----
            pR->cd();
            TGraphErrors* grR = new TGraphErrors(n, &ang[0], &ratioV[0], 0, &ratioE[0]);
            grR->SetTitle(";Position angle [degree];Rep.2 / Rep.1");
            grR->SetMarkerStyle(20); grR->SetMarkerColor(kBlack); grR->SetLineColor(kBlack); grR->SetMarkerSize(0.6);
            grR->Draw("APZ");
            grR->GetXaxis()->SetTitleSize(0.125); grR->GetXaxis()->SetLabelSize(0.105); grR->GetXaxis()->SetTitleOffset(1.15);
            grR->GetYaxis()->SetTitleSize(0.105); grR->GetYaxis()->SetLabelSize(0.095); grR->GetYaxis()->SetTitleOffset(0.62);
            grR->GetYaxis()->SetNdivisions(505);

            double xLo = ang.front() - 3, xHi = ang.back() + 3;
            TLine* l1 = new TLine(xLo, 1.0, xHi, 1.0);
            l1->SetLineColor(kRed + 1); l1->SetLineWidth(2); l1->Draw();
            // +-observed scatter, so points beyond the band are visible at a glance
            for (int sgn : {-1, 1}) {
                TLine* l = new TLine(xLo, 1.0 + sgn * st.obsSpread, xHi, 1.0 + sgn * st.obsSpread);
                l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
            }
        }

        if (drewAny) c->Print(outPath);
    }

    // ---- final page: pull distribution over every axis/PMT ----
    if (allPulls.size() >= 2) {
        c->Clear(); c->cd();
        gPad->SetGrid(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.13);

        TH1F* hPull = new TH1F("hPullAll",
                                "Pull distribution   (Rep.2 - Rep.1) / #sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}};Pull;Entries",
                                25, -6, 6);
        double sum2 = 0;
        for (double p : allPulls) { hPull->Fill(p); sum2 += p * p; }
        double rms = std::sqrt(sum2 / allPulls.size());

        hPull->SetLineColor(kBlue + 2); hPull->SetLineWidth(2); hPull->SetFillColorAlpha(kBlue + 2, 0.25);
        hPull->GetXaxis()->SetTitleSize(0.045); hPull->GetXaxis()->SetLabelSize(0.038);
        hPull->GetYaxis()->SetTitleSize(0.045); hPull->GetYaxis()->SetLabelSize(0.038);
        // Headroom so the stats box sits above the tallest bin rather than on it.
        hPull->SetMaximum(hPull->GetMaximum() * 1.45);
        hPull->Draw("hist");

        // Unit-Gaussian reference scaled to this histogram: "RMS(pull)=1.05"
        // is far easier to judge against a drawn curve than as a bare number.
        TF1* ref = new TF1("pullRef", "[0]*TMath::Gaus(x,0,1,kTRUE)", -6, 6);
        ref->SetParameter(0, allPulls.size() * hPull->GetBinWidth(1));
        ref->SetLineColor(kGray + 2); ref->SetLineStyle(2); ref->SetLineWidth(2);
        ref->Draw("same");

        TLegend* legP = new TLegend(0.16, 0.80, 0.44, 0.88);
        legP->SetBorderSize(1); legP->SetFillColor(kWhite); legP->SetFillStyle(1001);
        legP->SetTextFont(132); legP->SetTextSize(0.032);
        legP->AddEntry(hPull, "Pull", "f");
        legP->AddEntry(ref, "Gaussian", "l");
        legP->Draw();

        TPaveText* pb = new TPaveText(0.62, 0.70, 0.93, 0.87, "NDC");
        pb->SetBorderSize(2); pb->SetFillColor(kWhite); pb->SetFillStyle(1001);
        pb->SetTextFont(132); pb->SetTextSize(0.035); pb->SetTextAlign(12);
        pb->AddText(Form("N = %d", (int)allPulls.size()));
        pb->AddText(Form("RMS(pull) = %.2f #pm %.2f",
                         rms, rms / std::sqrt(2.0 * allPulls.size())));
        pb->Draw();
        c->Print(outPath);

        std::cout << Form("\n\033[1;33m[ All axes/PMTs combined ]  N = %d,  RMS(pull) = %.2f\033[0m",
                           (int)allPulls.size(), rms) << std::endl;
    }

    c->Print(outPath + "]");
    f1->Close(); f2->Close();
    std::cout << "\n[INFO] Saved: " << outPath << std::endl;
}
