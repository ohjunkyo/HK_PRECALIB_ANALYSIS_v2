// Draws Ratio_i vs Position angle for all 4 (PMT, axis) groups, with the
// linear-trend fit overlaid, plus a bootstrap histogram of sigma_observed.
// Standalone visualization companion to CheckTrend.C -- does not modify any
// existing macro.

#include <TFile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

struct MatchedPoint { double angle, v1, e1, v2, e2, ratio; };

static std::vector<MatchedPoint> LoadMatched(TFile* f1, TFile* f2, TString serial,
                                             TString axis, TString metric, double angleTol) {
    std::vector<MatchedPoint> pts;
    TGraphErrors* g1 = (TGraphErrors*)f1->Get(Form("gr_%s_%s_%s", serial.Data(), axis.Data(), metric.Data()));
    TGraphErrors* g2 = (TGraphErrors*)f2->Get(Form("gr_%s_%s_%s", serial.Data(), axis.Data(), metric.Data()));
    if (!g1 || !g2) return pts;
    for (int i = 0; i < g1->GetN(); ++i) {
        double an1 = g1->GetX()[i], v1 = g1->GetY()[i], e1 = g1->GetEY()[i];
        for (int j = 0; j < g2->GetN(); ++j) {
            double an2 = g2->GetX()[j], v2 = g2->GetY()[j], e2 = g2->GetEY()[j];
            if (std::abs(an1 - an2) > angleTol) continue;
            MatchedPoint mp; mp.angle = an1; mp.v1 = v1; mp.e1 = e1; mp.v2 = v2; mp.e2 = e2;
            mp.ratio = v2 / v1;
            pts.push_back(mp);
            break;
        }
    }
    std::sort(pts.begin(), pts.end(), [](const MatchedPoint& x, const MatchedPoint& y){ return x.angle < y.angle; });
    return pts;
}

void PlotTrendCheck(TString tag1 = "20260803_100_145", TString tag2 = "20260803_200_245",
                     TString metric = "MonNorm_CorrRawQE", double angleTol = 1.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile* f1 = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag1.Data()));
    TFile* f2 = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag2.Data()));
    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
        std::cout << "[ERROR] cannot open one of the tag files." << std::endl; return;
    }

    TString sns[2] = {"EM5370", "EL9590"};
    TString axes[2] = {"X", "Y"};

    // ---- Page 1: Ratio vs angle, 2x2 grid, with linear fit ----
    TCanvas* c1 = new TCanvas("c_trend", "Ratio vs Angle", 1400, 1000);
    c1->Divide(2, 2, 0.01, 0.01);

    int pad = 1;
    for (int a = 0; a < 2; ++a) {
        for (int s = 0; s < 2; ++s) {
            c1->cd(pad++);
            gPad->SetGrid();
            auto pts = LoadMatched(f1, f2, sns[s], axes[a], metric, angleTol);
            int n = (int)pts.size();
            if (n < 2) continue;

            std::vector<double> ang(n), rat(n), zero(n, 0);
            for (int i = 0; i < n; ++i) { ang[i] = pts[i].angle; rat[i] = (pts[i].ratio - 1.0) * 100.0; }

            TGraph* gr = new TGraph(n, &ang[0], &rat[0]);
            gr->SetTitle(Form("%s, %s-axis;Position angle [deg];(Ratio - 1) [%%]", sns[s].Data(), axes[a].Data()));
            gr->SetMarkerStyle(20);
            gr->SetMarkerColor(a == 0 ? kBlue + 1 : kRed + 1);
            gr->SetMarkerSize(1.1);
            gr->Draw("AP");

            TF1* lin = new TF1(Form("lin_%d", pad), "pol1", ang.front(), ang.back());
            gr->Fit(lin, "Q");
            double slope = lin->GetParameter(1), slopeErr = lin->GetParError(1);
            double sig = (slopeErr > 0) ? slope / slopeErr : 0;
            lin->SetLineColor(kBlack);
            lin->SetLineStyle(2);
            lin->Draw("SAME");

            TLine* zeroLine = new TLine(ang.front(), 0, ang.back(), 0);
            zeroLine->SetLineColor(kGray + 1);
            zeroLine->SetLineStyle(3);
            zeroLine->Draw("SAME");

            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.045);
            lat.DrawLatex(0.15, 0.87, Form("slope significance: %.2f#sigma", sig));
            lat.DrawLatex(0.15, 0.81, sig < 2.0 && sig > -2.0 ? "No significant trend" : "Trend present");
        }
    }
    c1->SaveAs("./Data/UNIFORMITY/TrendCheck_RatioVsAngle.png");

    // ---- Page 2: Bootstrap distribution of sigma_observed, all 4 groups ----
    TCanvas* c2 = new TCanvas("c_boot", "Bootstrap sigma_obs", 1400, 1000);
    c2->Divide(2, 2, 0.01, 0.01);

    pad = 1;
    TRandom3 rng(42);
    for (int a = 0; a < 2; ++a) {
        for (int s = 0; s < 2; ++s) {
            c2->cd(pad++);
            auto pts = LoadMatched(f1, f2, sns[s], axes[a], metric, angleTol);
            int n = (int)pts.size();
            if (n < 2) continue;
            std::vector<double> rat(n);
            for (int i = 0; i < n; ++i) rat[i] = pts[i].ratio;

            double meanRelErr = 0;
            for (auto& p : pts) meanRelErr += 0.5 * (p.e1 / p.v1 + p.e2 / p.v2);
            meanRelErr /= n;
            double sigma_exp = std::sqrt(2) * meanRelErr;

            int nboot = 5000;
            TH1F* h = new TH1F(Form("hboot_%d", pad), Form("%s, %s-axis;#sigma_{observed} [%%];Bootstrap samples", sns[s].Data(), axes[a].Data()), 60, 1.5, 6.0);
            double obsCentral;
            {
                double mean = 0; for (auto& r : rat) mean += r; mean /= n;
                double obs2 = 0; for (auto& r : rat) obs2 += (r - mean) * (r - mean);
                obsCentral = std::sqrt(obs2 / (n - 1));
            }
            for (int b = 0; b < nboot; ++b) {
                std::vector<double> sample(n);
                for (int i = 0; i < n; ++i) sample[i] = rat[rng.Integer(n)];
                double m = 0; for (auto& r : sample) m += r; m /= n;
                double s2 = 0; for (auto& r : sample) s2 += (r - m) * (r - m);
                h->Fill(std::sqrt(s2 / (n - 1)) * 100.0);
            }
            h->SetFillColorAlpha(a == 0 ? kBlue - 9 : kRed - 9, 0.7);
            h->SetLineColor(kBlack);
            h->Draw("HIST");

            TLine* expLine = new TLine(sigma_exp * 100, 0, sigma_exp * 100, h->GetMaximum());
            expLine->SetLineColor(kGreen + 2);
            expLine->SetLineWidth(2);
            expLine->Draw("SAME");
            TLine* obsLine = new TLine(obsCentral * 100, 0, obsCentral * 100, h->GetMaximum());
            obsLine->SetLineColor(kBlack);
            obsLine->SetLineWidth(2);
            obsLine->SetLineStyle(2);
            obsLine->Draw("SAME");

            TLegend* leg = new TLegend(0.55, 0.7, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.035);
            leg->AddEntry(expLine, Form("Expected = %.2f%%", sigma_exp * 100), "l");
            leg->AddEntry(obsLine, Form("Observed = %.2f%%", obsCentral * 100), "l");
            leg->Draw();
        }
    }
    c2->SaveAs("./Data/UNIFORMITY/TrendCheck_Bootstrap.png");

    std::cout << "\n[INFO] Saved:" << std::endl;
    std::cout << "  ./Data/UNIFORMITY/TrendCheck_RatioVsAngle.png" << std::endl;
    std::cout << "  ./Data/UNIFORMITY/TrendCheck_Bootstrap.png" << std::endl;
}
