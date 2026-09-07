// Draws the Pull distribution -- pull_i = (ratio_i - mean) / sigma_expected,i
// -- as a histogram for each (PMT, axis) group, overlaid with the standard
// normal N(0,1) that pure counting statistics would predict. RMS(pull) is
// printed on each pad; RMS > 1 is the same systematic signal already reported
// as sigma_syst in Draw_Reproducibility.C, shown here as "how many sigma did
// this point move" instead of "% excess variance".
//
// Standalone, does not modify Draw_Reproducibility.C.

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

struct MatchedPoint { double angle, v1, e1, v2, e2, ratio, ratioErr, pull; };

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
            // Relative error propagation for a ratio of two independent values.
            mp.ratioErr = mp.ratio * std::sqrt((e1/v1)*(e1/v1) + (e2/v2)*(e2/v2));
            pts.push_back(mp);
            break;
        }
    }
    std::sort(pts.begin(), pts.end(), [](const MatchedPoint& x, const MatchedPoint& y){ return x.angle < y.angle; });
    return pts;
}

void PlotPullDist(TString tag1 = "20260803_100_145", TString tag2 = "20260803_200_245",
                   TString metric = "MonNorm_CorrRawQE", double angleTol = 1.0) {
    gStyle->SetOptStat(0);
    // Times New Roman everywhere (ROOT font code 132 = Times, bold, precision 3
    // -- scales with pad size rather than staying a fixed pixel size).
    gStyle->SetTextFont(132);
    gStyle->SetLabelFont(132, "XYZ");
    gStyle->SetTitleFont(132, "XYZ");
    gStyle->SetLegendFont(132);
    gStyle->SetTitleFont(132, "");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetFrameLineWidth(2);

    TFile* f1 = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag1.Data()));
    TFile* f2 = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag2.Data()));
    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
        std::cout << "[ERROR] cannot open one of the tag files." << std::endl; return;
    }

    TString sns[2] = {"EM5370", "EL9590"};
    TString axes[2] = {"X", "Y"};

    TCanvas* c = new TCanvas("cPull", "Pull distribution", 1400, 1000);
    c->Divide(2, 2, 0.01, 0.01);

    int pad = 1;
    std::vector<double> allPulls;
    for (int a = 0; a < 2; ++a) {
        for (int s = 0; s < 2; ++s) {
            c->cd(pad++);
            gPad->SetGrid();
            auto pts = LoadMatched(f1, f2, sns[s], axes[a], metric, angleTol);
            int n = (int)pts.size();
            if (n < 2) continue;

            double mean = 0; for (auto& p : pts) mean += p.ratio; mean /= n;
            std::vector<double> pulls(n);
            for (int i = 0; i < n; ++i) {
                pulls[i] = (pts[i].ratio - mean) / pts[i].ratioErr;
                allPulls.push_back(pulls[i]);
            }
            double rms = 0; for (auto& v : pulls) rms += v * v; rms = std::sqrt(rms / n);

            // Bin width fixed at 0.5 (pull units -- dimensionless, no "%":
            // numerator and denominator are both in the same units, e.g. %,
            // so they cancel), range -4 to 4 -> 16 bins. 0.25 (32 bins) left
            // visible gaps with only ~23 points per group.
            TH1F* h = new TH1F(Form("hpull_%d", pad), Form("%s, %s-axis;Pull;Entries / 0.5", sns[s].Data(), axes[a].Data()), 16, -4, 4);
            for (auto& v : pulls) h->Fill(v);
            h->SetFillStyle(0);          // unfilled -- outline only
            h->SetLineColor(kBlack);
            h->SetLineWidth(2);
            // Headroom so the stat box (top-left) never overlaps the peak
            // bin, whatever its height turns out to be.
            h->SetMaximum(h->GetMaximum() * 1.45);
            h->Draw("HIST");

            // Standard normal N(0,1), scaled to the histogram's area -- what
            // pure counting statistics (no systematic) predicts.
            TF1* gaus = new TF1(Form("gaus_%d", pad), "[0]*TMath::Gaus(x,0,1)", -4, 4);
            double binw = h->GetBinWidth(1);
            gaus->SetParameter(0, n * binw / std::sqrt(2 * TMath::Pi()));
            gaus->SetLineColor(kRed);
            gaus->SetLineWidth(2);
            gaus->SetLineStyle(1);
            gaus->Draw("SAME");

            // Stat-box style: solid white background + visible border, not a
            // legend with a symbol column -- same corner/size on every pad so
            // the eye doesn't have to re-find it per plot. Top-left is clear
            // of the histogram peak (centered near Pull=0) on every group.
            TLegend* leg = new TLegend(0.16, 0.81, 0.40, 0.88);
            leg->SetBorderSize(2);
            leg->SetFillColor(kWhite);
            leg->SetFillStyle(1001);
            leg->SetTextFont(132);
            leg->SetTextSize(0.04);
            leg->SetMargin(0.08);
            leg->AddEntry((TObject*)0, Form("RMS Pull: %.2f", rms), "");
            leg->Draw();
        }
    }
    c->SaveAs("./Data/UNIFORMITY/PullDistribution_byGroup.png");

    // ---- Combined 92-point pull distribution on its own canvas ----
    TCanvas* c2 = new TCanvas("cPullAll", "Pull distribution, all groups", 900, 700);
    c2->SetGrid();
    int nAll = (int)allPulls.size();
    double rmsAll = 0; for (auto& v : allPulls) rmsAll += v * v; rmsAll = std::sqrt(rmsAll / nAll);

    TH1F* hAll = new TH1F("hPullAll", Form("All %d points (X+Y, both PMTs);Pull;Entries / 0.5", nAll), 16, -4, 4);
    for (auto& v : allPulls) hAll->Fill(v);
    hAll->SetFillStyle(0);
    hAll->SetLineColor(kBlack);
    hAll->SetLineWidth(2);
    hAll->SetMaximum(hAll->GetMaximum() * 1.45);
    hAll->Draw("HIST");

    TF1* gausAll = new TF1("gausAll", "[0]*TMath::Gaus(x,0,1)", -4, 4);
    double binwAll = hAll->GetBinWidth(1);
    gausAll->SetParameter(0, nAll * binwAll / std::sqrt(2 * TMath::Pi()));
    gausAll->SetLineColor(kRed);
    gausAll->SetLineWidth(2);
    gausAll->SetLineStyle(1);
    gausAll->Draw("SAME");

    // Same corner and box size as the 4-pad version above, kept in the same
    // NDC coordinates so both plots read identically.
    TLegend* legAll = new TLegend(0.16, 0.81, 0.40, 0.88);
    legAll->SetBorderSize(2);
    legAll->SetFillColor(kWhite);
    legAll->SetFillStyle(1001);
    legAll->SetTextFont(132);
    legAll->SetTextSize(0.04);
    legAll->SetMargin(0.08);
    legAll->AddEntry((TObject*)0, Form("RMS Pull: %.2f", rmsAll), "");
    legAll->Draw();

    c2->SaveAs("./Data/UNIFORMITY/PullDistribution_combined.png");

    std::cout << "\n[INFO] Saved:" << std::endl;
    std::cout << "  ./Data/UNIFORMITY/PullDistribution_byGroup.png" << std::endl;
    std::cout << "  ./Data/UNIFORMITY/PullDistribution_combined.png" << std::endl;
}
