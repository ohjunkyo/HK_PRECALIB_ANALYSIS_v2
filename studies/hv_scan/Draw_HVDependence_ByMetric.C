// Usage:  root -l -b -q Draw_HVDependence_ByMetric.C
//
// How much each per-angle metric shifts with test-PMT HV, measured against the
// nominal (Recom.) point. Reads the same Graphs_Uniformity_<tag>.root files the
// Overlay report draws from, so the numbers match that report exactly.
//
// Two views per metric:
//   * page plot : mean deviation from Recom. vs HV, |angle| < kAngleCut,
//                 Rot1/Rot2 x X/Y axis
//   * stdout    : the same deviation split into |angle| bands, which is where
//                 the low-HV points fall apart (>70 deg roughly doubles the
//                 deviation -- the low-light region where P/V collapses).
//
// Block/label mapping is verified against the real RunInfo HV2/HV3 -- note
// 20260829_100_145 is -200V and 20260829_0_45 is -100V (they read the other way
// round in a couple of older scripts; see project memory
// "systematic-error-summary-20260901", item 9).

#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include <vector>
#include <cmath>
#include <cstdio>

namespace {

struct HVPoint { const char* tag; const char* label; double dHV; };

// Ordered low -> high HV. The Recom. block is the reference.
const std::vector<HVPoint> kHV = {
    {"20260829_100_145", "-200V", -200},
    {"20260829_200_245", "-150V", -150},
    {"20260829_0_45",    "-100V", -100},
    {"20260827_0_45",    "-50V",   -50},
    {"20260815_0_45",    "Recom.",   0},
    {"20260827_100_145", "+50V",    50},
    {"20260830_0_45",    "+100V",  100},
    {"20260828_0_45",    "+150V",  150},
    {"20260828_100_145", "+200V",  200},
};
const char* kRefTag = "20260815_0_45";

struct Metric { const char* graph; const char* title; };
const std::vector<Metric> kMetrics = {
    {"MonNorm_CorrRawQE", "Relative QE (Monitor-normalised)"},
    {"TTS",               "TTS"},
    {"RawSPE",            "Gain (SPE mean charge)"},
    {"ChargeResolution",  "Charge Resolution"},
};

const double kAngleCut = 70.0;   // deviations beyond this blow up (low-light fits)

// spe_mean is stored in pC; gains are quoted in units of 1e7 electrons.
const double kPcToGain1e7 = 1e-12 / 1.602176634e-19 / 1e7;

struct Series { std::vector<double> x, y; };

Series ReadGraph(const char* tag, const char* pmt, const char* axis, const char* metric) {
    Series s;
    TString path = Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag);
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { if (f) f->Close(); return s; }
    TGraph* g = (TGraph*)f->Get(Form("gr_%s_%s_%s", pmt, axis, metric));
    if (g) for (int i = 0; i < g->GetN(); ++i) { s.x.push_back(g->GetX()[i]); s.y.push_back(g->GetY()[i]); }
    f->Close();
    return s;
}

// Mean of (value/reference - 1) over the angles both series share, restricted
// to |angle| in [lo, hi). Angles are matched by nearest x within 1 deg because
// the two runs' stage angles agree to well under that.
double MeanDev(const Series& s, const Series& ref, double lo, double hi, int* nOut = nullptr) {
    double sum = 0; int n = 0;
    for (size_t i = 0; i < s.x.size(); ++i) {
        double a = std::fabs(s.x[i]);
        if (a < lo || a >= hi) continue;
        int best = -1; double bestD = 1e9;
        for (size_t j = 0; j < ref.x.size(); ++j) {
            double d = std::fabs(ref.x[j] - s.x[i]);
            if (d < bestD) { bestD = d; best = (int)j; }
        }
        if (best < 0 || bestD > 1.0 || ref.y[best] == 0) continue;
        sum += (s.y[i] / ref.y[best] - 1.0) * 100.0;
        n++;
    }
    if (nOut) *nOut = n;
    return n ? sum / n : std::nan("");
}

// Mean gain (in 1e7 electrons) over |angle| < kAngleCut, both scan axes pooled.
// Used to label each HV point with the gain it actually corresponds to -- the
// two PMTs differ, so this goes on as two text rows rather than a second axis
// (2026-09-02, user: "x축이 Gain도 같이 나타낼 수 있을까?").
double MeanGain1e7(const char* tag, const char* pmt) {
    double sum = 0; int n = 0;
    for (const char* ax : {"X", "Y"}) {
        Series s = ReadGraph(tag, pmt, ax, "RawSPE");
        for (size_t i = 0; i < s.x.size(); ++i)
            if (std::fabs(s.x[i]) < kAngleCut && s.y[i] > 0) { sum += s.y[i]; n++; }
    }
    return n ? sum / n * kPcToGain1e7 : std::nan("");
}

} // namespace

void Draw_HVDependence_ByMetric() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(1);

    const char* pmts[2]  = {"EM6400", "EL5150"};
    const char* pmtLab[2] = {"Rot1 (EM6400)", "Rot2 (EL5150)"};
    const char* axes[2]  = {"X", "Y"};
    int          cols[2]  = {kAzure + 2, kRed + 1};
    int          mkFill[2] = {20, 21};
    int          mkOpen[2] = {24, 25};

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + "HVDependence_ByMetric.pdf";

    TCanvas* c = new TCanvas("cHVDep", "HV dependence by metric", 900, 800);
    bool firstPage = true;

    for (const auto& m : kMetrics) {
        c->Clear();
        c->SetGridx(); c->SetGridy();
        c->SetLeftMargin(0.11); c->SetBottomMargin(0.13);

        TMultiGraph* mg = new TMultiGraph();
        // Legend sits in a reserved strip across the top rather than in a
        // corner: the curves run bottom-left -> top-right for QE/Gain but the
        // other way for TTS, so no fixed corner is free for every metric and a
        // corner legend ended up covering points (2026-09-02, user: "모든
        // 그림이 Legend에 의해서 Data가 가려져"). Headroom is added below.
        TLegend* leg = new TLegend(0.12, 0.85, 0.92, 0.955);
        leg->SetNColumns(2);
        leg->SetTextFont(132); leg->SetTextSize(0.030);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        double yLo = 1e30, yHi = -1e30;

        printf("\n=== %s : deviation from Recom. [%%] by |angle| band ===\n", m.title);

        bool anyPoint = false;
        for (int p = 0; p < 2; ++p) {
            for (int a = 0; a < 2; ++a) {
                Series ref = ReadGraph(kRefTag, pmts[p], axes[a], m.graph);
                if (ref.x.empty()) continue;

                std::vector<double> hv, dev;
                printf("%-6s %s-axis | %8s %8s %8s %8s %8s\n",
                       pmts[p], axes[a], "|a|<40", "40-60", "60-70", "70-80", ">80");
                for (const auto& h : kHV) {
                    Series s = ReadGraph(h.tag, pmts[p], axes[a], m.graph);
                    if (s.x.empty()) continue;
                    hv.push_back(h.dHV);
                    dev.push_back(MeanDev(s, ref, 0.0, kAngleCut));
                    printf("  %-6s      | %+8.1f %+8.1f %+8.1f %+8.1f %+8.1f\n", h.label,
                           MeanDev(s, ref,  0, 40), MeanDev(s, ref, 40, 60),
                           MeanDev(s, ref, 60, 70), MeanDev(s, ref, 70, 80),
                           MeanDev(s, ref, 80, 95));
                }
                if (hv.empty()) continue;
                anyPoint = true;
                for (double v : dev) if (std::isfinite(v)) { yLo = std::min(yLo, v); yHi = std::max(yHi, v); }

                TGraph* g = new TGraph((int)hv.size(), hv.data(), dev.data());
                g->SetLineColor(cols[p]); g->SetMarkerColor(cols[p]);
                g->SetMarkerStyle(a == 0 ? mkFill[p] : mkOpen[p]);
                g->SetLineStyle(a == 0 ? 1 : 2);
                g->SetMarkerSize(1.3); g->SetLineWidth(2);
                mg->Add(g, "LP");
                leg->AddEntry(g, Form("%s  %s-axis", pmtLab[p], axes[a]), "lp");
            }
        }
        if (!anyPoint) continue;

        mg->SetTitle(Form("%s vs HV  --  deviation from Recom., |angle| < %.0f#circ;"
                          "#DeltaHV from Recom. [V];Deviation from Recom. [%%]",
                          m.title, kAngleCut));
        mg->Draw("A");
        mg->GetXaxis()->SetLimits(-230, 230);
        mg->GetXaxis()->SetTitleFont(132); mg->GetYaxis()->SetTitleFont(132);
        mg->GetXaxis()->SetLabelFont(132); mg->GetYaxis()->SetLabelFont(132);
        // Reserve headroom above the data for the legend strip + the two gain
        // annotation rows.
        double axLo = yLo, axHi = yHi;
        if (yHi > yLo) {
            double span = yHi - yLo;
            axLo = yLo - 0.10 * span;
            axHi = yHi + 0.95 * span;
            mg->GetHistogram()->GetYaxis()->SetRangeUser(axLo, axHi);
        }

        TLine* z = new TLine(-230, 0, 230, 0);
        z->SetLineStyle(2); z->SetLineWidth(2); z->Draw();
        leg->Draw();

        // Gain that each HV point actually corresponds to, one row per PMT.
        if (axHi > axLo) {
            double rowR1 = axLo + 0.700 * (axHi - axLo);
            double rowR2 = axLo + 0.625 * (axHi - axLo);
            TLatex tg; tg.SetTextFont(132); tg.SetTextSize(0.026); tg.SetTextAlign(22);
            TLatex tl; tl.SetTextFont(132); tl.SetTextSize(0.026); tl.SetTextAlign(32);
            tl.SetTextColor(cols[0]); tl.DrawLatex(-165, rowR1, "Rot1 gain [#times10^{7}]:");
            tl.SetTextColor(cols[1]); tl.DrawLatex(-165, rowR2, "Rot2 gain [#times10^{7}]:");
            for (const auto& h : kHV) {
                if (h.dHV < -100) continue;   // left of the label text
                double g1 = MeanGain1e7(h.tag, pmts[0]);
                double g2 = MeanGain1e7(h.tag, pmts[1]);
                tg.SetTextColor(cols[0]);
                if (std::isfinite(g1)) tg.DrawLatex(h.dHV, rowR1, Form("%.2f", g1));
                tg.SetTextColor(cols[1]);
                if (std::isfinite(g2)) tg.DrawLatex(h.dHV, rowR2, Form("%.2f", g2));
            }
        }

        if (firstPage) { c->Print(outPath + "("); firstPage = false; }
        else           { c->Print(outPath); }
    }
    if (!firstPage) c->Print(outPath + ")");
    printf("\n[INFO] Saved: %s\n", outPath.Data());
}
