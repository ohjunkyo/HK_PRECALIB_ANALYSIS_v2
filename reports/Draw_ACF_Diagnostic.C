// Draw_ACF_Diagnostic.C
// ----------------------------------------------------------------------------
// Autocorrelation diagnostic for a fixed-angle Stability run's QE time series:
// per channel (Mon./Rot#1/Rot#2), 4 panels --
//   1) raw QE(t) with its slow trend (wide moving average) overlaid
//   2) residual after subtracting that trend
//   3) ACF(lag) computed on the RAW (mean-only-removed) series
//   4) ACF(lag) computed on the DETRENDED residual
// Panels 3/4 both show the +-95% white-noise bound (Bartlett, 1.96/sqrt(N))
// and a marker at the candidate period, so the printed ACF numbers from the
// Python analysis can be read directly off the plot -- same numbers, ROOT
// rendering.
//
// Input CSV columns: run,epoch,qmon,q1,q2  (as produced for the Stability ACF
// work -- one row per point, epoch = unix seconds of acquisition).
//
// Usage:
//   root -l -b -q 'Draw_ACF_Diagnostic.C("set1_qe.csv", "Set1_20260813", 40.2)'
// ----------------------------------------------------------------------------
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

namespace {

struct Row { double epoch, qmon, q1, q2; };

std::vector<Row> LoadCSV(const TString& path) {
    std::vector<Row> rows;
    std::ifstream in(path.Data());
    if (!in.is_open()) { std::cout << "[ERROR] cannot open " << path << std::endl; return rows; }
    std::string line;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string f[5];
        for (int i = 0; i < 5 && std::getline(ss, f[i], ','); ++i) {}
        try {
            rows.push_back({std::stod(f[1]), std::stod(f[2]), std::stod(f[3]), std::stod(f[4])});
        } catch (...) { continue; }
    }
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b){ return a.epoch < b.epoch; });
    return rows;
}

// Resample an irregular (t,y) series onto a uniform grid via linear
// interpolation -- ACF assumes even spacing, but real points drift a little
// (watchdog waits etc.), same reasoning as Draw_Stability_v1.C.
std::vector<double> ResampleUniform(const std::vector<double>& t, const std::vector<double>& y,
                                    const std::vector<double>& grid) {
    TGraph g(t.size(), t.data(), y.data());
    std::vector<double> out(grid.size());
    for (size_t i = 0; i < grid.size(); ++i) out[i] = g.Eval(grid[i]);
    return out;
}

// Centered moving-average trend over `windowSteps` grid points -- much wider
// than the candidate period, so it tracks the slow drift only.
std::vector<double> MovingAvgTrend(const std::vector<double>& y, int windowSteps) {
    int n = (int)y.size(), half = windowSteps / 2;
    std::vector<double> trend(n);
    for (int i = 0; i < n; ++i) {
        int lo = std::max(0, i - half), hi = std::min(n, i + half + 1);
        double s = 0; for (int k = lo; k < hi; ++k) s += y[k];
        trend[i] = s / (hi - lo);
    }
    return trend;
}

// ACF(lag) on an already-mean-removed series, normalized so ACF(0)=1.
void Acf(const std::vector<double>& ym, double dt, int maxLagSteps,
        std::vector<double>& lags, std::vector<double>& rk) {
    int n = (int)ym.size();
    double denom = 0; for (double v : ym) denom += v * v;
    lags.push_back(0); rk.push_back(1.0);
    for (int k = 1; k <= maxLagSteps && k < n; ++k) {
        double num = 0; for (int i = 0; i < n - k; ++i) num += ym[i] * ym[i + k];
        lags.push_back(k * dt); rk.push_back(num / denom);
    }
}

}  // namespace

void Draw_ACF_Diagnostic(TString csvPath = "set1_qe.csv", TString tag = "Set1",
                         double candidatePeriodMin = 40.2) {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132);
    gStyle->SetLabelFont(132, "XYZ");
    gStyle->SetTitleFont(132, "XYZ");
    gStyle->SetFrameLineWidth(2);

    std::vector<Row> rows = LoadCSV(csvPath);
    std::cout << "[INFO] " << rows.size() << " points loaded from " << csvPath << std::endl;
    if (rows.size() < 10) { std::cout << "[ERROR] not enough points." << std::endl; return; }

    std::vector<double> t, qmon, q1, q2;
    double t0 = rows[0].epoch;
    for (auto& r : rows) {
        t.push_back((r.epoch - t0) / 60.0);
        qmon.push_back(r.qmon); q1.push_back(r.q1); q2.push_back(r.q2);
    }

    // Median spacing -> uniform resample grid (same logic as prep in
    // Draw_Stability_v1.C's own time handling).
    std::vector<double> dts;
    for (size_t i = 1; i < t.size(); ++i) dts.push_back(t[i] - t[i-1]);
    std::sort(dts.begin(), dts.end());
    double dtMed = dts[dts.size()/2];
    std::vector<double> grid;
    for (double g = 0; g < t.back(); g += dtMed) grid.push_back(g);
    int N = (int)grid.size();
    double bartlett = 1.96 / std::sqrt((double)N);
    int windowSteps = std::max(3, (int)std::round(180.0 / dtMed));   // 3h trend window
    int maxLagSteps = std::min(20, N / 3);

    std::cout << Form("[INFO] N=%d grid pts, dt=%.2fmin, Bartlett bound=+-%.3f, trend window=%.0fmin",
                      N, dtMed, bartlett, windowSteps * dtMed) << std::endl;

    struct Ch { const char* name; int color; std::vector<double>* y; };
    Ch chans[3] = {{"Mon.", kBlack, &qmon}, {"Rot#1", kAzure+2, &q1}, {"Rot#2", kRed+1, &q2}};

    // One page per analysis step, three PMTs stacked on each -- the earlier
    // 4x3 grid squeezed twelve pads onto a single canvas, which made the
    // trend/ACF curves too small to read.
    struct Series { std::vector<double> yi, trend, resid, lagsRaw, rkRaw, lagsDt, rkDt; };
    Series S[3];
    for (int row = 0; row < 3; ++row) {
        S[row].yi = ResampleUniform(t, *chans[row].y, grid);
        double mean = 0; for (double v : S[row].yi) mean += v; mean /= N;
        std::vector<double> ymRaw(N); for (int i = 0; i < N; ++i) ymRaw[i] = S[row].yi[i] - mean;
        S[row].trend = MovingAvgTrend(S[row].yi, windowSteps);
        S[row].resid.resize(N);
        for (int i = 0; i < N; ++i) S[row].resid[i] = S[row].yi[i] - S[row].trend[i];
        Acf(ymRaw,        dtMed, maxLagSteps, S[row].lagsRaw, S[row].rkRaw);
        Acf(S[row].resid, dtMed, maxLagSteps, S[row].lagsDt,  S[row].rkDt);
    }

    gSystem->mkdir("./Data/image/Stability", kTRUE);
    TString pdf = Form("./Data/image/Stability/ACF_Steps_%s.pdf", tag.Data());
    TCanvas* c = new TCanvas("cACF", "ACF steps", 950, 1300);
    c->Print(pdf + "[");

    auto stylePad = [&]() {
        gPad->SetGrid();
        gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.09);  gPad->SetBottomMargin(0.14);
    };
    auto styleAxes = [&](TGraph* g) {
        g->GetXaxis()->SetLabelSize(0.045); g->GetXaxis()->SetTitleSize(0.050);
        g->GetYaxis()->SetLabelSize(0.045); g->GetYaxis()->SetTitleSize(0.050);
        g->GetYaxis()->SetTitleOffset(0.95);
    };

    // ---- Page 1: raw QE(t) with the moving-average trend on top -------------
    c->Clear(); c->Divide(1, 3, 0.002, 0.015);
    for (int row = 0; row < 3; ++row) {
        c->cd(row + 1); stylePad();
        TGraph* gRaw = new TGraph(N, grid.data(), S[row].yi.data());
        gRaw->SetMarkerStyle(20); gRaw->SetMarkerSize(0.6); gRaw->SetMarkerColor(chans[row].color);
        gRaw->SetLineColor(chans[row].color); gRaw->SetLineWidth(1);
        gRaw->SetTitle(Form("%s: Raw QE + trend;Elapsed [min];QE [%%]", chans[row].name));
        styleAxes(gRaw); gRaw->Draw("APL");
        TGraph* gTrend = new TGraph(N, grid.data(), S[row].trend.data());
        gTrend->SetLineColor(kRed+1); gTrend->SetLineWidth(3);
        gTrend->Draw("L same");
    }
    c->Print(pdf);

    // ---- Page 2: what is left once that trend is subtracted ----------------
    c->Clear(); c->Divide(1, 3, 0.002, 0.015);
    for (int row = 0; row < 3; ++row) {
        c->cd(row + 1); stylePad();
        TGraph* gRes = new TGraph(N, grid.data(), S[row].resid.data());
        gRes->SetMarkerStyle(20); gRes->SetMarkerSize(0.6); gRes->SetMarkerColor(chans[row].color);
        gRes->SetLineColor(chans[row].color); gRes->SetLineWidth(1);
        gRes->SetTitle(Form("%s: Residual after detrend;Elapsed [min];QE - trend [%%]", chans[row].name));
        styleAxes(gRes); gRes->Draw("APL");
        TLine* zero = new TLine(grid.front(), 0, grid.back(), 0);
        zero->SetLineColor(kGray+1); zero->SetLineStyle(2); zero->Draw();
    }
    c->Print(pdf);

    // ---- Pages 3/4: ACF before and after the detrend ------------------------
    // Same axes on both so the two pages can be compared directly: if the
    // candidate-period peak survives page 4, it is a real oscillation rather
    // than the slow drift leaking into the autocorrelation.
    for (int pass = 0; pass < 2; ++pass) {
        c->Clear(); c->Divide(1, 3, 0.002, 0.015);
        for (int row = 0; row < 3; ++row) {
            c->cd(row + 1); stylePad();
            std::vector<double>& lags = pass == 0 ? S[row].lagsRaw : S[row].lagsDt;
            std::vector<double>& rk   = pass == 0 ? S[row].rkRaw   : S[row].rkDt;
            int col = pass == 0 ? kGray + 2 : chans[row].color;
            TGraph* gA = new TGraph(lags.size(), lags.data(), rk.data());
            gA->SetMarkerStyle(20); gA->SetMarkerSize(0.8); gA->SetMarkerColor(col);
            gA->SetLineColor(col); gA->SetLineWidth(2);
            gA->SetTitle(Form("%s: ACF %s detrend;Lag [min];ACF",
                              chans[row].name, pass == 0 ? "BEFORE" : "AFTER"));
            gA->SetMinimum(-0.5); gA->SetMaximum(1.05);
            styleAxes(gA); gA->Draw("APL");
            for (double b : {bartlett, -bartlett}) {
                TLine* l = new TLine(0, b, lags.back(), b);
                l->SetLineColor(kRed+1); l->SetLineStyle(2); l->Draw();
            }
            TLine* mark = new TLine(candidatePeriodMin, -0.5, candidatePeriodMin, 1.05);
            mark->SetLineColor(kGreen+2); mark->SetLineStyle(3); mark->SetLineWidth(2); mark->Draw();
        }
        c->Print(pdf);
    }

    c->Print(pdf + "]");
    std::cout << "[INFO] Saved: " << pdf << std::endl;
}
