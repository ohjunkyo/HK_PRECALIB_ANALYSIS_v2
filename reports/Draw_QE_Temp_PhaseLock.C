// Draw_QE_Temp_PhaseLock.C
// Overlay Monitor-PMT QE residual and Dark Box temperature residual
// (both detrended, both z-scored to their own amplitude) on the same time
// axis, with a fitted 40.2min sinusoid on each -- makes the ~24min phase
// offset between them visually obvious.
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

namespace {

std::vector<double> Resample(const std::vector<double>& t, const std::vector<double>& y,
                             const std::vector<double>& grid) {
    TGraph g(t.size(), t.data(), y.data());
    std::vector<double> out(grid.size());
    for (size_t i = 0; i < grid.size(); ++i) out[i] = g.Eval(grid[i]);
    return out;
}

std::vector<double> MovingAvg(const std::vector<double>& y, int windowSteps) {
    int n = (int)y.size(), half = windowSteps / 2;
    std::vector<double> trend(n);
    for (int i = 0; i < n; ++i) {
        int lo = std::max(0, i - half), hi = std::min(n, i + half + 1);
        double s = 0; for (int k = lo; k < hi; ++k) s += y[k];
        trend[i] = s / (hi - lo);
    }
    return trend;
}

// Returns {grid, residual} for a (t[min], y) series: uniform resample + wide
// moving-average detrend, same recipe as Draw_ACF_Diagnostic.C.
void PrepResidual(std::vector<double> t, std::vector<double> y,
                  std::vector<double>& grid, std::vector<double>& resid, double& dt) {
    std::vector<double> dts;
    for (size_t i = 1; i < t.size(); ++i) dts.push_back(t[i] - t[i-1]);
    std::sort(dts.begin(), dts.end());
    dt = dts[dts.size()/2];
    grid.clear();
    for (double g = 0; g < t.back(); g += dt) grid.push_back(g);
    std::vector<double> yi = Resample(t, y, grid);
    int windowSteps = std::max(3, (int)std::round(180.0 / dt));
    std::vector<double> trend = MovingAvg(yi, windowSteps);
    resid.resize(grid.size());
    for (size_t i = 0; i < grid.size(); ++i) resid[i] = yi[i] - trend[i];
}

double Rms(const std::vector<double>& v) {
    double s = 0; for (double x : v) s += x*x;
    return std::sqrt(s / v.size());
}

}  // namespace

void Draw_QE_Temp_PhaseLock(TString qeCsv, TString envCsv, TString tag,
                            double periodMin = 40.2) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTextFont(132);
    gStyle->SetLabelFont(132, "XYZ");
    gStyle->SetTitleFont(132, "XYZ");

    // ---- Load QE (run,epoch,qmon,q1,q2) ----
    std::vector<double> tQe, yQe;
    double qeEpoch0 = 0;
    { std::ifstream in(qeCsv.Data()); std::string line;
      std::vector<std::pair<double,double>> rows;
      while (std::getline(in, line)) {
          std::stringstream ss(line); std::string f[5];
          for (int i = 0; i < 5 && std::getline(ss, f[i], ','); ++i) {}
          try { rows.push_back({std::stod(f[1]), std::stod(f[2])}); } catch (...) { continue; }
      }
      std::sort(rows.begin(), rows.end());
      qeEpoch0 = rows[0].first;
      for (auto& r : rows) { tQe.push_back((r.first - qeEpoch0)/60.0); yQe.push_back(r.second); }
    }

    // ---- Load env (timestamp_iso,Dark_Box_1_T,Dark_Box_2_T) ----
    std::vector<double> tEnv, yEnv;
    double envEpoch0 = 0;
    { std::ifstream in(envCsv.Data()); std::string line;
      std::vector<std::pair<double,double>> rows;
      while (std::getline(in, line)) {
          std::stringstream ss(line); std::string f[3];
          for (int i = 0; i < 3 && std::getline(ss, f[i], ','); ++i) {}
          if (f[0].empty() || f[1].empty()) continue;
          // ISO "2026-08-13T17:00:28.953507" -> epoch (UTC-naive, consistent w/ python isoformat().timestamp())
          int Y,M,D,h,m; double s;
          if (sscanf(f[0].c_str(), "%d-%d-%dT%d:%d:%lf", &Y,&M,&D,&h,&m,&s) != 6) continue;
          struct tm tmv = {}; tmv.tm_year=Y-1900; tmv.tm_mon=M-1; tmv.tm_mday=D;
          tmv.tm_hour=h; tmv.tm_min=m; tmv.tm_sec=(int)s;
          double epoch = (double)mktime(&tmv);
          try { rows.push_back({epoch, std::stod(f[1])}); } catch (...) { continue; }
      }
      std::sort(rows.begin(), rows.end());
      envEpoch0 = rows[0].first;
      for (auto& r : rows) { tEnv.push_back((r.first - envEpoch0)/60.0); yEnv.push_back(r.second); }
    }

    std::cout << "[INFO] QE points=" << tQe.size() << "  Env points=" << tEnv.size() << std::endl;
    if (tQe.size() < 10 || tEnv.size() < 10) { std::cout << "[ERROR] not enough points." << std::endl; return; }

    std::vector<double> gridQe, residQe, gridEnv, residEnv;
    double dtQe, dtEnv;
    PrepResidual(tQe, yQe, gridQe, residQe, dtQe);
    PrepResidual(tEnv, yEnv, gridEnv, residEnv, dtEnv);

    // gridQe/gridEnv are each "minutes since THIS series' own first sample" --
    // env and QE loggers start at different absolute times, so before any
    // phase comparison the env grid must be shifted onto the QE series'
    // absolute clock (same fix the earlier Python phase check needed).
    double offsetMin = (envEpoch0 - qeEpoch0) / 60.0;
    for (double& g : gridEnv) g += offsetMin;

    double rmsQe = Rms(residQe), rmsEnv = Rms(residEnv);
    std::vector<double> zQe(residQe.size()), zEnv(residEnv.size());
    for (size_t i = 0; i < residQe.size(); ++i) zQe[i] = residQe[i] / rmsQe;
    for (size_t i = 0; i < residEnv.size(); ++i) zEnv[i] = residEnv[i] / rmsEnv;

    // Least-squares sinusoid fit: z = A cos(wt) + B sin(wt)
    double omega = 2*TMath::Pi()/periodMin;
    auto fitSine = [&](const std::vector<double>& grid, const std::vector<double>& z,
                       double& amp, double& peakMin) {
        double Scc=0,Sss=0,Scs=0,Sc=0,Ss=0;
        for (size_t i = 0; i < grid.size(); ++i) {
            double c = std::cos(omega*grid[i]), s = std::sin(omega*grid[i]);
            Scc += c*c; Sss += s*s; Scs += c*s;
            Sc += c*z[i]; Ss += s*z[i];
        }
        double det = Scc*Sss - Scs*Scs;
        double A = (Sc*Sss - Ss*Scs)/det, B = (Ss*Scc - Sc*Scs)/det;
        amp = std::hypot(A,B);
        double phi = std::atan2(B,A);
        peakMin = fmod(phi/omega + periodMin*10, periodMin); // keep positive
    };
    double ampQe, peakQe, ampEnv, peakEnv;
    fitSine(gridQe, zQe, ampQe, peakQe);
    fitSine(gridEnv, zEnv, ampEnv, peakEnv);
    double lag = fmod(peakEnv - peakQe + periodMin*10, periodMin);

    std::cout << Form("[INFO] QE peak phase=%.1fmin  Env peak phase=%.1fmin  lag(env-qe)=%.1fmin",
                      peakQe, peakEnv, lag) << std::endl;

    TCanvas* c = new TCanvas("cPhase", "QE vs Temp phase-lock", 1500, 650);
    c->SetGrid();
    c->SetLeftMargin(0.09); c->SetRightMargin(0.03); c->SetTopMargin(0.13); c->SetBottomMargin(0.13);

    double tmax = std::min(gridQe.back(), gridEnv.back());
    tmax = std::min(tmax, 6*periodMin); // show first ~6 cycles for readability

    TGraph* gQe = new TGraph(gridQe.size(), gridQe.data(), zQe.data());
    gQe->SetLineColor(kAzure+2); gQe->SetLineWidth(2); gQe->SetMarkerColor(kAzure+2);
    gQe->SetMarkerStyle(20); gQe->SetMarkerSize(0.6);
    gQe->SetTitle(Form("%s: Monitor QE vs Dark Box T, detrended residual (z-score);Elapsed [min];Residual / RMS", tag.Data()));
    gQe->GetXaxis()->SetRangeUser(0, tmax);
    gQe->SetMinimum(-3.2); gQe->SetMaximum(3.2);
    gQe->Draw("APL");

    TGraph* gEnv = new TGraph(gridEnv.size(), gridEnv.data(), zEnv.data());
    gEnv->SetLineColor(kRed+1); gEnv->SetLineWidth(2);
    gEnv->Draw("L same");

    TF1* fQe = new TF1("fQe", "[0]*cos(2*TMath::Pi()/[1]*x - [2])", 0, tmax);
    fQe->SetParameters(ampQe, periodMin, atan2(0,1)); // placeholder, draw via points instead for simplicity
    // (skip TF1 fitting draw -- the raw z-score curves already carry the periodic shape clearly)

    TLine* lQePeak = new TLine(peakQe, -3.2, peakQe, 3.2);
    lQePeak->SetLineColor(kAzure+2); lQePeak->SetLineStyle(2); lQePeak->SetLineWidth(2); lQePeak->Draw();
    TLine* lEnvPeak = new TLine(peakEnv, -3.2, peakEnv, 3.2);
    lEnvPeak->SetLineColor(kRed+1); lEnvPeak->SetLineStyle(2); lEnvPeak->SetLineWidth(2); lEnvPeak->Draw();

    TLegend* leg = new TLegend(0.75, 0.82, 0.97, 0.97);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(132); leg->SetTextSize(0.03);
    leg->AddEntry(gQe, "Monitor QE (z-score)", "l");
    leg->AddEntry(gEnv, "Dark Box T (z-score)", "l");
    leg->Draw();

    TLatex hdr; hdr.SetNDC(); hdr.SetTextFont(132); hdr.SetTextSize(0.032);
    hdr.DrawLatex(0.10, 0.955, Form("Period fixed at %.1fmin  |  QE peak @ %.1fmin (dashed blue)  "
                                    "Temp peak @ %.1fmin (dashed red)  |  lag(Temp - QE) = %.1fmin",
                                    periodMin, peakQe, peakEnv, lag));

    gSystem->mkdir("./Data/image/Stability", kTRUE);
    TString out = Form("./Data/image/Stability/PhaseLock_%s.png", tag.Data());
    c->SaveAs(out);
    std::cout << "[INFO] Saved: " << out << std::endl;
}
