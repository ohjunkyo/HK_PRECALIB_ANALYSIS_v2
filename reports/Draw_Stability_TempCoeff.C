// Usage:  root -l -b -q Draw_Stability_TempCoeff.C
//
// Temperature coefficient (dQE/dT): linear fit of Raw QE vs Dark Box 1
// temperature for each PMT, over the 2026-08-13/14 Stability run
// (20260813_800~899). Reuses the same file-loading/env-log logic as
// Draw_Stability_v1.C.

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

struct QEPoint { double epoch, qe, qeErr; };
struct EnvPoint { double epoch; double db1t; };

TString EpochToDbTimestamp(double epoch) {
    time_t t = (time_t)epoch;
    struct tm tmv;
    localtime_r(&t, &tmv);
    char buf[32];
    strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tmv);
    return TString(buf);
}

std::vector<EnvPoint> ReadEnvLog(double t0, double t1) {
    std::vector<EnvPoint> out;
    const char* dbPath = "/home/precalkor/Integrated_Control_SW/HV_Control_SW/monitoring_log.db";
    if (gSystem->AccessPathName(dbPath)) { std::cout << "[WARNING] monitoring_log.db not found." << std::endl; return out; }
    TString lo = EpochToDbTimestamp(t0), hi = EpochToDbTimestamp(t1);
    TString cmd = Form("sqlite3 -separator , \"%s\" "
                       "\"SELECT timestamp, Dark_Box_1_T FROM monitoring_data "
                       "WHERE timestamp BETWEEN '%s' AND '%s' ORDER BY timestamp;\"",
                       dbPath, lo.Data(), hi.Data());
    FILE* p = gSystem->OpenPipe(cmd, "r");
    if (!p) return out;
    char line[256];
    while (fgets(line, sizeof(line), p)) {
        struct tm tmv = {};
        double t1v;
        char tsbuf[40];
        if (sscanf(line, "%39[^,],%lf", tsbuf, &t1v) != 2) continue;
        if (sscanf(tsbuf, "%d-%d-%dT%d:%d:%d", &tmv.tm_year, &tmv.tm_mon, &tmv.tm_mday,
                   &tmv.tm_hour, &tmv.tm_min, &tmv.tm_sec) != 6) continue;
        tmv.tm_year -= 1900; tmv.tm_mon -= 1; tmv.tm_isdst = -1;
        double epoch = (double)mktime(&tmv);
        out.push_back({epoch, t1v});
    }
    gSystem->ClosePipe(p);
    return out;
}

double NearestTemp(const std::vector<EnvPoint>& env, double epoch) {
    if (env.empty()) return -999;
    double best = 1e18, bestT = -999;
    for (auto& e : env) { double d = std::abs(e.epoch - epoch); if (d < best) { best = d; bestT = e.db1t; } }
    return bestT;
}

void DrawChannelFit(TVirtualPad* pad, const char* label, int ch, const std::vector<QEPoint>& qePts,
                     const std::vector<EnvPoint>& env, int color) {
    std::vector<double> temp, qe, qeErr, tempErr;
    for (auto& p : qePts) {
        double t = NearestTemp(env, p.epoch);
        if (t < -500) continue;
        temp.push_back(t); qe.push_back(p.qe); qeErr.push_back(p.qeErr); tempErr.push_back(0.05);
    }
    pad->cd();
    pad->SetGridx(); pad->SetGridy();
    if (temp.size() < 3) { std::cout << label << ": not enough points\n"; return; }

    TGraphErrors* g = new TGraphErrors(temp.size(), temp.data(), qe.data(), tempErr.data(), qeErr.data());
    g->SetMarkerStyle(20); g->SetMarkerColor(color); g->SetLineColor(color); g->SetMarkerSize(0.7);
    g->SetTitle(Form("%s: Raw QE vs Dark Box 1 Temp;Temperature [#circC];Raw QE [%%]", label));
    g->Draw("AP");

    TF1* fit = new TF1(Form("fit_%s", label), "pol1", *std::min_element(temp.begin(), temp.end()), *std::max_element(temp.begin(), temp.end()));
    g->Fit(fit, "RQ");
    fit->SetLineColor(kRed+1); fit->SetLineWidth(2); fit->Draw("SAME");

    double slope = fit->GetParameter(1), slopeErr = fit->GetParError(1);
    double meanQE = 0; for (double v : qe) meanQE += v; meanQE /= qe.size();
    double relPct = 100.0 * slope / meanQE;

    TLatex lat; lat.SetTextFont(132); lat.SetTextSize(0.04); lat.SetNDC();
    lat.DrawLatex(0.15, 0.85, Form("dQE/dT = %.4f #pm %.4f %%/#circC", slope, slopeErr));
    lat.DrawLatex(0.15, 0.79, Form("(%.2f%%/#circC relative)", relPct));

    std::cout << Form("%-6s  dQE/dT = %.4f +/- %.4f %%/degC  (%.2f%%rel/degC)  n=%zu",
                      label, slope, slopeErr, relPct, temp.size()) << std::endl;
}

void Draw_Stability_TempCoeff(TString tag = "20260813", int run_start = 800, int run_end = 899) {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");

    std::vector<QEPoint> qe[3];
    double t0 = -1, tLast = -1;
    static const char* kRawDirs[] = { "./Data/RAW/Laser", "/media/precalkor/HD-EDS-E/Data_Backup/RAW/Laser" };

    for (int run = run_start; run <= run_end; ++run) {
        TString path = Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", tag.Data(), run);
        if (gSystem->AccessPathName(path)) continue;

        double mtime = 0; struct stat st;
        for (const char* dir : kRawDirs) {
            TString rawPath = Form("%s/precal_raw_kor_run_%s_%03d.root", dir, tag.Data(), run);
            if (stat(rawPath.Data(), &st) == 0) { mtime = (double)st.st_mtime; break; }
        }
        if (mtime == 0) { if (stat(path.Data(), &st) == 0) mtime = (double)st.st_mtime; else continue; }
        if (t0 < 0) t0 = mtime;
        tLast = std::max(tLast, mtime);

        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        for (int c = 0; c < 3; ++c) {
            TTree* tr = (TTree*)f->Get(Form("tree_ch%d", c));
            if (!tr) continue;
            double qeRaw = 0, qeRawErr = 0;
            bool hasRaw = tr->GetBranch("relativeQE_raw") != nullptr;
            if (hasRaw) { tr->SetBranchAddress("relativeQE_raw", &qeRaw); tr->SetBranchAddress("relativeQE_raw_err", &qeRawErr); }
            tr->GetEntry(0);
            if (qeRaw > 0) qe[c].push_back({mtime, qeRaw, qeRawErr});
        }
        f->Close();
    }
    if (t0 < 0) { std::cout << "No files found.\n"; return; }

    std::vector<EnvPoint> env = ReadEnvLog(t0 - 300, tLast + 300);
    std::cout << "Env points: " << env.size() << std::endl;

    TCanvas* c = new TCanvas("cTempCoeff", "dQE/dT", 1500, 500);
    c->Divide(3, 1);
    const char* labels[3] = {"Mon.", "Rot1", "Rot2"};
    int colors[3] = {kBlack, kAzure+2, kRed+1};
    for (int ch = 0; ch < 3; ++ch) DrawChannelFit(c->cd(ch+1), labels[ch], ch, qe[ch], env, colors[ch]);

    TString outPath = "./Data/image/ScanReport/Stability_TempCoeff.pdf";
    c->SaveAs(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
