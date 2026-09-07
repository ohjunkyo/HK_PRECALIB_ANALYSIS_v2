// Usage:  root -l -b -q Draw_GainCurve_Center.C
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

// One HV point = one run (not a lo..hi angle-sweep range).
struct HVRun { const char* date; int run; };
static const std::vector<HVRun> kRuns = {
    {"20260903", 800}, {"20260903", 801}, {"20260903", 802}, {"20260903", 803},
    {"20260903", 804}, {"20260903", 805}, {"20260903", 806}, {"20260903", 807},
    {"20260903", 808},
};

// 0 = no cut. Was 1.1 (matching Draw_HVScan_Report.C's campaign-QA
// threshold), but for a quick center-only Gain-curve test the low-HV points
// are exactly the ones worth seeing even with a barely-resolved SPE peak --
// 2026-09-03, user: "P/V이거 하한선 없애면..?".
const double MIN_PV_RATIO = 0.0;
const double MIN_REL_ERROR = 0.003;
const double PC_TO_GAIN_1E7 = 1e-12 / 1.602176634e-19 / 1e7;

TString RunPath(const HVRun& r) {
    return Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", r.date, r.run);
}

bool RunGain(const HVRun& r, int ch, double& meanOut, double& errOut) {
    TString path = RunPath(r);
    if (gSystem->AccessPathName(path)) return false;
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { if (f) f->Close(); return false; }
    TTree* tr = (TTree*)f->Get(Form("tree_ch%d", ch));
    bool ok = false;
    if (tr && tr->GetEntries() > 0) {
        double spe = 0, speErr = 0, pv = 0;
        tr->SetBranchAddress("spe_mean", &spe);
        tr->SetBranchAddress("spe_mean_error", &speErr);
        if (tr->GetBranch("peak_to_valley")) tr->SetBranchAddress("peak_to_valley", &pv);
        tr->GetEntry(0);
        if (pv > 0 && pv < MIN_PV_RATIO) {
            std::cout << "[WARNING] " << r.date << "_" << Form("%03d", r.run) << " ch" << ch
                      << ": excluded, P/V=" << pv << " < " << MIN_PV_RATIO << std::endl;
            f->Close();
            return false;
        }
        if (std::isfinite(spe) && spe > 0 && std::isfinite(speErr) && speErr >= 0) {
            double relErr = speErr / spe;
            if (relErr < MIN_REL_ERROR) speErr = spe * MIN_REL_ERROR;
            meanOut = spe; errOut = speErr; ok = true;
        }
    }
    f->Close();
    return ok;
}

bool RunHV(const HVRun& r, int& hv2, int& hv3) {
    TString path = RunPath(r);
    if (gSystem->AccessPathName(path)) return false;
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { if (f) f->Close(); return false; }
    TTree* info = (TTree*)f->Get("RunInfo");
    if (!info) { f->Close(); return false; }
    info->SetBranchAddress("HV2", &hv2);
    info->SetBranchAddress("HV3", &hv3);
    info->GetEntry(0);
    f->Close();
    return true;
}

// PMT serial number, straight from RunInfo -- was hardcoded as "EM6400
// (Rot1)"/"EL5150 (Rot2)" in the plot label, which goes stale the moment a
// PMT is swapped -- 2026-09-03, user: "PMT 이름은 자동으로 불러와야지".
TString RunSN(const HVRun& r, int ch) {
    TString path = RunPath(r);
    if (gSystem->AccessPathName(path)) return "";
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { if (f) f->Close(); return ""; }
    TTree* info = (TTree*)f->Get("RunInfo");
    if (!info) { f->Close(); return ""; }
    char sn[64] = "";
    TString branch = Form("SN%d", ch + 1);   // ch1 -> SN2, ch2 -> SN3
    if (info->GetBranch(branch)) info->SetBranchAddress(branch, sn);
    info->GetEntry(0);
    f->Close();
    return TString(sn);
}

double CaenVoltageError(double V) { return 2.0 + 0.0002 * std::abs(V); }

void FitAndDrawGain(TVirtualPad* outerPad, TGraphErrors* g, const char* pmtLabel) {
    outerPad->cd();
    TPad* padTop = new TPad(Form("padTop_%s", pmtLabel), "", 0.0, 0.33, 1.0, 1.0);
    TPad* padBot = new TPad(Form("padBot_%s", pmtLabel), "", 0.0, 0.0, 1.0, 0.33);
    padTop->SetBottomMargin(0.02);
    padBot->SetTopMargin(0.03);
    padBot->SetBottomMargin(0.38);
    for (TPad* p : {padTop, padBot}) {
        p->SetLeftMargin(0.16);
        p->SetGridx(); p->SetGridy();
        p->Draw();
    }
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(1);
    gStyle->SetGridWidth(1);
    gStyle->SetEndErrorSize(6);

    padTop->cd();
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.7);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
    g->SetLineWidth(1);
    g->SetTitle(Form("%s : Gain vs HV;;Gain [#times10^{7}]", pmtLabel));
    g->GetXaxis()->SetLabelSize(0);
    g->GetYaxis()->SetLabelFont(132); g->GetYaxis()->SetTitleFont(132);
    g->GetYaxis()->SetLabelSize(0.060); g->GetYaxis()->SetTitleSize(0.065);
    g->GetYaxis()->SetTitleOffset(1.05);

    double vLo = g->GetX()[0], vHi = g->GetX()[0];
    for (int k = 1; k < g->GetN(); ++k) {
        vLo = std::min(vLo, g->GetX()[k]);
        vHi = std::max(vHi, g->GetX()[k]);
    }
    double axLo = vLo - (vHi - vLo) * 0.05, axHi = vHi + (vHi - vLo) * 0.05;
    g->Draw("AP");
    g->GetXaxis()->SetLimits(axLo, axHi);

    const double V0 = 1000.0;
    TF1* fit = new TF1(Form("fit_%s", pmtLabel), "[0]*pow(x/[1],[2])", vLo, vHi);
    fit->FixParameter(1, V0);
    fit->SetParameters(g->GetY()[0] * std::pow(V0 / g->GetX()[0], -7.0), V0, 7.0);
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(1);
    g->Fit(fit, "RQ");

    std::vector<int> order(g->GetN());
    for (int k = 0; k < g->GetN(); ++k) order[k] = k;
    std::sort(order.begin(), order.end(),
             [&](int a, int b) { return g->GetX()[a] < g->GetX()[b]; });
    std::cout << "\n[INFO] " << pmtLabel << " -- Gain vs Voltage:" << std::endl;
    for (int idx : order) {
        std::cout << Form("    V = %6.1f +- %.1f V   Gain = %.4f +- %.4f (x10^7)",
                          g->GetX()[idx], g->GetEX()[idx], g->GetY()[idx], g->GetEY()[idx])
                  << std::endl;
    }

    // Fit inverted (V = V0*(G/a)^(1/b)) at round 0.5-step Gain values, so the
    // report also reads directly as "what HV gives Gain=1.5x10^7" instead of
    // only the as-measured points -- 2026-09-03, user: "0.5~3.5 * 10^7도
    // 보여주면 안되나? 데이터 구간 내에서 해당하는 Gain의 딱 정수값이라
    // 해야하나?" (turned out to mean the printed table, not the plot axis).
    // Restricted to the DATA's own gain range -- extrapolating the fit
    // beyond measured points isn't a real prediction.
    {
        double gLo = g->GetY()[0], gHi = g->GetY()[0];
        for (int k = 1; k < g->GetN(); ++k) {
            gLo = std::min(gLo, g->GetY()[k]);
            gHi = std::max(gHi, g->GetY()[k]);
        }
        double a = fit->GetParameter(0), b = fit->GetParameter(2);
        std::cout << "[INFO] " << pmtLabel << " -- Voltage for round Gain values"
                      " (within measured range " << Form("%.2f", gLo) << "-"
                  << Form("%.2f", gHi) << " x10^7):" << std::endl;
        for (double gTarget = 0.5; gTarget <= gHi + 1e-9; gTarget += 0.5) {
            if (gTarget < gLo) continue;
            double v = V0 * std::pow(gTarget / a, 1.0 / b);
            std::cout << Form("    Gain = %.1f x10^7  ->  V = %6.1f V", gTarget, v) << std::endl;
        }
    }

    TLegend* box = new TLegend(0.16, 0.68, 0.55, 0.90);
    box->SetTextFont(132);
    box->SetTextSize(0.032);
    box->SetBorderSize(1);
    box->SetFillColor(kWhite);
    box->AddEntry(fit, pmtLabel, "l");
    box->AddEntry((TObject*)nullptr, "Fit: G = a #times V^{b}", "");
    box->AddEntry((TObject*)nullptr,
                 Form("a_{1kV} = %.3f #pm %.3f", fit->GetParameter(0), fit->GetParError(0)), "");
    box->AddEntry((TObject*)nullptr,
                 Form("b = %.3f #pm %.3f", fit->GetParameter(2), fit->GetParError(2)), "");
    box->AddEntry((TObject*)nullptr,
                 Form("#chi^{2}/ndf = %.2f / %d", fit->GetChisquare(), fit->GetNDF()), "");
    box->Draw();

    padBot->cd();
    std::vector<double> pullX(g->GetN()), pullY(g->GetN()), pullXE(g->GetN());
    for (int k = 0; k < g->GetN(); ++k) {
        pullX[k] = g->GetX()[k];
        pullXE[k] = g->GetEX()[k];
        double pred = fit->Eval(pullX[k]);
        double err = g->GetEY()[k];
        pullY[k] = (err > 0) ? (g->GetY()[k] - pred) / err : 0.0;
    }
    TGraphErrors* gPull = new TGraphErrors(g->GetN(), pullX.data(), pullY.data(), pullXE.data(), nullptr);
    gPull->SetMarkerStyle(20);
    gPull->SetMarkerSize(0.7);
    gPull->SetMarkerColor(kBlack);
    gPull->SetLineColor(kBlack);
    gPull->SetLineWidth(1);
    gPull->SetTitle(";Voltage [V];Pull");
    gPull->GetXaxis()->SetLabelFont(132); gPull->GetXaxis()->SetTitleFont(132);
    gPull->GetYaxis()->SetLabelFont(132); gPull->GetYaxis()->SetTitleFont(132);
    gPull->GetXaxis()->SetLabelSize(0.11); gPull->GetXaxis()->SetTitleSize(0.13);
    gPull->GetYaxis()->SetLabelSize(0.11); gPull->GetYaxis()->SetTitleSize(0.13);
    gPull->GetYaxis()->SetTitleOffset(0.45);
    gPull->GetYaxis()->SetNdivisions(505);
    gPull->GetXaxis()->SetLimits(axLo, axHi);
    gPull->Draw("AP");
    TLine* zero = new TLine(axLo, 0, axHi, 0);
    zero->SetLineStyle(2); zero->SetLineColor(kRed + 1); zero->SetLineWidth(1);
    zero->Draw();
}

void Draw_GainCurve_Center() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFont(132, "");
    gStyle->SetTitleFontSize(0.040);

    std::vector<double> hvR1, gR1, gR1e, hvR2, gR2, gR2e;
    for (const auto& r : kRuns) {
        int hv2 = 0, hv3 = 0;
        if (!RunHV(r, hv2, hv3)) {
            std::cerr << "[WARN] No RunInfo found for " << r.date << "_" << r.run << ", skipping.\n";
            continue;
        }
        double m1, e1, m2, e2;
        if (RunGain(r, 1, m1, e1)) {
            hvR1.push_back(hv2); gR1.push_back(m1 * PC_TO_GAIN_1E7); gR1e.push_back(e1 * PC_TO_GAIN_1E7);
        }
        if (RunGain(r, 2, m2, e2)) {
            hvR2.push_back(hv3); gR2.push_back(m2 * PC_TO_GAIN_1E7); gR2e.push_back(e2 * PC_TO_GAIN_1E7);
        }
        std::cout << "[INFO] " << r.date << "_" << r.run << "  HV2=" << hv2 << " HV3=" << hv3 << std::endl;
    }

    // First run that actually has a readable SN wins -- a PMT's SN doesn't
    // change run to run within one campaign, so any hit is representative.
    TString sn1, sn2;
    for (const auto& r : kRuns) {
        if (sn1.IsNull()) sn1 = RunSN(r, 1);
        if (sn2.IsNull()) sn2 = RunSN(r, 2);
        if (!sn1.IsNull() && !sn2.IsNull()) break;
    }
    TString label1 = sn1.IsNull() ? "Rot1" : (sn1 + " (Rot1)");
    TString label2 = sn2.IsNull() ? "Rot2" : (sn2 + " (Rot2)");

    TCanvas* c = new TCanvas("cGainCurveCenter", "Gain Curve (center)", 1400, 850);
    c->Divide(2, 1);
    if (!hvR1.empty()) {
        std::vector<double> hvR1e(hvR1.size());
        for (size_t k = 0; k < hvR1.size(); ++k) hvR1e[k] = CaenVoltageError(hvR1[k]);
        TGraphErrors* g1 = new TGraphErrors(hvR1.size(), hvR1.data(), gR1.data(), hvR1e.data(), gR1e.data());
        FitAndDrawGain(c->cd(1), g1, label1);
    }
    if (!hvR2.empty()) {
        std::vector<double> hvR2e(hvR2.size());
        for (size_t k = 0; k < hvR2.size(); ++k) hvR2e[k] = CaenVoltageError(hvR2[k]);
        TGraphErrors* g2 = new TGraphErrors(hvR2.size(), hvR2.data(), gR2.data(), hvR2e.data(), gR2e.data());
        FitAndDrawGain(c->cd(2), g2, label2);
    }

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + "GainCurve_Center.pdf";
    c->Print(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
