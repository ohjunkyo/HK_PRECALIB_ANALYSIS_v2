// Usage:  root -l -b -q Draw_HVScan_Report.C

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

// ============================================================ Shared block list
struct HVBlock {
    const char* date;
    int lo, hi;         // inclusive run range actually on disk for this block
    const char* label;  // matches the Overlay Uniformity report's tag labels
};

static const std::vector<HVBlock> kBlocks = {
    {"20260829", 100, 145, "-200V"},
    {"20260829", 200, 245, "-150V"},
    {"20260829", 0,   45,  "-100V"},
    {"20260827", 0,   45,  "-50V"},
    {"20260815", 0,   45,  "Recom."},
    {"20260827", 100, 145, "+50V"},
    {"20260830", 0,   45,  "+100V"},
    {"20260828", 0,   45,  "+150V"},
    {"20260828", 100, 145, "+200V"},
};

int kColors[] = {kRed, kOrange+1, kPink+7, kGreen+2, kBlack,
                 kCyan+2, kAzure+2, kViolet+1, kBlue+1};

const double ADC_TO_MV = 0.1220703125;
const double THRESHOLD_MV = 3.0;

const double MIN_PV_RATIO = 1.1;
const double MIN_REL_ERROR = 0.003;

TString BlockPath(const HVBlock& b, int run) {
    return Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", b.date, run);
}
int CenterRun(const HVBlock& b) { return b.lo + 11; }

// ============================================================ Page 1: Gain Curve
bool BlockMeanGain(const HVBlock& b, int ch, double& meanOut, double& errOut) {
    int centerRun = CenterRun(b);
    if (centerRun > b.hi) return false;
    TString path = BlockPath(b, centerRun);
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
            std::cout << "[WARNING] " << b.date << "_" << Form("%03d", centerRun) << " ch" << ch
                      << ": excluded, P/V=" << pv << " < " << MIN_PV_RATIO
                      << " (pedestal/SPE peaks not resolved)" << std::endl;
            f->Close();
            return false;
        }
        if (std::isfinite(spe) && spe > 0 && std::isfinite(speErr) && speErr >= 0) {
            double relErr = speErr / spe;
            if (relErr < MIN_REL_ERROR) {
                std::cout << "[INFO] " << b.date << "_" << Form("%03d", centerRun) << " ch" << ch
                          << ": error floor applied (raw relErr=" << Form("%.4f", relErr * 100)
                          << "% -> " << MIN_REL_ERROR * 100 << "%)" << std::endl;
                speErr = spe * MIN_REL_ERROR;
            }
            meanOut = spe; errOut = speErr; ok = true;
        }
    }
    f->Close();
    return ok;
}

bool BlockHV(const HVBlock& b, int& hv2, int& hv3) {
    for (int run = b.lo; run <= b.hi; ++run) {
        TString path = BlockPath(b, run);
        if (gSystem->AccessPathName(path)) continue;
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        TTree* info = (TTree*)f->Get("RunInfo");
        if (!info) { f->Close(); continue; }
        info->SetBranchAddress("HV2", &hv2);
        info->SetBranchAddress("HV3", &hv3);
        info->GetEntry(0);
        f->Close();
        return true;
    }
    return false;
}

double CaenVoltageError(double V) { return 2.0 + 0.0002 * std::abs(V); }

const double PC_TO_GAIN_1E7 = 1e-12 / 1.602176634e-19 / 1e7;

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

TCanvas* BuildGainCurvePage() {
    std::vector<double> hvR1, gR1, gR1e, hvR2, gR2, gR2e;
    for (const auto& b : kBlocks) {
        int hv2 = 0, hv3 = 0;
        if (!BlockHV(b, hv2, hv3)) {
            std::cerr << "[WARN] No RunInfo found for " << b.date << "_" << b.lo << "~" << b.hi << ", skipping.\n";
            continue;
        }
        double m1, e1, m2, e2;
        if (BlockMeanGain(b, 1, m1, e1)) {
            hvR1.push_back(hv2); gR1.push_back(m1 * PC_TO_GAIN_1E7); gR1e.push_back(e1 * PC_TO_GAIN_1E7);
        }
        if (BlockMeanGain(b, 2, m2, e2)) {
            hvR2.push_back(hv3); gR2.push_back(m2 * PC_TO_GAIN_1E7); gR2e.push_back(e2 * PC_TO_GAIN_1E7);
        }
        std::cout << "[INFO] " << b.date << "_" << b.lo << "~" << b.hi
                  << "  HV2=" << hv2 << " HV3=" << hv3 << std::endl;
    }

    TCanvas* c = new TCanvas("cGainCurve", "Gain Curve", 1400, 850);
    c->Divide(2, 1);
    if (!hvR1.empty()) {
        std::vector<double> hvR1e(hvR1.size());
        for (size_t k = 0; k < hvR1.size(); ++k) hvR1e[k] = CaenVoltageError(hvR1[k]);
        TGraphErrors* g1 = new TGraphErrors(hvR1.size(), hvR1.data(), gR1.data(), hvR1e.data(), gR1e.data());
        FitAndDrawGain(c->cd(1), g1, "EM6400 (Rot1)");
    }
    if (!hvR2.empty()) {
        std::vector<double> hvR2e(hvR2.size());
        for (size_t k = 0; k < hvR2.size(); ++k) hvR2e[k] = CaenVoltageError(hvR2[k]);
        TGraphErrors* g2 = new TGraphErrors(hvR2.size(), hvR2.data(), gR2.data(), hvR2e.data(), gR2e.data());
        FitAndDrawGain(c->cd(2), g2, "EL5150 (Rot2)");
    }
    return c;
}

// ============================================================ Page 1b: TTS Curve
bool BlockMeanTTS(const HVBlock& b, int ch, double& tnsOut, double& errOut) {
    int centerRun = CenterRun(b);
    if (centerRun > b.hi) return false;
    TString path = BlockPath(b, centerRun);
    if (gSystem->AccessPathName(path)) return false;
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { if (f) f->Close(); return false; }
    TTree* tr = (TTree*)f->Get(Form("tree_ch%d", ch));
    bool ok = false;
    if (tr && tr->GetEntries() > 0) {
        double tts = 0, ttsRelErr = 0;   // rms_exG in samples; rms_exG_err is RELATIVE
        tr->SetBranchAddress("rms_exG", &tts);
        tr->SetBranchAddress("rms_exG_err", &ttsRelErr);
        tr->GetEntry(0);
        if (std::isfinite(tts) && tts > 0) {
            tnsOut = tts * 2.0;   // 2ns/sample
            errOut = (ttsRelErr > 0) ? tnsOut * ttsRelErr : 0.0;
            ok = true;
        }
    }
    f->Close();
    return ok;
}

void DrawTTSPanel(TVirtualPad* pad, TGraphErrors* g, const char* pmtLabel) {
    pad->cd();
    pad->SetLeftMargin(0.16);
    pad->SetGridx(); pad->SetGridy();
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(1);
    gStyle->SetEndErrorSize(6);

    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
    g->SetLineWidth(1);
    g->SetTitle(Form("%s : TTS vs HV;Voltage [V];TTS [ns]", pmtLabel));
    g->GetXaxis()->SetLabelFont(132); g->GetXaxis()->SetTitleFont(132);
    g->GetYaxis()->SetLabelFont(132); g->GetYaxis()->SetTitleFont(132);
    g->GetXaxis()->SetLabelSize(0.038); g->GetXaxis()->SetTitleSize(0.042);
    g->GetYaxis()->SetLabelSize(0.038); g->GetYaxis()->SetTitleSize(0.042);
    g->Draw("AP");
}

TCanvas* BuildTTSCurvePage() {
    std::vector<double> hvR1, tR1, tR1e, hvR2, tR2, tR2e;
    for (const auto& b : kBlocks) {
        int hv2 = 0, hv3 = 0;
        if (!BlockHV(b, hv2, hv3)) continue;
        double t1, e1, t2, e2;
        if (BlockMeanTTS(b, 1, t1, e1)) { hvR1.push_back(hv2); tR1.push_back(t1); tR1e.push_back(e1); }
        if (BlockMeanTTS(b, 2, t2, e2)) { hvR2.push_back(hv3); tR2.push_back(t2); tR2e.push_back(e2); }
    }

    TCanvas* c = new TCanvas("cTTSCurve", "TTS Curve", 1400, 650);
    c->Divide(2, 1);
    if (!hvR1.empty()) {
        std::vector<double> hvR1e(hvR1.size());
        for (size_t k = 0; k < hvR1.size(); ++k) hvR1e[k] = CaenVoltageError(hvR1[k]);
        TGraphErrors* g1 = new TGraphErrors(hvR1.size(), hvR1.data(), tR1.data(), hvR1e.data(), tR1e.data());
        DrawTTSPanel(c->cd(1), g1, "EM6400 (Rot1)");
    }
    if (!hvR2.empty()) {
        std::vector<double> hvR2e(hvR2.size());
        for (size_t k = 0; k < hvR2.size(); ++k) hvR2e[k] = CaenVoltageError(hvR2[k]);
        TGraphErrors* g2 = new TGraphErrors(hvR2.size(), hvR2.data(), tR2.data(), hvR2e.data(), tR2e.data());
        DrawTTSPanel(c->cd(2), g2, "EL5150 (Rot2)");
    }
    return c;
}

// ============================================================ Page 2: Pulse/Time/Charge
struct DistDef { const char* branch; const char* label; double xlo, xhi; bool showThreshold; bool normalizePeak; };
static const std::vector<DistDef> kDists = {
    {"Max_ch",  "Pulse Height [ADC Counts]", -20, 400, true,  false},
    {"Diff_ch", "Time [ADC Samples] (peak-norm.)", 170, 220, false, true},
    {"Pico_ch", "Charge [pC]",               -2,  10,  false, false},
};

void DrawDistPad(TVirtualPad* pad, const DistDef& d, int ch, const char* chLabel) {
    pad->cd();
    if (!d.normalizePeak) pad->SetLogy();
    pad->SetGridx(); pad->SetGridy();
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(1);

    TLegend* leg = new TLegend(0.58, 0.55, 0.94, 0.90);
    leg->SetTextFont(132);
    leg->SetTextSize(0.028);
    leg->SetBorderSize(1);

    struct Loaded { TH1D* h; int blockIdx; };
    std::vector<Loaded> loaded;
    double autoLo = 1e18, autoHi = -1e18;

    for (size_t i = 0; i < kBlocks.size(); ++i) {
        const auto& b = kBlocks[i];
        int centerRun = CenterRun(b);
        if (centerRun > b.hi) continue;
        TString path = BlockPath(b, centerRun);
        if (gSystem->AccessPathName(path)) continue;
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        TH1D* h = (TH1D*)f->Get(Form("%s%d", d.branch, ch));
        if (!h) { f->Close(); continue; }
        h->SetDirectory(0);
        f->Close();

        h->GetYaxis()->UnZoom();

        h->SetLineColor(kColors[i % 9]);
        h->SetLineWidth(1);
        h->SetTitle(Form("%s : %s;%s;Entries", chLabel, d.label, d.label));
        h->GetXaxis()->SetLabelFont(132); h->GetXaxis()->SetTitleFont(132);
        h->GetYaxis()->SetLabelFont(132); h->GetYaxis()->SetTitleFont(132);
        h->GetXaxis()->SetLabelSize(0.038); h->GetXaxis()->SetTitleSize(0.042);
        h->GetYaxis()->SetLabelSize(0.038); h->GetYaxis()->SetTitleSize(0.042);

        if (d.xlo < d.xhi) {
            h->GetXaxis()->SetRangeUser(d.xlo, d.xhi);
        } else {
            int peakBin = h->GetMaximumBin();
            double peakX = h->GetXaxis()->GetBinCenter(peakBin);
            autoLo = std::min(autoLo, peakX - 30);
            autoHi = std::max(autoHi, peakX + 30);
            h->GetYaxis()->SetTitle("Normalized [a.u.]");
        }

        if (d.normalizePeak && h->GetMaximum() > 0) h->Scale(1.0 / h->GetMaximum());

        loaded.push_back({h, (int)i});
    }

    if (d.xlo >= d.xhi && autoHi > autoLo) {
        for (auto& L : loaded) L.h->GetXaxis()->SetRangeUser(autoLo, autoHi);
    }

    double drawMax = 0;
    for (auto& L : loaded) if (L.h->GetMaximum() > drawMax) drawMax = L.h->GetMaximum();

    std::vector<Loaded> drawOrder(loaded.rbegin(), loaded.rend());
    for (size_t k = 0; k < drawOrder.size(); ++k) {
        if (k == 0 && drawMax > 0) {
            drawOrder[k].h->SetMaximum(drawMax * (d.normalizePeak ? 1.25 : 3.0));
            if (d.normalizePeak) drawOrder[k].h->SetMinimum(0);
        }
        drawOrder[k].h->DrawCopy(k == 0 ? "HIST" : "HIST SAME");
    }

    std::vector<TH1D*> keep;
    for (auto& L : loaded) {
        leg->AddEntry(L.h, kBlocks[L.blockIdx].label, "l");
        keep.push_back(L.h);
    }

    if (d.showThreshold) {
        double thrADC = THRESHOLD_MV / ADC_TO_MV;
        TLine* thrLine = new TLine(thrADC, 0.5, thrADC, drawMax > 0 ? drawMax * 3.0 : 1e5);
        thrLine->SetLineColor(kGray + 2);
        thrLine->SetLineStyle(2);
        thrLine->SetLineWidth(2);
        //thrLine->Draw();
       // leg->AddEntry(thrLine, "3mV threshold", "l");
        
    }

    leg->Draw();
}

TCanvas* BuildDistributionPage() {
    TCanvas* c = new TCanvas("cDistHVScan", "Pulse/Time/Charge vs HV", 2000, 1900);
    c->Divide(3, 3);
    const char* chLabels[3] = {"Mon (EM2740)", "Rot1 (EM6400)", "Rot2 (EL5150)"};
    int padIdx = 1;
    for (const auto& d : kDists) {
        for (int ch = 0; ch < 3; ++ch) {
            DrawDistPad(c->cd(padIdx++), d, ch, chLabels[ch]);
        }
    }
    return c;
}

// ============================================================ Entry point
void Draw_HVScan_Report() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFont(132, "");
    gStyle->SetTitleFontSize(0.040);

    TCanvas* cGain = BuildGainCurvePage();
    TCanvas* cTTS  = BuildTTSCurvePage();
    TCanvas* cDist = BuildDistributionPage();

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString firstDate = kBlocks.front().date, lastDate = kBlocks.front().date;
    for (const auto& b : kBlocks) {
        if (TString(b.date) < firstDate) firstDate = b.date;
        if (TString(b.date) > lastDate) lastDate = b.date;
    }
    TString outPath = outDir + "HVScan_Report_" + firstDate + "_" + lastDate + ".pdf";

    cGain->Print(outPath + "(");
    cTTS->Print(outPath);
    cDist->Print(outPath + ")");
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
