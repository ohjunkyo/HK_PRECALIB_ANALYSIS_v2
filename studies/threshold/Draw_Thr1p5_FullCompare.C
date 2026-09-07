// Usage:  root -l -b -q Draw_Thr1p5_FullCompare.C
//
// Full comparison, 8/15 Recom. block (46 angles), production 3mV vs the
// diagnostic 1.5mV reprocess (Analysis_Threshold_mV changed in BOTH the
// amplitude cut AND the diff/timing-crossing calculation together -- see
// project memory "pending-adaptive-threshold-fix", 2026-09-01 decision to
// keep them as one shared constant). Pages: QE vs angle, Gain vs angle,
// TTS vs angle, Charge distribution (center run).
//
// DIAGNOSTIC ONLY -- reads Data/FinalResult (existing 3mV production) and
// Data/FinalResult_thr1p5 (new 1.5mV test run), does not touch either.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include <vector>

const double PC_TO_GAIN_1E7 = 1e-12 / 1.602176634e-19 / 1e7;
const char* kDate = "20260815";
const int kLo = 0, kHi = 45;

struct RunVals { double angle=0, qe=0, qeErr=0, gain=0, gainErr=0, tts=0, ttsErr=0; };

std::vector<RunVals> LoadSeries(const char* dir, int ch) {
    std::vector<RunVals> out;
    for (int run = kLo; run <= kHi; ++run) {
        TString path = Form("%s/precal_result_kor_run_%s_%03d.root", dir, kDate, run);
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        TTree* info = (TTree*)f->Get("RunInfo");
        TTree* tr = (TTree*)f->Get(Form("tree_ch%d", ch));
        int tilt2 = 0, tilt3 = 0;
        if (info) { info->SetBranchAddress("RawTiltAngle2", &tilt2); info->SetBranchAddress("RawTiltAngle3", &tilt3); info->GetEntry(0); }
        RunVals v;
        v.angle = (ch == 1) ? tilt2 : tilt3;
        if (tr) {
            double spe=0, speErr=0, qe=0, qeErr=0, tts=0, ttsErr=0;
            tr->SetBranchAddress("spe_mean", &spe); tr->SetBranchAddress("spe_mean_error", &speErr);
            tr->SetBranchAddress("relativeQE", &qe); tr->SetBranchAddress("relativeQE_err", &qeErr);
            tr->SetBranchAddress("rms_exG", &tts); tr->SetBranchAddress("rms_exG_err", &ttsErr);
            tr->GetEntry(0);
            v.gain = spe * PC_TO_GAIN_1E7; v.gainErr = speErr * PC_TO_GAIN_1E7;
            v.qe = qe; v.qeErr = qeErr;
            v.tts = tts * 2.0; v.ttsErr = (tts > 0) ? ttsErr * tts * 2.0 : 0;  // rms_exG_err stored relative
        }
        f->Close();
        out.push_back(v);
    }
    return out;
}

void DrawVsAngle(TVirtualPad* pad, int ch, const char* chLabel, const char* what,
                  std::vector<RunVals>& v3, std::vector<RunVals>& v15) {
    pad->cd();
    gPad->SetGridx(); gPad->SetGridy();
    int n = v3.size();
    std::vector<double> ang(n), y3(n), y3e(n), y15(n), y15e(n), ang0(n, 0);
    for (int i = 0; i < n; ++i) {
        ang[i] = v3[i].angle;
        if (TString(what) == "QE")   { y3[i]=v3[i].qe;   y3e[i]=v3[i].qeErr;   y15[i]=v15[i].qe;   y15e[i]=v15[i].qeErr; }
        if (TString(what) == "Gain") { y3[i]=v3[i].gain; y3e[i]=v3[i].gainErr; y15[i]=v15[i].gain; y15e[i]=v15[i].gainErr; }
        if (TString(what) == "TTS")  { y3[i]=v3[i].tts;  y3e[i]=v3[i].ttsErr;  y15[i]=v15[i].tts;  y15e[i]=v15[i].ttsErr; }
    }
    TGraphErrors* g3 = new TGraphErrors(n, ang.data(), y3.data(), ang0.data(), y3e.data());
    TGraphErrors* g15 = new TGraphErrors(n, ang.data(), y15.data(), ang0.data(), y15e.data());
    g3->SetMarkerStyle(20); g3->SetMarkerColor(kRed+1); g3->SetLineColor(kRed+1); g3->SetMarkerSize(0.7);
    g15->SetMarkerStyle(21); g15->SetMarkerColor(kBlue+1); g15->SetLineColor(kBlue+1); g15->SetMarkerSize(0.7);
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(g3, "P"); mg->Add(g15, "P");
    const char* yTitle = (TString(what)=="QE") ? "QE [%]" : (TString(what)=="Gain") ? "Gain [#times10^{7}]" : "TTS [ns]";
    mg->SetTitle(Form("%s : %s vs Tilt Angle;Tilt Angle [deg];%s", chLabel, what, yTitle));
    mg->Draw("AP");
    TLegend* leg = new TLegend(0.6, 0.15, 0.9, 0.3);
    leg->SetBorderSize(1); leg->SetFillColor(kWhite); leg->SetTextFont(132); leg->SetTextSize(0.03);
    leg->AddEntry(g3, "3mV (production)", "p");
    leg->AddEntry(g15, "1.5mV (test)", "p");
    leg->Draw();
}

void DrawChargeDist(TVirtualPad* pad, int ch, const char* chLabel) {
    pad->cd();
    pad->SetLogy();
    gPad->SetGridx(); gPad->SetGridy();
    TFile* f3 = TFile::Open(Form("./Data/FinalResult/precal_result_kor_run_%s_011.root", kDate), "READ");
    TFile* f15 = TFile::Open(Form("./Data/FinalResult_thr1p5/precal_result_kor_run_%s_011.root", kDate), "READ");
    TH1D* h3 = f3 ? (TH1D*)f3->Get(Form("Pico_ch%d", ch)) : nullptr;
    TH1D* h15 = f15 ? (TH1D*)f15->Get(Form("Pico_ch%d", ch)) : nullptr;
    if (!h3 || !h15) return;
    h3->SetDirectory(0); h15->SetDirectory(0);
    h3->GetYaxis()->UnZoom(); h15->GetYaxis()->UnZoom();
    f3->Close(); f15->Close();
    h3->SetLineColor(kRed+1); h3->SetLineWidth(2);
    h15->SetLineColor(kBlue+1); h15->SetLineWidth(2);
    h3->SetTitle(Form("%s : Charge Distribution (center run, Recom.);Charge [pC];Entries", chLabel));
    h3->GetXaxis()->SetRangeUser(-2, 10);
    double ymax = std::max(h3->GetMaximum(), h15->GetMaximum());
    h3->SetMaximum(ymax * 2);
    h3->DrawCopy("HIST");
    h15->DrawCopy("HIST SAME");
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.88);
    leg->SetBorderSize(1); leg->SetFillColor(kWhite); leg->SetTextFont(132); leg->SetTextSize(0.03);
    leg->AddEntry(h3, "3mV (production)", "l");
    leg->AddEntry(h15, "1.5mV (test)", "l");
    leg->Draw();
}

void Draw_Thr1p5_FullCompare() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");

    TCanvas* c = new TCanvas("cThr1p5", "3mV vs 1.5mV full compare", 1500, 750);
    TString outPath = "./Data/image/ScanReport/Thr1p5_FullCompare_20260815.pdf";

    bool first = true;
    for (int ch : {1, 2}) {
        const char* chLabel = (ch == 1) ? "Rot1 (EM6400)" : "Rot2 (EL5150)";
        auto v3  = LoadSeries("./Data/FinalResult", ch);
        auto v15 = LoadSeries("./Data/FinalResult_thr1p5", ch);

        for (const char* what : {"QE", "Gain", "TTS"}) {
            c->Clear();
            DrawVsAngle(c, ch, chLabel, what, v3, v15);
            if (first) { c->Print(outPath + "("); first = false; }
            else c->Print(outPath);
        }
        c->Clear();
        DrawChargeDist(c, ch, chLabel);
        c->Print(outPath);
    }
    c->Print(outPath + ")");
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
