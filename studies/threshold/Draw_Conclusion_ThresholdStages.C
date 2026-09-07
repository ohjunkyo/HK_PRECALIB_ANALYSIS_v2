//
// Stage definitions (recovered from the stored branches):
//   Raw QE        = relativeQE_raw                       (amplitude cut only)
//   Timing cut    = relativeQE / (1 - dark_frac_phc/100) (amplitude AND timing)
//   Timing+Dark   = relativeQE                           (accidentals removed)
//   Monitor-Norm  = Timing+Dark(ch) / Timing+Dark(ch0)   (dimensionless)

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TRegexp.h"
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

struct HVBlock { const char* date; int blk; const char* label; };

static const std::vector<HVBlock> kHV = {
    {"20260829", 1, "-200V"}, {"20260829", 2, "-150V"}, {"20260829", 0, "-100V"},
    {"20260827", 0, "-50V"},  {"20260815", 0, "Recom."}, {"20260827", 1, "+50V"},
    {"20260830", 0, "+100V"}, {"20260828", 0, "+150V"},  {"20260828", 1, "+200V"},
};

static const int kNAngle = 46;

// Pulls |theta_Ham| out of the Max_ch1 histogram title (e.g.
// "...#theta_{Ham}: 81.5#circ)"), which read_ntp_v7.C writes from the SAME
// angle-conversion the DAQ used to move the stage -- avoids recomputing (and
// possibly mis-recomputing) the rot-fold mapping here. ch1/ch2 have equal
// |theta_Ham| for a given run by construction, so one read covers both.
// Returns -1 if the file/histogram/pattern can't be found (caller must skip).
static double GetAbsThetaHam(const char* path) {
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { if (f) f->Close(); return -1; }
    TH1D* h = (TH1D*)f->Get("Max_ch1");
    if (!h) { f->Close(); return -1; }
    TString title = h->GetTitle();
    f->Close();

    TString tag = "#theta_{Ham}: ";
    int i = title.Index(tag);
    if (i < 0) return -1;
    i += tag.Length();
    int j = title.Index("#circ", i);
    if (j < 0) return -1;
    TString num = title(i, j - i);
    return std::fabs(num.Atof());
}

// Per-channel, per-stage accumulator.
struct Acc {
    double sum[4] = {0, 0, 0, 0};   // Raw, Timing, Timing+Dark, MonNorm
    long   n[4]   = {0, 0, 0, 0};
    void add(int s, double v) { if (std::isfinite(v)) { sum[s] += v; n[s]++; } }
    double mean(int s) const { return n[s] ? sum[s] / n[s] : 0.0; }
};

// Reads one FinalResult file and appends its four stage values for ch1/ch2.
static void Harvest(const char* path, Acc acc[3]) {
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { if (f) f->Close(); return; }

    double qraw[3] = {0}, qcorr[3] = {0}, dfrac[3] = {0};
    bool ok[3] = {false, false, false};
    for (int ch = 0; ch < 3; ++ch) {
        TTree* t = (TTree*)f->Get(Form("tree_ch%d", ch));
        if (!t || t->GetEntries() < 1) continue;
        t->SetBranchAddress("relativeQE_raw", &qraw[ch]);
        t->SetBranchAddress("relativeQE",     &qcorr[ch]);
        t->SetBranchAddress("dark_frac_phc",  &dfrac[ch]);
        t->GetEntry(0);
        ok[ch] = true;
    }
    // The monitor normalisation is only meaningful when ch0 itself is valid.
    for (int ch = 1; ch <= 2; ++ch) {
        if (!ok[ch]) continue;
        double denom = 1.0 - dfrac[ch] / 100.0;
        acc[ch].add(0, qraw[ch]);
        if (denom > 0) acc[ch].add(1, qcorr[ch] / denom);
        acc[ch].add(2, qcorr[ch]);
        if (ok[0] && qcorr[0] > 0) acc[ch].add(3, qcorr[ch] / qcorr[0]);
    }
    f->Close();
}

// One 3-group bar panel (Raw / Timing / Timing+Dark) for a single channel.
static void DrawStagePanel(TVirtualPad* pad, const Acc& a15, const Acc& a3,
                           const char* title, const char* yTitle,
                           int nStage, int stageOffset, const char* const* names) {
    pad->cd();
    pad->SetGridy();
    pad->SetBottomMargin(0.14);
    pad->SetLeftMargin(0.16);
    gStyle->SetGridColor(kGray);

    TH1D* h3  = new TH1D(Form("h3_%s",  title), "", nStage, 0, nStage);
    TH1D* h15 = new TH1D(Form("h15_%s", title), "", nStage, 0, nStage);
    double top = 0;
    for (int i = 0; i < nStage; ++i) {
        int s = stageOffset + i;
        h3->SetBinContent(i + 1, a3.mean(s));
        h15->SetBinContent(i + 1, a15.mean(s));
        h3->GetXaxis()->SetBinLabel(i + 1, names[i]);
        h15->GetXaxis()->SetBinLabel(i + 1, names[i]);
        top = std::max(top, std::max(a3.mean(s), a15.mean(s)));
    }

    h3->SetFillColor(kGray + 1);   h3->SetLineColor(kGray + 2);   h3->SetBarWidth(0.36); h3->SetBarOffset(0.10);
    h15->SetFillColor(kRed - 7);   h15->SetLineColor(kRed + 1);   h15->SetBarWidth(0.36); h15->SetBarOffset(0.52);

    h3->SetTitle(Form("%s;;%s", title, yTitle));
    h3->SetMaximum(top * 1.70);
    h3->SetMinimum(0);
    h3->SetStats(0);
    h3->GetXaxis()->SetLabelFont(132); h3->GetXaxis()->SetLabelSize(0.058);
    h3->GetYaxis()->SetLabelFont(132); h3->GetYaxis()->SetTitleFont(132);
    h3->GetYaxis()->SetTitleSize(0.052); h3->GetYaxis()->SetTitleOffset(1.45);
    h3->SetTitleFont(132, "");
    h3->Draw("bar");
    h15->Draw("bar same");

    TLatex tx; tx.SetTextFont(132); tx.SetTextSize(0.044); tx.SetTextAlign(21);
    for (int i = 0; i < nStage; ++i) {
        int s = stageOffset + i;
        double v3 = a3.mean(s), v15 = a15.mean(s);
        tx.SetTextColor(kBlack);   tx.DrawLatex(i + 0.28, v3  + top * 0.035, Form("%.3f", v3));
        tx.SetTextColor(kRed + 2); tx.DrawLatex(i + 0.70, v15 + top * 0.035, Form("%.3f", v15));
        if (v3 > 0) {
            tx.SetTextColor(kBlue + 1);
            tx.DrawLatex(i + 0.5, top * 1.22, Form("%+.1f%%", (v15 - v3) / v3 * 100.0));
        }
    }

    TLegend* leg = new TLegend(0.30, 0.815, 0.925, 0.88);
    leg->SetNColumns(2);
    leg->SetTextFont(132); leg->SetTextSize(0.042);
    leg->SetBorderSize(0); leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
    leg->AddEntry(h3,  "Fixed 3mV",            "f");
    leg->AddEntry(h15, "Pedestal #mu+5#sigma (1.5mV)", "f");
    leg->SetMargin(0.25);
    leg->Draw();
}

// Grouped bar panel across the 9 HV points, one stage, one channel. Shows
// whether the 3mV->1.5mV gap (top labels) is flat across HV or drifts with it
// -- the averaged panel alone cannot answer that.
static void DrawByHVPanel(TVirtualPad* pad, const Acc a15[9], const Acc a3[9],
                          int nHVUsed, int stage, const char* title, const char* yTitle) {
    pad->cd();
    pad->SetGridy();
    pad->SetBottomMargin(0.16);
    pad->SetLeftMargin(0.14);
    gStyle->SetGridColor(kGray);

    TH1D* h3  = new TH1D(Form("hv3_%s_%d",  title, stage), "", nHVUsed, 0, nHVUsed);
    TH1D* h15 = new TH1D(Form("hv15_%s_%d", title, stage), "", nHVUsed, 0, nHVUsed);
    double top = 0;
    for (int i = 0; i < nHVUsed; ++i) {
        h3->SetBinContent(i + 1, a3[i].mean(stage));
        h15->SetBinContent(i + 1, a15[i].mean(stage));
        h3->GetXaxis()->SetBinLabel(i + 1, kHV[i].label);
        h15->GetXaxis()->SetBinLabel(i + 1, kHV[i].label);
        top = std::max(top, std::max(a3[i].mean(stage), a15[i].mean(stage)));
    }

    h3->SetFillColor(kGray + 1);  h3->SetLineColor(kGray + 2);  h3->SetBarWidth(0.36); h3->SetBarOffset(0.10);
    h15->SetFillColor(kRed - 7);  h15->SetLineColor(kRed + 1);  h15->SetBarWidth(0.36); h15->SetBarOffset(0.52);

    h3->SetTitle(Form("%s;;%s", title, yTitle));
    h3->SetMaximum(top * 1.42);
    h3->SetMinimum(0);
    h3->SetStats(0);
    h3->GetXaxis()->SetLabelFont(132); h3->GetXaxis()->SetLabelSize(0.048);
    h3->GetYaxis()->SetLabelFont(132); h3->GetYaxis()->SetTitleFont(132);
    h3->GetYaxis()->SetTitleSize(0.045); h3->GetYaxis()->SetTitleOffset(1.35);
    h3->SetTitleFont(132, "");
    h3->Draw("bar");
    h15->Draw("bar same");

    TLatex tx; tx.SetTextFont(132); tx.SetTextSize(0.032); tx.SetTextAlign(21); tx.SetTextColor(kBlue + 1);
    for (int i = 0; i < nHVUsed; ++i) {
        double v3 = a3[i].mean(stage), v15 = a15[i].mean(stage);
        if (v3 > 0) tx.DrawLatex(i + 0.5, top * 1.15, Form("%+.1f%%", (v15 - v3) / v3 * 100.0));
    }

    TLegend* leg = new TLegend(0.14, 0.88, 0.60, 0.955);
    leg->SetNColumns(2);
    leg->SetTextFont(132); leg->SetTextSize(0.034);
    leg->SetBorderSize(0); leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
    leg->AddEntry(h3,  "Fixed 3mV", "f");
    leg->AddEntry(h15, "Pedestal #mu+5#sigma (1.5mV)", "f");
    leg->Draw();
}

// Runs the whole harvest+draw+report pipeline for one angle selection.
//   maxAbsTheta = 70.0  -> "Inner" (|theta_Ham| < 70 deg only)
//   maxAbsTheta = 999.0 -> "AllAngles" (no cut)
// tag is used in filenames/titles/logs so the two passes never collide.
static void RunOnePass(double maxAbsTheta, const char* tag, const char* titleSuffix,
                       TCanvas*& outC1, TCanvas*& outC2) {
    const int nHV = (int)kHV.size();
    Acc a15[3], a3[3];                 // pooled over all HV (averaged figure)
    Acc a15_hv[9][3], a3_hv[9][3];     // per-HV, so the gap's HV-dependence can be shown
    int nFile15 = 0, nFile3 = 0, nSkippedAngle = 0;

    for (int h = 0; h < nHV; ++h) {
        const auto& hv = kHV[h];
        for (int a = 0; a < kNAngle; ++a) {
            int run = hv.blk * 100 + a;
            TString base = Form("precal_result_kor_run_%s_%03d.root", hv.date, run);
            TString p15  = "./Data/FinalResult/" + base;
            TString p3   = "./Data/FinalResult_thr3compare/" + base;
            // Only count a point when BOTH thresholds have it, so the two
            // averages are over an identical (angle, HV) set.
            if (gSystem->AccessPathName(p15) || gSystem->AccessPathName(p3)) continue;

            double theta = GetAbsThetaHam(p15);
            if (theta < 0) continue;               // couldn't read the angle -- skip, don't guess
            if (theta >= maxAbsTheta) { nSkippedAngle++; continue; }

            Harvest(p15, a15); Harvest(p15, a15_hv[h]); nFile15++;
            Harvest(p3,  a3);  Harvest(p3,  a3_hv[h]);  nFile3++;
        }
    }

    printf("\n[INFO][%s] paired (angle,HV) points used: %d of %d (excluded by angle cut: %d)\n",
           tag, nFile15, (int)kHV.size() * kNAngle, nSkippedAngle);
    if (nFile15 == 0) { printf("[ERROR][%s] nothing to draw.\n", tag); return; }

    static const char* stageNames[3] = {"Raw QE", "Timing cut", "Timing + Dark"};
    static const char* normName[1]   = {"Normalized by monitor"};

    TCanvas* c = new TCanvas(Form("cConcl_%s", tag), Form("Threshold stages, %s", tag), 1800, 1273);
    c->Divide(2, 2);
    DrawStagePanel(c->cd(1), a15[1], a3[1], Form("Rotation 1 (EM6400)%s", titleSuffix), "QE [%]", 3, 0, stageNames);
    DrawStagePanel(c->cd(2), a15[1], a3[1], Form("Rot 1: Normalized by monitor%s", titleSuffix), "QE_{Rot}/QE_{Mon} [a.u.]", 1, 3, normName);
    DrawStagePanel(c->cd(3), a15[2], a3[2], Form("Rotation 2 (EL5150)%s", titleSuffix), "QE [%]", 3, 0, stageNames);
    DrawStagePanel(c->cd(4), a15[2], a3[2], Form("Rot 2: Normalized by monitor%s", titleSuffix), "QE_{Rot}/QE_{Mon} [a.u.]", 1, 3, normName);

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
   // c->SaveAs(outDir + Form("Conclusion_ThresholdStages_AllHV_%s.png", tag));   // quick-look only

    // ---- Second figure: same comparison broken out per HV point, so the gap's
    // HV-dependence is visible instead of hidden inside one pooled average.
    // Stage shown: Timing+Dark (Corr QE) and Monitor-Norm -- the two numbers the
    // slide's conclusion actually rests on.
    // a15_hv/a3_hv are indexed [HV][channel]; DrawByHVPanel wants a length-9
    // array indexed by HV for ONE fixed channel, so repack per channel first.
    Acc a15_ch1[9], a3_ch1[9], a15_ch2[9], a3_ch2[9];
    for (int h = 0; h < nHV; ++h) {
        a15_ch1[h] = a15_hv[h][1]; a3_ch1[h] = a3_hv[h][1];
        a15_ch2[h] = a15_hv[h][2]; a3_ch2[h] = a3_hv[h][2];
    }

    TCanvas* c2 = new TCanvas(Form("cConclByHV_%s", tag), Form("Threshold stages by HV, %s", tag), 1800, 1273);
    c2->Divide(2, 2);
    DrawByHVPanel(c2->cd(1), a15_ch1, a3_ch1, nHV, 2, Form("Rot1 (EM6400) Corr QE by HV%s", titleSuffix), "QE [%]");
    DrawByHVPanel(c2->cd(2), a15_ch1, a3_ch1, nHV, 3, Form("Rot1 Monitor-Norm by HV%s", titleSuffix),     "QE_{Rot}/QE_{Mon} [a.u.]");
    DrawByHVPanel(c2->cd(3), a15_ch2, a3_ch2, nHV, 2, Form("Rot2 (EL5150) Corr QE by HV%s", titleSuffix), "QE [%]");
    DrawByHVPanel(c2->cd(4), a15_ch2, a3_ch2, nHV, 3, Form("Rot2 Monitor-Norm by HV%s", titleSuffix),     "QE_{Rot}/QE_{Mon} [a.u.]");

//    c2->SaveAs(outDir + Form("Conclusion_ThresholdStages_ByHV_%s.png", tag));   // quick-look only

    // Console table, same shape as the slide's summary table.
    printf("\n===== [%s] Threshold comparison, averaged over all HV x all angles (n=%d) =====\n", tag, nFile15);
    printf("%-22s %18s %18s %10s\n", "QE Calculation process", "3mV (Rot1/Rot2)", "1.5mV (Rot1/Rot2)", "Diff");
    const char* rows[4] = {"Raw", "Timing cut", "Timing+Dark (Corr)", "Monitor-Norm"};
    for (int s = 0; s < 4; ++s) {
        double r1_3 = a3[1].mean(s), r2_3 = a3[2].mean(s);
        double r1_15 = a15[1].mean(s), r2_15 = a15[2].mean(s);
        double d1 = r1_3 > 0 ? (r1_15 - r1_3) / r1_3 * 100.0 : 0;
        double d2 = r2_3 > 0 ? (r2_15 - r2_3) / r2_3 * 100.0 : 0;
        printf("%-22s %8.3f /%8.3f %8.3f /%8.3f %6.1f%% /%5.1f%%\n",
               rows[s], r1_3, r2_3, r1_15, r2_15, d1, d2);
    }

    printf("\n===== [%s] Corr QE (Timing+Dark) and Monitor-Norm, by HV point =====\n", tag);
    printf("%-8s | %-30s | %-30s\n", "HV", "Corr QE  3mV / 1.5mV / Diff", "Mon-Norm 3mV / 1.5mV / Diff");
    for (int ch = 1; ch <= 2; ++ch) {
        printf("--- ch%d (%s) ---\n", ch, ch == 1 ? "Rot1 EM6400" : "Rot2 EL5150");
        for (int h = 0; h < nHV; ++h) {
            double c3 = a3_hv[h][ch].mean(2), c15 = a15_hv[h][ch].mean(2);
            double m3 = a3_hv[h][ch].mean(3), m15 = a15_hv[h][ch].mean(3);
            double dc = c3 > 0 ? (c15 - c3) / c3 * 100.0 : 0;
            double dm = m3 > 0 ? (m15 - m3) / m3 * 100.0 : 0;
            printf("%-8s | %6.3f /%6.3f /%+6.1f%%   | %6.3f /%6.3f /%+6.1f%%\n",
                   kHV[h].label, c3, c15, dc, m3, m15, dm);
        }
    }
    outC1 = c;
    outC2 = c2;
}

void Draw_Conclusion_ThresholdStages() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");
    gStyle->SetPaperSize(26.0, 14.4);   // cm, matches the 1800x1000 canvas

    // Per the Gain-curve study: the Gain fit itself is only reliable within
    // |theta_Ham| < 70deg (beyond that, -20% to -40% fit bias from low-light
    // Peak/Valley degradation -- not real physics). Report both so the
    // Inner-only conclusion can't be confused with one diluted by that
    // extreme-angle artefact, and the difference between the two IS itself a
    // check for whether the artefact is leaking into the threshold comparison.
    TCanvas *cInner1, *cInner2, *cAll1, *cAll2;
    RunOnePass(70.0, "Inner70", "  [|#theta_{Ham}| < 70#circ]", cInner1, cInner2);
    RunOnePass(999.0, "AllAngles", "  [all angles]", cAll1, cAll2);

    // Single combined deliverable: one multi-page PDF, in the order
    // Inner-averaged -> Inner-by-HV -> AllAngles-averaged -> AllAngles-by-HV.
    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString combined = outDir + "Conclusion_ThresholdStages.pdf";
    cInner1->Print(combined + "(");
    cInner2->Print(combined);
    cAll1->Print(combined);
    cAll2->Print(combined + ")");
    printf("\n[INFO] Saved combined PDF (4 pages): %s\n", combined.Data());
}
