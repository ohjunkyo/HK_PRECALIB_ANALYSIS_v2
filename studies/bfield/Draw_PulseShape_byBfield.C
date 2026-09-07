// Usage:  root -l -b -q Draw_PulseShape_byBfield.C
//
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TF1.h"
#include "TVectorD.h"

const char* kCachePath = "./Data/image/ScanReport/PulseShape_Cache.root";
#include <cstring>
#include <cmath>

const double ADC_TO_MV = 0.1220703125;
const double NS_PER_SAMPLE = 2.0;
const double kThrMV = 1.5;    // same analysis threshold as prod_ntp_v7.C
const int    kFitLo = 550;    // display window (samples)
const int    kFitHi = 700;
const int    kHitLo = 585;
const int    kHitHi = 625;

const int    TriggerCh      = 3;        // config3.h
const double kTrigThrADC    = 14000.0;  // absolute raw-ADC threshold, as in prod_ntp_v7.C
const int    kValidDuration = 3;        // consecutive samples below threshold (Analysis.cpp)

inline int TimeAtThreshold(const UShort_t* w, double thr, int lo, int hi) {
    int consec = 0, first = -1;
    for (int s = lo; s < hi; ++s) {
        if ((double)w[s] < thr) {
            if (consec == 0) first = s;
            if (++consec >= kValidDuration) return first;
        } else consec = 0;
    }
    return -1;
}

inline double CrossAt(const std::vector<double>& v, double lev, int from, int to) {
    int step = (to > from) ? 1 : -1;
    for (int s = from; s != to; s += step) {
        double a = v[s], b = v[s + step];
        if ((a > lev && b <= lev) || (a < lev && b >= lev)) {
            if (a == b) return s;
            return s + step * (a - lev) / (a - b);
        }
    }
    return std::nan("");
}
const long long kMaxEvents = 3000;   // enough for a stable persistence plot, fast

struct HVPoint { const char* date; int run; const char* label; };
// B-field campaign (2026-08-31 ~ 09-01), center run (tilt=0, rot=135) of each
// block -- Rot1's B=0 label is per calibration intent, not the (mis-oriented)
// sensor readout, per [[usb-bus001-contention-20260904]]-style convention:
// see project memory / 2026-09-04 chat with user on this.
static const std::vector<HVPoint> kPoints = {
    {"20260831", 11,  "Nominal"}, {"20260831", 111, "B=0 (Rot1 target)"},
    {"20260901", 111, "Bx=-70mG"}, {"20260901", 211, "Bx=-150mG"},
};

void DrawPulseGrid(int ch, const char* chLabel, TTree* avgTree,
                    char labelBuf[32], Int_t& runNum, Int_t& chOut, Int_t& sampleOut, Double_t& voltOut) {
    gStyle->SetOptStat(0);
    static const char* dirs[2] = {"./Data/RAW/Laser", "/media/precalkor/HD-EDS-E/Data_Backup/RAW/Laser"};

    TString cname = Form("cPulse_ch%d", ch);
    TCanvas* c = new TCanvas(cname, "pulse", 1200, 1000);
    c->Divide(2, 2);   // 4 B-field points, not the HV campaign's 9

    std::vector<TGraph*> overlayGraphs;
    std::vector<double> overlayThrT, overlayTts;
    std::vector<TString> overlayLabel;
    std::vector<int> overlayColor;

    std::vector<double> thrCmpHV[3], thrCmpT[3], thrCmpTTS[3];

    for (size_t i = 0; i < kPoints.size(); ++i) {
        const auto& p = kPoints[i];
        TVirtualPad* pad = c->cd(i + 1);
        pad->SetRightMargin(0.14);
        pad->SetLeftMargin(0.14);

        TString path;
        for (int di = 0; di < 2; ++di) {
            TString cand = Form("%s/precal_raw_kor_run_%s_%03d.root", dirs[di], p.date, p.run);
            if (!gSystem->AccessPathName(cand)) { path = cand; break; }
        }
        if (path.IsNull()) { std::cout << p.label << " ch" << ch << ": RAW not found\n"; continue; }

        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        TTree* tr = (TTree*)f->Get("T");
        unsigned int RecordLength = 1024;
        if (tr->GetBranch("RecordLength")) tr->SetBranchAddress("RecordLength", &RecordLength);
        tr->GetEntry(0);
        int nSamples = RecordLength;

        UShort_t* buf = new UShort_t[8 * nSamples];
        tr->SetBranchStatus("*", 0); tr->SetBranchStatus("ADC", 1);
        tr->SetBranchAddress("ADC", buf);

        TString key = Form("%s_%d_ch%d", p.date, p.run, ch);
        TString hname = Form("h_pulse_%s", key.Data());
        TString htitle = Form("%s (%s_%03d);Sample [2ns];Voltage [mV]", p.label, p.date, p.run);

        TH2D* h = nullptr;
        TH1D* hDiff = nullptr;
        std::vector<double> sumV(nSamples, 0.0);
        long long tN = 0, n = 0;
        double phSum = 0, phSum2 = 0;

        TFile* fCache = TFile::Open(kCachePath, "UPDATE");
        bool cacheHit = false;
        if (fCache && !fCache->IsZombie()) {
            TH2D* hc = (TH2D*)fCache->Get("Waveform_h2d_" + key);
            TH1D* hdc = (TH1D*)fCache->Get("Waveform_hdiff_" + key);
            TVectorD* sv = (TVectorD*)fCache->Get("Waveform_sumV_" + key);
            TVectorD* meta = (TVectorD*)fCache->Get("Waveform_meta_" + key);
            if (hc && hdc && sv && meta && sv->GetNoElements() == nSamples) {
                h = (TH2D*)hc->Clone(hname); h->SetDirectory(0);
                hDiff = (TH1D*)hdc->Clone("hdiff_" + key); hDiff->SetDirectory(0);
                for (int s = 0; s < nSamples; ++s) sumV[s] = (*sv)[s];
                n = (long long)(*meta)[0]; tN = (long long)(*meta)[1];
                phSum = (*meta)[2]; phSum2 = (*meta)[3];
                cacheHit = true;
            }
        }

        if (!cacheHit) {
            h = new TH2D(hname, htitle, nSamples, 0, nSamples, 600, -50, 10);
            h->SetDirectory(0);
            hDiff = new TH1D("hdiff_" + key, "diff", 600, -100, 500);
            hDiff->SetDirectory(0);
            n = tr->GetEntries();
            for (long long e = 0; e < n; ++e) {
                tr->GetEntry(e);
                double ped = 0;
                for (int s = 50; s < 400; ++s) ped += buf[ch * nSamples + s];
                ped /= 350.0;
                for (int s = 0; s < nSamples; ++s)
                    h->Fill(s, (buf[ch * nSamples + s] - ped) * ADC_TO_MV);

                double thrADC = ped - kThrMV / ADC_TO_MV;
                int tFall = TimeAtThreshold(buf + ch * nSamples, thrADC, kHitLo, kHitHi);
                if (tFall < 0) continue;                // no pulse -- excluded below
                int tTrig = TimeAtThreshold(buf + TriggerCh * nSamples, kTrigThrADC, 0, nSamples);
                if (tTrig >= 0) hDiff->Fill(tFall - tTrig);
                tN++;

                double evPeak = 0;
                for (int s = kHitLo; s < kHitHi && s < nSamples; ++s) {
                    double v = (buf[ch * nSamples + s] - ped) * ADC_TO_MV;
                    if (v < evPeak) evPeak = v;
                }
                phSum += evPeak; phSum2 += evPeak * evPeak;

                for (int s = 0; s < nSamples; ++s)
                    sumV[s] += (buf[ch * nSamples + s] - ped) * ADC_TO_MV;
            }

            if (fCache && !fCache->IsZombie()) {
                fCache->cd();
                h->Clone("Waveform_h2d_" + key)->Write("", TObject::kOverwrite);
                hDiff->Clone("Waveform_hdiff_" + key)->Write("", TObject::kOverwrite);
                TVectorD sv(nSamples); for (int s = 0; s < nSamples; ++s) sv[s] = sumV[s];
                sv.Write("Waveform_sumV_" + key, TObject::kOverwrite);
                TVectorD meta(4); meta[0] = n; meta[1] = tN; meta[2] = phSum; meta[3] = phSum2;
                meta.Write("Waveform_meta_" + key, TObject::kOverwrite);
            }
        }
        // Threshold comparison (1/2/3mV) -- 2026-09-02, user: "Threshold
        // 1,2,3mV에 대한 값들도 비교 및 그림 그려줄 수 있나". Kept separate from
        // the main 1.5mV analysis pass above (which stays the reference/
        // official-threshold numbers): a different threshold changes which
        // events count as hits, so each needs its own mean waveform, hit
        // count, and diff/TTS histogram -- not just a different crossing
        // level on the SAME mean waveform. Reuses the same RAW read (same
        // event loop) when a fresh read is needed at all, rather than
        // re-opening the file per threshold.
        static const double kThrCompareList[3] = {1.0, 2.0, 3.0};
        std::vector<double> cmpSumV[3];
        long long cmpTN[3] = {0, 0, 0};
        double cmpPhSum[3] = {0, 0, 0}, cmpPhSum2[3] = {0, 0, 0};
        TH1D* cmpHDiff[3] = {nullptr, nullptr, nullptr};
        for (int t = 0; t < 3; ++t) cmpSumV[t].assign(nSamples, 0.0);

        bool cmpCacheHit = true;
        if (fCache && !fCache->IsZombie()) {
            for (int t = 0; t < 3; ++t) {
                TString tkey = Form("%s_thr%dmV", key.Data(), (int)kThrCompareList[t]);
                TH1D* hdc = (TH1D*)fCache->Get("Waveform_thrcmp_hdiff_" + tkey);
                TVectorD* sv = (TVectorD*)fCache->Get("Waveform_thrcmp_sumV_" + tkey);
                TVectorD* meta = (TVectorD*)fCache->Get("Waveform_thrcmp_meta_" + tkey);
                if (!hdc || !sv || !meta || sv->GetNoElements() != nSamples) { cmpCacheHit = false; break; }
                cmpHDiff[t] = (TH1D*)hdc->Clone("cmp_hdiff_" + tkey); cmpHDiff[t]->SetDirectory(0);
                for (int s = 0; s < nSamples; ++s) cmpSumV[t][s] = (*sv)[s];
                cmpTN[t] = (long long)(*meta)[0]; cmpPhSum[t] = (*meta)[1]; cmpPhSum2[t] = (*meta)[2];
            }
        } else cmpCacheHit = false;

        if (!cmpCacheHit) {
            for (int t = 0; t < 3; ++t) {
                cmpHDiff[t] = new TH1D(Form("cmp_hdiff_%s_thr%dmV", key.Data(), (int)kThrCompareList[t]),
                                       "diff", 600, -100, 500);
                cmpHDiff[t]->SetDirectory(0);
            }
            long long nCmp = tr->GetEntries();
            for (long long e = 0; e < nCmp; ++e) {
                tr->GetEntry(e);
                double ped = 0;
                for (int s = 50; s < 400; ++s) ped += buf[ch * nSamples + s];
                ped /= 350.0;

                for (int t = 0; t < 3; ++t) {
                    double thrADC = ped - kThrCompareList[t] / ADC_TO_MV;
                    int tFall = TimeAtThreshold(buf + ch * nSamples, thrADC, kHitLo, kHitHi);
                    if (tFall < 0) continue;
                    int tTrig = TimeAtThreshold(buf + TriggerCh * nSamples, kTrigThrADC, 0, nSamples);
                    if (tTrig >= 0) cmpHDiff[t]->Fill(tFall - tTrig);
                    cmpTN[t]++;

                    double evPeak = 0;
                    for (int s = kHitLo; s < kHitHi && s < nSamples; ++s) {
                        double v = (buf[ch * nSamples + s] - ped) * ADC_TO_MV;
                        if (v < evPeak) evPeak = v;
                    }
                    cmpPhSum[t] += evPeak; cmpPhSum2[t] += evPeak * evPeak;

                    for (int s = 0; s < nSamples; ++s)
                        cmpSumV[t][s] += (buf[ch * nSamples + s] - ped) * ADC_TO_MV;
                }
            }

            if (fCache && !fCache->IsZombie()) {
                fCache->cd();
                for (int t = 0; t < 3; ++t) {
                    TString tkey = Form("%s_thr%dmV", key.Data(), (int)kThrCompareList[t]);
                    cmpHDiff[t]->Clone("Waveform_thrcmp_hdiff_" + tkey)->Write("", TObject::kOverwrite);
                    TVectorD sv(nSamples); for (int s = 0; s < nSamples; ++s) sv[s] = cmpSumV[t][s];
                    sv.Write("Waveform_thrcmp_sumV_" + tkey, TObject::kOverwrite);
                    TVectorD meta(3); meta[0] = cmpTN[t]; meta[1] = cmpPhSum[t]; meta[2] = cmpPhSum2[t];
                    meta.Write("Waveform_thrcmp_meta_" + tkey, TObject::kOverwrite);
                }
            }
        }

        for (int t = 0; t < 3; ++t) {
            std::vector<double> cmpMeanV(nSamples);
            for (int s = 0; s < nSamples; ++s) cmpMeanV[s] = cmpTN[t] ? cmpSumV[t][s] / cmpTN[t] : 0.0;

            int cpk = kFitLo;
            for (int s = kFitLo; s <= kFitHi && s < nSamples; ++s)
                if (cmpMeanV[s] < cmpMeanV[cpk]) cpk = s;
            double cLThr = -kThrCompareList[t];
            double cTThrR = CrossAt(cmpMeanV, cLThr, cpk, kFitLo);

            double cTts = std::nan("");
            if (cmpHDiff[t]->GetEntries() > 100) {
                double cDiffPk = cmpHDiff[t]->GetBinCenter(cmpHDiff[t]->GetMaximumBin());
                TF1 cfg("cfg", "gaus", cDiffPk - 8, cDiffPk + 8);
                if (cmpHDiff[t]->Fit(&cfg, "RQ0N") == 0) cTts = cfg.GetParameter(2) * NS_PER_SAMPLE;
            }

            thrCmpHV[t].push_back(i);
            thrCmpT[t].push_back(std::isfinite(cTThrR) ? cTThrR * NS_PER_SAMPLE : std::nan(""));
            thrCmpTTS[t].push_back(cTts);
        }
        for (int t = 0; t < 3; ++t) delete cmpHDiff[t];

        if (fCache) fCache->Close();
        f->Close();
        delete[] buf;

        std::cout << "[" << p.label << "] ch" << ch << (cacheHit ? " (cache)" : " (RAW)") << "\n";

        std::vector<double> tAxis(nSamples), meanV(nSamples);
        for (int s = 0; s < nSamples; ++s) { tAxis[s] = s; meanV[s] = tN ? sumV[s] / tN : 0.0; }

        strncpy(labelBuf, p.label, 31); labelBuf[31] = '\0';
        runNum = p.run;
        chOut = ch;
        for (int s = 0; s < nSamples; ++s) {
            sampleOut = s; voltOut = meanV[s];
            avgTree->Fill();
        }

        pad->SetLogz();
        /// h->GetXaxis()->SetRangeUser(0, nSamples);
        h->GetXaxis()->SetRangeUser(kFitLo, kFitHi);

        h->GetYaxis()->SetRangeUser(-30, 3);

//        h->GetXaxis()->SetLabelSize(0.045); //h->GetXaxis()->SetTitleFont(132);
  //      h->GetXaxis()->SetTitleSize(0.045); 
    //    h->GetYaxis()->SetTitleSize(0.045); 
      //  h->GetYaxis()->SetLabelSize(0.045); //h->GetYaxis()->SetTitleFont(132);

        h->Draw("COLZ");

        TGraph* gMeanHalo = new TGraph(nSamples, tAxis.data(), meanV.data());
        gMeanHalo->SetLineColor(kBlack); gMeanHalo->SetLineWidth(4);
        gMeanHalo->Draw("L SAME");
        TGraph* gMean = new TGraph(nSamples, tAxis.data(), meanV.data());
        gMean->SetLineColor(kMagenta); gMean->SetLineWidth(1);
        gMean->Draw("L SAME");

        int pk = kFitLo;
        for (int s = kFitLo; s <= kFitHi && s < nSamples; ++s)
            if (meanV[s] < meanV[pk]) pk = s;
        double amp = meanV[pk];                       // negative peak [mV]
        double l10 = 0.10 * amp, l50 = 0.50 * amp, l90 = 0.90 * amp;

        double lThr = -kThrMV;
        double tThrR = CrossAt(meanV, lThr, pk, kFitLo);   // enters cut (leading)
        double tThrF = CrossAt(meanV, lThr, pk, kFitHi);   // exits cut (falling)

        double t10r = CrossAt(meanV, l10, pk, kFitLo);   // leading edge
        double t90r = CrossAt(meanV, l90, pk, kFitLo);
        double t90f = CrossAt(meanV, l90, pk, kFitHi);   // falling edge
        double t10f = CrossAt(meanV, l10, pk, kFitHi);
        double t50r = CrossAt(meanV, l50, pk, kFitLo);
        double t50f = CrossAt(meanV, l50, pk, kFitHi);

        auto edge = [&](double tLo, double tHi, int col) {
            if (!std::isfinite(tLo) || !std::isfinite(tHi)) return;
            TLine* ln = new TLine(tLo, l10, tHi, l90);
            ln->SetLineColor(col); ln->SetLineWidth(1);
            ln->Draw();
            TMarker* m1 = new TMarker(tLo, l10, 20); m1->SetMarkerColor(col); m1->SetMarkerSize(0.7); m1->Draw();
            TMarker* m2 = new TMarker(tHi, l90, 20); m2->SetMarkerColor(col); m2->SetMarkerSize(0.7); m2->Draw();
        };
        //edge(t10r, t90r, kCyan + 2);
        //edge(t10f, t90f, kRed + 1);

        if (std::isfinite(t50r) && std::isfinite(t50f)) {
            TLine* fw = new TLine(t50r, l50, t50f, l50);
            fw->SetLineColor(kOrange + 7); fw->SetLineWidth(1);
            fw->Draw();
        }

        if (std::isfinite(tThrR)) {
            TMarker* mThr = new TMarker(tThrR, lThr, kFullCircle);
            mThr->SetMarkerColor(kSpring); mThr->SetMarkerSize(1.3);
            mThr->Draw();
        }

        double rise  = (std::isfinite(t90r) && std::isfinite(t10r)) ? (t90r - t10r) * NS_PER_SAMPLE : std::nan("");
        double fall  = (std::isfinite(t10f) && std::isfinite(t90f)) ? (t10f - t90f) * NS_PER_SAMPLE : std::nan("");
        double fwhm  = (std::isfinite(t50f) && std::isfinite(t50r)) ? (t50f - t50r) * NS_PER_SAMPLE : std::nan("");
        double tts = std::nan(""), diffPk = std::nan("");
        if (hDiff->GetEntries() > 100) {
            diffPk = hDiff->GetBinCenter(hDiff->GetMaximumBin());
            TF1 fg("fg", "gaus", diffPk - 8, diffPk + 8);
            if (hDiff->Fit(&fg, "RQ0N") == 0) tts = fg.GetParameter(2) * NS_PER_SAMPLE;
        }

        double phMean = tN ? phSum / tN : std::nan("");
        double phRms  = tN > 1 ? std::sqrt(std::max(0.0, phSum2 / tN - phMean * phMean)) : std::nan("");
        double hitFrac = n ? 100.0 * tN / n : 0.0;

        TPaveText* pt = new TPaveText(0.46, 0.148, 0.855, 0.42, "NDC");
        pt->SetFillColor(kWhite); pt->SetFillStyle(1001);
        pt->SetBorderSize(1); pt->SetLineColor(kGray + 2);
        pt->SetTextFont(132); pt->SetTextSize(0.038); pt->SetTextAlign(12);
        pt->SetMargin(0.03);
        auto row = [&](const char* s, int col) { pt->AddText(s)->SetTextColor(col); };
        row(Form("Peak = %.2f mV", phMean), kBlack);
        row(Form("Leading = %.1f ns", rise), kBlack);
        row(Form("Falling = %.1f ns", fall), kBlack);
        row(Form("FWHM = %.1f ns", fwhm), kBlack);
        row(Form("TTS(#sigma) = %.2f ns", tts), kBlack);
        if (std::isfinite(tThrR))
            row(Form("Thr(1.5mV) @ %.1f ns", tThrR * NS_PER_SAMPLE), kRed);
        pt->Draw();

        printf("[%-6s] ch%d  <PH>=%6.2f+-%5.2f mV  mean-wf peak=%6.2f  rise=%5.1f  fall=%5.1f  "
                "FWHM=%5.1f  TTS_sigma=%5.2f ns  thr@%7.1f ns  fall90@%7.1f ns  fall10@%7.1f ns  hit=%5.1f%% (n=%lld/%lld)\n",
                p.label, ch, phMean, phRms, amp, rise, fall, fwhm, tts,
                tThrR * NS_PER_SAMPLE, t90f * NS_PER_SAMPLE, t10f * NS_PER_SAMPLE, hitFrac, tN, n);

        // Blue (-200V) -> red (+200V) gradient, ordered by kPoints position.
        double frac = kPoints.size() > 1 ? (double)i / (double)(kPoints.size() - 1) : 0.0;
        int col = TColor::GetColor((int)(255 * frac), 0, (int)(255 * (1 - frac)));
        // Color map
        TGraph* gOv = new TGraph(nSamples, tAxis.data(), meanV.data());
        gOv->SetLineColor(col); gOv->SetLineWidth(2);
        overlayGraphs.push_back(gOv);
        overlayThrT.push_back(tThrR);
        overlayTts.push_back(tts);
        overlayLabel.push_back(p.label);
        overlayColor.push_back(col);
    }

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + Form("PulseShape_byBfield_ch%d.pdf", ch);
    c->Print(outPath + "(", "pdf");   // page 1: the 3x3 grid


    TCanvas* c2 = new TCanvas(Form("cOverlay_ch%d", ch), "overlay", 1200, 800);
    c2->SetLeftMargin(0.10); c2->SetRightMargin(0.04); c2->SetGrid();
    TH1F* frame = c2->DrawFrame(kFitLo, -20, kFitHi, 3);
    frame->SetTitle(Form("%s;Sample [2ns];Voltage [mV]", chLabel));
    TLegend* leg = new TLegend(0.72, 0.20, 0.95, 0.55);
    leg->SetBorderSize(3); leg->SetFillStyle(1); leg->SetTextFont(132); leg->SetTextSize(0.028);
    for (size_t i = 0; i < overlayGraphs.size(); ++i) {
        overlayGraphs[i]->Draw("L SAME");
        leg->AddEntry(overlayGraphs[i], overlayLabel[i], "l");
        if (std::isfinite(overlayThrT[i])) {
            TMarker* m = new TMarker(overlayThrT[i], -kThrMV, kFullCircle);
            m->SetMarkerColor(overlayColor[i]); m->SetMarkerSize(1.0);
            m->Draw();
            if (std::isfinite(overlayTts[i])) {
                double halfWidthSamples = overlayTts[i] / NS_PER_SAMPLE;   // TTS(sigma) in samples
                TLine* ttsBar = new TLine(overlayThrT[i] - halfWidthSamples, -kThrMV,
                                          overlayThrT[i] + halfWidthSamples, -kThrMV);
                ttsBar->SetLineColor(overlayColor[i]); ttsBar->SetLineWidth(2);
                ttsBar->Draw();
            }
        }
    }
    leg->Draw();
 //   TLatex note; note.SetNDC(); note.SetTextFont(132); note.SetTextSize(0.026); note.SetTextColor(kGray + 2);
 //   note.DrawLatex(0.10, 0.02, "dot = Thr(1.5mV) crossing on mean waveform;  bar width = #pm TTS(#sigma)");
    c2->Print(outPath, "pdf");   // page 2

    // Page 3: threshold comparison (1/2/3mV) -- how much the crossing time
    // and TTS themselves depend on the choice of analysis threshold, across
    // the same 9 HV points.
    static const double kThrCompareList[3] = {1.0, 2.0, 3.0};
    static const int    kThrCompareColor[3] = {kBlue + 1, kGreen + 2, kRed + 1};
    TCanvas* c3 = new TCanvas(Form("cThrCmp_ch%d", ch), "thr compare", 1200, 800);
    c3->Divide(1, 2);

    c3->cd(1); gPad->SetGrid(); gPad->SetLeftMargin(0.10);
    TH1F* fr1 = (TH1F*)gPad->DrawFrame(-0.5, 1175, kPoints.size() - 0.5, 1215);
    fr1->SetTitle(Form("%s -- Thr crossing time vs HV, by threshold;HV point;Thr crossing [ns]", chLabel));
    for (int i = 0; i < (int)kPoints.size(); ++i)
        fr1->GetXaxis()->SetBinLabel(fr1->GetXaxis()->FindBin(i), kPoints[i].label);
    fr1->GetXaxis()->LabelsOption("h");
    TLegend* leg1 = new TLegend(0.75, 0.72, 0.94, 0.90);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextFont(132); leg1->SetTextSize(0.03);
    for (int t = 0; t < 3; ++t) {
        TGraph* g = new TGraph((int)thrCmpHV[t].size(), thrCmpHV[t].data(), thrCmpT[t].data());
        g->SetLineColor(kThrCompareColor[t]); g->SetMarkerColor(kThrCompareColor[t]);
        g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetLineWidth(2);
        g->Draw("LP SAME");
        leg1->AddEntry(g, Form("%.0f mV", kThrCompareList[t]), "lp");
    }
    leg1->Draw();

    c3->cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.10);
    TH1F* fr2 = (TH1F*)gPad->DrawFrame(-0.5, 0, kPoints.size() - 0.5, 3.0);
    fr2->SetTitle(";HV point;TTS(#sigma) [ns]");
    for (int i = 0; i < (int)kPoints.size(); ++i)
        fr2->GetXaxis()->SetBinLabel(fr2->GetXaxis()->FindBin(i), kPoints[i].label);
    fr2->GetXaxis()->LabelsOption("h");
    for (int t = 0; t < 3; ++t) {
        TGraph* g = new TGraph((int)thrCmpHV[t].size(), thrCmpHV[t].data(), thrCmpTTS[t].data());
        g->SetLineColor(kThrCompareColor[t]); g->SetMarkerColor(kThrCompareColor[t]);
        g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetLineWidth(2);
        g->Draw("LP SAME");
    }
    c3->Print(outPath + ")", "pdf");   // page 3, closes the file
    delete c2; delete c3;

    std::cout << "[INFO] Saved: " << outPath << std::endl;
}

void Draw_PulseShape_byBfield() {
    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString ntuplePath = outDir + "PulseShape_MeanWaveform_Bfield.root";
    TFile* fOut = new TFile(ntuplePath, "RECREATE");
    TTree* avgTree = new TTree("avg_waveform", "Mean pulse shape, one row per (HV point, channel, sample)");
    char labelBuf[32] = "";
    Int_t runNum = 0, chOut = 0, sampleOut = 0;
    Double_t voltOut = 0;
    avgTree->Branch("label", labelBuf, "label/C");
    avgTree->Branch("run", &runNum, "run/I");
    avgTree->Branch("ch", &chOut, "ch/I");
    avgTree->Branch("sample", &sampleOut, "sample/I");
    avgTree->Branch("voltage_mV", &voltOut, "voltage_mV/D");

    DrawPulseGrid(1, "Rot1 (EM6400)", avgTree, labelBuf, runNum, chOut, sampleOut, voltOut);
    DrawPulseGrid(2, "Rot2 (EL5150)", avgTree, labelBuf, runNum, chOut, sampleOut, voltOut);

    fOut->cd();
    avgTree->Write();
    fOut->Close();
    std::cout << "[INFO] Saved: " << ntuplePath << std::endl;
}
