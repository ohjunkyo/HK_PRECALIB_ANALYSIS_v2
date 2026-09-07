#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TLegend.h>
#include <TParameter.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <set>
#include <algorithm>
#include <TNamed.h>

void ApplyYMargin(TMultiGraph* mg, double frac = 0.12, double fracTop = 0.32) {/*{{{*/
    if (!mg || !mg->GetListOfGraphs()) return;
    double lo = 1e30, hi = -1e30;
    for (TObject* o : *mg->GetListOfGraphs()) {
        TGraphErrors* g = dynamic_cast<TGraphErrors*>(o);
        if (!g) continue;
        for (int k = 0; k < g->GetN(); ++k) {
            double x, y; g->GetPoint(k, x, y);
            double ey = g->GetErrorY(k); if (ey < 0) ey = 0;
            lo = std::min(lo, y - ey); hi = std::max(hi, y + ey);
        }
    }
    if (lo > hi) return;
    double range = hi - lo;
    double mBottom = range * frac;
    double mTop = range * fracTop;
    if (mBottom <= 0) { mBottom = std::abs(hi) * 0.1 + 1e-6; mTop = mBottom * (fracTop / frac); }
    mg->SetMinimum(lo - mBottom); mg->SetMaximum(hi + mTop);
}
/*}}}*/

void DrawErrorBands(TMultiGraph* mg) {/*{{{*/
    if (!mg || !mg->GetListOfGraphs()) return;
    for (TObject* o : *mg->GetListOfGraphs()) {
        TGraphErrors* g = dynamic_cast<TGraphErrors*>(o);
        if (!g) continue;
        g->Sort();   
                     
        g->SetFillColorAlpha(g->GetLineColor(), 0.18);
        g->Draw("3 same");
    }
    mg->Draw("P same");
}
/*}}}*/
// ==========================================
const bool drawMonitor = false;
const bool drawRaw = false;
const bool drawCorrected = true;   // "Count Corr" = PHC + TimingCut (dark-subtracted) QE, shown alongside Poisson
const bool drawPoissonRaw= false;
const bool drawPoisson = false;
// ==========================================
void ApplyGlobalStyle() {/*{{{*/
    gStyle->SetTextFont(132);
    gStyle->SetLabelFont(132, "xyz");
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetTitleFont(22, ""); 
    gStyle->SetLegendFont(132);
    gStyle->SetLabelSize(0.06, "xyz");
    gStyle->SetTitleSize(0.06, "xyz");
    gStyle->SetTitleSize(0.08, "t");
    gStyle->SetTitleOffset(0.95, "x");
    gStyle->SetTitleOffset(1.10, "y");
    gStyle->SetFrameLineWidth(2);
    gStyle->SetOptTitle(1);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
}/*}}}*/
#include "angle_convert.h"
void NormalizeData(const TString& label, const std::vector<double>& angles, std::vector<double>& values, std::vector<double>& errors) {/*{{{*/
    if (values.empty() || angles.size() != values.size()) return;

    double sum = 0;
    double sumErrSq = 0;   // for the standard error of centralAvg, treating the
                            // |tilt|<48deg points as independent measurements
    int count = 0;
    for (size_t i = 0; i < values.size(); ++i) {
        if (std::abs(angles[i]) < 48.0) { sum += values[i]; sumErrSq += errors[i] * errors[i]; count++; }
    }
    if (count == 0) return;

    double centralAvg = sum / count;
    // Standard error of the mean: centralAvg = (1/N) sum v_i, so
    // Var(centralAvg) = (1/N^2) sum e_i^2 for independent v_i.
    double centralAvgErr = std::sqrt(sumErrSq) / count;
    if (centralAvg != 0) {
        double relAvgErr = centralAvgErr / centralAvg;
        for (size_t i = 0; i < values.size(); ++i) {
            double relPointErr = (values[i] != 0) ? errors[i] / values[i] : 0;
            values[i] /= centralAvg;
            errors[i] = values[i] * std::sqrt(relPointErr * relPointErr + relAvgErr * relAvgErr);
        }
    }
}/*}}}*/

/*{{{*/
void MonitorNormalize(const std::vector<double>& tiltRaw, const std::vector<double>& values, const std::vector<double>& errors,
                      const std::vector<double>& monTiltRaw, const std::vector<double>& monValues, const std::vector<double>& monErrors,
                      std::vector<double>& outValues, std::vector<double>& outErrors) {
    outValues.assign(tiltRaw.size(), 0); outErrors.assign(tiltRaw.size(), 0);
    for (size_t i = 0; i < tiltRaw.size(); ++i) {
        for (size_t j = 0; j < monTiltRaw.size(); ++j) {
            if (std::abs(tiltRaw[i] - monTiltRaw[j]) < 0.5) {
                if (monValues[j] != 0) {
                    outValues[i] = values[i] / monValues[j];
                    double relErrTest = (values[i] != 0) ? errors[i] / values[i] : 0;
                    double relErrMon = (monValues[j] != 0) ? monErrors[j] / monValues[j] : 0;
                    outErrors[i] = outValues[i] * std::sqrt(relErrTest * relErrTest + relErrMon * relErrMon);
                }
                break;
            }
        }
    }
}/*}}}*/

// Naming convention used throughout this file
// -------------------------------------------
//   ...Abs      : absolute value, NOT divided by the |tilt|<48deg average
//                 (i.e. NormalizeData() was never applied to it)
//   ...Norm     : angle-normalized -- divided by the |tilt|<48deg average,
//                 so the central region sits at 1.0
//   ...DarkCorr : dark-hit-subtracted  (ntuple branch "relativeQE")
//                 no DarkCorr in the name  ==  NOT dark-subtracted
//                 (ntuple branch "relativeQE_raw")
//   monNorm...  : additionally divided by the monitor PMT (CH0) at the same
//                 raw tilt, correcting shot-to-shot laser-intensity drift
//
// NOTE: the old code used the prefix "raw_" for BOTH "not angle-normalized"
// (members) and "not dark-subtracted" (ntuple branches). Those two unrelated
// meanings are now spelled Abs and <absence of DarkCorr> respectively.
struct PMTData {/*{{{*/
    std::vector<TString> fileName;
    std::vector<double> rotAngle, tiltRaw, tiltFlipped, tiltHamamatsu;
    std::vector<double> speNorm, speNormErr, speAbs, speAbsErr;                          // SPE charge (gain)
    std::vector<double> chargeRes, chargeResErr;                                          // SPE_Sigma/SPE_Gain [%]
    std::vector<double> qeNorm, qeNormErr, qeAbs, qeAbsErr;                              // PHC-counting QE, no dark subtraction
    std::vector<double> qeDarkCorrNorm, qeDarkCorrNormErr, qeDarkCorrAbs, qeDarkCorrAbsErr;  // PHC-counting QE, dark-subtracted
    std::vector<double> darkRate, darkRateErr, darkCount, darkCountErr;
    std::vector<double> poissonQeAbs, poissonQeAbsErr, poissonQeNorm, poissonQeNormErr;  // QE from the Poisson-mu fit
    std::vector<double> tts, ttsErr;
    std::vector<double> fwhm, fwhmErr, sigma, sigmaErr;   // same exGaus fit as TTS
    // -- 장지승 박사님 (2026-09-03, via chat): "FWHM으로 같은 방식으로 봐보세"
    std::vector<double> pedMean, pedSigma;
    // Monitor(CH0)-normalized: Test_absolute / Monitor_absolute at the matching raw
    // tilt, computed BEFORE NormalizeData() mutates the source vectors in place --
    // corrects for shot-to-shot laser-intensity drift using the non-rotating
    // monitor PMT as reference, independent of the existing |tilt|<48deg centering.
    std::vector<double> monNormQe, monNormQeErr, monNormQeDarkCorr, monNormQeDarkCorrErr;
    std::vector<double> monNormPoissonQeAbs, monNormPoissonQeAbsErr, monNormPoissonQe, monNormPoissonQeErr;
    std::vector<double> monNormSpe, monNormSpeErr;

    TGraphErrors *grSpeNorm = nullptr, *grSpeAbs = nullptr;
    TGraphErrors *grChargeRes = nullptr;
    TGraphErrors *grQeNorm = nullptr, *grQeAbs = nullptr;
    TGraphErrors *grQeDarkCorrNorm = nullptr, *grQeDarkCorrAbs = nullptr;
    TGraphErrors *grDarkRate = nullptr, *grDarkCount = nullptr;
    TGraphErrors *grTts = nullptr;
    TGraphErrors *grFwhm = nullptr, *grSigma = nullptr;
    TGraphErrors *grPedStability = nullptr; // Y = Pedestal Mean [mV], error = Pedestal Sigma [mV]
    TGraphErrors *grPoissonQeNorm = nullptr, *grPoissonQeAbs = nullptr;
    TGraphErrors *grRawStageAngle = nullptr;  // y = raw stage tilt at each point
    TGraphErrors *grMonNormQe = nullptr, *grMonNormQeDarkCorr = nullptr;
    TGraphErrors *grMonNormPoissonQeAbs = nullptr, *grMonNormPoissonQe = nullptr, *grMonNormSpe = nullptr;
};/*}}}*/

struct PMTUnit {/*{{{*/
    TString serial; char direction; TString hv;
    int channel;                 // digitizer channel: 0 = monitor PMT, 1/2 = PMTs under test
    int xScanRotAngle, yScanRotAngle;  // stage rotation that selects the X- / Y-axis scan
    PMTData X, Y;
};/*}}}*/

void FillData(PMTUnit& p, TString dirPath, TString tag, int run_start, int run_end, bool isX, std::set<TString>& processed) {/*{{{*/
    auto& d = isX ? p.X : p.Y;

    TSystemDirectory sysdir(dirPath, dirPath);
    TList *files = sysdir.GetListOfFiles();
    if (!files) return;

    for (int i = 0; i < files->GetSize(); ++i) {
        TSystemFile *file = (TSystemFile*)files->At(i);
        TString fname = file->GetName();

        if (file->IsDirectory() || !fname.EndsWith(".root") || !fname.Contains("result")) continue;
        if (!fname.Contains(tag)) continue;

        TObjArray* tokens = fname.Tokenize("_.");
        int run_num = -1;
        if (tokens->GetEntries() >= 2) {
            TString numStr = ((TObjString*)tokens->At(tokens->GetEntries() - 2))->GetString();
            run_num = numStr.Atoi();
        }
        delete tokens;

        if (run_num < run_start || run_num > run_end) continue;

        TString path = dirPath + fname;
        TFile* fResult = TFile::Open(path, "READ");
        if (!fResult || fResult->IsZombie()) continue;

        TTree* info = (TTree*)fResult->Get("RunInfo");
        if (!info) { fResult->Close(); continue; }

        int rawRotateAngle2=0, rawTiltAngle2=0, rawRotateAngle3=0, rawTiltAngle3=0;
        char runModeBuf[20]="";

        info->SetBranchStatus("*", 0);
        info->SetBranchStatus("RawRotateAngle2", 1); info->SetBranchStatus("RawTiltAngle2", 1);
        info->SetBranchStatus("RawRotateAngle3", 1); info->SetBranchStatus("RawTiltAngle3", 1);
        if (info->GetBranch("RunMode")) info->SetBranchStatus("RunMode", 1);

        info->SetBranchAddress("RawRotateAngle2", &rawRotateAngle2); info->SetBranchAddress("RawTiltAngle2", &rawTiltAngle2);
        info->SetBranchAddress("RawRotateAngle3", &rawRotateAngle3); info->SetBranchAddress("RawTiltAngle3", &rawTiltAngle3);
        if (info->GetBranch("RunMode")) info->SetBranchAddress("RunMode", runModeBuf);
        info->GetEntry(0);

        if (TString(runModeBuf) == "Dark") { fResult->Close(); continue; }

        int currentRotAngle = (p.channel == 1 || p.channel == 0) ? rawRotateAngle2 : rawRotateAngle3;
        // For ch0 this is NOT the monitor's own orientation -- the monitor PMT
        // does not rotate at all. RunInfo only stores stage angles for Devices
        // 2/3, so ch0 borrows Device 2's rotation purely as a key for "was this
        // run taken during the X-axis block or the Y-axis block".
        //

        int targetRotAngle = isX ? p.xScanRotAngle : p.yScanRotAngle;
        if (currentRotAngle != targetRotAngle) { fResult->Close(); continue; }

        if (p.channel == 0 && processed.count(path)) { fResult->Close(); continue; }
        if (p.channel == 0) processed.insert(path);

        double tiltRawVal = (p.channel == 1 || p.channel == 0) ? rawTiltAngle2 : rawTiltAngle3;

        // Incidence angle (and its sign) comes from angle_convert.h

        const char* axisLabel = "?-axis";
        double tiltHamamatsuVal =
            GetHamamatsuAngle(p.direction, tiltRawVal, (double)currentRotAngle, axisLabel);
        // Raw stage magnitude carrying the corrected sign (kept for the
        // per-point table below, which reports raw vs flipped vs Hamamatsu).
        double tiltFlippedVal = (tiltHamamatsuVal < 0 ? -1.0 : 1.0) * std::abs(tiltRawVal);

        double darkRateVal = 0, darkRateErrVal = 0;
        long long darkCountVal = 0;
        double pedMeanVal = 0, pedSigmaVal = 0;
        TString prd = path;
        prd.ReplaceAll("FinalResult", "production"); prd.ReplaceAll("result", "prd");
        TFile* fPrd = TFile::Open(prd, "READ");
        if (fPrd) {
            TParameter<double>* pRate = (TParameter<double>*)fPrd->Get(Form("NoiseCountRate_ch%d", p.channel));
            TParameter<Long64_t>* pCount = (TParameter<Long64_t>*)fPrd->Get(Form("NoiseCount_ch%d", p.channel));
            if (pRate && pCount) {
                darkRateVal = pRate->GetVal(); darkCountVal = pCount->GetVal();
                darkRateErrVal = (darkCountVal > 0) ? darkRateVal / sqrt((double)darkCountVal) : 0;
            }
            // Baseline stability (core-window mean/sigma, mV) -- see prod_ntp_v7.C's
            // cPedStab canvas for how these are computed.
            TParameter<double>* pPedMean = (TParameter<double>*)fPrd->Get(Form("PedMean_mV_ch%d", p.channel));
            TParameter<double>* pPedSigma = (TParameter<double>*)fPrd->Get(Form("PedSigma_mV_ch%d", p.channel));
            if (pPedMean) pedMeanVal = pPedMean->GetVal();
            if (pPedSigma) pedSigmaVal = pPedSigma->GetVal();
            fPrd->Close();
        }

        TTree* tr = (TTree*)fResult->Get(Form("tree_ch%d", p.channel));
        if (tr) {
            double speMeanVal, speMeanErrVal, qeDarkCorrVal, qeDarkCorrErrVal, qeAbsVal, qeAbsErrVal, ttsVal = 0, ttsErrVal = 0;
            tr->SetBranchAddress("spe_mean", &speMeanVal); tr->SetBranchAddress("spe_mean_error", &speMeanErrVal);
            tr->SetBranchAddress("relativeQE", &qeDarkCorrVal); tr->SetBranchAddress("relativeQE_err", &qeDarkCorrErrVal);
            if (tr->GetBranch("rms_exG")) tr->SetBranchAddress("rms_exG", &ttsVal); // Read TTS (in Sample)
            if (tr->GetBranch("rms_exG_err")) tr->SetBranchAddress("rms_exG_err", &ttsErrVal); // relative error, or -1 if the fit was invalid
            double fwhmVal = 0, fwhmErrVal = 0, sigmaVal = 0, sigmaErrVal = 0;
            if (tr->GetBranch("fwhm_exG")) tr->SetBranchAddress("fwhm_exG", &fwhmVal);
            if (tr->GetBranch("fwhm_exG_err")) tr->SetBranchAddress("fwhm_exG_err", &fwhmErrVal); // relative error, same convention as rms_exG_err
            if (tr->GetBranch("sigma_exG")) tr->SetBranchAddress("sigma_exG", &sigmaVal);
            if (tr->GetBranch("sigma_exG_err")) tr->SetBranchAddress("sigma_exG_err", &sigmaErrVal); // ABSOLUTE error (direct fit GetParError), not relative like TTS/FWHM

            double chargeResVal = 0, chargeResErrVal = 0;
            if (tr->GetBranch("charge_resolution")) {
                tr->SetBranchAddress("charge_resolution", &chargeResVal);
                tr->SetBranchAddress("charge_resolution_err", &chargeResErrVal);
            }
            bool hasRawQe = false;
            if (tr->GetBranch("relativeQE_raw")) {
                tr->SetBranchAddress("relativeQE_raw", &qeAbsVal); tr->SetBranchAddress("relativeQE_raw_err", &qeAbsErrVal);
                hasRawQe = true;
            }

            double poissonQeVal = 0;
            double poissonQeAbsVal = 0;
            bool hasPoisson = false;

            if (tr->GetBranch("poisson_qe")) {
                tr->SetBranchAddress("poisson_qe", &poissonQeVal);
                if (tr->GetBranch("poisson_qe_raw")) {
                    tr->SetBranchAddress("poisson_qe_raw", &poissonQeAbsVal);
                } else {
                    poissonQeAbsVal = poissonQeVal;
                }
                hasPoisson = true;
            }

            tr->GetEntry(0);
            d.fileName.push_back(gSystem->BaseName(path)); d.tiltRaw.push_back(tiltRawVal); d.tiltFlipped.push_back(tiltFlippedVal); d.tiltHamamatsu.push_back(tiltHamamatsuVal);
            d.rotAngle.push_back(currentRotAngle);
            d.speNorm.push_back(speMeanVal); d.speNormErr.push_back(speMeanErrVal); d.speAbs.push_back(speMeanVal); d.speAbsErr.push_back(speMeanErrVal);
            d.chargeRes.push_back(chargeResVal); d.chargeResErr.push_back(chargeResErrVal);

            d.darkRate.push_back(darkRateVal); d.darkRateErr.push_back(darkRateErrVal);
            d.darkCount.push_back((double)darkCountVal); d.darkCountErr.push_back(darkCountVal > 0 ? sqrt((double)darkCountVal) : 0);
            d.pedMean.push_back(pedMeanVal);
            d.pedSigma.push_back(pedSigmaVal);
            d.tts.push_back(ttsVal * 2.0);
            // ttsErrVal read above is a RELATIVE error (sigma_TTS/TTS) from the
            // ExGaussian fit's sigma/tau covariance; -1 means "fit invalid", same
            // sentinel convention as ttsVal itself.
            d.ttsErr.push_back((ttsVal > 0 && ttsErrVal >= 0) ? ttsErrVal * ttsVal * 2.0 : 0.0);
            d.fwhm.push_back(fwhmVal * 2.0);
            d.fwhmErr.push_back((fwhmVal > 0 && fwhmErrVal >= 0) ? fwhmErrVal * fwhmVal * 2.0 : 0.0);
            d.sigma.push_back(sigmaVal * 2.0);
            d.sigmaErr.push_back((sigmaVal > 0 && sigmaErrVal >= 0) ? sigmaErrVal * 2.0 : 0.0);   // absolute error, only needs the ns conversion

            if (hasRawQe) {
                d.qeAbs.push_back(qeAbsVal); d.qeAbsErr.push_back(qeAbsErrVal); d.qeNorm.push_back(qeAbsVal); d.qeNormErr.push_back(qeAbsErrVal);
                d.qeDarkCorrAbs.push_back(qeDarkCorrVal); d.qeDarkCorrAbsErr.push_back(qeDarkCorrErrVal); d.qeDarkCorrNorm.push_back(qeDarkCorrVal); d.qeDarkCorrNormErr.push_back(qeDarkCorrErrVal);
            } else {
                d.qeAbs.push_back(qeDarkCorrVal); d.qeAbsErr.push_back(qeDarkCorrErrVal); d.qeNorm.push_back(qeDarkCorrVal); d.qeNormErr.push_back(qeDarkCorrErrVal);
                double legacyDarkCorrVal = qeDarkCorrVal - (darkRateVal * 1e-5);
                double legacyDarkCorrErrVal = sqrt(qeDarkCorrErrVal*qeDarkCorrErrVal + (darkRateErrVal * 1e-5)*(darkRateErrVal * 1e-5));
                d.qeDarkCorrAbs.push_back(legacyDarkCorrVal); d.qeDarkCorrAbsErr.push_back(legacyDarkCorrErrVal); d.qeDarkCorrNorm.push_back(legacyDarkCorrVal); d.qeDarkCorrNormErr.push_back(legacyDarkCorrErrVal);
            }
            if (hasPoisson) {
                d.poissonQeNorm.push_back(poissonQeVal);
                d.poissonQeAbs.push_back(poissonQeAbsVal);

                double poissonQeErrVal = (qeDarkCorrVal > 0) ? (qeDarkCorrErrVal / qeDarkCorrVal) * poissonQeVal : 0;
                // Raw 값에 대한 에러는 Raw Counting의 에러 비율을 연동하여 계산
                double poissonQeAbsErrVal = (qeAbsVal > 0) ? (qeAbsErrVal / qeAbsVal) * poissonQeAbsVal : 0;

                d.poissonQeNormErr.push_back(poissonQeErrVal);
                d.poissonQeAbsErr.push_back(poissonQeAbsErrVal);
            } else {
                d.poissonQeNorm.push_back(0);
                d.poissonQeNormErr.push_back(0);
                d.poissonQeAbs.push_back(0);
                d.poissonQeAbsErr.push_back(0);
            }
        }
        fResult->Close();
    }
}/*}}}*/

void BuildScanDistributionReport(TString dirPath, TString tag, int run_start, int run_end, TString outSuffix) {/*{{{*/
    TSystemDirectory sysdir(dirPath, dirPath);
    TList *files = sysdir.GetListOfFiles();
    if (!files) return;

    std::vector<std::pair<int, TString>> runFiles;
    for (int i = 0; i < files->GetSize(); ++i) {
        TSystemFile *file = (TSystemFile*)files->At(i);
        TString fname = file->GetName();
        if (file->IsDirectory() || !fname.EndsWith(".root") || !fname.Contains("result")) continue;
        if (!fname.Contains(tag)) continue;

        TObjArray* tokens = fname.Tokenize("_.");
        int run_num = -1;
        if (tokens->GetEntries() >= 2) {
            TString numStr = ((TObjString*)tokens->At(tokens->GetEntries() - 2))->GetString();
            run_num = numStr.Atoi();
        }
        delete tokens;
        if (run_num < run_start || run_num > run_end) continue;

        runFiles.push_back({run_num, dirPath + fname});
    }
    if (runFiles.empty()) {
        std::cout << "[WARNING] BuildScanDistributionReport: no result files found for " << tag << std::endl;
        return;
    }
    std::sort(runFiles.begin(), runFiles.end());

    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + "ScanReport_" + outSuffix + ".pdf";

    for (size_t idx = 0; idx < runFiles.size(); ++idx) {
        int run_num = runFiles[idx].first;
        TFile* f = TFile::Open(runFiles[idx].second, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }

        char b_RunMode[20] = "";
        int b_RawRotateAngle2 = 0, b_RawTiltAngle2 = 0;
        TTree* info = (TTree*)f->Get("RunInfo");
        if (info) {
            if (info->GetBranch("RunMode")) info->SetBranchAddress("RunMode", b_RunMode);
            if (info->GetBranch("RawRotateAngle2")) info->SetBranchAddress("RawRotateAngle2", &b_RawRotateAngle2);
            if (info->GetBranch("RawTiltAngle2")) info->SetBranchAddress("RawTiltAngle2", &b_RawTiltAngle2);
            info->GetEntry(0);
        }

        std::vector<int> activeChannels;
        for (int ch = 0; ch < 8; ++ch) {
            if (f->Get(Form("Diff_ch%d", ch)) || f->Get(Form("Pico_ch%d", ch))) activeChannels.push_back(ch);
        }
        if (activeChannels.empty()) { f->Close(); continue; }

        int nCh = (int)activeChannels.size();
        TCanvas *c = new TCanvas(Form("cScanReport_run%d", run_num), "Time + Charge", 700 * nCh, 900);
        c->Divide(nCh, 2);

        for (int k = 0; k < nCh; ++k) {
            int ch = activeChannels[k];
            c->cd(k + 1);
            TH1D* hTime = (TH1D*)f->Get(Form("Diff_ch%d", ch));
            if (hTime) { hTime->SetLineWidth(2); hTime->SetLineColor(kBlack); hTime->DrawCopy("Hist"); }

            c->cd(nCh + k + 1);
            gPad->SetLogy();
            TH1D* hCharge = (TH1D*)f->Get(Form("Pico_ch%d", ch));
            if (hCharge) { hCharge->SetLineWidth(2); hCharge->SetLineColor(kBlack); hCharge->DrawCopy("Hist"); }
        }

        TString pageTitle = Form("Run %03d  (%s)  Rot(R,T)=(%d,%d)", run_num, b_RunMode, b_RawRotateAngle2, b_RawTiltAngle2);
        TString printPath = outPath;
        if (runFiles.size() > 1) {
            if (idx == 0) printPath += "(";
            else if (idx == runFiles.size() - 1) printPath += ")";
        }
        c->Print(printPath, Form("pdf Title:%s", pageTitle.Data()));

        delete c;
        f->Close();
    }
    std::cout << "[INFO] Scan distribution report written: " << outPath << std::endl;
}/*}}}*/

static bool kUseRawStageAxis = false;
static const char* AngleAxisTitle() {
    return kUseRawStageAxis ? "Raw Stage Angle [degree]" : "Position Angle [degree]";
}

void Draw_Uniformity_Norm_v7(TString tag = "20260331", int run_start = 0, int run_end = 45,
        TString chsel = "0,1,2", TString xaxis = "hamamatsu") {
    ApplyGlobalStyle();

    // "hamamatsu" (default) = incidence angle in the Hamamatsu convention, as
    // every report before this used. "rawstage" = the stage tilt the scan
    // actually commanded, matching the DAQ GUI's live view / scan matrix.
    kUseRawStageAxis = xaxis.Contains("raw", TString::kIgnoreCase);
    std::cout << "[INFO] X axis: " << AngleAxisTitle()
        << (kUseRawStageAxis ? "  (raw stage angle mode)" : "") << std::endl;

    // Parse the channel-selection string ("0,1,2") into a list of channel indices.
    // Only these channels are drawn in the combined plots. Defaults to all 3.
    std::vector<int> chs;
    {
        TObjArray* toks = chsel.Tokenize(",");
        for (int i = 0; i < toks->GetEntries(); ++i) {
            int c = ((TObjString*)toks->At(i))->GetString().Atoi();
            if (c >= 0 && c < 3 &&
                    std::find(chs.begin(), chs.end(), c) == chs.end()) chs.push_back(c);
        }
        delete toks;
        if (chs.empty()) chs = {0, 1, 2};
        std::sort(chs.begin(), chs.end());
    }
    std::cout << "[INFO] Drawing channels: ";
    for (int c : chs) std::cout << "CH" << c << " ";
    std::cout << std::endl;

    TString dir = "./Data/FinalResult/";             
    std::set<TString> procX, procY;
    TString info_Bfield = "FullCurrent";

    TString info_Expert = "Unknown";
    TString info_Shifter = "Unknown";
    int info_Laser = 0;
    TString sn1_str = "CH0_PMT", sn2_str = "CH1_PMT", sn3_str = "CH2_PMT";
    TString hv1_str = "0", hv2_str = "0", hv3_str = "0";
    // Cable direction (A~H) -- was hardcoded to 'B' for both PMT slots below,
    // which silently mis-scanned/mis-analyzed any unit wired with a different
    // direction (e.g. 'H'). Read it from RunInfo like SN/HV already are.
    char dir2_c = 'B', dir3_c = 'B';

    TSystemDirectory sysdir(dir, dir);
    TList *fileList = sysdir.GetListOfFiles();
    if (fileList) {
        for (int i = 0; i < fileList->GetSize(); ++i) {
            TSystemFile *f = (TSystemFile*)fileList->At(i);
            TString fname = f->GetName();

            if (fname.Contains(tag) && fname.EndsWith(".root") && fname.Contains("result")) {
                TFile *tmpFile = TFile::Open(dir + fname, "READ");
                if (tmpFile && !tmpFile->IsZombie()) {
                    TTree* info = (TTree*)tmpFile->Get("RunInfo");
                    if (info) {
                        char expert_name[100] = "";
                        char shifter_name[100] = "";
                        int laser_current = 0;
                        char sn1[100] = "", sn2[100] = "", sn3[100] = "";
                        int hv1 = 0, hv2 = 0, hv3 = 0;
                        char dir2[10] = "", dir3[10] = "";

                        if (info->GetBranch("Expert")) info->SetBranchAddress("Expert", expert_name);
                        if (info->GetBranch("Shifter")) info->SetBranchAddress("Shifter", shifter_name);
                        if (info->GetBranch("Laser_mA")) info->SetBranchAddress("Laser_mA", &laser_current);
                        if (info->GetBranch("SN1")) info->SetBranchAddress("SN1", sn1);
                        if (info->GetBranch("SN2")) info->SetBranchAddress("SN2", sn2);
                        if (info->GetBranch("SN3")) info->SetBranchAddress("SN3", sn3);
                        if (info->GetBranch("HV1")) info->SetBranchAddress("HV1", &hv1);
                        if (info->GetBranch("HV2")) info->SetBranchAddress("HV2", &hv2);
                        if (info->GetBranch("HV3")) info->SetBranchAddress("HV3", &hv3);
                        if (info->GetBranch("Direction2")) info->SetBranchAddress("Direction2", dir2);
                        if (info->GetBranch("Direction3")) info->SetBranchAddress("Direction3", dir3);

                        info->GetEntry(0);

                        if (TString(expert_name) != "") info_Expert = expert_name;
                        if (TString(shifter_name) != "") info_Shifter = shifter_name;
                        info_Laser = laser_current;
                        if (TString(sn1) != "") sn1_str = sn1;
                        if (TString(sn2) != "") sn2_str = sn2;
                        if (TString(sn3) != "") sn3_str = sn3;
                        hv1_str = Form("%d", hv1);
                        hv2_str = Form("%d", hv2);
                        hv3_str = Form("%d", hv3);
                        if (dir2[0] != '\0') dir2_c = toupper(dir2[0]);
                        if (dir3[0] != '\0') dir3_c = toupper(dir3[0]);
                    }
                    tmpFile->Close();
                    break; 
                }
            }
        }
    }
    // =========================================================================

    // xScanRotAngle/yScanRotAngle derived from each PMT's actual cable
    // direction (angle_convert.h) instead of assuming everything is wired
    // 'B' -- a unit wired 'H' (or any other letter) now gets the correct
    // rotation targets instead of being silently mis-scanned.
    int x_rot2, y_rot2, x_rot3, y_rot3;
    GetXYRotForDirection(dir2_c, x_rot2, y_rot2);
    GetXYRotForDirection(dir3_c, x_rot3, y_rot3);

    // ch0 is the monitor PMT: it never rotates, so it has no scan rotation of
    // its own. It takes Device 2's targets because RunInfo records only Devices
    // 2/3, making Device 2's angle the de-facto "which axis block" tag on every
    // run (see FillData). For a 'B'-wired Device 2 these are 45/135 -- identical
    // to the value ch0 hardcoded before, so older data is unaffected.
    PMTUnit pmt[3] = {
        // ch0's `direction` is Device 2's cable, NOT the monitor's own ('A'):
        // every angle ch0 is plotted against is Device 2's stage angle (see
        // FillData), so the sign convention has to be read with Device 2's
        // cable too or the monitor's points land mirrored against the very
        // scan they annotate.
        {sn1_str, dir2_c, hv1_str, 0, x_rot2, y_rot2},
        {sn2_str, dir2_c, hv2_str, 1, x_rot2, y_rot2},
        {sn3_str, dir3_c, hv3_str, 2, x_rot3, y_rot3}
    };

    TString file_suffix = Form("%s_%d_%d", tag.Data(), run_start, run_end);
    // Every result canvas below goes into ONE multi-page PDF (same convention
    // as BuildScanDistributionReport's own PDF and Draw_Overlay_Uniformity_v7.C):
    // first Print() call opens with "(", last one closes with ")".
    TString uniPdfPath = Form("./Data/image/Uniformity/Uniformity_Report_%s.pdf", file_suffix.Data());

    BuildScanDistributionReport(dir, tag, run_start, run_end, file_suffix);

    for (int i = 0; i < 3; ++i) {
        FillData(pmt[i], dir, tag,  run_start, run_end, true, procX);
        FillData(pmt[i], dir, tag,  run_start, run_end, false, procY);
    }

    // Monitor(CH0)-normalization: Test_absolute / Monitor_absolute, matched by raw
    // tilt. Must run BEFORE NormalizeData() below (which mutates raw_pqe/pqe/q/cq
    // in place), so both sides of the ratio are still absolute values.
    for (int i = 1; i < 3; ++i) {
        auto monNorm = [&](PMTData& d, PMTData& mon) {
            MonitorNormalize(d.tiltRaw, d.qeAbs, d.qeAbsErr, mon.tiltRaw, mon.qeAbs, mon.qeAbsErr, d.monNormQe, d.monNormQeErr);
            MonitorNormalize(d.tiltRaw, d.qeDarkCorrAbs, d.qeDarkCorrAbsErr, mon.tiltRaw, mon.qeDarkCorrAbs, mon.qeDarkCorrAbsErr, d.monNormQeDarkCorr, d.monNormQeDarkCorrErr);
            MonitorNormalize(d.tiltRaw, d.poissonQeAbs, d.poissonQeAbsErr, mon.tiltRaw, mon.poissonQeAbs, mon.poissonQeAbsErr, d.monNormPoissonQeAbs, d.monNormPoissonQeAbsErr);
            MonitorNormalize(d.tiltRaw, d.poissonQeNorm, d.poissonQeNormErr, mon.tiltRaw, mon.poissonQeNorm, mon.poissonQeNormErr, d.monNormPoissonQe, d.monNormPoissonQeErr);
            MonitorNormalize(d.tiltRaw, d.speAbs, d.speAbsErr, mon.tiltRaw, mon.speAbs, mon.speAbsErr, d.monNormSpe, d.monNormSpeErr);
        };
        monNorm(pmt[i].X, pmt[0].X);
        monNorm(pmt[i].Y, pmt[0].Y);
    }

    for (int i = 0; i < 3; ++i) {
        NormalizeData(pmt[i].serial + " X-SPE", pmt[i].X.tiltRaw, pmt[i].X.speNorm, pmt[i].X.speNormErr);
        NormalizeData(pmt[i].serial + " Y-SPE", pmt[i].Y.tiltRaw, pmt[i].Y.speNorm, pmt[i].Y.speNormErr);
        NormalizeData(pmt[i].serial + " X-QE",  pmt[i].X.tiltRaw, pmt[i].X.qeNorm, pmt[i].X.qeNormErr);
        NormalizeData(pmt[i].serial + " Y-QE",  pmt[i].Y.tiltRaw, pmt[i].Y.qeNorm, pmt[i].Y.qeNormErr);
        NormalizeData(pmt[i].serial + " X-CorrQE", pmt[i].X.tiltRaw, pmt[i].X.qeDarkCorrNorm, pmt[i].X.qeDarkCorrNormErr);
        NormalizeData(pmt[i].serial + " Y-CorrQE", pmt[i].Y.tiltRaw, pmt[i].Y.qeDarkCorrNorm, pmt[i].Y.qeDarkCorrNormErr);
        NormalizeData(pmt[i].serial + " X-RawPoissonQE", pmt[i].X.tiltRaw, pmt[i].X.poissonQeAbs, pmt[i].X.poissonQeAbsErr);
        NormalizeData(pmt[i].serial + " Y-RawPoissonQE", pmt[i].Y.tiltRaw, pmt[i].Y.poissonQeAbs, pmt[i].Y.poissonQeAbsErr);
        NormalizeData(pmt[i].serial + " X-PoissonQE", pmt[i].X.tiltRaw, pmt[i].X.poissonQeNorm, pmt[i].X.poissonQeNormErr);
        NormalizeData(pmt[i].serial + " Y-PoissonQE", pmt[i].Y.tiltRaw, pmt[i].Y.poissonQeNorm, pmt[i].Y.poissonQeNormErr);
    }

    // Console Output
    for (int p = 0; p < 3; ++p) {
        if (std::find(chs.begin(), chs.end(), p) == chs.end()) continue;  // only selected channels
        std::cout << "\n" << std::string(90, '=') << "\n [TABLE] " << pmt[p].serial << "\n" << std::left
            << std::setw(15) << "Serial" << std::setw(8) << "Rot" << std::setw(8) << "Tilt" << std::setw(8) << "Flip" << std::setw(8) << "Ham" 
            << std::setw(12) << "SPE(Raw)" << std::setw(12) << "Rate(Hz)" << std::endl;

        auto print_row = [&](PMTData& d) {
            if (d.fileName.empty()) return;
            std::vector<size_t> indices(d.fileName.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](size_t i1, size_t i2) {
                    if (d.rotAngle[i1] != d.rotAngle[i2]) return d.rotAngle[i1] < d.rotAngle[i2];
                    return d.tiltRaw[i1] < d.tiltRaw[i2];
                    });
            for (size_t i : indices) {
                std::cout << std::left << std::setw(15) << pmt[p].serial 
                    << std::fixed << std::setprecision(1) << std::setw(8) << d.rotAngle[i] << std::setw(8) << d.tiltRaw[i] 
                    << std::setw(8) << d.tiltFlipped[i] << std::setprecision(3) << std::setw(8) << d.tiltHamamatsu[i] 
                    << std::setprecision(2) << std::setw(12) << d.speAbs[i] << std::setw(12) << d.darkRate[i] << std::endl;
            }
        };
        print_row(pmt[p].X); print_row(pmt[p].Y);
    }

    int clr[] = {kBlack, kRed, kBlue};
    auto styleGraph = [&](TGraphErrors* g, int pmt_idx, int method) {
        if (!g) return;
        int color;
        int marker = 20;

        if(pmt_idx == 0) color = kBlack;
        else if (pmt_idx == 1) color = kRed;
        else if (pmt_idx == 2) color = kBlue;

        //   if (pmt_idx == 1)      color = (method == 2 || method == 3) ? kRed : kBlue;
        //   else if (pmt_idx == 2) color = (method == 2 || method == 3) ? kRed : kBlue;
        //   else if (pmt_idx == 0) color = (method == 2 || method == 3) ? kRed : kBlue;

        // Each PMT gets its own marker family so CH0 (monitor) and CH2 are
        // distinguishable by shape, not only by colour. Within a family the
        // open/filled shape still encodes the QE method (Counting vs Poisson)
        // -- Counting (PHC + Timing cut) is now the PRIMARY QE method, so it
        // gets the FILLED (solid) marker; Poisson is the secondary/reference
        // method and gets the OPEN marker. (Previously the other way around,
        // from when Poisson was primary.)
        if (pmt_idx == 0) { // Monitor CH0
            if      (method == 0 || method == 1) marker = 20; // Counting: Filled Circle
            else if (method == 2 || method == 3) marker = 25; // Poisson : Open Square
            else marker = 20;
        }
        else if (pmt_idx == 1) { // CH1
            if      (method == 0 || method == 1) marker = 21; // Counting: Filled Square
            else if (method == 2 || method == 3) marker = 24; // Poisson : Open Circle
            else marker = 20;
        }
        else { // CH2
            if      (method == 0 || method == 1) marker = 33; // Counting: Filled Diamond
            else if (method == 2 || method == 3) marker = 26; // Poisson : Open Triangle
            else marker = 22; // Default
        }

        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->SetMarkerStyle(marker);
        g->SetMarkerSize(3.5); 
    };

    for (int i = 0; i < 3; ++i) {
        auto prep = [&](PMTData& d) {
            if(d.tiltHamamatsu.empty()) return;
            // X axis source. Default stays the Hamamatsu incidence angle (the
            // convention every previous report used); xaxis="rawstage" plots
            // against the stage angle actually commanded during the scan
            // instead, which is what the DAQ GUI's live view and the scan
            // matrix are keyed on -- so a point can be matched up between a
            // finished report and the run that produced it.
            const std::vector<double>& X = kUseRawStageAxis ? d.tiltRaw : d.tiltHamamatsu;
            d.grSpeNorm = new TGraphErrors(X.size(), &X[0], &d.speNorm[0], 0, &d.speNormErr[0]);
            d.grQeNorm = new TGraphErrors(X.size(), &X[0], &d.qeNorm[0], 0, &d.qeNormErr[0]);
            d.grDarkRate = new TGraphErrors(X.size(), &X[0], &d.darkRate[0], 0, &d.darkRateErr[0]);
            d.grDarkCount = new TGraphErrors(X.size(), &X[0], &d.darkCount[0], 0, &d.darkCountErr[0]);
            d.grSpeAbs = new TGraphErrors(X.size(), &X[0], &d.speAbs[0], 0, &d.speAbsErr[0]);
            d.grChargeRes = new TGraphErrors(X.size(), &X[0], &d.chargeRes[0], 0, &d.chargeResErr[0]);
            d.grQeAbs = new TGraphErrors(X.size(), &X[0], &d.qeAbs[0], 0, &d.qeAbsErr[0]);
            d.grQeDarkCorrNorm = new TGraphErrors(X.size(), &X[0], &d.qeDarkCorrNorm[0], 0, &d.qeDarkCorrNormErr[0]);
            d.grQeDarkCorrAbs = new TGraphErrors(X.size(), &X[0], &d.qeDarkCorrAbs[0], 0, &d.qeDarkCorrAbsErr[0]);
            d.grTts = new TGraphErrors(X.size(), &X[0], &d.tts[0], 0, &d.ttsErr[0]); // TTS Graph
            d.grFwhm = new TGraphErrors(X.size(), &X[0], &d.fwhm[0], 0, &d.fwhmErr[0]);
            d.grSigma = new TGraphErrors(X.size(), &X[0], &d.sigma[0], 0, &d.sigmaErr[0]);
            d.grPedStability = new TGraphErrors(X.size(), &X[0], &d.pedMean[0], 0, &d.pedSigma[0]);
            d.grPoissonQeNorm = new TGraphErrors(X.size(), &X[0], &d.poissonQeNorm[0], 0, &d.poissonQeNormErr[0]);
            d.grPoissonQeAbs = new TGraphErrors(X.size(), &X[0], &d.poissonQeAbs[0], 0, &d.poissonQeAbsErr[0]);

            // Raw-stage-angle companion: y = raw stage tilt at the same point,
            // x = whatever the report is plotted against. Downstream tools
            // (Overlay, the live GUI view) read this to label a point with the
            // stage angle that produced it without having to re-derive the
            // sign convention themselves.
            d.grRawStageAngle = new TGraphErrors(X.size(), &X[0], &d.tiltRaw[0], 0, 0);

            // Monitor-normalized graphs: empty for the monitor itself (i==0, never
            // filled above), so only build these when the vectors are populated.
            if (!d.monNormQe.empty())      d.grMonNormQe      = new TGraphErrors(X.size(), &X[0], &d.monNormQe[0], 0, &d.monNormQeErr[0]);
            if (!d.monNormQeDarkCorr.empty())     d.grMonNormQeDarkCorr     = new TGraphErrors(X.size(), &X[0], &d.monNormQeDarkCorr[0], 0, &d.monNormQeDarkCorrErr[0]);
            if (!d.monNormPoissonQeAbs.empty()) d.grMonNormPoissonQeAbs = new TGraphErrors(X.size(), &X[0], &d.monNormPoissonQeAbs[0], 0, &d.monNormPoissonQeAbsErr[0]);
            if (!d.monNormPoissonQe.empty())    d.grMonNormPoissonQe    = new TGraphErrors(X.size(), &X[0], &d.monNormPoissonQe[0], 0, &d.monNormPoissonQeErr[0]);
            if (!d.monNormSpe.empty())    d.grMonNormSpe    = new TGraphErrors(X.size(), &X[0], &d.monNormSpe[0], 0, &d.monNormSpeErr[0]);

            std::vector<TGraphErrors*> graphs = {d.grSpeNorm, d.grQeNorm, d.grDarkRate, d.grDarkCount, d.grSpeAbs, d.grQeAbs, d.grQeDarkCorrNorm, d.grQeDarkCorrAbs};
            for(auto* g : graphs) { 
                //  if(g) { g->SetMarkerStyle(20+i); g->SetMarkerColor(clr[i]); g->SetLineColor(clr[i]); g->GetXaxis()->SetTitle("Position Angle [degree]"); }
            }
            std::vector<TGraphErrors*> raw_graphs = {d.grSpeNorm, d.grQeNorm, d.grDarkRate, d.grDarkCount, d.grSpeAbs, d.grQeAbs};
            for(auto* g : raw_graphs) {
                //  if(g) { g->SetMarkerStyle(24+i); g->SetMarkerColor(clr[i]); g->SetLineColor(clr[i]); g->GetXaxis()->SetTitle("Position Angle [degree]"); }
            }
            std::vector<TGraphErrors*> corr_graphs = {d.grQeDarkCorrNorm, d.grQeDarkCorrAbs};
            for(auto* g : corr_graphs) {
                // if(g) { g->SetMarkerStyle(20+i); g->SetMarkerColor(clr[i]); g->SetLineColor(clr[i]); g->GetXaxis()->SetTitle("Position Angle [degree]"); }
            }
        };
        prep(pmt[i].X); prep(pmt[i].Y);
    } 

    const char* single_metrics[] = {"Raw SPE", "Dark Rate", "Raw Count", "Raw QE"}; 
    const char* y_titles[] = {"Charge [pC]", "Dark Rate [Hz]", "Dark Count [entries]", "Raw Rel.QE [%]"};

    for (int m = 0; m < 4; ++m) { 
        TCanvas *c = new TCanvas(Form("c_%s", single_metrics[m]), single_metrics[m], 4800,2400);
        c->Divide(2, 1); 

        // Top band, inside the extra headroom ApplyYMargin() reserves above
        // the data -- never overlaps a curve regardless of its shape.
        double legX1 = 0.62, legY1 = 0.78, legX2 = 0.90, legY2 = 0.92;

        for (int a = 0; a < 2; ++a) {
            c->cd(a + 1); gPad->SetGrid();
            gPad->SetLeftMargin(m == 1 || m == 2 ? 0.20 : 0.16); gPad->SetBottomMargin(0.12);

            TMultiGraph *mg = new TMultiGraph();
            TLegend *leg = new TLegend(legX1, legY1, legX2, legY2);
            leg->SetBorderSize(2); leg->SetFillStyle(1001); leg->SetFillColor(kWhite);

            for (int i : chs) {
                TGraphErrors* g = nullptr;
                if      (m == 0) g = (a==0 ? pmt[i].X.grSpeAbs : pmt[i].Y.grSpeAbs);
                else if (m == 1) g = (a==0 ? pmt[i].X.grDarkRate      : pmt[i].Y.grDarkRate);
                else if (m == 2) g = (a==0 ? pmt[i].X.grDarkCount      : pmt[i].Y.grDarkCount);
                else if (m == 3) g = (a==0 ? pmt[i].X.grQeAbs  : pmt[i].Y.grQeAbs);

                if (m < 3) {
                    if (g) {
                        styleGraph(g, i, 4);
                        mg->Add(g);
                        TString label = pmt[i].serial; label.ToUpper();
                        leg->AddEntry(g, label, "lp");
                    }
                } else { // m == 3 (Raw QE)
                    if (drawRaw && g) {
                        styleGraph(g, i, 0);
                        mg->Add(g); leg->AddEntry(g, pmt[i].serial + " (Count Raw)", "p");
                    }
                    TGraphErrors* g_corr = (a==0 ? pmt[i].X.grQeDarkCorrAbs : pmt[i].Y.grQeDarkCorrAbs);
                    if (drawCorrected && g_corr) {
                        styleGraph(g_corr, i, 1);
                        mg->Add(g_corr); leg->AddEntry(g_corr, pmt[i].serial + " (Count Corr)", "p");
                    }
                    TGraphErrors* g_raw_pqe = (a==0 ? pmt[i].X.grPoissonQeAbs : pmt[i].Y.grPoissonQeAbs);
                    if (drawPoissonRaw && g_raw_pqe) {
                        styleGraph(g_raw_pqe, i, 2);
                        mg->Add(g_raw_pqe); leg->AddEntry(g_raw_pqe, pmt[i].serial + " (Poi Raw)", "p");
                    }
                    TGraphErrors* g_pqe = (a==0 ? pmt[i].X.grPoissonQeNorm : pmt[i].Y.grPoissonQeNorm);
                    if (drawPoisson && g_pqe) {
                        styleGraph(g_pqe, i, 3);
                        mg->Add(g_pqe); leg->AddEntry(g_pqe, pmt[i].serial + " (Poi Corr)", "p");
                    }
                }
            }

            mg->SetTitle(Form("%s (%s-Axis);%s;%s", single_metrics[m], (a==0?"X":"Y"), AngleAxisTitle(), y_titles[m]));
            ApplyYMargin(mg);
            mg->Draw("AP");
            mg->GetYaxis()->SetTitleOffset(m == 1 || m == 2 ? 2.0 : 1.2);
            leg->Draw();
        }
        c->Print(m == 0 ? uniPdfPath + "(" : uniPdfPath);
    }

    TCanvas *c_comb = new TCanvas("c_combined_qe_gain", "Rel QE & Gain", 4800, 3600);
    c_comb->Divide(2, 2); 
    const char* combined_names[] = {"Rel QE", "Rel Gain"};

    for (int row = 0; row < 2; ++row) {     
        for (int col = 0; col < 2; ++col) { 
            c_comb->cd(row * 2 + col + 1); gPad->SetGrid();
            TMultiGraph *mg = new TMultiGraph();
            TLegend *leg = new TLegend(0.50, 0.76, 0.90, 0.92);
            leg->SetBorderSize(2); leg->SetFillStyle(1001); leg->SetFillColor(kWhite);
            leg->SetNColumns(2); leg->SetTextSize(0.028);

            for (int i : chs) {
                if (row == 0) { // Rel QE
                    TGraphErrors* g_raw = (col == 0) ? pmt[i].X.grQeNorm : pmt[i].Y.grQeNorm;
                    TGraphErrors* g_corr = (col == 0) ? pmt[i].X.grQeDarkCorrNorm : pmt[i].Y.grQeDarkCorrNorm;
                    TGraphErrors* g_raw_pqe = (col == 0) ? pmt[i].X.grPoissonQeAbs : pmt[i].Y.grPoissonQeAbs;
                    TGraphErrors* g_pqe = (col == 0) ? pmt[i].X.grPoissonQeNorm : pmt[i].Y.grPoissonQeNorm;

                    if (drawRaw && g_raw) {
                        styleGraph(g_raw, i, 0);
                        mg->Add(g_raw); leg->AddEntry(g_raw, pmt[i].serial + " (Count Raw)", "p");
                    }
                    if (drawCorrected && g_corr) {
                        styleGraph(g_corr, i, 1);
                        mg->Add(g_corr); leg->AddEntry(g_corr, pmt[i].serial + " (Count Corr)", "p");
                    }
                    if (drawPoissonRaw && g_raw_pqe) {
                        styleGraph(g_raw_pqe, i, 2);
                        mg->Add(g_raw_pqe); leg->AddEntry(g_raw_pqe, pmt[i].serial + " (Poi Raw)", "p");
                    }
                    if (drawPoisson && g_pqe) {
                        styleGraph(g_pqe, i, 3);
                        mg->Add(g_pqe); leg->AddEntry(g_pqe, pmt[i].serial + " (Poi Corr)", "p");
                    }

                } else { // Rel Gain
                    TGraphErrors* g_spe = (col == 0) ? pmt[i].X.grSpeNorm : pmt[i].Y.grSpeNorm;
                    if (g_spe) {
                        styleGraph(g_spe, i, 4);
                        mg->Add(g_spe); leg->AddEntry(g_spe, pmt[i].serial, "p");
                    }
                }
            }
            TString yAxisTitle = (row == 0) ? "Relative QE" : "Relative Gain";
            mg->SetTitle(Form("%s (%s-Axis);%s;%s", combined_names[row], (col==0?"X":"Y"), AngleAxisTitle(), yAxisTitle.Data()));
            ApplyYMargin(mg);
            mg->Draw("AP"); leg->Draw();
        }
    }
    c_comb->Print(uniPdfPath);

    TCanvas *c_qe_dark = new TCanvas("c_qe_dark", "Raw QE vs Dark Count", 4800, 3600);
    TPad *pad_info = new TPad("pad_info", "Information", 0.0, 0.85, 1.0, 1.0);
    pad_info->SetBottomMargin(0); pad_info->Draw();
    TPad *pad_graphs = new TPad("pad_graphs", "Graphs", 0.0, 0.0, 1.0, 0.85);
    pad_graphs->SetTopMargin(0); pad_graphs->Draw(); pad_graphs->Divide(2, 1); 

    pad_info->cd();
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
    pt->SetBorderSize(2); pt->SetFillColor(kWhite);
    pt->SetTextAlign(22); pt->SetTextFont(42); pt->SetTextSize(0.25); 
    pt->AddText("[ Measurement Settings & Info ]");
    pt->AddText(Form("Block: %d ~ %d", run_start, run_end));
    pt->Draw();

    for (int col = 0; col < 2; ++col) { 
        pad_graphs->cd(col + 1);
        TPad *pad_qe = new TPad(Form("pad_qe_%d", col), "Raw QE", 0.0, 0.4, 1.0, 1.0);
        TPad *pad_dark = new TPad(Form("pad_dark_%d", col), "Dark Count", 0.0, 0.0, 1.0, 0.4);
        pad_qe->SetBottomMargin(0.005); pad_dark->SetTopMargin(0.005); pad_dark->SetBottomMargin(0.25);
        pad_qe->SetLeftMargin(0.15); pad_dark->SetLeftMargin(0.15);
        pad_qe->Draw(); pad_dark->Draw();

        pad_qe->cd(); gPad->SetGrid();
        TMultiGraph *mg_qe = new TMultiGraph();
        TLegend *leg_qe = new TLegend(0.30, 0.77, 0.70, 0.90); leg_qe->SetBorderSize(2); leg_qe->SetNColumns(2);
        leg_qe->SetTextSize(0.03);

        for (int i : chs) {
            TGraphErrors* g_raw  = (col == 0) ? pmt[i].X.grQeAbs : pmt[i].Y.grQeAbs;
            TGraphErrors* g_corr = (col == 0) ? pmt[i].X.grQeDarkCorrAbs : pmt[i].Y.grQeDarkCorrAbs;
            TGraphErrors* g_raw_pqe  = (col == 0) ? pmt[i].X.grPoissonQeAbs : pmt[i].Y.grPoissonQeAbs;
            TGraphErrors* g_pqe  = (col == 0) ? pmt[i].X.grPoissonQeNorm : pmt[i].Y.grPoissonQeNorm;

            if (drawRaw && g_raw) {
                styleGraph(g_raw, i, 0);
                mg_qe->Add(g_raw); leg_qe->AddEntry(g_raw, pmt[i].serial + " (Count Raw)", "p");
            }
            if (drawCorrected && g_corr) {
                styleGraph(g_corr, i, 1);
                mg_qe->Add(g_corr); leg_qe->AddEntry(g_corr, pmt[i].serial + " (Count Corr)", "p");
            }
            if (drawPoissonRaw && g_raw_pqe) {
                styleGraph(g_raw_pqe, i, 2);
                mg_qe->Add(g_raw_pqe); leg_qe->AddEntry(g_raw_pqe, pmt[i].serial + " (Poi Raw)", "p");
            }
            if (drawPoisson && g_pqe) {
                styleGraph(g_pqe, i, 3);
                mg_qe->Add(g_pqe); leg_qe->AddEntry(g_pqe, pmt[i].serial + " (Poi Corr)", "p");
            }
        }

        mg_qe->SetTitle(Form("Raw QE (%s-Axis); ;Relative QE [%%]", (col==0?"X":"Y")));
        ApplyYMargin(mg_qe);
        mg_qe->Draw("AP"); mg_qe->GetXaxis()->SetLabelSize(0);
        leg_qe->Draw();

        pad_dark->cd(); gPad->SetGrid();
        TMultiGraph *mg_dark = new TMultiGraph();
        for (int i : chs) {
            TGraphErrors* g = (col == 0) ? pmt[i].X.grDarkCount : pmt[i].Y.grDarkCount; 
            if (g) {
                styleGraph(g, i, 4);
                mg_dark->Add(g);
            }
        }

        mg_dark->SetTitle(Form(";%s;Dark Count [entries]", AngleAxisTitle()));
        ApplyYMargin(mg_dark);
        mg_dark->Draw("AP"); mg_dark->GetYaxis()->SetTitleOffset(0.9);
        mg_dark->GetXaxis()->SetTitleSize(0.10); mg_dark->GetXaxis()->SetLabelSize(0.08);
        mg_dark->GetXaxis()->SetTitleOffset(1.05);
        mg_dark->GetYaxis()->SetTitleSize(0.08); mg_dark->GetYaxis()->SetLabelSize(0.08);
    }
    c_qe_dark->Print(uniPdfPath);

    const char* row_names[] = {"Raw SPE", "Dark Rate", "Raw Count", "Rel QE", "Rel Gain"};
    const char* row_units[] = {"Charge [pC]", "Dark Rate [Hz]", "Dark Count [entries]", "Relative QE", "Relative Gain"};

    for (int i = 0; i < 3; ++i) {
        if (std::find(chs.begin(), chs.end(), i) == chs.end()) continue;  // skip unselected channel
        TCanvas *c_ind = new TCanvas(Form("c_summary_%s", pmt[i].serial.Data()), pmt[i].serial, 3600, 6000);
        c_ind->Divide(2, 5); 
        for (int r = 0; r < 5; ++r) { 
            for (int c = 0; c < 2; ++c) { 
                c_ind->cd(r * 2 + c + 1); gPad->SetGrid();
                gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.12);
                TGraphErrors* g = nullptr;
                auto& d = (c == 0) ? pmt[i].X : pmt[i].Y; 

                if      (r == 0) g = d.grSpeAbs; else if (r == 1) g = d.grDarkRate;      
                else if (r == 2) g = d.grDarkCount;      else if (r == 3) g = d.grQeNorm;     
                else if (r == 4) g = d.grSpeNorm;     

                if (g || r == 3) {
                    TMultiGraph *mg_temp = new TMultiGraph();
                    TLegend *leg = new TLegend(0.55, 0.80, 0.90, 0.95); leg->SetBorderSize(2); leg->SetFillColor(kWhite);
                    leg->SetNColumns(2); leg->SetTextSize(0.05);

                    if (r == 3) { // Rel QE
                        TGraphErrors* g_raw = d.grQeNorm;
                        TGraphErrors* g_corr = d.grQeDarkCorrNorm;
                        TGraphErrors* g_raw_pqe = d.grPoissonQeAbs;
                        TGraphErrors* g_pqe = d.grPoissonQeNorm;

                        if (drawRaw && g_raw) {
                            styleGraph(g_raw, i, 0);
                            mg_temp->Add(g_raw); leg->AddEntry(g_raw, "Count Raw", "lp");
                        }
                        if (drawCorrected && g_corr) {
                            styleGraph(g_corr, i, 1);
                            mg_temp->Add(g_corr); leg->AddEntry(g_corr, "Count Corr", "lp");
                        }
                        if (drawPoissonRaw && g_raw_pqe) {
                            styleGraph(g_raw_pqe, i, 2);
                            mg_temp->Add(g_raw_pqe); leg->AddEntry(g_raw_pqe, "Poi Raw", "lp");
                        }
                        if (drawPoisson && g_pqe) {
                            styleGraph(g_pqe, i, 3);
                            mg_temp->Add(g_pqe); leg->AddEntry(g_pqe, "Poi Corr", "lp");
                        }
                    } else if (g) { // Other Metrics (Raw SPE, Dark Rate, etc.)
                        styleGraph(g, i, 4);
                        mg_temp->Add(g);
                    } 

                    if (mg_temp->GetListOfGraphs()) {
                        mg_temp->SetTitle(Form("%s - %s (%s-Axis);%s;%s", pmt[i].serial.Data(),row_names[r], (c==0?"X":"Y"), AngleAxisTitle(), row_units[r]));
                        ApplyYMargin(mg_temp);
                        mg_temp->Draw("AP");
                        mg_temp->GetYaxis()->SetTitleOffset(r == 1 || r == 2 ? 1.6 : 1.3);
                        if(r==3) leg->Draw();
                    }
                }


            }
        }
        c_ind->Print(uniPdfPath);
        delete c_ind; 
    }

    TCanvas *c_info = new TCanvas("c_info", "Run Information", 4800, 4000);
    TPaveText *pt_info = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
    pt_info->SetBorderSize(2); pt_info->SetFillColor(kWhite);
    pt_info->SetTextAlign(12); pt_info->SetTextFont(42); pt_info->SetTextSize(0.035);

    pt_info->AddText("#bf{[ Analysis Run Summary ]}");
    pt_info->AddText("----------------------------------------------------------------------");
    pt_info->AddText(Form("Data Taken Date : %s", tag.Data()));
    pt_info->AddText(Form("Run Block       : %d ~ %d", run_start, run_end));

    pt_info->AddText(Form("Laser Intensity : %d mA", info_Laser));
    pt_info->AddText(Form("Expert          : %s", info_Expert.Data()));
    pt_info->AddText(Form("Shifter         : %s", info_Shifter.Data()));
    pt_info->AddText(Form("Comp Coil       : %s", info_Bfield.Data()));
    // Cable-position labels matching HV_Control_SW's Dark Box naming (Ch0=Mon.,
    // Ch1=Rot#1, Ch2=Rot#2), so the report identifies PMTs by physical cable
    // position instead of the internal channel index.
    const char* chDisplayName[] = {"Mon.", "Rot#1", "Rot#2"};
    for (int k = 0; k < 3; ++k) {
        if (std::find(chs.begin(), chs.end(), k) == chs.end()) continue;  // only selected channels
        pt_info->AddText(Form("  - %s: #color[2]{%s} | High Voltage: #bf{%s V}", chDisplayName[k], pmt[k].serial.Data(), pmt[k].hv.Data()));
    }
    pt_info->AddText("");
    pt_info->AddText("#bf{2. Angle Transformation}");
    pt_info->AddText("  Formula: y = -0.0049x^{2} + 1.7515x - 0.0402"); 
    pt_info->AddText("----------------------------------------------------------------------");

    pt_info->Draw();
    c_info->Print(uniPdfPath + ")");
    std::cout << "[INFO] Saved combined report: " << uniPdfPath << std::endl;

    TFile *fOut = new TFile(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", file_suffix.Data()), "RECREATE");
    for (int i = 0; i < 3; ++i) {
        if (std::find(chs.begin(), chs.end(), i) == chs.end()) continue;  // skip unselected channel
        auto save_graphs = [&](PMTData& d, TString axis) {
            if (!d.grSpeNorm) return;
            d.grSpeNorm->Write(Form("gr_%s_%s_RelGain", pmt[i].serial.Data(), axis.Data()));
            d.grQeNorm->Write(Form("gr_%s_%s_RelQE", pmt[i].serial.Data(), axis.Data()));
            d.grDarkRate->Write(Form("gr_%s_%s_DarkRate", pmt[i].serial.Data(), axis.Data()));
            d.grSpeAbs->Write(Form("gr_%s_%s_RawSPE", pmt[i].serial.Data(), axis.Data()));
            if(d.grChargeRes) d.grChargeRes->Write(Form("gr_%s_%s_ChargeResolution", pmt[i].serial.Data(), axis.Data()));
            d.grQeAbs->Write(Form("gr_%s_%s_RawQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grQeDarkCorrNorm) d.grQeDarkCorrNorm->Write(Form("gr_%s_%s_CorrRelQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grQeDarkCorrAbs) d.grQeDarkCorrAbs->Write(Form("gr_%s_%s_CorrRawQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grDarkCount) d.grDarkCount->Write(Form("gr_%s_%s_DarkCount", pmt[i].serial.Data(), axis.Data()));
            if(d.grTts) d.grTts->Write(Form("gr_%s_%s_TTS", pmt[i].serial.Data(), axis.Data())); // Save TTS
            if(d.grFwhm) d.grFwhm->Write(Form("gr_%s_%s_FWHM", pmt[i].serial.Data(), axis.Data()));
            if(d.grSigma) d.grSigma->Write(Form("gr_%s_%s_Sigma", pmt[i].serial.Data(), axis.Data()));
            if(d.grPedStability) d.grPedStability->Write(Form("gr_%s_%s_PedStability", pmt[i].serial.Data(), axis.Data()));
            if(d.grPoissonQeAbs) d.grPoissonQeAbs->Write(Form("gr_%s_%s_RawPoissonQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grPoissonQeNorm) d.grPoissonQeNorm->Write(Form("gr_%s_%s_PoissonQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grRawStageAngle) d.grRawStageAngle->Write(Form("gr_%s_%s_RawStageAngle", pmt[i].serial.Data(), axis.Data()));
            // Monitor(CH0)-normalized: Test_absolute / Monitor_absolute. Not
            // written for the monitor itself (no graph built, i==0 skipped above).
            if(d.grMonNormQe)      d.grMonNormQe->Write(Form("gr_%s_%s_MonNorm_RawQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grMonNormQeDarkCorr)     d.grMonNormQeDarkCorr->Write(Form("gr_%s_%s_MonNorm_CorrRawQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grMonNormPoissonQeAbs) d.grMonNormPoissonQeAbs->Write(Form("gr_%s_%s_MonNorm_RawPoissonQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grMonNormPoissonQe)    d.grMonNormPoissonQe->Write(Form("gr_%s_%s_MonNorm_PoissonQE", pmt[i].serial.Data(), axis.Data()));
            if(d.grMonNormSpe)    d.grMonNormSpe->Write(Form("gr_%s_%s_MonNorm_RelGain", pmt[i].serial.Data(), axis.Data()));
        };
        save_graphs(pmt[i].X, "X");
        save_graphs(pmt[i].Y, "Y");
    }
    TString chmap = "";
    for (int i = 0; i < 3; ++i) {
        if (std::find(chs.begin(), chs.end(), i) == chs.end()) continue; 
        TString s = pmt[i].serial; s.ToUpper();
        chmap += Form("%d:%s;", i, s.Data());
    }
    TNamed chmapObj("ChannelMap", chmap.Data());
    chmapObj.Write();
    TNamed axisObj("AngleAxis", kUseRawStageAxis ? "rawstage" : "hamamatsu");
    axisObj.Write();
    fOut->Close();
}
