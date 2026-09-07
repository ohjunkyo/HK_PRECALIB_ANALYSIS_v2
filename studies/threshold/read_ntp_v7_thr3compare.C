#include <vector>
#include <string>
#include <iostream>
#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TLine.h>
#include <TMath.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TNamed.h>
#include <TLegend.h>
#include <TLatex.h>

#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"
#include "path_builder2.h"
#include "./Base/analysisCode/DrawFormat.cpp"
#include "./Base/analysisCode/Analysis.cpp"

const int N_max = 3; // max p.e. component in the Poisson-Gaussian model
const double PE_DRAW_FRAC = 0.02; // n p.e. shown if P(n;mu) >= 2% of P(1;mu)

void HistFont(TH1D *hist){ /*{{{*/
    hist->GetXaxis()->SetTitleSize(0.06); 
    hist->GetXaxis()->SetLabelSize(0.06); 
    hist->GetYaxis()->SetTitleSize(0.06); 
    hist->GetYaxis()->SetLabelSize(0.06); 
}/*}}}*/

// Timing PDF
Double_t exGausPDF(Double_t *x, Double_t *par) {/*{{{*/
    Double_t xx = x[0]; Double_t A = par[0]; Double_t mu = par[1]; Double_t sigma = par[2]; Double_t tau = par[3];
    Double_t bg = par[4];
    Double_t invTau = 1.0 / tau;
    Double_t z = (sigma * invTau - (xx - mu) / sigma) / TMath::Sqrt(2.0);
    Double_t expo = (sigma*sigma*invTau*invTau/2.0) - (xx - mu)*invTau;
    return A * invTau/2.0 * TMath::Exp(expo) * TMath::Erfc(z) + bg;
}/*}}}*/
// SPE Fitting Function with Erf-based Backscattering Model{{{
double spe_fitting_function(double *x, double *par) {
    double amplitude = par[0];
    double mu        = par[1];
    double pedestal  = par[2];
    double speMean  = par[3]; // SPE Gain
    double pedSigma = par[4];
    double speSigma = par[5];
    double w         = par[6]; // Backscattering weight fraction

    double totalSum = 0;
    double xx = x[0];

    totalSum += TMath::Poisson(0, mu) * TMath::Gaus(xx, pedestal, pedSigma, true);

    for (int n = 1; n <= N_max; ++n) {
        double prob = TMath::Poisson(n, mu);
        double meanN = pedestal + n * speMean;
        double sigmaN = TMath::Sqrt(pedSigma*pedSigma + n * speSigma*speSigma);

        double gausComponent = TMath::Gaus(xx, meanN, sigmaN, true);

        double erf1 = TMath::Erf((xx - pedestal) / par[4]);
        double erf2 = TMath::Erf((xx - meanN) / sigmaN);
        double backComponent = (erf1 - erf2) / (2.0 * n * speMean); 

        totalSum += prob * ((1.0 - w) * gausComponent + w * backComponent);
    }
    return amplitude * totalSum;
}/*}}}*/
double single_pe_comp(double *x, double *par) {/*{{{*/
    double xx = x[0];
    int n = TMath::Nint(par[7]); 
    double prob = TMath::Poisson(n, par[1]);

    if (n == 0) return par[0] * prob * TMath::Gaus(xx, par[2], par[4], true);

    double meanN = par[2] + n * par[3];
    double sigmaN = TMath::Sqrt(par[4]*par[4] + n * par[5]*par[5]);

    double gausComponent = TMath::Gaus(xx, meanN, sigmaN, true);
    double erf1 = TMath::Erf((xx - par[2]) / par[4]);
    double erf2 = TMath::Erf((xx - meanN) / sigmaN);
    double backComponent = (erf1 - erf2) / (2.0 * n * par[3]);

    return par[0] * prob * ((1.0 - par[6]) * gausComponent + par[6] * backComponent);
}/*}}}*/
TF1* PreFit(TH1D* hist, double rangeMin, double rangeMax, const std::string& name) {/*{{{*/
    int binMin = hist->GetXaxis()->FindBin(rangeMin);
    int binMax = hist->GetXaxis()->FindBin(rangeMax);
    hist->GetXaxis()->SetRange(binMin, binMax);
    double initialMean = hist->GetMean(); double initialSigma = hist->GetStdDev();
    double initialAmp = hist->GetBinContent(hist->FindBin(initialMean));
    hist->GetXaxis()->SetRange(); 
    TF1* preFit = new TF1(name.c_str(), "gaus", rangeMin, rangeMax);
    preFit->SetParameters(initialAmp, initialMean, initialSigma);
    hist->Fit(preFit, "RQN+"); return preFit;
}/*}}}*/
TF1* FitGaussian(TH1D* hist, double pedestalMin, double pedestalMax, double signalMin, double signalMax,/*{{{*/
        double& signalMean, double& signalError, double& resolutionOut, double& resolutionErrOut,
        TH1D *hist2, std::string unit = "pC", double* peakToValleyOut = nullptr) {

    double histMin = hist ? hist->GetXaxis()->GetXmin() : -5.0;
    double histMax = hist ? hist->GetXaxis()->GetXmax() : 10.0;

    // No data to fit (e.g., a disconnected channel) → skip fitting cleanly instead
    // of letting ROOT spam errors. Caller-visible result is flagged "Fit Empty".
    if (!hist || !hist2 || hist->GetEntries() < 1 || hist2->GetEntries() < 1) {
        signalMean = 0; signalError = 0; resolutionOut = 0; resolutionErrOut = 0;
        if (peakToValleyOut) *peakToValleyOut = 0;
        std::cout << "  (C) SPE Fitting Results : \033[1;33mFit Empty\033[0m"
                  << "  (no signal events)" << std::endl;
        TF1* dummy = new TF1("total_fit", "gaus", histMin, histMax);
        dummy->SetParameters(0.0, 0.0, 1.0);   // flat: Eval()=0, GetParameter(n)=0
        return dummy;
    }
    TF1* prePedFit = PreFit(hist, pedestalMin, pedestalMax, "pre_ped_gaus");
    double pMean = prePedFit->GetParameter(1); 
    double pSigma = prePedFit->GetParameter(2);

    double sigPeakX = hist2->GetXaxis()->GetBinCenter(hist2->GetMaximumBin());
    double fitWindow = sigPeakX * 0.5; 
    if (fitWindow < 2.5) fitWindow = 2.5; 

    TF1* preSigFit = PreFit(hist2, sigPeakX - fitWindow, sigPeakX + fitWindow, "pre_sig_gaus");
    double roughMean = preSigFit->GetParameter(1); 
    double roughSigma = preSigFit->GetParameter(2);

    double window = roughSigma * 0.5; 
    if(window < 2.0) window = 2.5; 
    TF1* finalSigFit = new TF1("final_sig_fit", "gaus", roughMean - window, roughMean + window);
    finalSigFit->SetParameters(preSigFit->GetParameters()); 
    hist2->Fit(finalSigFit, "RQ0");

    double initGain = finalSigFit->GetParameter(1) - pMean; 
    double initSigma = finalSigFit->GetParameter(2);
    if (initSigma > initGain * 0.4) initSigma = initGain * 0.35; 

    TF1* total_fit = new TF1("total_fit", spe_fitting_function, histMin, histMax, 7);
    total_fit->SetParNames("Amplitude", "Mu", "Ped_Mean", "SPE_Gain", "Ped_Sigma", "SPE_Sigma", "Back_Weight");
    double totalArea = hist->GetEntries() * hist->GetBinWidth(1);
    total_fit->SetParameters(totalArea, 0.1, pMean, initGain, pSigma, initSigma); total_fit->SetParameter(6, 0.05);

    total_fit->SetParLimits(1, 0.001, 2.5); // Mu range 
    total_fit->SetParLimits(2, pMean - pSigma, pMean + pSigma); // Ped_Mean
    total_fit->SetParLimits(4, pSigma * 0.5, pSigma * 1.3);        // Ped_Sigma

    if (unit == "pC") {
        if (initGain < 1.6) { //
            total_fit->SetParLimits(3, initGain * 0.85, initGain * 1.15); // SPE_Gain
            total_fit->SetParLimits(5, pSigma * 0.8, initGain * 0.45);   // SPE_Sigma
            total_fit->SetParameter(6, 0.03);
            total_fit->SetParLimits(6, 0.0, 0.15);                         // Back_Weight
        } else { // Normally
            total_fit->SetParLimits(3, initGain * 0.75, initGain * 1.25);
            total_fit->SetParLimits(5, pSigma * 0.5, initGain * 0.55);
            total_fit->SetParameter(6, 0.08);
            total_fit->SetParLimits(6, 0.0, 0.35);
        }

    } else {
        if (initGain < 40.0) {
            total_fit->SetParLimits(3, initGain * 0.80, initGain * 1.20); 
            total_fit->SetParLimits(5, pSigma * 0.4, initGain * 0.32);   
            total_fit->SetParameter(6, 0.01);
            total_fit->SetParLimits(6, 0.0, 0.08);                         
        } else {
            total_fit->SetParLimits(3, initGain * 0.65, initGain * 1.35);
            total_fit->SetParLimits(5, pSigma * 0.3, initGain * 0.45);
            total_fit->SetParameter(6, 0.02);
            total_fit->SetParLimits(6, 0.0, 0.20);
        }
    }

    hist->Fit(total_fit, "RMQ+");

    // signalMean is a SUM of two fit parameters (Ped_Mean + SPE_Gain), so its
    // error needs both parameters' own uncertainty AND their covariance --
    // using GetParError(3) alone silently dropped Ped_Mean's contribution
    // entirely, which is usually small but not negligible, and let the
    // reported error be smaller than it should be even on healthy fits
    // (2026-08-30/31 investigation into the Gain Curve's implausibly tight
    // errors at low HV; see project memory "gain-curve-spe-fit-decision").
    // Note this does NOT fully fix the low-HV (P/V~1.0) case -- there the
    // Hesse matrix itself is unreliable regardless of how its entries are
    // combined, which is why Draw_GainCurve_v1.C separately excludes those
    // points by P/V ratio. This is the correct fix for every other point.
    signalMean = total_fit->GetParameter(2) + total_fit->GetParameter(3);
    {
        double errPed = total_fit->GetParError(2);
        double errSpe = total_fit->GetParError(3);
        TVirtualFitter* vf = TVirtualFitter::GetFitter();
        double covPedSpe = vf ? vf->GetCovarianceMatrixElement(2, 3) : 0.0;
        double combinedVar = errPed * errPed + errSpe * errSpe + 2.0 * covPedSpe;
        signalError = (combinedVar > 0) ? std::sqrt(combinedVar) : errSpe;
    }
    double finalMu = total_fit->GetParameter(1);
    double resolution = total_fit->GetParameter(5) / total_fit->GetParameter(3) * 100;
    // Charge resolution = SPE_Sigma/SPE_Gain, propagated the same
    // independent-parameter-error way as every other ratio in this file
    // (e.g. MonitorNormalize) -- SPE_Sigma and SPE_Gain are technically
    // correlated fit parameters, but this codebase's convention elsewhere
    // is the simpler independent approximation, not the full covariance.
    {
        double speSigma = total_fit->GetParameter(5), speSigmaErr = total_fit->GetParError(5);
        double speGain = total_fit->GetParameter(3), speGainErr = total_fit->GetParError(3);
        resolutionOut = (speGain != 0) ? (speSigma / speGain) * 100.0 : 0.0;
        double relErrSigma = (speSigma != 0) ? speSigmaErr / speSigma : 0.0;
        double relErrGain = (speGain != 0) ? speGainErr / speGain : 0.0;
        resolutionErrOut = resolutionOut * std::sqrt(relErrSigma * relErrSigma + relErrGain * relErrGain);
    }

    std::cout << "  (C) SPE Fitting Results (Poisson-Gaussian)" << std::endl;
    std::cout << "  --------------------------------------------------------------" << std::endl;
    if (signalMean > 2.0) {
        std::cout << Form("  %-25s : \033[1;31m%4.3f +/- %.3f [%s]\033[0m", "Signal Mean", signalMean, signalError, unit.c_str()) << std::endl;
    } else {
        std::cout << Form("  %-25s : \033[1;34m%4.3f +/- %.3f [%s]\033[0m", "Signal Mean", signalMean, signalError, unit.c_str()) << std::endl;
    }
    std::cout << Form("  %-25s : %4.2f%%", "Resolution", resolution) << std::endl;
    std::cout << Form("  %-25s : %4.4f +/- %.4f", "Pedestal Sigma", total_fit->GetParameter(4), total_fit->GetParError(4)) << std::endl;
    std::cout << Form("  %-25s : %4.4f +/- %.4f", "SPE Sigma", total_fit->GetParameter(5), total_fit->GetParError(5)) << std::endl;
    std::cout << Form("  %-25s : %4.4f +/- %.4f", "Fit Mu", finalMu, total_fit->GetParError(1)) << std::endl;
    std::cout << "  --------------------------------------------------------------" << std::endl;

    int colors[] = {kOrange, kBlack, kBlue, kMagenta}; // ped, 1pe, 2pe, 3pe
    double fitMuDraw = total_fit->GetParameter(1);   // occupancy, for legend gating
    TLegend *leg = new TLegend(0.55, 0.55, 0.93, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(132); leg->SetTextSize(0.04);
    leg->AddEntry(hist, "Data", "l");
    leg->AddEntry(total_fit, "Total Fit", "l");

    // n = 0 Component (Pure Pedestal)
    TF1* f_ped = new TF1("f_comp_0", single_pe_comp, histMin, histMax, 8);
    for(int p=0; p<7; p++) f_ped->SetParameter(p, total_fit->GetParameter(p));
    f_ped->SetParameter(7, 0.0); 
    f_ped->SetLineColor(colors[0]); 
    f_ped->SetLineStyle(1); f_ped->Draw("SAME");
    leg->AddEntry(f_ped, "Pedestal", "l");

    // Draw order: Backscattering (thin) → Gaussian (thick) → Total p.e. (medium) → Total Fit (last/on top)
    for (int n = 1; n <= N_max; ++n) {
        // Only draw/label p.e. peaks that actually carry weight. 1 p.e. is always
        // shown (original format); higher peaks appear only if the fitted occupancy
        // makes them non-negligible vs the 1 p.e. peak.
        bool showN = (n == 1) ||
                      (TMath::Poisson(n, fitMuDraw) >= PE_DRAW_FRAC * TMath::Poisson(1, fitMuDraw));
        if (!showN) continue;

        // Backscattering only — thin solid
        TF1* f_back_only = new TF1(Form("f_back_%d", n),
                "[0] * TMath::Poisson([7], [1]) * [6] * (TMath::Erf((x-[2])/[4]) - TMath::Erf((x-([2]+[7]*[3]))/TMath::Sqrt([4]*[4]+[7]*[5]*[5]))) / (2.0 * [7]*[3])",
                histMin, histMax);
        for(int p=0; p<7; p++) f_back_only->SetParameter(p, total_fit->GetParameter(p));
        f_back_only->SetParameter(7, (double)n);
        f_back_only->SetLineColor(colors[n]);
        f_back_only->SetLineStyle(1);
        f_back_only->SetLineWidth(1);
        f_back_only->Draw("SAME");

        // Gaussian only — thick solid
        TF1* f_gaus_only = new TF1(Form("f_gaus_%d", n),
                "[0] * TMath::Poisson([7], [1]) * (1.0 - [6]) * TMath::Gaus(x, [2] + [7]*[3], TMath::Sqrt([4]*[4] + [7]*[5]*[5]), 1)",
                histMin, histMax);
        for(int p=0; p<7; p++) f_gaus_only->SetParameter(p, total_fit->GetParameter(p));
        f_gaus_only->SetParameter(7, (double)n);
        f_gaus_only->SetLineColor(colors[n]);
        f_gaus_only->SetLineStyle(1);
        f_gaus_only->SetLineWidth(3);
        f_gaus_only->Draw("SAME");

        // Total p.e. component (Gaus + Back) — medium solid
        TF1* f_total = new TF1(Form("f_total_%d", n), single_pe_comp, histMin, histMax, 8);
        for(int p=0; p<7; p++) f_total->SetParameter(p, total_fit->GetParameter(p));
        f_total->SetParameter(7, (double)n);
        f_total->SetLineColor(colors[n]);
        f_total->SetLineStyle(1);
        f_total->SetLineWidth(2);
        f_total->Draw("SAME");
        leg->AddEntry(f_total, Form("%d p.e.", n), "l");
    }
    leg->Draw("same");

    total_fit->Draw("SAME");

    int startBin = hist->GetXaxis()->FindBin(pMean);
    int endBin   = hist->GetXaxis()->FindBin(signalMean);
    double valleyY = 1e9; int valleyBin = -1;
    for (int b = startBin; b <= endBin; ++b) {
        double currentBin = hist->GetBinContent(b);
        if (currentBin > 0 && currentBin < valleyY) { 
            valleyY = currentBin; valleyBin = b; 
        }
    }
    double valleyX = (valleyBin != -1) ? hist->GetXaxis()->GetBinCenter(valleyBin) : (pMean + roughMean) / 2.0;
    // Peak and valley must come from the SAME curve/histogram that's actually
    // drawn (hist, and its red Total Fit curve), not finalSigFit -- that's
    // only a throwaway single-Gaussian pre-fit used to seed total_fit's
    // initial parameters, fit on hist2 (the amplitude-thresholded SUBSET),
    // while valleyY above is read from hist (the full, unfiltered spectrum).
    // Mixing a subset-histogram pre-fit peak with a full-histogram valley
    // made P/V inconsistent with what the plot shows.
    double peakY = total_fit->Eval(signalMean);
    double peakToValleyRatio = (valleyY > 0) ? peakY / valleyY : 0.0;
    if (peakToValleyOut) *peakToValleyOut = peakToValleyRatio;

    TLine* valleyLine = new TLine(valleyX, gPad->GetUymin(), valleyX, valleyY);
    valleyLine->SetLineColor(kBlack); 
    valleyLine->SetLineStyle(2); 
    valleyLine->SetLineWidth(2); 
    //valleyLine->Draw("same");

    TLine* speMeanLine = new TLine(signalMean, hist->GetMinimum(), signalMean, hist->GetMaximum());
    speMeanLine->SetLineColor(kRed); speMeanLine->SetLineWidth(1); speMeanLine->SetLineStyle(1);
    //speMeanLine->Draw("same");

    TLatex* speMeanText = new TLatex(0.55, 0.52, Form("#mu: %.2f #pm %.3f [%s]", signalMean, signalError, unit.c_str()));
    speMeanText->SetTextFont(132); 
    speMeanText->SetTextSize(0.04); 
    speMeanText->SetNDC(true); 
    speMeanText->Draw();

    TLatex* resText = new TLatex(0.55, 0.47, Form("#sigma/#mu: %.2f %%", resolution));
    resText->SetTextFont(132); 
    resText->SetTextSize(0.04); 
    resText->SetNDC(true); 
    resText->Draw();

    TLatex* pvText = new TLatex(0.55, 0.42, Form("P/V: %.2f", peakToValleyRatio));
    pvText->SetTextFont(132); 
    pvText->SetTextSize(0.04); 
    pvText->SetNDC(true); 
    pvText->Draw();

    TLatex* pedSigText = new TLatex(0.55, 0.37, Form("#sigma_{ped}: %.3f #pm %.3f", total_fit->GetParameter(4), total_fit->GetParError(4)));
    pedSigText->SetTextFont(132); pedSigText->SetTextSize(0.035); pedSigText->SetNDC(true); pedSigText->Draw();

    TLatex* speSigText = new TLatex(0.55, 0.33, Form("#sigma_{spe}: %.3f #pm %.3f", total_fit->GetParameter(5), total_fit->GetParError(5)));
    speSigText->SetTextFont(132); speSigText->SetTextSize(0.035); speSigText->SetNDC(true); speSigText->Draw();

    double redChi2 = (total_fit->GetNDF() > 0) ? total_fit->GetChisquare() / total_fit->GetNDF() : 0.0;
    TLatex* chi2Text = new TLatex(0.55, 0.29, Form("#chi^{2}/ndf: %.2f", redChi2));
    chi2Text->SetTextFont(132); chi2Text->SetTextSize(0.035); chi2Text->SetNDC(true); chi2Text->Draw();

    delete prePedFit; delete finalSigFit; delete preSigFit;
    return total_fit;
}/*}}}*/

void read_ntp_v7_thr3compare(int run, const char* processedFilePath = "") {
    gErrorIgnoreLevel = kError;   // suppress ROOT fit warnings (empty/failed fits)
    gStyle->SetTitleFont(22,""); gStyle->SetTitleSize(0.06); gStyle->SetFrameLineWidth(3); gStyle->SetLineWidth(3);

    std::string openPath; std::string savePath; 
    if (strlen(processedFilePath) > 0) {
        openPath = processedFilePath;
        TString baseName = gSystem->BaseName(openPath.c_str());
        baseName.ReplaceAll("prd", "result"); 
        baseName.ReplaceAll("raw", "result"); 
        if (!baseName.Contains("result")) baseName.Prepend("result_");
        savePath = "./Data/FinalResult_thr3compare/" + std::string(baseName.Data());
    } else { 
        openPath = find_filepath(run, PROCESSED_DATA);
        savePath = find_filepath(run, FINAL_RESULT);
    }

    int runDateTag = 0;
    {
        TString bn = gSystem->BaseName(openPath.c_str());
        TObjArray* toks = bn.Tokenize("_");
        for (int i = 0; i < toks->GetEntries(); ++i) {
            TString t = ((TObjString*)toks->At(i))->GetString();
            if (t.Length() == 8 && t.IsDigit()) { runDateTag = t.Atoi(); break; }
        }
        delete toks;
    }
    const int CABLE_SHORTEN_DATE = 20260720;
    const int SECOND_SHIFT_DATE = 20260728;
    bool useShortCable = (runDateTag == 0) || (runDateTag >= CABLE_SHORTEN_DATE);
    bool useSecondShift = (runDateTag == 0) || (runDateTag >= SECOND_SHIFT_DATE);

    TFile *file = TFile::Open(openPath.c_str(), "READ");
    if (!file || file->IsZombie()) { std::cerr << "[ERROR] Could not open processed file: " << openPath << std::endl; return; }

    int hv1=0, hv2=0, hv3=0, rot2=0, tilt2=0, rot3=0, tilt3=0;
    char sn1[50]="", sn2[50]="", sn3[50]="", runMode[20]="";
    TTree* infoTree = (TTree*)file->Get("RunInfo");
    if (infoTree) {
        infoTree->SetBranchAddress("SN1", sn1); infoTree->SetBranchAddress("SN2", sn2); infoTree->SetBranchAddress("SN3", sn3);
        infoTree->SetBranchAddress("HV1", &hv1); infoTree->SetBranchAddress("HV2", &hv2); infoTree->SetBranchAddress("HV3", &hv3);
        infoTree->SetBranchAddress("RawRotateAngle2", &rot2); infoTree->SetBranchAddress("RawTiltAngle2", &tilt2);
        infoTree->SetBranchAddress("RawRotateAngle3", &rot3); infoTree->SetBranchAddress("RawTiltAngle3", &tilt3);
        if (infoTree->GetBranch("RunMode")) infoTree->SetBranchAddress("RunMode", runMode);
        infoTree->GetEntry(0);

        std::cout << "\n=======================================================" << std::endl;
        std::cout << "[INFO] Current Run Metadata:" << std::endl;
        std::cout << "  - Run Mode    : " << runMode << std::endl;
        std::cout << "  - Mon.  : " << sn1 << " | HV: " << hv1 << " V" << std::endl;
        std::cout << "  - Rot#1 : " << sn2 << " | Angle(R,T): (" << rot2 << ", " << tilt2 << ")" << std::endl;
        std::cout << "  - Rot#2 : " << sn3 << " | Angle(R,T): (" << rot3 << ", " << tilt3 << ")" << std::endl;
        std::cout << "=======================================================\n" << std::endl;
    }

    std::vector<int> activeChannels;
    for (int i = 0; i < 8; ++i) { if (file->Get(TString::Format("tree_ch%d", i))) activeChannels.push_back(i); }
    int nAnalysisCh = 0; for (int ch : activeChannels) { if (ch != TriggerCh) nAnalysisCh++; }

    gSystem->mkdir(FinalResultPath.c_str(), kTRUE);
    TFile *outFile = TFile::Open(savePath.c_str(), "RECREATE");
    if (infoTree) { outFile->cd(); infoTree->ResetBranchAddresses(); TTree *outRunInfo = infoTree->CloneTree(); outRunInfo->Write(); }

    std::vector<TTree*> outTrees(activeChannels.size());
    TCanvas *c1 = new TCanvas(Form("c1_read_run%d", run), "Time Analysis", 1800, (nAnalysisCh > 0) ? 400 * nAnalysisCh : 800);
    TCanvas *c2 = new TCanvas(Form("c2_read_run%d", run), "Charge Analysis", (nAnalysisCh > 0) ? 1200 * nAnalysisCh : 1200, 1000);
    TCanvas *c3 = new TCanvas(Form("c3_read_run%d", run), "Pulse Height Analysis", (nAnalysisCh > 0) ? 800 * nAnalysisCh : 800, 1000);

    if(nAnalysisCh > 0) { c1->Divide(2, nAnalysisCh); c2->Divide(nAnalysisCh, 1); c3->Divide(nAnalysisCh, 1); }

    double qeDarkCorrB, qeDarkCorrErrB, ttsB, ttsErrB, speMeanB, speMeanErrorB, qeAbsB, qeAbsErrB;
    double poissonMuB, poissonQeB, poissonQeRawB, poissonMuMaxB, poissonQeMaxB, poissonQeRawMaxB;
    double chargeResolutionB, chargeResolutionErrB;
    double pvRatioB;   // peak-to-valley -- previously only drawn as TLatex text, never saved
    double darkFracB;  // dark fraction of PHC-passing counts (%) -- ditto, only printed before

    int pad_time = 1; int pad_charge = 1;
    // Record length: read from the produced file (written by prod_ntp_v7). Fall
    // back to the Time_ch histogram's bin count, then to 1024, if absent.
    int nSamples = 1024;
    if (TParameter<int>* pRL = (TParameter<int>*)file->Get("Config_RecordLength")) {
        nSamples = pRL->GetVal();
    } else {
        for (int ch : activeChannels) {
            TH1D* hT = (TH1D*)file->Get(Form("Time_ch%d", ch));
            if (hT) { nSamples = hT->GetNbinsX(); break; }
        }
    }

    // Cable-position labels matching HV_Control_SW's Dark Box naming
    // (Ch0=Mon., Ch1=Rot#1, Ch2=Rot#2), used in place of the raw channel
    // index everywhere a plot title or console summary names a PMT.
    const char* chDisplayName[] = {"Mon.", "Rot#1", "Rot#2"};

    for (size_t i = 0; i < activeChannels.size(); ++i) {
        int ch = activeChannels[i];
        TString pmtSN = (ch == 0) ? sn1 : (ch == 1) ? sn2 : (ch == 2) ? sn3 : Form("Ch%d", ch);
        TString chLabel = (ch >= 0 && ch <= 2) ? chDisplayName[ch] : Form("Ch%d", ch);
        outTrees[i] = new TTree(Form("tree_ch%d", ch), Form("Channel %d result", ch));
        outTrees[i]->Branch("rms_exG", &ttsB, "rms_exG/D");
        outTrees[i]->Branch("rms_exG_err", &ttsErrB, "rms_exG_err/D");
        outTrees[i]->Branch("spe_mean", &speMeanB, "spe_mean/D");
        outTrees[i]->Branch("spe_mean_error", &speMeanErrorB, "spe_mean_error/D");
        outTrees[i]->Branch("charge_resolution", &chargeResolutionB, "charge_resolution/D");
        outTrees[i]->Branch("charge_resolution_err", &chargeResolutionErrB, "charge_resolution_err/D");
        outTrees[i]->Branch("relativeQE", &qeDarkCorrB, "relativeQE/D");
        outTrees[i]->Branch("relativeQE_err", &qeDarkCorrErrB, "relativeQE_err/D"); 
        outTrees[i]->Branch("relativeQE_raw", &qeAbsB, "relativeQE_raw/D");
        outTrees[i]->Branch("relativeQE_raw_err", &qeAbsErrB, "relativeQE_raw_err/D"); 
        outTrees[i]->Branch("poisson_mu", &poissonMuB, "poisson_mu/D");
        outTrees[i]->Branch("poisson_qe", &poissonQeB, "poisson_qe/D"); 
        outTrees[i]->Branch("poisson_qe_raw", &poissonQeRawB, "poisson_qe_raw/D");
        outTrees[i]->Branch("poisson_mu_max", &poissonMuMaxB, "poisson_mu_max/D"); 
        outTrees[i]->Branch("poisson_qe_max", &poissonQeMaxB, "poisson_qe_max/D");
        outTrees[i]->Branch("poisson_qe_raw_max", &poissonQeRawMaxB, "poisson_qe_raw_max/D");
        outTrees[i]->Branch("peak_to_valley", &pvRatioB, "peak_to_valley/D");
        outTrees[i]->Branch("dark_frac_phc", &darkFracB, "dark_frac_phc/D");

        if (ch == TriggerCh) {
            ttsB=0; ttsErrB=0; speMeanB=0; speMeanErrorB=0; chargeResolutionB=0; chargeResolutionErrB=0; qeDarkCorrB=0; qeDarkCorrErrB=0; qeAbsB=0; qeAbsErrB=0;
            poissonMuB=0; poissonQeB=0; poissonQeRawB=0; poissonMuMaxB=0; poissonQeMaxB=0; poissonQeRawMaxB=0; pvRatioB=0; darkFracB=0;
            outFile->cd(); 
            outTrees[i]->Fill(); 
            continue; 
        }

        TTree* rawTree = (TTree*)file->Get(Form("tree_ch%d", ch)); 
        if (!rawTree) continue;
        double tPico, tMax, tDiff, tLive, tPed;
        rawTree->SetBranchAddress("pico", &tPico);
        rawTree->SetBranchAddress("max", &tMax);
        rawTree->SetBranchAddress("diff", &tDiff);
        rawTree->SetBranchAddress("LiveTime", &tLive);
        rawTree->SetBranchAddress("pedestal", &tPed);

        TParameter<double>*  pThr = (TParameter<double>*)file->Get(Form("Config_Threshold_mV_ch%d", ch));
        TParameter<double>*  pSig = (TParameter<double>*)file->Get(Form("Config_SigWindow_ns_ch%d", ch));
        TParameter<double>*  pNR  = (TParameter<double>*)file->Get(Form("NoiseCountRate_ch%d", ch));
        TParameter<Long64_t>* pNC = (TParameter<Long64_t>*)file->Get(Form("NoiseCount_ch%d", ch));
        if (!pThr || !pSig || !pNR || !pNC) {
            std::cerr << "[ERROR] ch" << ch << ": missing analysis parameters "
                      << "(Config_Threshold_mV / NoiseCountRate ...). The input is probably a RAW file; "
                      << "read_ntp_v7 needs the PRODUCED (prd) file from prod_ntp_v7. Skipping channel." << std::endl;
            continue;
        }
        double Config_Threshold_mV = pThr->GetVal();
        double Config_SigWindow_ns  = pSig->GetVal();
        double NoiseRate            = pNR->GetVal();
        double NoiseCount           = (double)pNC->GetVal();
        // Carry NoiseCountRate through to the result file -- prod_ntp_v7.C
        // writes it into the prd file but nothing here re-wrote it into the
        // result file, so downstream consumers (e.g. Draw_Stability_v1.C's
        // Dark Rate page) that only ever open result files couldn't see it.
        outFile->cd();
        (new TParameter<double>(Form("NoiseCountRate_ch%d", ch), NoiseRate))->Write();

        double TotalNoiseTime = 0.0;
        if (file->Get(Form("TotalNoiseLiveTime_ch%d", ch))) {
            TotalNoiseTime = ((TParameter<double>*)file->Get(Form("TotalNoiseLiveTime_ch%d", ch)))->GetVal();
        } else {
            TotalNoiseTime = (NoiseRate > 0) ? NoiseCount / NoiseRate : 0.0;
        }

        double tiltVal = (ch == 1) ? (double)tilt2 : ((ch == 2) ? (double)tilt3 : 0.0);
        double rotVal  = (ch == 1) ? (double)rot2  : ((ch == 2) ? (double)rot3  : 0.0);

        // Converted PMT incidence (Hamamatsu) angle + scan axis, saved by prod_ntp_v7.C.
        double hamVal = tiltVal;   
        const char* axisLabel = "?-axis";
        if (TParameter<double>* pHam = (TParameter<double>*)file->Get(Form("HamamatsuAngle_ch%d", ch)))
            hamVal = pHam->GetVal();
        if (TNamed* pAx = (TNamed*)file->Get(Form("ScanAxis_ch%d", ch)))
            axisLabel = pAx->GetTitle();
        // Common title suffix: (theta, phi | axis, Hamamatsu theta).
        TString angleInfo = Form("#theta: %.1f#circ, #phi: %.1f#circ | %s, #theta_{Ham}: %.1f#circ",
                                 tiltVal, rotVal, axisLabel, hamVal);

        // Pad 1-1: PMT Time Response
        std::string fallName = "Time_ch" + std::to_string(ch);
        TH1D *histFall = (TH1D*)file->Get(fallName.c_str());
        if (histFall) {
            c1->cd(pad_time); 
            gStyle->SetOptStat("emrou"); gPad->SetLeftMargin(0.2); gPad->SetBottomMargin(0.15); gPad->SetGrid();
            histFall->GetYaxis()->SetRangeUser(0, histFall->GetMaximum() * 1.2);
            histFall->SetTitle(Form("PMT Time Response %s - %s (%s);ADC Samples [2ns/sample];Entries", chLabel.Data(), pmtSN.Data(), angleInfo.Data()));
            HistFont(histFall); histFall->Draw("HIST");
        }

        // Pad 1-2: Timing Distribution Scan & Draw
        TH1D *histDiff = new TH1D(Form("Diff_ch%d", ch), Form("Trigger-Response %s - %s (%s);ADC Samples [2ns/sample];Entries", chLabel.Data(), pmtSN.Data(), angleInfo.Data()), nSamples, 0, nSamples);
        long long nEntries = rawTree->GetEntries();
        for (long long ev = 0; ev < nEntries; ++ev) { rawTree->GetEntry(ev); if (tDiff > -999) histDiff->Fill(tDiff); }

        c1->cd(pad_time + 1); gPad->SetLeftMargin(0.24); gPad->SetBottomMargin(0.20); gPad->SetGrid();
        histDiff->GetYaxis()->SetRangeUser(0, histDiff->GetMaximum() * 1.2); SetHistStyle(histDiff); HistFont(histDiff);
        histDiff->GetXaxis()->SetTitleSize(0.05); histDiff->GetXaxis()->SetLabelSize(0.05);
        histDiff->GetYaxis()->SetTitleSize(0.05); histDiff->GetYaxis()->SetLabelSize(0.05);
        histDiff->Draw("HIST");

        double meanFit = histDiff->GetMean(); double dChi2 = 0;
        double meanMax = histDiff->GetXaxis()->GetBinCenter(histDiff->GetMaximumBin());

        histDiff->GetXaxis()->SetRangeUser(meanMax - 20, meanMax + 10);
        double rmsWindowed = histDiff->GetRMS();
        histDiff->GetXaxis()->UnZoom();
        ttsB = -1;

        TF1* preTimeFit = new TF1(Form("pre_time_fit_ch%d",ch), "gaus", meanMax - 15, meanMax + 15); histDiff->Fit(preTimeFit, "RQN+");

        double bgSeed = 0;
        {
            int loBin = histDiff->FindBin(meanMax - 20), hiBin = histDiff->FindBin(meanMax - 15);
            int nBins = 0;
            for (int b = loBin; b <= hiBin; ++b) { bgSeed += histDiff->GetBinContent(b); nBins++; }
            if (nBins > 0) bgSeed /= nBins;
        }

        double peakHmax = histDiff->GetMaximum();
        double preMean = preTimeFit->GetParameter(1), preSig = preTimeFit->GetParameter(2), preA = preTimeFit->GetParameter(0);
        double muSeed  = (TMath::Abs(preMean - meanMax) < 5.0) ? preMean : meanMax;
        double sigSeed = (preSig > 0.3 && preSig < 5.0) ? preSig : 1.0;
        double aSeed   = (preA > 0) ? preA : peakHmax;

        bool hasPeak = (peakHmax > bgSeed + 3.0 * TMath::Sqrt(bgSeed + 1.0));

        TF1 *exGausFit = new TF1(Form("exGausFit_ch%d",ch), exGausPDF, meanMax - 20, meanMax + 10, 5);
        exGausFit->SetParameters(aSeed, muSeed, sigSeed, sigSeed, bgSeed);
        exGausFit->SetParLimits(1, meanMax - 5, meanMax + 5);   // mu: pinned near the real peak, can't drift onto the noise floor
        exGausFit->SetParLimits(2, 0.05, 10.0);   // sigma: keep physical, cannot rail
        exGausFit->SetParLimits(3, 0.01, 10.0);   // tau:   keep physical, cannot rail
        exGausFit->SetParLimits(4, 0.0, aSeed);   // bg: non-negative, can't outgrow the peak
        int fitStat = -1;
        if (hasPeak) fitStat = histDiff->Fit(exGausFit, "RQ");
        (void)fitStat;   
        double sigmaFit = exGausFit->GetParameter(2), tauFit = exGausFit->GetParameter(3), muFit = exGausFit->GetParameter(1);
        double ttsFit = TMath::Sqrt(TMath::Power(sigmaFit,2) + TMath::Power(tauFit,2));
        bool fitOK = hasPeak && ttsFit > 0 && ttsFit < 10.0 && TMath::Abs(muFit - meanMax) < 3.0;
        if (fitOK) {
            double hiParam = TMath::Max(sigmaFit, tauFit);
            if (hiParam > 3.0) fitOK = false;
        }
        if (fitOK) {
            meanFit = exGausFit->GetParameter(1); ttsB = ttsFit;
            // TTS = sqrt(sigma^2 + tau^2) -> standard error propagation, including
            // the sigma/tau covariance (they are strongly anti-correlated in this
            // model, so omitting it overstates the error substantially).
            double sigmaFitErr = exGausFit->GetParError(2), tauFitErr = exGausFit->GetParError(3);
            TVirtualFitter* lastFitter = TVirtualFitter::GetFitter();
            double sigmaTauCov = lastFitter ? lastFitter->GetCovarianceMatrixElement(2, 3) : 0;
            double ttsVariance = sigmaFit*sigmaFit*sigmaFitErr*sigmaFitErr
                                + tauFit*tauFit*tauFitErr*tauFitErr
                                + 2.0*sigmaFit*tauFit*sigmaTauCov;
            ttsErrB = (ttsVariance > 0 && ttsFit > 0) ? std::sqrt(ttsVariance) / ttsFit : 0;
            if (exGausFit->GetNDF() > 0) dChi2 = exGausFit->GetChisquare() / exGausFit->GetNDF();
        } else {
            ttsB = -1; ttsErrB = -1; dChi2 = -1;
            (void)rmsWindowed;
        }

        double tLow  = useSecondShift ? 170 : (useShortCable ? 180 : 195);
        double tHigh = useSecondShift ? 190 : (useShortCable ? 200 : 215);
        // Re-center on this run's own fitted timing peak (meanFit), same
        // width as before -- used to be ch0(Monitor)-only, leaving ch1/ch2
        // (Rot1/Rot2) on the old fixed absolute window even though the real
        // peak walks with gain (measured ~8 samples/16ns across the HV
        // campaign, HV-driven and monotonic, vs only ~3 samples of
        // angle-to-angle jitter at fixed HV -- see project memory
        // "pending-adaptive-threshold-fix"). Safely inside the fixed window
        // for every point tested so far, but margin had shrunk to ~3.5
        // samples at the high-HV end; this removes that fragility for all
        // three channels uniformly instead of just Monitor (2026-08-31,
        // user: "Timing cut은 Peak 기준으로 너비가 동일해야할 것 같은데,
        // 모든 PMT가").
        {
            double halfWidth = (tHigh - tLow) / 2.0;
            tLow  = meanFit - halfWidth;
            tHigh = meanFit + halfWidth;
        }

        double peakHeight = histDiff->GetMaximum();
        TLine *cutL = new TLine(tLow, 0, tLow, peakHeight); TLine *cutR = new TLine(tHigh, 0, tHigh, peakHeight);
        cutL->SetLineColor(kRed+2);
        cutL->SetLineWidth(3);
        cutR->SetLineColor(kRed+2);
        cutR->SetLineWidth(3);
        cutL->Draw("same");
        cutR->Draw("same");

        double fwhm = -1;
        TLegend* legTime = new TLegend(0.60, 0.40, 0.92, 0.68);
        legTime->SetBorderSize(0);
        legTime->SetFillStyle(0);
        legTime->SetTextFont(132);
        legTime->SetTextSize(0.05);

        if (fitOK) {
            double sigmaP = exGausFit->GetParameter(2), tauP = exGausFit->GetParameter(3);
            double bg_p = exGausFit->GetParameter(4);
            double peakX = exGausFit->GetMaximumX(meanFit - 5, meanFit + 20);
            double halfMax = bg_p + (exGausFit->Eval(peakX) - bg_p) / 2.0;
            double fwhmLow  = exGausFit->GetX(halfMax, peakX - 5*sigmaP - 2, peakX);
            double fwhmHigh = exGausFit->GetX(halfMax, peakX, peakX + 8*tauP + 5);
            fwhm = fwhmHigh - fwhmLow;

            histDiff->GetXaxis()->SetRangeUser(tLow - 5, tHigh + 5);

            exGausFit->SetLineColor(kRed); exGausFit->SetLineWidth(2); exGausFit->Draw("SAME");

            TLine *meanLine = new TLine(meanFit, 0, meanFit, peakHeight);
            meanLine->SetLineColor(kBlue+1); meanLine->SetLineStyle(2); meanLine->SetLineWidth(2); meanLine->Draw("same");

            TLine *fwhmLine = new TLine(fwhmLow, halfMax, fwhmHigh, halfMax);
            fwhmLine->SetLineColor(kGreen+2); fwhmLine->SetLineWidth(2); fwhmLine->Draw("same");
            TLine *fwhmTickL = new TLine(fwhmLow, 0, fwhmLow, halfMax); TLine *fwhmTickR = new TLine(fwhmHigh, 0, fwhmHigh, halfMax);
            fwhmTickL->SetLineColor(kGreen+2); fwhmTickL->SetLineStyle(3);
            fwhmTickR->SetLineColor(kGreen+2); fwhmTickR->SetLineStyle(3);
            fwhmTickL->Draw("same"); fwhmTickR->Draw("same");

            double ttsHeight = peakHeight * 0.15;
            double ttsSample = ttsB; // ttsB is in samples; *2.0 gives ns (2ns/sample)
            TLine *ttsLine = new TLine(meanFit - ttsSample / 2.0, ttsHeight, meanFit + ttsSample / 2.0, ttsHeight);
            ttsLine->SetLineColor(kMagenta+1); ttsLine->SetLineWidth(2); ttsLine->Draw("same");

            legTime->AddEntry(exGausFit, "Ex-Gaussian Fit", "l")->SetTextColor(kRed);
            legTime->AddEntry(meanLine, Form("Mean: %.2f", meanFit), "l")->SetTextColor(kBlue+1);
            legTime->AddEntry(fwhmLine, Form("FWHM: %.2f ns", fwhm * 2.0), "l")->SetTextColor(kGreen+2);
            legTime->AddEntry(ttsLine, Form("TTS: %.2f ns", ttsB * 2.0), "l")->SetTextColor(kMagenta+1);
        } else {
            histDiff->GetXaxis()->SetRangeUser(tLow - 5, tHigh + 5);   // fit failed: still show the full Timing Cut window
        }

        legTime->AddEntry(cutL, "Timing Cut", "l");
        legTime->AddEntry((TObject*)0,Form("Chi2/ndf: %.2f", dChi2),""); legTime->Draw("same");

        std::cout << "\n\033[1;36m==============================================================\033[0m" << std::endl;
        std::cout << "\033[1;36m   " << chLabel << " (" << pmtSN << ") SUMMARY\033[0m" << std::endl;
        std::cout << "\033[1;36m==============================================================\033[0m" << std::endl;
        std::cout << Form("  %-25s : %4.2f [Sample]", "Mean Position", meanFit) << std::endl;
        std::cout << Form("  %-25s : \033[1;36m%4.2f [Sample]\033[0m", "TTS (Imported)", ttsB) << std::endl;

        // Pad 2: Charge Distribution & Fit
        TH1D *histPico = new TH1D(Form("Pico_ch%d", ch), Form("Charge %s - %s (%s);Charge [pC];Entries", chLabel.Data(), pmtSN.Data(), angleInfo.Data()), 300, -5, 10);
        TH1D *histPicoAbove = new TH1D(Form("PicoAbove_ch%d", ch), "PHC Cut Only", 300, -5, 10);
        // Pad 3: Pulse Height Data Structure
        TH1D *histMax = new TH1D(Form("Max_ch%d", ch), Form("Pulse Height %s - %s (%s);ADC Counts;Entries", chLabel.Data(), pmtSN.Data(), angleInfo.Data()), 220, -20, 200);
        TH1D *histMaxAbove = new TH1D(Form("MaxAbove_ch%d", ch), "PHC Cut Active", 220, -20, 200);

        for (long long ev = 0; ev < nEntries; ++ev) {
            rawTree->GetEntry(ev); histPico->Fill(tPico); histMax->Fill(tMax);
        }

        double ADCthreshold = Config_Threshold_mV / Config::ADC_to_mV;

        double N_total = nEntries, N_phc = 0, N_all = 0;
        for (long long ev = 0; ev < nEntries; ++ev) {
            rawTree->GetEntry(ev);
            if (tMax > ADCthreshold) { histPicoAbove->Fill(tPico); N_phc++; }
            if (tMax > ADCthreshold && (tDiff > tLow && tDiff < tHigh)) N_all++;
            if (tMax > ADCthreshold) histMaxAbove->Fill(tMax);
        }

        double winSeqS = (tHigh - tLow) * 2.0 * 1e-9; 
        double kFact = (N_total * winSeqS) / TotalNoiseTime; 
        double realSig = N_all - (kFact * NoiseCount); if (realSig < 0) realSig = 0;
        double countingErr = sqrt(N_all + (kFact * kFact * NoiseCount)); 

        qeAbsB = (N_total > 0) ? (N_phc / N_total) * 100.0 : 0.0; 
        qeAbsErrB = (N_total > 0) ? (sqrt(N_phc) / N_total) * 100.0 : 0.0;
        qeDarkCorrB = (N_total > 0) ? (realSig / N_total) * 100.0 : 0.0; 
        qeDarkCorrErrB = (N_total > 0) ? (countingErr / N_total) * 100.0 : 0.0;

        // QE numbers are consolidated into a single table after the Poisson fit below.
        double expectedNoiseCounts = N_total * winSeqS * NoiseRate;

        c2->cd(pad_charge); TPad* mainPad = (TPad*)gPad;
        TPad *pad1 = new TPad(Form("pad1_ch%d",ch), "pad1", 0, 0.3, 1, 1.0); 
        pad1->SetBottomMargin(0); pad1->SetLeftMargin(0.15); 
        pad1->SetGrid(); 
        pad1->Draw(); pad1->cd(); gPad->SetLogy();

        histPico->SetLineWidth(2); 
        histPico->GetYaxis()->SetTitle("Entries"); 
        histPico->SetLineColor(kBlack);
        SetHistStyle(histPico); 
        HistFont(histPico); 
        histPico->Draw("Hist");
        histPicoAbove->SetFillColor(30);
        histPicoAbove->SetFillStyle(3004);
        //histPicoAbove->Draw("same");
        histPico->GetYaxis()->SetTitleOffset(1.2);
        gPad->SetTickx(); gPad->SetTicky();

        double signalMean = 0, signalError = 0, chargeResolution = 0, chargeResolutionErr = 0;
        double peakToValleyRatio = 0;
        TF1* total_fit_func = FitGaussian(histPico, -2., 0.9, 1.1, 2.4, signalMean, signalError,
                                           chargeResolution, chargeResolutionErr, histPicoAbove,
                                           "pC", &peakToValleyRatio);
        // Fixed x-range: default [-4, 10]; widen to [-4, 14] for high-gain PMTs (1 p.e. > 1.8 pC)
        double xChargeMax = (signalMean > 1.8) ? 14.0 : 10.0;
        histPico->GetXaxis()->SetRangeUser(-4, xChargeMax);

        poissonMuB = total_fit_func->GetParameter(1); poissonQeRawB = poissonMuB * 100.0;
        double muDarkIntegration = NoiseRate * (Config_SigWindow_ns * 1e-9); double muSigPoisson = poissonMuB - muDarkIntegration;
        poissonQeB = (muSigPoisson > 0 ? muSigPoisson : 0) * 100.0;

        // ── Consolidated QE summary (one clean table per channel) ──────────
        bool fitEmpty = (signalMean == 0 && poissonMuB == 0);
        std::cout << "\n  \033[1;33m[ Charge & QE Summary ]\033[0m" << std::endl;
        std::cout << "  --------------------------------------------------------------" << std::endl;
        std::cout << Form("    %-18s : %.0f events", "Total", N_total) << std::endl;
        double chargeRedChi2 = (total_fit_func->GetNDF() > 0)
                                 ? total_fit_func->GetChisquare() / total_fit_func->GetNDF() : 0.0;
        // Timing is always meaningful (trigger-response), shown even when SPE fit is empty.
        std::cout << Form("    %-18s : mean %.2f Sample,  TTS %.2f Sample (%.2f ns)",
                          "Timing", meanFit, ttsB, ttsB * 2.0) << std::endl;
        // Dark fraction within the PHC-passing (all-timing-cut) counting window --
        // how much of N_all is estimated to be dark hits rather than real photon
        // hits, i.e. expectedNoiseCounts / N_all. Not meaningful once N_all == 0.
        double darkFracOfPhc = (N_all > 0) ? (expectedNoiseCounts / N_all) * 100.0 : 0.0;
        darkFracB = darkFracOfPhc;
        if (fitEmpty) {
            std::cout << Form("    %-18s : \033[1;33mFit Empty\033[0m", "SPE Mean") << std::endl;
            std::cout << Form("    %-18s : \033[1;33mFit Empty\033[0m", "QE (Counting/Poi)") << std::endl;
        } else {
            std::cout << Form("    %-18s : %.3f +/- %.3f pC", "SPE Mean", signalMean, signalError) << std::endl;
            std::cout << Form("    %-18s : %.2f", "Charge Fit chi2/ndf", chargeRedChi2) << std::endl;
            std::cout << "    " << std::string(54, '-') << std::endl;
            std::cout << Form("    %-12s %12s %16s", "Method", "Raw QE", "Final (dark-sub)") << std::endl;
            std::cout << Form("    %-12s %11.3f%% %15.3f%%", "Counting", qeAbsB, qeDarkCorrB) << std::endl;
            std::cout << Form("    %-12s %11.3f%% %15.3f%%", "Poisson", poissonQeRawB, poissonQeB) << std::endl;
            std::cout << "    " << std::string(54, '-') << std::endl;
            std::cout << Form("    %-18s : %.5f  (dark %.5f -> %.5f)", "Fit Mu",
                              poissonMuB, muDarkIntegration, (muSigPoisson > 0 ? muSigPoisson : 0)) << std::endl;
            std::cout << Form("    %-18s : %.2f counts  (%.2f%% of N_all=%.0f)", "Subtracted dark", expectedNoiseCounts, darkFracOfPhc, N_all) << std::endl;
        }
        std::cout << "  --------------------------------------------------------------" << std::endl;

        // QE label on the charge plot, below the fit-parameter stack.
        //   QE_PHC  = counting (pulse-height cut) method,  Final (raw)
        // QE_Poisson dropped from this plot (charge-distribution page is PHC-only now);
        // the Poisson numbers are still computed/printed in the console table and
        // stored in the poisson_qe* branches for anyone who wants them downstream.
        pad1->cd();
        TLatex* qePhcText = new TLatex(0.55, 0.24, fitEmpty
            ? "QE_{PHC}: Fit Empty"
            : Form("QE_{PHC}: %.2f%% (raw %.2f%%)", qeDarkCorrB, qeAbsB));
        qePhcText->SetTextFont(132); qePhcText->SetTextSize(0.035);
        qePhcText->SetNDC(true); qePhcText->SetTextColor(kGreen+2); qePhcText->Draw();

        mainPad->cd();
        TPad *pad2 = new TPad(Form("pad2_ch%d",ch), "pad2", 0, 0.05, 1, 0.3); 
        pad2->SetTopMargin(0); pad2->SetBottomMargin(0.3); pad2->SetLeftMargin(0.15); 
        pad2->SetGrid(); pad2->Draw(); pad2->cd(); gPad->SetTickx(); gPad->SetTicky();
        TH1D *h_residual = (TH1D*)histPico->Clone(Form("h_residual_ch%d",ch));
        h_residual->SetTitle("");
        h_residual->Reset();
        for (int k = 1; k <= histPico->GetNbinsX(); ++k) {
            double data = histPico->GetBinContent(k); if (data == 0) continue;
            double fitVal = total_fit_func->Eval(histPico->GetBinCenter(k)); 
            h_residual->SetBinContent(k, (data - fitVal) / histPico->GetBinError(k)); h_residual->SetBinError(k, 1);
        }
        SetHistStyle(h_residual);
        h_residual->SetMarkerStyle(20);
        h_residual->SetMarkerSize(0.9);
        h_residual->GetYaxis()->SetTitle("Data-Fit/Err");
        h_residual->GetYaxis()->SetNdivisions(505);
        h_residual->GetYaxis()->SetTitleSize(0.14);
        h_residual->GetYaxis()->SetLabelSize(0.13);
        h_residual->GetYaxis()->SetTitleOffset(0.35);
        h_residual->GetYaxis()->CenterTitle();
        h_residual->GetXaxis()->SetTitle("Charge [pC]");
        h_residual->GetXaxis()->SetTitleSize(0.14);
        h_residual->GetXaxis()->SetLabelSize(0.13);
        h_residual->Draw("P"); 
        h_residual->GetXaxis()->SetRangeUser(-4, xChargeMax);
        TLine *zero_line = new TLine(-2, 0, 9, 0); 
        zero_line->SetLineColor(kRed); 
        zero_line->SetLineStyle(2); 
        zero_line->Draw("same");

        speMeanB = signalMean;
        speMeanErrorB = signalError;
        chargeResolutionB = chargeResolution;
        chargeResolutionErrB = chargeResolutionErr;
        pvRatioB = peakToValleyRatio;

        // Pad 3: Pulse Height Matrix Processing
        c3->cd(pad_charge); 
        gPad->SetLogy(); gPad->SetGrid(); histMax->SetLineWidth(2); 
        histMax->SetLineColor(kBlack); SetHistStyle(histMax); HistFont(histMax);
        histMax->Draw("Hist");
        histMaxAbove->SetFillColor(kBlue-9); histMaxAbove->SetFillStyle(3004); 
        histMaxAbove->Draw("HIST SAME");

        std::cout << "\n  \033[1;35m[ Pulse Height Poisson Fitting ]\033[0m" << std::endl;
        double maxSignalMean = 0, maxSignalError = 0, maxResolutionUnused = 0, maxResolutionErrUnused = 0;
        TF1* fitMax = FitGaussian(histMax, -2.0, 2.0, ADCthreshold, 30.0, maxSignalMean, maxSignalError,
                                   maxResolutionUnused, maxResolutionErrUnused, histMaxAbove, "ADC");
        histMax->GetXaxis()->SetRangeUser(-5, 200);

        poissonMuMaxB = fitMax->GetParameter(1); poissonQeRawMaxB = poissonMuMaxB * 100.0; double muSigPoissonMax = poissonMuMaxB - muDarkIntegration;
        poissonQeMaxB = (muSigPoissonMax > 0 ? muSigPoissonMax : 0) * 100.0;

        std::cout << Form("  PH Fit Mu (Raw)       : %6.5f", poissonMuMaxB) << std::endl;
        std::cout << Form("  PH Raw QE (No Dark)   : %4.3f %%", poissonQeRawMaxB) << std::endl;
        std::cout << Form("  PH Corrected Mu       : %6.5f", (muSigPoissonMax > 0 ? muSigPoissonMax : 0)) << std::endl;
        std::cout << Form("  PH Poisson P(0) QE    : \033[1;35m%4.3f %%\033[0m", poissonQeMaxB) << std::endl;
        double qeDiffMax = (poissonQeB > 0) ? ((poissonQeMaxB - poissonQeB) / poissonQeB) * 100.0 : 0;
        std::cout << Form("  Diff vs Charge QE     : %+4.2f %%", qeDiffMax) << std::endl;
        std::cout << "  --------------------------------------------------------------" << std::endl;

        outFile->cd(); histDiff->Write(); histPico->Write(); histMax->Write(); outTrees[i]->Fill(); pad_time += 2; pad_charge++;
    }

    outFile->cd(); for (auto *tree : outTrees) if(tree) tree->Write();
    std::string imageDir = ImagePath + "ByAnalysis/"; gSystem->mkdir(imageDir.c_str(), kTRUE); 
    std::string baseName = gSystem->BaseName(savePath.c_str()); baseName.replace(baseName.find(".root"), 5, ""); 

    //    TString cutSuffix = Form("_cut_%.1f_%.1f", tLowOffset, tHighOffset);
    c1->SaveAs((imageDir + baseName/* + cutSuffix*/ + "_time.png").c_str());
    c2->SaveAs((imageDir + baseName/* + cutSuffix*/ + "_charge.png").c_str());
    //c3->SaveAs((imageDir + baseName/* + cutSuffix*/ + "_pulseheight.png").c_str());
    outFile->Close(); file->Close();
}











