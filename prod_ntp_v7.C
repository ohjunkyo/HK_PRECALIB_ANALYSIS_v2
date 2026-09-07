#include <vector>
#include <string>
#include <iostream>
#include <cctype>
#include <algorithm>
#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLine.h>
#include <TMath.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TNamed.h>
#include <TLatex.h>
#include <TLegend.h>

#include "./Base/analysisCode/DrawFormat.cpp"
#include "./Base/analysisCode/Analysis.cpp"
#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"
#include "path_builder2.h"
#include "angle_convert.h"

// Shared threshold: amplitude cut, timing crossing and dark count all use it.
const double Analysis_Threshold_mV = 1.5;

Double_t exGausPDF(Double_t *x, Double_t *par) {/*{{{*/
    Double_t xx = x[0]; Double_t A = par[0]; Double_t mu = par[1]; Double_t sigma = par[2]; Double_t tau = par[3];
    Double_t invTau = 1.0 / tau;
    Double_t z = (sigma * invTau - (xx - mu) / sigma) / TMath::Sqrt(2.0);
    Double_t expo = (sigma*sigma*invTau*invTau/2.0) - (xx - mu)*invTau;
    return A * invTau/2.0 * TMath::Exp(expo) * TMath::Erfc(z);
}/*}}}*/
TF1 *PreFit(TH1D* hist, double histMin, double histMax, const std::string& name) {/*{{{*/
    int binMin = hist->GetXaxis()->FindBin(histMin); int binMax = hist->GetXaxis()->FindBin(histMax);
    double localMaxContent = -1.0; int localMaxBin = -1;
    for (int b = binMin; b <= binMax; ++b) {
        if (hist->GetBinContent(b) > localMaxContent) { localMaxContent = hist->GetBinContent(b); localMaxBin = b; }
    }
    TF1* preFit = new TF1(name.c_str(), "gaus", histMin, histMax);
    preFit->SetParameter(0, localMaxContent); preFit->SetParameter(1, hist->GetBinCenter(localMaxBin));
    hist->Fit(preFit, "RQMN+"); return preFit;
}/*}}}*/
void HistFont(TH1D *hist){ /*{{{*/
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
}/*}}}*/

void prod_ntp_v7(int run, const char* rawFilePath = "") {
    gErrorIgnoreLevel = kWarning;
    gStyle->SetTitleFont(22,""); gStyle->SetTitleSize(0.08); gStyle->SetFrameLineWidth(3);

    std::string OpenPath;
    std::string outName;
    bool darkmode = false;
    char runMode[100]="Unknown";

    if (strlen(rawFilePath) > 0) {
        OpenPath = rawFilePath;
        TString baseName = gSystem->BaseName(OpenPath.c_str());
        baseName.ReplaceAll("raw", "prd");
        if (!baseName.Contains("prd")) baseName.Prepend("prd_");
        outName = std::string(ProcessedDataPath) + "/" + baseName.Data();
    } else {
        OpenPath = find_filepath(run, RAW_DATA);
        outName = find_filepath(run, PROCESSED_DATA);
    }

    int runDateTag = 0;
    {
        TString bn = gSystem->BaseName(OpenPath.c_str());
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

    TFile *file = TFile::Open(OpenPath.c_str(), "READ");
    if (!file || file->IsZombie()) { std::cerr << "[ERROR] Could not open raw data file: " << OpenPath << std::endl; return; }

    TTree* infoTree = (TTree*)file->Get("RunInfo");
    if (infoTree) {
        infoTree->SetBranchStatus("*", 1);
        if (infoTree->GetBranch("RunMode")) infoTree->SetBranchAddress("RunMode", runMode);
        infoTree->GetEntry(0);
        if (TString(runMode) == "Dark" || TString(runMode) == "DARK") darkmode = true;
    }

    TTree *T = (TTree*)file->Get("T");
    long long NEntry = T->GetEntries();
    unsigned int RecordLength = 1024, PostTrigger = 60, TriggerTimeTag = 0;

    bool constsOnT = (T->GetBranch("RecordLength") != nullptr);

    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("TriggerTimeTag", 1);
    T->SetBranchAddress("TriggerTimeTag", &TriggerTimeTag);
    if (constsOnT) {
        T->SetBranchStatus("RecordLength", 1); T->SetBranchStatus("PostTrigger", 1);
        T->SetBranchAddress("RecordLength", &RecordLength);
        T->SetBranchAddress("PostTrigger", &PostTrigger);
    } else if (infoTree && infoTree->GetBranch("RecordLength")) {
        infoTree->SetBranchAddress("RecordLength", &RecordLength);
        infoTree->SetBranchAddress("PostTrigger", &PostTrigger);
        infoTree->GetEntry(0);
        infoTree->ResetBranchAddresses();   
    }
    T->GetEntry(0);

    int nSamples = RecordLength;

    UShort_t* ADC_buffer = new UShort_t[8 * nSamples];
    if (T->GetBranch("ADC")) {
        T->SetBranchStatus("ADC", 1);
        T->SetBranchAddress("ADC", ADC_buffer);
    } else {
        std::cerr << "[ERROR] 'ADC' branch not found!" << std::endl;
        delete[] ADC_buffer; return;
    }

    std::vector<int> activeChannels;
    for (int i = 0; i < 8; ++i) {
        if (T->GetBranch(Form("OffsetValue%d", i))) activeChannels.push_back(i);
    }
    if (activeChannels.empty() && infoTree) {
        if (infoTree->GetBranch("ChannelMask")) {
            char chMask[16] = {0};
            infoTree->SetBranchAddress("ChannelMask", chMask);
            infoTree->GetEntry(0);
            infoTree->ResetBranchAddresses();
            // The mask string is written MSB-first (leftmost char = ch7).
            int len = (int)strlen(chMask);
            for (int i = 0; i < len && i < 8; ++i)
                if (chMask[len - 1 - i] == '1') activeChannels.push_back(i);
        }
        if (activeChannels.empty()) {          // fall back to the moved branches
            for (int i = 0; i < 8; ++i)
                if (infoTree->GetBranch(Form("OffsetValue%d", i))) activeChannels.push_back(i);
        }
    }
    int channelCount = activeChannels.size();
    if (channelCount == 0) { delete[] ADC_buffer; return; }

    int trgPoint = (int)(nSamples * (1.0 - (PostTrigger / 100.0)));

    int darkSigCenter = trgPoint;
    if (darkmode) {
        long long preScanEventCount = std::min((long long)5000, NEntry);
        std::vector<double> avgWaveform(nSamples, 0.0);
        for (long long preScanEvt = 0; preScanEvt < preScanEventCount; ++preScanEvt) {
            T->GetEntry(preScanEvt);
            for (int sample = 0; sample < nSamples; ++sample) {
                avgWaveform[sample] += (double)ADC_buffer[TriggerCh * nSamples + sample];
            }
        }

        int troughSample = -1;
        double troughAdc = 1e18;
        for (int sample = 0; sample < nSamples; ++sample) {
            double avgAdc = avgWaveform[sample] / preScanEventCount;
            if (avgAdc < troughAdc) {
                troughAdc = avgAdc;
                troughSample = sample;
            }
        }

        double avgPedestalAdc = 0;
        for (int sample = 10; sample < 150; ++sample) avgPedestalAdc += avgWaveform[sample] / preScanEventCount;
        avgPedestalAdc /= 140.0;
        double troughDepthMv = (avgPedestalAdc - troughAdc) * Config::ADC_to_mV;

        if (troughSample > 0 && troughDepthMv > 5.0) {
            darkSigCenter = troughSample;
            std::cout << Form("[INFO] Dark mode: auto-detected self-trigger pulse at sample %d (depth %.1f mV, avg of %lld events)",
                    troughSample, troughDepthMv, preScanEventCount) << std::endl;
        } else {
            std::cout << Form("[WARNING] Dark mode: could not auto-detect pulse position (trough=%d, depth=%.1f mV) "
                    "-- falling back to PostTrigger-derived trgPoint=%d", troughSample, troughDepthMv, trgPoint) << std::endl;
        }
    }

    // Laser mode: signal is fixed at 595, so a wide pedestal baseline (50-550) is safe.
    // Dark mode: pedestal stays narrow to avoid overlapping the signal window.
    int pedStart = darkmode ? 100 : 50;
    int pedEnd   = darkmode ? 300 : 550;
    int sigStart = darkmode ? std::max(0, darkSigCenter - 15) : (useSecondShift ? 580 : (useShortCable ? 595 : 610));
    int sigEnd   = (sigStart + 50 > nSamples) ? nSamples : sigStart + 50;

    int noiseStart = 50;
    int noiseEnd = darkmode ? std::max(noiseStart + 10, sigStart - 20) : (trgPoint - 50);
    double Config_DarkWindow_ns = (double)(noiseEnd - noiseStart) * 2.0;
    double Config_SigWindow_ns  = (double)(sigEnd - sigStart) * 2.0;

    TFile *fOut = new TFile(outName.c_str(), "RECREATE");
    if (infoTree) {
        fOut->cd(); infoTree->ResetBranchAddresses(); TTree* cloneInfo = infoTree->CloneTree(); cloneInfo->Write();
    }

    char sn1[100], sn2[100], sn3[100];
    char dir1[8] = "", dir2[8] = "", dir3[8] = "";
    int tilt2 = 0, rot2 = 0;
    int tilt3 = 0, rot3 = 0;

    TTree *runInfo = (TTree*)file->Get("RunInfo");
    if (runInfo) {

        if (runInfo->GetBranch("SN1"))   runInfo->SetBranchAddress("SN1", &sn1);
        if (runInfo->GetBranch("SN2"))   runInfo->SetBranchAddress("SN2", &sn2);
        if (runInfo->GetBranch("SN3"))   runInfo->SetBranchAddress("SN3", &sn3);
        if (runInfo->GetBranch("Direction1")) runInfo->SetBranchAddress("Direction1", &dir1);
        if (runInfo->GetBranch("Direction2")) runInfo->SetBranchAddress("Direction2", &dir2);
        if (runInfo->GetBranch("Direction3")) runInfo->SetBranchAddress("Direction3", &dir3);
        if (runInfo->GetBranch("RawTiltAngle2"))   runInfo->SetBranchAddress("RawTiltAngle2", &tilt2);
        if (runInfo->GetBranch("RawRotateAngle2")) runInfo->SetBranchAddress("RawRotateAngle2", &rot2);
        if (runInfo->GetBranch("RawTiltAngle3"))   runInfo->SetBranchAddress("RawTiltAngle3", &tilt3);
        if (runInfo->GetBranch("RawRotateAngle3")) runInfo->SetBranchAddress("RawRotateAngle3", &rot3);
        runInfo->GetEntry(0);
    }

    auto dirOf = [&](int ch) -> char {
        const char* d = (ch == 0) ? dir1 : (ch == 1) ? dir2 : dir3;
        return (d && d[0] != '\0') ? d[0] : DirForCh(ch);
    };

    std::vector<TTree*> outTrees(channelCount);
    std::vector<TH1D*> histPed(channelCount), histMax(channelCount), histFall(channelCount), histDiff(channelCount), histPico(channelCount);
    std::vector<TH1D*> histPicoAbove(channelCount), histPicoAboveTCut(channelCount), histPicoAboveAllCut(channelCount);
    std::vector<std::vector<double>> maxVec(channelCount), picoVec(channelCount), diffVec(channelCount);
    // std::vector<TH2D*> histQT(channelCount);   // Charge vs Time -- replaced by Pedestal Stability below
    std::vector<TH1D*> histNoiseTime(channelCount);
    std::vector<TH2D*> histNoiseMT(channelCount);
    std::vector<TH1D*> histNoiseCharge(channelCount);   // Charge of each dark/noise pulse found in the noise-search window
    std::vector<TH1D*> histDarkWinCharge(channelCount);   // One integration per event over the WHOLE dark window (noiseStart..noiseEnd), like histPico is for sigStart..sigEnd
    std::vector<TProfile*> profPedStab(channelCount);   // Pedestal mean+RMS vs event index -- baseline stability over the run


    std::vector<long long> noiseHits(channelCount, 0LL);

    std::vector<long long> selfTrigCount(channelCount, 0LL);

    std::vector<long long> sigWinDarkHits(channelCount, 0LL);
    const double ADCthreshold_amp = Analysis_Threshold_mV / Config::ADC_to_mV;
    double LiveTime = 0.0, totalNoiseLiveTime = 0.0;
    double LiveTime_b, pedestalB, falltimeB, maxB, diffB, picoB, darkWinChargeB;


    for (int i = 0; i < channelCount; ++i) {
        int ch = activeChannels[i];
        TString pmtSN;

        if (ch == 0)      pmtSN = sn1;
        else if (ch == 1) pmtSN = sn2;
        else if (ch == 2) pmtSN = sn3;
        else              pmtSN = Form("CH%d", ch); 
        std::cout << Form("  > PMT %s:", pmtSN.Data()) << std::endl;

        double tiltVal = 0.0; double rotVal = 0.0;
        if (ch == 1) { tiltVal = (double)tilt2; rotVal = (double)rot2; }
        else if (ch == 2) { tiltVal = (double)tilt3; rotVal = (double)rot3; }

        // Convert rotation-stage angle to PMT incidence (Hamamatsu) angle, and
        // tag the scan axis. Reused as a common title suffix for every histogram.
        const char* axisLabel = "?-axis";
        double hamVal = GetHamamatsuAngle(dirOf(ch), tiltVal, rotVal, axisLabel);
        TString angleInfo = Form("#theta: %.1f#circ, #phi: %.1f#circ | %s, #theta_{Ham}: %.1f#circ",
                tiltVal, rotVal, axisLabel, hamVal);

        outTrees[i] = new TTree(Form("tree_ch%d", ch), Form("CH%d result", ch));
        outTrees[i]->Branch("LiveTime", &LiveTime_b, "LiveTime/D");
        outTrees[i]->Branch("pico", &picoB, "pico/D");
        outTrees[i]->Branch("darkWinCharge", &darkWinChargeB, "darkWinCharge/D");
        outTrees[i]->Branch("falltime", &falltimeB, "falltime/D");
        outTrees[i]->Branch("max", &maxB, "max/D");
        outTrees[i]->Branch("diff", &diffB, "diff/D");
        outTrees[i]->Branch("pedestal", &pedestalB, "pedestal/D");

        histPed[i]  = new TH1D(Form("Ped_ch%d", ch), Form("%s Pedestal (%s);ADC Counts;Entries", pmtSN.Data(), angleInfo.Data()), 500, 14000, 16000);
        histMax[i]  = new TH1D(Form("Max_ch%d", ch), Form("%s Pulse Height (%s);ADC Counts;Entries", pmtSN.Data(), angleInfo.Data()), 220, -20, 200);
        histFall[i] = new TH1D(Form("Time_ch%d", ch), Form("%s PMT Time Response (%s);ADC Samples [2ns/sample];Entries", pmtSN.Data(), angleInfo.Data()), nSamples, 0, nSamples);
        histDiff[i] = new TH1D(Form("Diff_ch%d", ch), Form("%s Trigger-Response (%s);ADC Samples [2ns/sample];Entries", pmtSN.Data(), angleInfo.Data()), nSamples, 0, nSamples);
        histPico[i] = new TH1D(Form("Pico_ch%d", ch), Form("%s Charge (%s);Charge [pC];Entries", pmtSN.Data(), angleInfo.Data()), 300, -5, 10);
        // histQT[i] = new TH2D(Form("QT_ch%d", ch), Form("%s Charge vs Time; Timing [ADC Samples]; Charge [pC]", pmtSN.Data()), nSamples, 0, nSamples, 300, -5, 10);
        histNoiseTime[i] = new TH1D(Form("NoiseTime_ch%d", ch), Form("%s Noise Arrival Time Distribution (%s);ADC Samples [2ns/sample];Entries", pmtSN.Data(), angleInfo.Data()), nSamples, 0, nSamples);
        histNoiseMT[i] = new TH2D(Form("NoiseMT_ch%d", ch), Form("%s Noise Amplitude vs Time (%s);Time [ADC Samples];Pulse Height [ADC Counts]", pmtSN.Data(), angleInfo.Data()), nSamples, 0, nSamples, 200, 0, 200);
        histNoiseCharge[i] = new TH1D(Form("NoiseCharge_ch%d", ch), Form("%s Noise/Dark Pulse Charge (%s);Charge [pC];Entries", pmtSN.Data(), angleInfo.Data()), 300, -5, 10);
        histDarkWinCharge[i] = new TH1D(Form("DarkWinCharge_ch%d", ch), Form("%s Dark Window Integrated Charge (%s);Charge [pC];Entries", pmtSN.Data(), angleInfo.Data()), 300, -5, 10);

        profPedStab[i] = new TProfile(Form("PedStability_ch%d", ch),
                Form("%s Pedestal Stability (%s);Event Index;Pedestal [ADC Counts]", pmtSN.Data(), angleInfo.Data()),
                200, 0, NEntry);
        profPedStab[i]->SetErrorOption("s");   // error bar = RMS per bin (baseline noise), not error-on-mean
    }

    unsigned int TriggerTimeTag_prev = 0;
    std::vector<unsigned int> ADC(nSamples, 0);

    std::cout << "===== Starting Event Loop =====" << std::endl;
    for (int iEntry = 0; iEntry < NEntry; ++iEntry) {
        T->GetEntry(iEntry);
        if (iEntry % 10000 == 0) std::cout << "Processing... " << (iEntry*100/NEntry) << "%\r" << std::flush;

        totalNoiseLiveTime += Config_DarkWindow_ns * 1e-9;
        if (iEntry > 0) {
            double diffTag = (double)TriggerTimeTag - (double)TriggerTimeTag_prev;
            if (diffTag < 0) diffTag += TMath::Power(2, 31);
            LiveTime += diffTag * 8.0 * 1e-9;
        }
        TriggerTimeTag_prev = TriggerTimeTag;

        for(int s=0; s<nSamples; ++s) ADC[s] = ADC_buffer[TriggerCh * nSamples + s];
        int triggerTime = GetTimeAtThreshold(ADC.data(), 14000.0, nSamples);

        for (int i = 0; i < channelCount; ++i) {
            int ch = activeChannels[i];
            for(int s=0; s<nSamples; ++s) ADC[s] = ADC_buffer[ch * nSamples + s];

            pedestalB = GetPedestal(ADC.data(), pedStart, pedEnd, nSamples);
            double thresholdB = pedestalB - (Analysis_Threshold_mV / Config::ADC_to_mV);

            // Dark hits counted over the same window totalNoiseLiveTime covers.
            std::vector<int> darkPulseTimes = GetTimesBelowThreshold(ADC.data(), thresholdB, nSamples, noiseStart, noiseEnd);
            if (ch != TriggerCh && darkPulseTimes.size() > 0) {
                noiseHits[i] += (long long)darkPulseTimes.size();
                for (int tNoise : darkPulseTimes) {
                    histNoiseTime[i]->Fill(tNoise);
                    double noiseAmp = pedestalB - static_cast<double>(ADC[tNoise]);
                    histNoiseMT[i]->Fill(tNoise, noiseAmp);

                    int noiseChargeStart = tNoise - 10;
                    int noiseChargeEnd = tNoise + 40;
                    if (noiseChargeStart < 0) noiseChargeStart = 0;
                    if (noiseChargeEnd > nSamples) noiseChargeEnd = nSamples;
                    double noisePico = GetCharge(ADC.data(), pedestalB, noiseChargeStart, noiseChargeEnd, nSamples);
                    histNoiseCharge[i]->Fill(noisePico);
                }
            }
            /*
               if (ch != TriggerCh && darkPulseTimes.size() > 0) {
               totalNoiseHits += darkPulseTimes.size();

               for (int tNoise : darkPulseTimes) {
               histNoiseTime[i]->Fill(tNoise);

               double noiseAmp = pedestalB - static_cast<double>(ADC[tNoise]);
               histNoiseMT[i]->Fill(tNoise, noiseAmp);
               }
               }*/

            falltimeB = GetTimeAtThreshold(ADC.data(), thresholdB, nSamples);
            diffB = (falltimeB >= 0) ? falltimeB - triggerTime : -1000;

            auto peak = GetPeakInfo(ADC.data(), pedestalB, sigStart, sigEnd, nSamples);
            maxB = peak.first;
            picoB = GetCharge(ADC.data(), pedestalB, sigStart, sigEnd, nSamples);

            darkWinChargeB = GetCharge(ADC.data(), pedestalB, noiseStart, noiseEnd, nSamples);

            if (darkmode && ch != TriggerCh && maxB > ADCthreshold_amp) selfTrigCount[i]++;

            if (darkmode && ch != TriggerCh) {
                sigWinDarkHits[i] += (long long)GetTimesBelowThreshold(ADC.data(), thresholdB, nSamples, sigStart, sigEnd).size();
            }

            // if (falltimeB >= 0) {
            //     int dynStart = falltimeB - 10; int dynEnd = falltimeB + 40;
            //     if (dynStart < 0) dynStart = 0; if (dynEnd > nSamples) dynEnd = nSamples;
            //     double dynamicPico = GetCharge(ADC.data(), pedestalB, dynStart, dynEnd, nSamples);
            //     histQT[i]->Fill(falltimeB, dynamicPico);
            // }

            if (pedestalB > 0) { histPed[i]->Fill(pedestalB); profPedStab[i]->Fill((double)iEntry, pedestalB); }
            histMax[i]->Fill(maxB);
            if (falltimeB >= 0) histFall[i]->Fill(falltimeB);
            if (diffB > -999) histDiff[i]->Fill(diffB);
            histPico[i]->Fill(picoB);
            histDarkWinCharge[i]->Fill(darkWinChargeB);

            maxVec[i].push_back(maxB); picoVec[i].push_back(picoB); diffVec[i].push_back(diffB);
            LiveTime_b = LiveTime;

            outTrees[i]->Fill();
        }
    }

    // TCanvas *cQT = new TCanvas("cQT", "Q-T Graphs", 1800, 600);
    // cQT->Divide(channelCount, 1);
    // for (int i = 0; i < channelCount; ++i) {
    //     cQT->cd(i + 1); gPad->SetLogz(); histQT[i]->Draw("COLZ");
    // }


    std::vector<double> pedMeanOut(channelCount, 0.0), pedSigmaOut(channelCount, 0.0);
    TCanvas *cPedStab = new TCanvas("cPedStab", "Pedestal Stability", 1800, 600);
    cPedStab->Divide(channelCount, 1);
    for (int i = 0; i < channelCount; ++i) {
        cPedStab->cd(i + 1); gPad->SetGrid();
        profPedStab[i]->SetMarkerStyle(20); profPedStab[i]->SetMarkerSize(0.6);
        profPedStab[i]->SetLineColor(kBlue + 1); profPedStab[i]->SetMarkerColor(kBlue + 1);
        profPedStab[i]->Draw("E1");

        if (histPed[i] && histPed[i]->GetEntries() > 0) {
            double coreMean = histPed[i]->GetMean();
            int b0 = histPed[i]->FindBin(coreMean - 100), b1 = histPed[i]->FindBin(coreMean + 50);
            double sumW = 0, sumWX = 0, sumWX2 = 0;
            for (int b = b0; b <= b1; ++b) {
                double x = histPed[i]->GetBinCenter(b), w = histPed[i]->GetBinContent(b);
                sumW += w; sumWX += w*x; sumWX2 += w*x*x;
            }
            double meanAdc = coreMean, sigmaAdc = 0;
            if (sumW > 0) { meanAdc = sumWX/sumW; sigmaAdc = TMath::Sqrt(TMath::Max(0.0, sumWX2/sumW - meanAdc*meanAdc)); }
            pedMeanOut[i] = meanAdc * Config::ADC_to_mV;
            pedSigmaOut[i] = sigmaAdc * Config::ADC_to_mV;
            TLatex *tSigma = new TLatex(0.15, 0.85, Form("#sigma = %.3f mV", pedSigmaOut[i]));
            tSigma->SetNDC(); tSigma->SetTextFont(132); tSigma->SetTextSize(0.06);
            tSigma->SetTextColor(kRed + 1); tSigma->Draw();
        }
    }

    // NoiseMap canvas removed; histNoise* still written to the output file.

    TCanvas *c1 = new TCanvas(Form("c1_prod_run%d", run), Form("Production for Run %d", run), 3200, 2000);
    int drawPads = 0;
    for (int ch : activeChannels) if (ch != TriggerCh) drawPads++;
    if (drawPads > 0) c1->Divide(5, drawPads);

    std::cout << "==============================================================" << std::endl;
    std::cout << "  - Raw File Path   : " << OpenPath << std::endl;
    std::cout << "  - Processed Path  : " << outName << std::endl;
    std::cout << "  - Run Mode        : " << runMode << std::endl;
    std::cout << "  - Total Events    : " << NEntry << " entries" << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;

    // Dark mode: rate = selfTrigCount/LiveTime (every event is a real trigger).
    // Laser mode: windowed count, since the laser sets the trigger rate.
    std::vector<long long> noiseCountOut(channelCount, 0LL);
    std::vector<double> noiseRateOut(channelCount, 0.0);
    std::vector<double> noiseLiveTimeOut(channelCount, 0.0);
    for (int i = 0; i < channelCount; ++i) {
        if (darkmode) {
            noiseCountOut[i] = selfTrigCount[i];
            noiseRateOut[i] = (LiveTime > 0) ? (double)selfTrigCount[i] / LiveTime : 0.0;
            noiseLiveTimeOut[i] = LiveTime;
        } else {
            noiseCountOut[i] = noiseHits[i];
            noiseRateOut[i] = (totalNoiseLiveTime > 0) ? (double)noiseHits[i] / totalNoiseLiveTime : 0.0;
            noiseLiveTimeOut[i] = totalNoiseLiveTime;
        }
    }

    // Dark rate measured inside the signal window. Dark mode only.
    std::vector<double> sigWinDarkRateOut(channelCount, 0.0);
    if (darkmode) {
        for (int i = 0; i < channelCount; ++i)
            sigWinDarkRateOut[i] = (LiveTime > 0) ? (double)sigWinDarkHits[i] / LiveTime : 0.0;
    }
    if (darkmode) std::cout << Form("  - Self-trig LiveTime : %.3f s (%.0f Hz combined trigger rate)",
            LiveTime, LiveTime > 0 ? NEntry / LiveTime : 0.0) << std::endl;

    int padIdx = 1;
    for (int i = 0; i < channelCount; ++i) {
        int ch = activeChannels[i]; if(ch == TriggerCh) continue;
        double ADCthreshold = Analysis_Threshold_mV / Config::ADC_to_mV;

        double tiltVal = (ch == 1) ? (double)tilt2 : ((ch == 2) ? (double)tilt3 : 0.0);
        double rotVal  = (ch == 1) ? (double)rot2  : ((ch == 2) ? (double)rot3  : 0.0);
        const char* axisLabel = "?-axis";
        double hamVal = GetHamamatsuAngle(dirOf(ch), tiltVal, rotVal, axisLabel);

        // Pad 1: Pedestal
        c1->cd(padIdx++); gPad->SetGrid(); gPad->SetLogy();
        histPed[i]->GetXaxis()->SetNdivisions(510); SetHistStyle(histPed[i]); HistFont(histPed[i]); 
        histPed[i]->Draw("HIST");
        if (histPed[i]->GetEntries() > 0) histPed[i]->GetXaxis()->SetRangeUser(histPed[i]->GetMean() - 100, histPed[i]->GetMean() + 50);

        // Pad 2: Pulse Height
        c1->cd(padIdx++); gPad->SetGrid(); gPad->SetLogy();
        SetHistStyle(histMax[i]);  HistFont(histMax[i]); 
        histMax[i]->Draw("HIST");
        TF1* Gaus = new TF1(Form("Gaus_ch%d_run%d", ch, run), "gaus", -2, 30);
        histMax[i]->Fit(Gaus, "RMQ"); 
        Gaus->Draw("same");
        TLine *lineThr = new TLine(ADCthreshold, gPad->GetUymin(), ADCthreshold, histMax[i]->GetMaximum());
        lineThr->SetLineColor(kRed); lineThr->SetLineStyle(2); 
        lineThr->Draw("same");

        // Pad 3: PMT Response
        c1->cd(padIdx++); SetHistStyle(histFall[i]); HistFont(histFall[i]); histFall[i]->Draw("HIST");
        double FallMax = histFall[i]->GetMaximum();
        int lines[] = {pedStart, pedEnd, sigStart, sigEnd}; 
        int colors[] = {kBlue, kBlue, kRed, kRed};
        for(int l=0; l<4; ++l) { 
            TLine *line = new TLine(lines[l], 0, lines[l], FallMax); 
            line->SetLineColor(colors[l]); 
            line->Draw("same"); 
        }

        // Pad 4: Timing difference Monitoring Fit
        c1->cd(padIdx++); 
        SetHistStyle(histDiff[i]); HistFont(histDiff[i]); histDiff[i]->Draw("HIST");
        double meanT = histDiff[i]->GetXaxis()->GetBinCenter(histDiff[i]->GetMaximumBin());

        TF1* preTimeFit = new TF1(Form("pre_time_fit_ch%d",ch), "gaus", meanT - 20, meanT + 20); histDiff[i]->Fit(preTimeFit, "RQ+N");
        TF1 *exGausFit = new TF1(Form("exGausFit_ch%d",ch), exGausPDF, meanT - 20, meanT + 50, 4);
        exGausFit->SetParameters(preTimeFit->GetParameter(0), preTimeFit->GetParameter(1), preTimeFit->GetParameter(2), preTimeFit->GetParameter(2));
        histDiff[i]->Fit(exGausFit, "RQ"); 
        exGausFit->SetLineColor(kRed); exGausFit->SetLineWidth(3); 
        exGausFit->Draw("SAME");

        double tLow  = useSecondShift ? 170 : (useShortCable ? 180 : 195);
        double tHigh = useSecondShift ? 190 : (useShortCable ? 200 : 215);
        TLine *L_cut = new TLine(tLow, 0, tLow, histDiff[i]->GetMaximum()); TLine *R_cut = new TLine(tHigh, 0, tHigh, histDiff[i]->GetMaximum());
        L_cut->SetLineColor(kGreen+2); L_cut->SetLineWidth(3); R_cut->SetLineColor(kGreen+2); R_cut->SetLineWidth(3);
        L_cut->Draw("same"); R_cut->Draw("same"); histDiff[i]->GetXaxis()->SetRangeUser(meanT - 20, meanT + 20);

        double tSigma = exGausFit->GetParameter(2);
        double tTau = exGausFit->GetParameter(3);
        double rms_exG = TMath::Sqrt(tSigma * tSigma + tTau * tTau);

        // Pad 5: Charge distribution Monitoring Fit
        c1->cd(padIdx++); gPad->SetLogy(); gPad->SetGrid();
        SetHistStyle(histPico[i]); HistFont(histPico[i]); 
        histPico[i]->Draw("HIST");
        histPicoAbove[i] = (TH1D*)histPico[i]->Clone(); histPicoAbove[i]->Reset();
        for (size_t idx = 0; idx < maxVec[i].size(); ++idx) { 
            if (maxVec[i][idx] > ADCthreshold) histPicoAbove[i]->Fill(picoVec[i][idx]); 
        }
        histPicoAbove[i]->SetLineColor(kGray+1); histPicoAbove[i]->SetFillColor(kGray); histPicoAbove[i]->SetFillStyle(3001);
        histPicoAbove[i]->Draw("HIST SAME");

        // Dark-pulse charge overlay -- COMMENTED OUT 2026-08-25 

        // TH1D* histNoiseChargeScaled = nullptr;
        // if (histNoiseCharge[i]->GetEntries() > 0 && histNoiseCharge[i]->GetMaximum() > 0) {
        //     histNoiseChargeScaled = (TH1D*)histNoiseCharge[i]->Clone(Form("NoiseChargeScaled_ch%d", ch));
        //     double scale = histPico[i]->GetMaximum() / histNoiseCharge[i]->GetMaximum();
        //     histNoiseChargeScaled->Scale(scale);
        //     histNoiseChargeScaled->SetLineColor(kRed + 1); histNoiseChargeScaled->SetLineWidth(2); histNoiseChargeScaled->SetLineStyle(2);
        //     histNoiseChargeScaled->Draw("HIST SAME");
        //     TLegend* legCharge = new TLegend(0.40, 0.78, 0.98, 0.90);
        //     legCharge->SetBorderSize(0); legCharge->SetFillStyle(0); legCharge->SetTextSize(0.028);
        //     legCharge->AddEntry(histPico[i], "Signal", "l");
        //     legCharge->AddEntry(histNoiseChargeScaled, "Dark (scaled)", "l");
        //     legCharge->Draw("same");
        // }

        bool fitEmpty = (histPicoAbove[i]->GetEntries() < 1);
        double Mean = 0, Err = 0;
        if (!fitEmpty) {
            TF1 *PreSigFit = PreFit(histPicoAbove[i], -5, 10, Form("presig_ch%d", ch));
            Mean = PreSigFit->GetParameter(1); Err = PreSigFit->GetParError(1);
        }
        TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(0.05);
        tex->DrawLatex(0.55, 0.55, fitEmpty ? "Fit Empty" : Form("%.2f +/- %.3f pC", Mean, Err));

        TString pmtSN = (ch == 0) ? sn1 : (ch == 1) ? sn2 : (ch == 2) ? sn3 : Form("Ch%d", ch);
        const char* chDisplayName[] = {"Mon.", "Rot#1", "Rot#2"};
        TString chLabel = (ch >= 0 && ch <= 2) ? chDisplayName[ch] : Form("Ch%d", ch);
        std::cout << Form("\n  [ %s  %s ]", chLabel.Data(), pmtSN.Data()) << std::endl;
        std::cout << "  --------------------------------------------------------------" << std::endl;
        std::cout << Form("    %-22s : (%.1f, %.1f) deg", "Angle (theta, phi)", tiltVal, rotVal) << std::endl;
        std::cout << Form("    %-22s : %s, theta_Ham = %.1f deg", "Scan Axis / PMT Angle", axisLabel, hamVal) << std::endl;
        if (fitEmpty)
            std::cout << Form("    %-22s : \033[1;33mFit Empty\033[0m", "Rough Charge") << std::endl;
        else
            std::cout << Form("    %-22s : %.3f +/- %.3f pC", "Rough Charge", Mean, Err) << std::endl;
        std::cout << Form("    %-22s : %.3f Sample (%.2f ns)", "Monitoring TTS", rms_exG, rms_exG * 2.0) << std::endl;
        std::cout << Form("    %-22s : %lld (rate %.1f Hz)%s", "Dark Count", noiseCountOut[i], noiseRateOut[i],
                darkmode ? "  [self-trigger: N events with this ch's own pulse / LiveTime]" : "") << std::endl;
        if (darkmode)
            std::cout << Form("    %-22s : %lld (rate %.1f Hz)  [measured in sigStart..sigEnd, real laser-off data]",
                    "Dark in Sig. Window", sigWinDarkHits[i], sigWinDarkRateOut[i]) << std::endl;
        std::cout << "  --------------------------------------------------------------" << std::endl;
        std::cout << "  --------------------------------------------------------------" << std::endl;
    }
    std::cout << "==============================================================\n" << std::endl;

    std::string imageDir = ImagePath + "ByProduce/"; gSystem->mkdir(imageDir.c_str(), kTRUE);
    std::string baseNameForImage = gSystem->BaseName(outName.c_str());
    baseNameForImage.replace(baseNameForImage.find(".root"), 5, "");
    c1->SaveAs((imageDir + baseNameForImage + ".png").c_str());
    cPedStab->SaveAs((imageDir + baseNameForImage + "_PedStability.png").c_str());

    fOut->cd();
    (new TParameter<int>("Config_RecordLength", nSamples))->Write();
    for (int i = 0; i < channelCount; ++i) {
        int ch = activeChannels[i]; if (ch == TriggerCh) continue;
        (new TParameter<Long64_t>(Form("NoiseCount_ch%d", ch), noiseCountOut[i]))->Write();
        (new TParameter<double>(Form("NoiseCountRate_ch%d", ch), noiseRateOut[i]))->Write();
        if (darkmode) {
            (new TParameter<Long64_t>(Form("SigWinDarkCount_ch%d", ch), sigWinDarkHits[i]))->Write();
            (new TParameter<double>(Form("SigWinDarkCountRate_ch%d", ch), sigWinDarkRateOut[i]))->Write();
        }
        (new TParameter<double>(Form("Config_Threshold_mV_ch%d", ch), Analysis_Threshold_mV))->Write();
        (new TParameter<double>(Form("Config_SigWindow_ns_ch%d", ch), Config_SigWindow_ns))->Write();
        (new TParameter<double>(Form("Config_DarkWindow_ns_ch%d", ch), Config_DarkWindow_ns))->Write();
        (new TParameter<double>(Form("TotalNoiseLiveTime_ch%d", ch), noiseLiveTimeOut[i]))->Write();
        (new TParameter<bool>(Form("IsSelfTrigRate_ch%d", ch), darkmode))->Write();

        (new TParameter<double>(Form("PedMean_mV_ch%d", ch), pedMeanOut[i]))->Write();
        (new TParameter<double>(Form("PedSigma_mV_ch%d", ch), pedSigmaOut[i]))->Write();

        double tv = (ch == 1) ? (double)tilt2 : ((ch == 2) ? (double)tilt3 : 0.0);
        double rv = (ch == 1) ? (double)rot2  : ((ch == 2) ? (double)rot3  : 0.0);
        const char* axLbl = "?-axis";
        double hvAng = GetHamamatsuAngle(dirOf(ch), tv, rv, axLbl);
        (new TParameter<double>(Form("HamamatsuAngle_ch%d", ch), hvAng))->Write();
        TNamed(Form("ScanAxis_ch%d", ch), axLbl).Write();
    }

    for (auto *tree : outTrees) if(tree) tree->Write();
    for (auto *hist : histPed) if(hist) hist->Write();
    for (auto *hist : histMax) if(hist) hist->Write();
    for (auto *hist : histFall) if(hist) hist->Write();
    for (auto *hist : histDiff) if(hist) hist->Write();
    for (auto *hist : histPico) if(hist) hist->Write();
    // for (auto *hist : histQT) if(hist) hist->Write();
    for (auto *hist : histNoiseCharge) if(hist) hist->Write();
    for (auto *prof : profPedStab) if(prof) prof->Write();
    for (auto *hist : histNoiseTime) if(hist) hist->Write();
    for (auto *hist : histNoiseMT) if(hist) hist->Write();

    fOut->Close(); file->Close(); delete[] ADC_buffer;
}
