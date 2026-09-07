#include "/home/precalkor/ADC/ADC_test/Base/analysisCode/DrawFormat.cpp"
#include "/home/precalkor/ADC/ADC_test/Base/analysisCode/Analysis.cpp"
#include "path_builder.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TTree.h"
#include "TMath.h"
#include "angle_convert.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cctype>
#include <stdint.h>

void Draw_Contour_v3(int run, const char* specific_file_path = "",
                     double y_lo =10.0, double y_hi = 3.0,
                     int x_start = -1, int x_end = -1) {
    gStyle->SetFrameLineWidth(4);
    gStyle->SetPalette(53);
    gStyle->SetNumberContours(100);
    gStyle->SetTitleFont(22,"");

    TString openPath;
    if (strlen(specific_file_path) > 0) openPath = specific_file_path;
    else openPath = find_filepath(run, RAW_DATA);

    if (gSystem->AccessPathName(openPath)) {
        std::cerr << "[ERROR] Could not find data file: " << openPath << std::endl;
        return;
    }
    std::cout << "[INFO] Opening file: " << openPath << std::endl;

    TChain *chIn = new TChain("T");
    chIn->Add(openPath);
    chIn->LoadTree(0);

    long long NEntry = chIn->GetEntries();
    if (NEntry == 0) {
        std::cerr << "[ERROR] No entries found." << std::endl;
        return;
    }
    std::cout << "[INFO] Total Entries: " << NEntry << std::endl;

    // RecordLength and the OffsetValue<N> active-channel markers moved from
    // "T" to "RunInfo" in the ADC_test8 schema (2026-08-08) -- see the same
    // fallback in prod_ntp_v7.C. Old raw files still have them on T; new ones
    // only have them on RunInfo (as ChannelMask, or as their own
    // OffsetValue<N> branches). Without this, RecordLength silently reads as
    // an uninitialized value and channelCount comes out 0 on any new-format
    // file (confirmed: "Error in <TChain::SetBranchAddress>: unknown branch
    // -> RecordLength" followed by "No active ADC channels found").
    TTree* runInfo = (TTree*)chIn->GetFile()->Get("RunInfo");

    unsigned int RecordLength = 1024;
    bool recordLengthOnT = (chIn->GetBranch("RecordLength") != nullptr);
    if (recordLengthOnT) {
        chIn->SetBranchAddress("RecordLength", &RecordLength);
        chIn->GetEntry(0);
    } else if (runInfo && runInfo->GetBranch("RecordLength")) {
        runInfo->SetBranchAddress("RecordLength", &RecordLength);
        runInfo->GetEntry(0);
        runInfo->ResetBranchAddresses();
    }

    const int nChannels = 8;
    const int nSamples = RecordLength;

    const int trigger_channel = 3;


    const double adc_to_mv_factor = 2000.0 / 16384.0;

    // Read 16-bit continuous buffer
    UShort_t* ADC_buffer = new UShort_t[nChannels * nSamples];
    chIn->SetBranchAddress("ADC", ADC_buffer);

    std::vector<int> channels;
    for (int i = 0; i < nChannels; ++i) {
        if (chIn->GetBranch(Form("OffsetValue%d", i))) channels.push_back(i);
    }
    if (channels.empty() && runInfo) {
        if (runInfo->GetBranch("ChannelMask")) {
            char chMask[16] = {0};
            runInfo->SetBranchAddress("ChannelMask", chMask);
            runInfo->GetEntry(0);
            runInfo->ResetBranchAddresses();
            int len = (int)strlen(chMask);
            for (int i = 0; i < len && i < 8; ++i)
                if (chMask[len - 1 - i] == '1') channels.push_back(i);
        }
        if (channels.empty()) {
            for (int i = 0; i < 8; ++i)
                if (runInfo->GetBranch(Form("OffsetValue%d", i))) channels.push_back(i);
        }
    }

    int channelCount = channels.size();
    if (channelCount == 0) {
        std::cout << "[WARNING] No active ADC channels found." << std::endl;
        delete[] ADC_buffer;
        return;
    }

    // Read stage angles + cable direction from RunInfo so titles can show the
    // converted PMT angle. Direction falls back to DirForCh for older files.
    int rt2 = 0, rr2 = 0, rt3 = 0, rr3 = 0;
    char cdir1[8] = "", cdir2[8] = "", cdir3[8] = "";
    if (TTree* runInfo = (TTree*)chIn->GetFile()->Get("RunInfo")) {
        if (runInfo->GetBranch("RawTiltAngle2"))   runInfo->SetBranchAddress("RawTiltAngle2", &rt2);
        if (runInfo->GetBranch("RawRotateAngle2")) runInfo->SetBranchAddress("RawRotateAngle2", &rr2);
        if (runInfo->GetBranch("RawTiltAngle3"))   runInfo->SetBranchAddress("RawTiltAngle3", &rt3);
        if (runInfo->GetBranch("RawRotateAngle3")) runInfo->SetBranchAddress("RawRotateAngle3", &rr3);
        if (runInfo->GetBranch("Direction1")) runInfo->SetBranchAddress("Direction1", &cdir1);
        if (runInfo->GetBranch("Direction2")) runInfo->SetBranchAddress("Direction2", &cdir2);
        if (runInfo->GetBranch("Direction3")) runInfo->SetBranchAddress("Direction3", &cdir3);
        runInfo->GetEntry(0);
    }
    auto dirOf = [&](int ch) -> char {
        const char* d = (ch == 0) ? cdir1 : (ch == 1) ? cdir2 : cdir3;
        return (d && d[0] != '\0') ? d[0] : DirForCh(ch);
    };

    std::vector<TH2D*> h_contours;
    for (int ch : channels) {
        double tilt_val = (ch == 1) ? (double)rt2 : ((ch == 2) ? (double)rt3 : 0.0);
        double rot_val  = (ch == 1) ? (double)rr2 : ((ch == 2) ? (double)rr3 : 0.0);
        const char* axisLabel = "?-axis";
        double ham_val = GetHamamatsuAngle(dirOf(ch), tilt_val, rot_val, axisLabel);

        // Y축을 아예 0 ~ 2000 mV 범위로 고정하여 선언합니다.
        // Bin의 개수는 해상도 확보를 위해 2000개(1mV/bin)로 설정했습니다.
        TString title = (ch == trigger_channel)
            ? Form("   Channel %d contour;Sample [2ns / Sample] @500 MHz;Voltage [mV]", ch)
            : Form("   Channel %d contour (#theta: %.1f#circ, #phi: %.1f#circ | %s, #theta_{Ham}: %.1f#circ);Sample [2ns / Sample] @500 MHz;Voltage [mV]",
                   ch, tilt_val, rot_val, axisLabel, ham_val);
        h_contours.push_back(new TH2D(
                    Form("h_contour_ch%d", ch),
                    title,
                    nSamples, 0, nSamples,
                    16384, 0, 2000));
    }

    // 데이터 채우기 루프
    for (long long entry = 0; entry < NEntry; ++entry) {
        chIn->GetEntry(entry);
        if ((entry % 1000) == 0) std::cout << "Processing... " << (entry * 100.0) / NEntry << "%\r" << std::flush;

        for (size_t i = 0; i < channels.size(); ++i) {
            int ch = channels[i];
            for (unsigned int sample = 0; sample < nSamples; ++sample) {
                double voltage_mV = ADC_buffer[ch * nSamples + sample] * adc_to_mv_factor;
                h_contours[i]->Fill(sample, voltage_mV);
            }
        }
    }
    std::cout << "\n[INFO] Processing Done." << std::endl;

    TCanvas *c1 = new TCanvas("c1", "ADC Contour Plots", 3600, 900);
    c1->Divide(channelCount, 1);

    for (size_t i = 0; i < h_contours.size(); ++i) {
        int actual_ch = channels[i];
        c1->cd(i + 1);
        gPad->SetGrid(); 
        gPad->SetLeftMargin(0.2); 
        gPad->SetRightMargin(0.15); 
        gPad->SetBottomMargin(0.15); 
        gPad->SetLogz();
        SetHist2DStyle(h_contours[i]);


        if (actual_ch == trigger_channel) {
            h_contours[i]->GetYaxis()->SetRangeUser(0, 2000);
            //    if(i==0)h_contours[i]->GetYaxis()->SetRangeUser(1830, 1840);
            //  else if(i==1)h_contours[i]->GetYaxis()->SetRangeUser(1835, 1840);
            //   else if(i==2)h_contours[i]->GetYaxis()->SetRangeUser(1848, 1855);
        } else {
            TH1D* projY = h_contours[i]->ProjectionY("_py", 1, 100); 
            int maxBin = projY->GetMaximumBin();
            double dynamic_pedestal_mV = projY->GetXaxis()->GetBinCenter(maxBin);
            delete projY;

            h_contours[i]->GetYaxis()->SetRangeUser(dynamic_pedestal_mV - y_lo, dynamic_pedestal_mV + y_hi);
        }
        if (x_start >= 0 && x_end > x_start)
            h_contours[i]->GetXaxis()->SetRangeUser(x_start, x_end);
        h_contours[i]->GetXaxis()->SetNdivisions(505);
        h_contours[i]->Draw("COLZ");
    }
    c1->Update();

    std::string imageDir = "./Data/image/Contour/";
    gSystem->mkdir(imageDir.c_str(), kTRUE);
    std::string baseName = gSystem->BaseName(openPath.Data());
    size_t pos = baseName.find(".root");
    if (pos != std::string::npos) baseName.replace(pos, 5, "");

    std::string finalImagePath = imageDir + baseName + "_Contour.png";
    c1->SaveAs(finalImagePath.c_str());
    std::cout << "[INFO] Contour image saved to: " << finalImagePath << std::endl;

    delete[] ADC_buffer;
    for(auto h : h_contours) delete h;
    delete c1; delete chIn;
}
