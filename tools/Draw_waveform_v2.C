#include "/home/precalkor/ADC/ADC_test/Base/analysisCode/DrawFormat.cpp"
#include "/home/precalkor/ADC/ADC_test/Base/analysisCode/Analysis.cpp"
#include "path_builder2.h"
#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"
#include <stdint.h>

#include <vector> // 이거 없으면 에러 날 수 있으니 꼭 확인해주세요.

int GetCountBelowThreshold(std::vector<unsigned int>& data, float threshold, int n, int start, int end) {
    int count = 0;
    for (int i = start; i < end && i < n; ++i) {
        if (data[i] < threshold) count++;
    }
    return count;
}

void textsize(TH1D* hist){
    gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15);
    hist->GetXaxis()->SetLabelSize(0.05); hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05); hist->GetYaxis()->SetTitleSize(0.05);
}

void Draw_waveform_v2(int run, const char* filePath = "") {    
    gStyle->SetTitleFont(22,""); gStyle->SetTitleSize(0.09); gStyle->SetFrameLineWidth(2);

    TString openPath;
    if (std::string(filePath) != "" && std::string(filePath) != "\"\"") openPath = filePath;
    else openPath = find_filepath(run, RAW_DATA).c_str();

    if (gSystem->AccessPathName(openPath)) { 
        std::cerr << "[ERROR] Could not find raw data file: " << openPath << std::endl; 
        return;
    }
    std::cout << "[INFO] Opening file: " << openPath << std::endl;

    TChain *chIn = new TChain("T");
    chIn->Add(openPath);
    chIn->LoadTree(0);
    int NEntry = chIn->GetEntries();
    std::cout << "[INFO] Total Entries: " << NEntry << std::endl;

    unsigned int RecordLength;
    chIn->SetBranchAddress("RecordLength", &RecordLength);
    chIn->GetEntry(0); 

    const int nChannels = 8;
    const int nSamples = RecordLength;

    // Read 16-bit continuous buffer
    UShort_t* ADC_buffer = new UShort_t[nChannels * nSamples];
    chIn->SetBranchAddress("ADC", ADC_buffer);

    std::vector<int> channels;
    for (int i = 0; i < nChannels; ++i) {
        if (chIn->GetBranch(Form("OffsetValue%d", i))) channels.push_back(i);
    }
    int channelCount = channels.size();
    std::cout << "[INFO] Active Channels: " << channelCount << std::endl;

    unsigned int triggerValues[nChannels] = {0}; 
    for (int ch : channels) {
        if (chIn->GetBranch(Form("TriggerValue%d", ch))) {
            chIn->SetBranchAddress(Form("TriggerValue%d", ch), &triggerValues[ch]);
        } else triggerValues[ch] = 15000; 
    }

    TCanvas *c1 = new TCanvas("c1", "Event Viewer", 1600, 900);
    c1->Divide(2, (channelCount + 1) / 2);

    string inputQuit;
    std::vector<TH1D*> histograms(channelCount);
    long long entry = 0;

    while (entry < NEntry) {
        chIn->GetEntry(entry);
        for (size_t i = 0; i < channels.size(); ++i) {
            int ch = channels[i];
            if (!histograms[i]) {
                histograms[i] = new TH1D(Form("Channel_%d", ch),
                        Form("Channel %d; ADC Sample [2ns / Sample] @500MHz; Voltage [mV] @Max 2V", ch), nSamples, 0, nSamples);
            } else histograms[i]->Reset();

            for (unsigned int sample = 0; sample < nSamples; ++sample) {
                histograms[i]->SetBinContent(sample + 1, ADC_buffer[ch * nSamples + sample] * Config::ADC_to_mV);
            }

            c1->cd(i + 1);
            gPad->SetLeftMargin(0.15); gPad->SetGrid();
            SetHistStyle(histograms[i]);
            histograms[i]->SetTitle(Form("CH %d #%lld", ch, entry));
            textsize(histograms[i]);
            histograms[i]->Draw("HIST");

            double ThrValue = triggerValues[ch] * Config::ADC_to_mV;
            TLine *thrLine = new TLine(0, ThrValue, nSamples, ThrValue);
            thrLine->SetLineColor(kRed); thrLine->SetLineWidth(2); thrLine->Draw("same");
        }
        c1->Update();

        std::cout << "[INFO] Entry #" << entry << " / " << NEntry - 1 
            << ". [n: next, s: search hit, q: quit, or #]: ";
        getline(cin, inputQuit);

        if (inputQuit == "q" || inputQuit == "Q") break; 
        else if (inputQuit == "s" || inputQuit == "S") {
            bool found = false;
            long long searchEntry = entry + 1; 
            std::cout << "[SEARCHING] Looking for next threshold crossing..." << std::endl;

            while (searchEntry < NEntry) {
                chIn->GetEntry(searchEntry);
                for (int ch : channels) {
                    if (ch == TriggerCh) continue;
                    std::vector<unsigned int> temp_adc(nSamples);
                    for(unsigned int s=0; s<nSamples; ++s) temp_adc[s] = ADC_buffer[ch * nSamples + s];
                    double pedestal = GetPedestal(temp_adc.data(), 0, 100, nSamples);
                    //double pedestal = GetPedestal(temp_adc, 0, 100, nSamples);
                    if (GetCountBelowThreshold(temp_adc, pedestal - 30, nSamples, 0, 500) > 0) {
                        found = true; break;
                    }
                }
                if (found) {
                    entry = searchEntry;
                    std::cout << "[FOUND] Hit detected at Entry #" << entry << std::endl;
                    break;
                }
                searchEntry++;
                if (searchEntry % 5000 == 0) std::cout << "  Searching... " << searchEntry << "\r" << std::flush;
            }
            if (!found) {
                std::cout << "[INFO] No more hit events found." << std::endl;
                entry = NEntry; 
            }
        }
        else if (inputQuit.empty() || inputQuit == "n" || inputQuit == "N") entry++;
        else {
            try {
                long long targetEntry = std::stoll(inputQuit);
                if (targetEntry >= 0 && targetEntry < NEntry) entry = targetEntry; 
            } catch (...) {
                std::cout << "[WARNING] Invalid input." << std::endl;
            }
        }
    }
    delete[] ADC_buffer;
}
