// Usage: root -l -b -q 'Draw_IntWindow.C(20260831,11,11)'   // (dateTag, startRun, endRun) -- all channels
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <iostream>

void Draw_IntWindow(int dateTag, int startRun, int endRun) {
    TString inDir  = "/home/precalkor/ADC/ADC_test/Data/FinalResult/IntWindowStudy";
    TString outDir = "/home/precalkor/ADC/ADC_test/Data/FinalResult/IntWindowStudy";
    gSystem->mkdir(outDir.Data(), kTRUE);

    std::vector<int> widths = {100, 200, 300};
    struct MetricDef { const char* histPrefix; const char* label; };
    std::vector<MetricDef> metrics = {
        {"Pico_ch", "Charge"},
        {"Diff_ch", "Timing"},
        {"Max_ch",  "Pulse Height"},
    };

    // must match read_ntp_v7_intWindow.C's output naming
    auto resPath = [&](int run, int w) {
        return TString(Form("%s/precal_result_kor_run_%d_%03d_int%dns.root", inDir.Data(), dateTag, run, w));
    };

    std::vector<int> runs;
    for (int r = startRun; r <= endRun; ++r)
        if (!gSystem->AccessPathName(resPath(r, widths[0]))) runs.push_back(r);

    if (runs.empty()) {
        std::cerr << "[ERROR] no files found, expected e.g. " << resPath(startRun, widths[0]) << std::endl;
        std::cerr << "[ERROR] run prod_ntp_v7_intWindow.C + read_ntp_v7_intWindow.C for 100/200/300 first." << std::endl;
        return;
    }

    // channels taken from the file itself, so Monitor + both Rot PMTs all land
    // in the one PDF without having to be listed here
    std::vector<int> chans;
    {
        TFile* f = TFile::Open(resPath(runs[0], widths[0]), "READ");
        if (f && !f->IsZombie())
            for (int c0 = 0; c0 < 8; ++c0)
                if (f->Get(Form("Pico_ch%d", c0))) chans.push_back(c0);
        if (f) { f->Close(); delete f; }
    }
    if (chans.empty()) { std::cerr << "[ERROR] no Pico_ch* histograms in the result files" << std::endl; return; }

    gStyle->SetOptStat(0);
    TString pdfPath = Form("%s/Grid_%d_run%d-%d.pdf", outDir.Data(), dateTag, startRun, endRun);
    TCanvas c("c", "c", 1500, 1400);
    bool first = true;

    for (int run : runs) {
      for (int ch : chans) {
        c.Clear();
        c.Divide(3, 3, 0.002, 0.002);

        for (size_t row = 0; row < metrics.size(); ++row) {
            for (size_t col = 0; col < widths.size(); ++col) {
                c.cd((int)(row * 3 + col + 1));
                gPad->SetGrid();
                if (TString(metrics[row].label) != "Timing") gPad->SetLogy();

                TFile* f = TFile::Open(resPath(run, widths[col]), "READ");
                if (!f || f->IsZombie()) {
                    std::cerr << "[WARN] missing " << resPath(run, widths[col]) << std::endl;
                    if (f) delete f;
                    continue;
                }
                TH1D* h = (TH1D*)f->Get(Form("%s%d", metrics[row].histPrefix, ch));
                if (h) {
                    h = (TH1D*)h->Clone(Form("clone_%zu_%zu", row, col));
                    h->SetDirectory(0);   // detach, else f->Close() deletes it
                    h->SetLineColor(kAzure + 2); h->SetLineWidth(2);
                    h->SetTitle(Form("%s @ %d ns;%s", metrics[row].label, widths[col],
                                      h->GetXaxis()->GetTitle()));
                    h->Draw("HIST");
                }
                f->Close(); delete f;
            }
        }

        // the histogram title already carries PMT SN + angle, so reuse it as the page header
        TString headerTitle = Form("run %d  ch%d", run, ch);
        TFile* f0 = TFile::Open(resPath(run, widths[0]), "READ");
        if (f0 && !f0->IsZombie()) {
            TH1D* h0 = (TH1D*)f0->Get(Form("Pico_ch%d", ch));
            if (h0) headerTitle = TString(h0->GetTitle());
            f0->Close();
        }
        if (f0) delete f0;

        c.cd(0);
        TLatex* pageTitle = new TLatex(0.5, 0.985, headerTitle);
        pageTitle->SetNDC(); pageTitle->SetTextAlign(22); pageTitle->SetTextFont(22);
        pageTitle->SetTextSize(0.018); pageTitle->Draw();

        c.Print(pdfPath + (first ? "(" : ""), "pdf");
        first = false;
      }
    }

    c.Print(pdfPath + ")", "pdf");
    std::cout << "[INFO] Saved: " << pdfPath << std::endl;
}
