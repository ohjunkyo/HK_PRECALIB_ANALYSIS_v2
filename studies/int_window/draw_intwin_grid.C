// Usage: root -l -b -q 'draw_intwin_grid.C("Rep1",20260831,11,11,1)'   // (repTag, dateTag, startRun, endRun, channel)
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <iostream>

void draw_intwin_grid(const char* repTag, int dateTag, int startRun, int endRun, int ch) {
    TString inDir  = "/home/precalkor/ADC/ADC_test/Data/production/IntWindowStudy";
    TString outDir = "/home/precalkor/ADC/ADC_test/Data/FinalResult/IntWindowStudy";
    gSystem->mkdir(outDir.Data(), kTRUE);

    std::vector<int> widths = {100, 200, 300};
    struct MetricDef { const char* histPrefix; const char* label; };
    std::vector<MetricDef> metrics = {
        {"Pico_ch", "Charge"},
        {"Diff_ch", "Timing"},
        {"Max_ch",  "Pulse Height"},
    };

    // must match the output naming of prod_ntp_v7_intwin.C
    auto prdPath = [&](int run, int w) {
        return TString(Form("%s/precal_prd_kor_run_%d_%03d_int%dns.root",
                            inDir.Data(), dateTag, run, w));
    };

    std::vector<int> runs;
    for (int r = startRun; r <= endRun; ++r)
        if (!gSystem->AccessPathName(prdPath(r, widths[0]))) runs.push_back(r);

    if (runs.empty()) {
        std::cerr << "[ERROR] no files found for date " << dateTag
                  << ", runs " << startRun << "-" << endRun << std::endl;
        std::cerr << "[ERROR] expected e.g. " << prdPath(startRun, widths[0]) << std::endl;
        std::cerr << "[ERROR] run prod_ntp_v7_intwin.C for each width (100/200/300) first." << std::endl;
        return;
    }

    TString pdfPath = Form("%s/Grid_ch%d_%s.pdf", outDir.Data(), ch, repTag);
    TCanvas c("c", "c", 1500, 1400);
    bool first = true;

    for (int run : runs) {
        c.Clear();
        c.Divide(3, 3, 0.002, 0.002);

        for (size_t row = 0; row < metrics.size(); ++row) {
            for (size_t col = 0; col < widths.size(); ++col) {
                c.cd((int)(row * 3 + col + 1));
                gPad->SetGrid();
                if (TString(metrics[row].label) != "Timing") gPad->SetLogy();

                TFile* file = TFile::Open(prdPath(run, widths[col]), "READ");
                if (!file || file->IsZombie()) {
                    std::cerr << "[WARN] missing " << prdPath(run, widths[col]) << std::endl;
                    if (file) delete file;
                    continue;
                }
                TH1D* h = (TH1D*)file->Get(Form("%s%d", metrics[row].histPrefix, ch));
                if (h) {
                    h = (TH1D*)h->Clone(Form("clone_%zu_%zu", row, col));
                    h->SetDirectory(0);   // detach, else file->Close() deletes it
                    h->SetLineColor(kAzure + 2); h->SetLineWidth(2);
                    h->SetTitle(Form("%s @ %d ns;%s", metrics[row].label, widths[col],
                                      h->GetXaxis()->GetTitle()));
                    h->Draw("HIST");
                }
                file->Close();
                delete file;
            }
        }

        // the histogram title already carries PMT SN + angle, so reuse it as the page header
        TString headerTitle = Form("run %d", run);
        TFile* f0 = TFile::Open(prdPath(run, widths[0]), "READ");
        if (f0 && !f0->IsZombie()) {
            TH1D* h0 = (TH1D*)f0->Get(Form("Pico_ch%d", ch));
            if (h0) headerTitle = TString(h0->GetTitle());
            f0->Close();
        }
        if (f0) delete f0;

        c.cd(0);
        TLatex* pageTitle = new TLatex(0.5, 0.985, Form("%s   [%s]", headerTitle.Data(), repTag));
        pageTitle->SetNDC(); pageTitle->SetTextAlign(22); pageTitle->SetTextFont(22);
        pageTitle->SetTextSize(0.018); pageTitle->Draw();

        c.Print(pdfPath + (first ? "(" : ""), "pdf");
        first = false;
    }

    c.Print(pdfPath + ")", "pdf");
    std::cout << "[INFO] Saved: " << pdfPath << std::endl;
}
