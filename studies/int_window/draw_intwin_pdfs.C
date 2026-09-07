// Builds one multi-page PDF per (metric, channel) for a charge-integration-
// window study variant -- e.g. PulseHeight_ch1_Rep1_200ns.pdf, one page per
// scan angle in that rep. Reads the histograms straight out of the
// prod_ntp_v7_intwin.C outputs under Data/production/IntWindowStudy/
// (Max_ch<N> = Pulse Height, Diff_ch<N> = Timing, Pico_ch<N> = Charge --
// titles already carry the angle info baked in by prod_ntp_v7_intwin.C).
#include <vector>
#include <string>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TObjString.h>

void draw_intwin_pdfs(const char* repTag, int startRun, int endRun, int intWindowNs) {
    TString inDir = "/home/precalkor/ADC/ADC_test/Data/production/IntWindowStudy";
    TString outDir = "/home/precalkor/ADC/ADC_test/Data/FinalResult/IntWindowStudy";
    gSystem->mkdir(outDir.Data(), kTRUE);

    // Collect matching files, sorted by run number, within [startRun, endRun].
    TSystemDirectory sd("d", inDir.Data());
    TList *files = sd.GetListOfFiles();
    std::vector<std::pair<int, TString>> runFiles;
    if (files) {
        TIter next(files);
        TSystemFile *f;
        TString suffix = Form("_int%dns.root", intWindowNs);
        while ((f = (TSystemFile*)next())) {
            TString name = f->GetName();
            if (f->IsDirectory() || !name.EndsWith(suffix)) continue;
            TObjArray *toks = name.Tokenize("_");
            int runNum = -1;
            for (int i = 0; i < toks->GetEntries(); ++i) {
                TString t = ((TObjString*)toks->At(i))->GetString();
                if (t.IsDigit() && t.Length() <= 3) runNum = t.Atoi();
            }
            delete toks;
            if (runNum >= startRun && runNum <= endRun) {
                runFiles.push_back({runNum, inDir + "/" + name});
            }
        }
    }
    std::sort(runFiles.begin(), runFiles.end());
    if (runFiles.empty()) {
        std::cerr << "[ERROR] No files found for run range " << startRun << "-" << endRun
                   << " with suffix _int" << intWindowNs << "ns.root in " << inDir << std::endl;
        return;
    }
    std::cout << "[INFO] " << repTag << " @ " << intWindowNs << "ns: " << runFiles.size() << " runs found." << std::endl;

    struct MetricDef { const char* histPrefix; const char* label; };
    std::vector<MetricDef> metrics = {
        {"Max_ch",  "PulseHeight"},
        {"Diff_ch", "Timing"},
        {"Pico_ch", "Charge"},
    };

    for (int ch = 0; ch <= 2; ++ch) {
        for (auto &m : metrics) {
            TString pdfPath = Form("%s/%s_ch%d_%s_%dns.pdf", outDir.Data(), m.label, ch, repTag, intWindowNs);
            TCanvas c("c", "c", 1000, 700);
            bool first = true;
            for (auto &rf : runFiles) {
                TFile *file = TFile::Open(rf.second, "READ");
                if (!file || file->IsZombie()) continue;
                TH1D *h = (TH1D*)file->Get(Form("%s%d", m.histPrefix, ch));
                if (h) {
                    h = (TH1D*)h->Clone();
                    h->SetDirectory(0);
                    c.cd(); c.Clear();
                    if (TString(m.label) != "Timing") gPad->SetLogy();
                    gPad->SetGrid();
                    h->SetLineColor(kAzure + 2); h->SetLineWidth(2);
                    h->Draw("HIST");
                    c.Print(pdfPath + (first ? "(" : ""), "pdf");
                    first = false;
                    delete h;
                }
                file->Close();
                delete file;
            }
            if (!first) {
                c.Print(pdfPath + ")", "pdf");   // close multi-page pdf
                std::cout << "  -> " << pdfPath << std::endl;
            }
        }
    }
}
