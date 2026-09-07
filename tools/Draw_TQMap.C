// TQ map: 2D histogram of Timing (diff = falltime - triggerTime) vs Charge
// (pico, fixed sigStart..sigEnd window), read directly from the already-
// reprocessed production/ ntuple -- no RAW re-read needed.
// Usage: root -l -b -q 'Draw_TQMap.C(20260815,11)'   // (dateTag, run)
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>

void Draw_TQMap(int dateTag, int run) {
    TString path = Form("/home/precalkor/ADC/ADC_test/Data/production/precal_prd_kor_run_%d_%03d.root", dateTag, run);
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { std::cerr << "[ERROR] cannot open " << path << std::endl; return; }

    gStyle->SetOptStat(0);
    TString outDir = "/home/precalkor/ADC/ADC_test/Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString pdfPath = Form("%s/TQMap_%d_%03d.pdf", outDir.Data(), dateTag, run);

    TCanvas c("c", "c", 1500, 500);
    c.Divide(3, 1);
    bool any = false;

    for (int ch = 0; ch < 8; ++ch) {
        TTree* t = (TTree*)f->Get(Form("tree_ch%d", ch));
        if (!t || !t->GetBranch("diff") || !t->GetBranch("pico")) continue;
        any = true;

        c.cd(ch % 3 + 1);
        gPad->SetLogz();
        TString hname = Form("hTQ_ch%d", ch);
        TString sel = Form("pico:diff>>%s(200,0,%d,300,-5,10)", hname.Data(), 1024);
        t->Draw(sel, "diff>-999", "COLZ");
        TH2D* h = (TH2D*)gDirectory->Get(hname);
        if (h) { h->SetTitle(Form("ch%d Q-T Map;Timing [samples];Charge [pC]", ch)); h->Draw("COLZ"); }
    }

    if (!any) { std::cerr << "[ERROR] no diff/pico branches found in " << path << std::endl; return; }
    c.Print(pdfPath);
    std::cout << "[INFO] Saved: " << pdfPath << std::endl;
    f->Close();
}
