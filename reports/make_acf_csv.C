#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>

void make_acf_csv(TString tag, int run_start, int run_end, TString outCsv) {
    static const char* kRawDirs[] = {
        "./Data/RAW/Laser",
        "/home/precalkor/external_HDD_1_4T/Data_Backup/RAW/Laser",
    };
    std::ofstream out(outCsv.Data());
    for (int run = run_start; run <= run_end; ++run) {
        TString path = Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", tag.Data(), run);
        if (gSystem->AccessPathName(path)) continue;
        double mtime = 0; struct stat st;
        for (const char* dir : kRawDirs) {
            TString rp = Form("%s/precal_raw_kor_run_%s_%03d.root", dir, tag.Data(), run);
            if (stat(rp.Data(), &st) == 0) { mtime = (double)st.st_mtime; break; }
        }
        if (mtime == 0) continue;
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        double q[3] = {0,0,0}; bool ok[3] = {false,false,false};
        for (int c = 0; c < 3; ++c) {
            TTree* tr = (TTree*)f->Get(Form("tree_ch%d", c));
            if (!tr) continue;
            double qeRaw = 0;
            if (!tr->GetBranch("relativeQE_raw")) continue;
            tr->SetBranchAddress("relativeQE_raw", &qeRaw);
            tr->GetEntry(0);
            q[c] = qeRaw; ok[c] = true;
        }
        f->Close();
        if (ok[0] && ok[1] && ok[2])
            out << run << "," << (long long)mtime << "," << q[0] << "," << q[1] << "," << q[2] << "\n";
    }
    out.close();
    std::cout << "[INFO] wrote " << outCsv << std::endl;
}
