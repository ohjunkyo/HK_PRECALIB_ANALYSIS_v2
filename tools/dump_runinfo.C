// dump_runinfo.C  —  prints the RunInfo metadata of a .root file as KEY=VALUE lines.
// Used by the GUI (Data Files tab) to show SN / HV / angles / shifter without uproot.
#include <TFile.h>
#include <TTree.h>
#include <cstdio>

void dump_runinfo(const char* path) {
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { printf("ERR=cannot open file\n"); return; }
    TTree* t = (TTree*)f->Get("RunInfo");
    if (!t) { printf("ERR=no RunInfo tree\n"); return; }

    char sn1[128]="", sn2[128]="", sn3[128]="", runmode[64]="", expert[128]="",
         shifter[128]="", note[256]="";
    int hv1=0, hv2=0, hv3=0, r2=0, t2=0, r3=0, t3=0, laser=0, wl=0;

    auto SB = [&](const char* n, void* p){ if (t->GetBranch(n)) t->SetBranchAddress(n, p); };
    SB("SN1", sn1); SB("SN2", sn2); SB("SN3", sn3);
    SB("HV1", &hv1); SB("HV2", &hv2); SB("HV3", &hv3);
    SB("RawRotateAngle2", &r2); SB("RawTiltAngle2", &t2);
    SB("RawRotateAngle3", &r3); SB("RawTiltAngle3", &t3);
    SB("RunMode", runmode); SB("Expert", expert); SB("Shifter", shifter);
    SB("Wavelength", &wl); SB("Laser_mA", &laser); SB("NOTE", note);

    t->GetEntry(0);

    printf("RunMode=%s\n", runmode);
    printf("Shifter=%s\n", shifter);
    printf("Expert=%s\n", expert);
    printf("SN1=%s\n", sn1); printf("SN2=%s\n", sn2); printf("SN3=%s\n", sn3);
    printf("HV1=%d\n", hv1); printf("HV2=%d\n", hv2); printf("HV3=%d\n", hv3);
    printf("Rot2=%d\n", r2); printf("Tilt2=%d\n", t2);
    printf("Rot3=%d\n", r3); printf("Tilt3=%d\n", t3);
    printf("Wavelength=%d\n", wl); printf("Laser_mA=%d\n", laser); printf("NOTE=%s\n", note);
    f->Close();
}
