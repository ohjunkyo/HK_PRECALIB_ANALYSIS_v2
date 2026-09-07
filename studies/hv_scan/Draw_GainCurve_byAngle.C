// Usage:  root -l -b -q Draw_GainCurve_byAngle.C
//


#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TColor.h"
#include "angle_convert.h"
#include <vector>
#include <map>

const double PC_TO_GAIN_1E7 = 1e-12 / 1.602176634e-19 / 1e7;

struct HVBlock { const char* date; int lo, hi; };
static const std::vector<HVBlock> kBlocks = {
    {"20260903", 800, 808},
/*    {"20260829", 100, 145}, {"20260829", 200, 245}, {"20260829", 0, 45},
    {"20260827", 0, 45},    {"20260815", 0, 45},    {"20260827", 100, 145},
    {"20260830", 0, 45},    {"20260828", 0, 45},    {"20260828", 100, 145},*/
};

static const int kAngleColors[5] = {kRed, kOrange+1, kGreen+2, kAzure+2, kViolet+1};

void DrawByAngle(TVirtualPad* pad, int ch, const char* chLabel, const char* wantAxis) {
    pad->cd();
    pad->SetGridx(); pad->SetGridy();
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(1);
    pad->SetLeftMargin(0.14);

    std::map<int, std::vector<std::pair<double,double>>> byAngle;   // tilt -> (HV, gain)
    std::map<int, std::vector<double>> errByAngle;                  // tilt -> gain err
    std::map<int, double> hamOfTilt;

    for (const auto& b : kBlocks) {
        for (int run = b.lo; run <= b.hi; ++run) {
            TString path = Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", b.date, run);
            TFile* f = TFile::Open(path, "READ");
            if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
            TTree* info = (TTree*)f->Get("RunInfo");
            TTree* tr = (TTree*)f->Get(Form("tree_ch%d", ch));
            int tilt2 = 0, tilt3 = 0, rot2 = 0, rot3 = 0, hv2 = 0, hv3 = 0;
            double spe = 0, speErr = 0;
            char dir2[8] = "", dir3[8] = "";
            if (info) {
                info->SetBranchAddress("RawTiltAngle2", &tilt2); info->SetBranchAddress("RawTiltAngle3", &tilt3);
                info->SetBranchAddress("RawRotateAngle2", &rot2); info->SetBranchAddress("RawRotateAngle3", &rot3);
                info->SetBranchAddress("HV2", &hv2); info->SetBranchAddress("HV3", &hv3);
                if (info->GetBranch("Direction2")) info->SetBranchAddress("Direction2", &dir2);
                if (info->GetBranch("Direction3")) info->SetBranchAddress("Direction3", &dir3);
                info->GetEntry(0);
            }
            if (tr) { tr->SetBranchAddress("spe_mean", &spe); tr->SetBranchAddress("spe_mean_error", &speErr); tr->GetEntry(0); }
            f->Close();
            int tilt = (ch == 1) ? tilt2 : tilt3;
            int rot  = (ch == 1) ? rot2  : rot3;
            int hv   = (ch == 1) ? hv2   : hv3;
            char dir = (ch == 1) ? (dir2[0] ? dir2[0] : DirForCh(1)) : (dir3[0] ? dir3[0] : DirForCh(2));
            const char* axisLabel = "";
            double ham = GetHamamatsuAngle(dir, (double)tilt, (double)rot, axisLabel);
            // At tilt==0 the X-axis and Y-axis scans are the SAME point
            // (zero incidence either way), so it must show up on both pages
            // even when the rotate angle doesn't happen to match this cable's
            // computed x_rot/y_rot ("?-axis") -- 2026-09-03, user: "Draw
            // GainCurve에서 Rot2만 그려지고 Rot1그림은 안그려지네" while only
            // taking center-only test runs (Direction2='H' cable's x_rot/
            // y_rot are 135/45, but the center runs sit at rot2=0, neither of
            // which is a real scan-axis mismatch -- it's just the center).
            bool wanted = (tilt == 0) || (TString(axisLabel) == TString(wantAxis));
            if (spe > 0 && wanted) {
                byAngle[tilt].push_back({(double)hv, spe * PC_TO_GAIN_1E7});
                errByAngle[tilt].push_back(speErr * PC_TO_GAIN_1E7);
                hamOfTilt[tilt] = ham;
            }
        }
    }

    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.15, 0.7, 0.5, 0.92);
    leg->SetTextFont(132); leg->SetTextSize(0.03); leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetFillColor(kWhite);
    leg->SetNColumns(2);

    static const std::vector<int> kShowAngles = {-55, -45, 0, 45, 55};

    int i = 0;
    for (int tilt : kShowAngles) {
        if (!byAngle.count(tilt)) continue;
        auto& pts = byAngle[tilt];
        auto& errs = errByAngle[tilt];
        std::vector<double> hv, g, hvErr;
        for (size_t k = 0; k < pts.size(); ++k) { hv.push_back(pts[k].first); g.push_back(pts[k].second); hvErr.push_back(0.0); }
        TGraphErrors* gr = new TGraphErrors(hv.size(), hv.data(), g.data(), hvErr.data(), errs.data());
        int color = kAngleColors[i % 5];
        gr->SetLineColor(color); gr->SetMarkerColor(color);
        gr->SetMarkerStyle(20); gr->SetMarkerSize(0.5); gr->SetLineWidth(1);
        mg->Add(gr, "LP");
        double ham = hamOfTilt.count(tilt) ? hamOfTilt[tilt] : (double)tilt;
        leg->AddEntry(gr, Form("#theta_{Ham}: %.1f#circ", ham), "p");
        ++i;
    }

    mg->SetTitle(Form("%s : Gain vs HV, per angle (#theta_{Ham}, %s only);Voltage [V];Gain [#times10^{7}]", chLabel, wantAxis));
    mg->Draw("AP");
    leg->Draw();
}

void Draw_GainCurve_byAngle() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(132, "");

    TCanvas* c = new TCanvas("cGainByAngle", "Gain vs HV per angle", 1600, 750);
    TString outDir = "./Data/image/ScanReport/";
    gSystem->mkdir(outDir, kTRUE);
    TString outPath = outDir + "GainCurve_byAngle.pdf";

    c->Clear(); c->Divide(2, 1);
    DrawByAngle(c->cd(1), 1, "Rot1 (EM6400)", "X-axis");
    DrawByAngle(c->cd(2), 2, "Rot2 (EL5150)", "X-axis");
    c->Print(outPath + "(");

    c->Clear(); c->Divide(2, 1);
    DrawByAngle(c->cd(1), 1, "Rot1 (EM6400)", "Y-axis");
    DrawByAngle(c->cd(2), 2, "Rot2 (EL5150)", "Y-axis");
    c->Print(outPath + ")");

    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
