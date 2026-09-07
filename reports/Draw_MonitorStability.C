// Monitor PMT (CH0) stability analysis: separates "the laser output drifted"
// from "the Monitor PMT itself drifted", and quantifies how much of the
// run-to-run scatter is plain counting statistics versus a real drift.
//
// Why this works
// --------------
// The Monitor PMT never rotates, so its geometry (and hence its intrinsic
// photon-detection probability) is fixed across every run of a scan. Two
// independent observables let us split the two possible drift sources:
//
//   spe_mean  (SPE charge, i.e. dynode gain)
//       Independent of how many photons arrive -- it is the charge produced by
//       ONE photoelectron. If this drifts, the PMT/HV is the culprit; the laser
//       cannot move it.
//
//   relativeQE_raw  (PHC counting QE = N_phc / N_total)
//       Directly proportional to how many photons actually landed in that run.
//       If the gain above is stable but this still scatters, the light level
//       itself moved -> the laser.
//
// Expected vs observed scatter
// ----------------------------
// Each run's QE already carries its own Poisson counting error
// (relativeQE_raw_err = sqrt(N_phc)/N_total x100%). So for M runs:
//
//   sigma_expected^2 = (1/M) * SUM_i (err_i)^2         <- pure counting stats
//   sigma_observed^2 = (1/(M-1)) * SUM_i (QE_i - mean)^2
//   sigma_excess     = sqrt(max(0, sigma_obs^2 - sigma_exp^2))
//
// The ratio sigma_obs/sigma_exp is the square root of the reduced chi2 against
// a constant model: ~1 means counting statistics explain everything (no
// evidence of drift), >>1 means there is genuine excess variance.
//
// NOTE on the laser telemetry: LOG/LASER/*.csv records the driver SETPOINT
// (pulse_ma / bias_ma), not a measured optical power, so it cannot reveal
// shot-to-shot or run-to-run light-level wobble. The Monitor PMT counting rate
// IS our only real power meter -- which is exactly what this macro exploits.

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

struct RunPoint {
    int run;
    double qe, qeErr;        // relativeQE_raw [%] (PHC counting, no dark subtraction)
    double gain, gainErr;    // spe_mean [pC]
};

// Mean / observed-sigma / expected-sigma / excess for one series.
struct ScatterStats {
    double mean = 0, obsSigma = 0, expSigma = 0, excess = 0, ratio = 0;
    int n = 0;
};

static ScatterStats ComputeScatter(const std::vector<double>& vals, const std::vector<double>& errs) {
    ScatterStats s;
    s.n = (int)vals.size();
    if (s.n < 2) return s;

    for (double v : vals) s.mean += v;
    s.mean /= s.n;

    double sumSq = 0;
    for (double v : vals) sumSq += (v - s.mean) * (v - s.mean);
    s.obsSigma = std::sqrt(sumSq / (s.n - 1));   // sample std dev

    double sumErrSq = 0;
    for (double e : errs) sumErrSq += e * e;
    s.expSigma = std::sqrt(sumErrSq / s.n);      // RMS of the per-run statistical errors

    double excessVar = s.obsSigma * s.obsSigma - s.expSigma * s.expSigma;
    s.excess = (excessVar > 0) ? std::sqrt(excessVar) : 0.0;
    s.ratio  = (s.expSigma > 0) ? s.obsSigma / s.expSigma : 0.0;
    return s;
}

void Draw_MonitorStability(TString tag = "20260803", int run_start = 200, int run_end = 245) {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(22, "");
    gStyle->SetTitleSize(0.05, "");
    gStyle->SetTitleFont(132, "XY"); gStyle->SetLabelFont(132, "XY");

    TString dir = "./Data/FinalResult/";
    std::vector<RunPoint> pts;

    TSystemDirectory sysdir(dir, dir);
    TList* files = sysdir.GetListOfFiles();
    if (!files) { std::cout << "[ERROR] cannot list " << dir << std::endl; return; }

    for (int i = 0; i < files->GetSize(); ++i) {
        TSystemFile* file = (TSystemFile*)files->At(i);
        TString fname = file->GetName();
        if (file->IsDirectory() || !fname.EndsWith(".root") || !fname.Contains("result")) continue;
        if (!fname.Contains(tag)) continue;

        TObjArray* tokens = fname.Tokenize("_.");
        int runNum = -1;
        if (tokens->GetEntries() >= 2)
            runNum = ((TObjString*)tokens->At(tokens->GetEntries() - 2))->GetString().Atoi();
        delete tokens;
        if (runNum < run_start || runNum > run_end) continue;

        TFile* f = TFile::Open(dir + fname, "READ");
        if (!f || f->IsZombie()) continue;

        // Skip dedicated dark runs -- no laser light, so they carry no
        // information about laser stability and would bias the QE scatter.
        TTree* info = (TTree*)f->Get("RunInfo");
        if (info && info->GetBranch("RunMode")) {
            char runModeBuf[20] = "";
            info->SetBranchAddress("RunMode", runModeBuf);
            info->GetEntry(0);
            if (TString(runModeBuf) == "Dark") { f->Close(); continue; }
        }

        TTree* tr = (TTree*)f->Get("tree_ch0");   // CH0 = Monitor PMT
        if (!tr) { f->Close(); continue; }

        double qe = 0, qeErr = 0, gain = 0, gainErr = 0;
        if (!tr->GetBranch("relativeQE_raw")) { f->Close(); continue; }
        tr->SetBranchAddress("relativeQE_raw", &qe);
        tr->SetBranchAddress("relativeQE_raw_err", &qeErr);
        tr->SetBranchAddress("spe_mean", &gain);
        tr->SetBranchAddress("spe_mean_error", &gainErr);
        tr->GetEntry(0);
        f->Close();

        // Guard against failed fits / sentinel values.
        if (qe <= 0 || qe > 100 || gain <= 0) continue;
        pts.push_back({runNum, qe, qeErr, gain, gainErr});
    }

    if (pts.size() < 2) { std::cout << "[ERROR] need >=2 usable runs, got " << pts.size() << std::endl; return; }
    std::sort(pts.begin(), pts.end(), [](const RunPoint& a, const RunPoint& b){ return a.run < b.run; });

    std::vector<double> runX, qeV, qeE, gainV, gainE;
    for (auto& p : pts) {
        runX.push_back(p.run);
        qeV.push_back(p.qe);     qeE.push_back(p.qeErr);
        gainV.push_back(p.gain); gainE.push_back(p.gainErr);
    }

    ScatterStats qeS   = ComputeScatter(qeV, qeE);
    ScatterStats gainS = ComputeScatter(gainV, gainE);

    // ---- console report ----
    std::cout << "\n\033[1;33m[ Monitor PMT (CH0) Stability -- " << tag
              << " run " << run_start << "-" << run_end << " ]\033[0m" << std::endl;
    std::cout << "  ------------------------------------------------------------" << std::endl;
    std::cout << Form("  Runs used                 : %d", qeS.n) << std::endl;
    std::cout << "\n  \033[1;36mSPE Gain (spe_mean) -- PMT hardware stability\033[0m" << std::endl;
    std::cout << Form("    mean                    : %.4f pC", gainS.mean) << std::endl;
    std::cout << Form("    observed sigma          : %.4f pC  (%.2f %% of mean)",
                       gainS.obsSigma, gainS.mean > 0 ? gainS.obsSigma / gainS.mean * 100 : 0) << std::endl;
    std::cout << Form("    expected (fit err)      : %.4f pC", gainS.expSigma) << std::endl;
    std::cout << Form("    obs/exp ratio           : %.2f", gainS.ratio) << std::endl;
    std::cout << Form("    excess sigma            : %.4f pC  (%.2f %% of mean)",
                       gainS.excess, gainS.mean > 0 ? gainS.excess / gainS.mean * 100 : 0) << std::endl;

    std::cout << "\n  \033[1;36mCounting QE (relativeQE_raw) -- light level + PMT\033[0m" << std::endl;
    std::cout << Form("    mean                    : %.4f %%", qeS.mean) << std::endl;
    std::cout << Form("    observed sigma          : %.4f %%  (%.2f %% of mean)",
                       qeS.obsSigma, qeS.mean > 0 ? qeS.obsSigma / qeS.mean * 100 : 0) << std::endl;
    std::cout << Form("    expected (Poisson)      : %.4f %%", qeS.expSigma) << std::endl;
    std::cout << Form("    obs/exp ratio           : %.2f", qeS.ratio) << std::endl;
    std::cout << Form("    excess sigma            : %.4f %%  (%.2f %% of mean)",
                       qeS.excess, qeS.mean > 0 ? qeS.excess / qeS.mean * 100 : 0) << std::endl;

    std::cout << "\n  \033[1;33mInterpretation\033[0m" << std::endl;
    if (qeS.ratio < 1.3) {
        std::cout << "    QE scatter is consistent with pure counting statistics" << std::endl;
        std::cout << "    -> no evidence of laser or Monitor-PMT drift at this precision." << std::endl;
    } else if (gainS.ratio > 2.0 && gainS.mean > 0 && gainS.excess / gainS.mean > 0.02) {
        std::cout << "    BOTH gain and QE show excess scatter" << std::endl;
        std::cout << "    -> Monitor PMT / HV instability is at least part of the cause." << std::endl;
    } else {
        std::cout << "    QE shows excess scatter while the gain does not" << std::endl;
        std::cout << "    -> points to the LASER output level moving run-to-run," << std::endl;
        std::cout << "       since PMT gain (photon-count independent) stayed put." << std::endl;
    }
    std::cout << "  ------------------------------------------------------------" << std::endl;

    // ---- plot: 2 stacked panels ----
    TCanvas* c = new TCanvas("c_monstab", "Monitor PMT Stability", 1500, 1000);
    c->Divide(1, 2, 0.001, 0.001);

    auto drawPanel = [&](int pad, const char* title, const char* yTitle,
                          std::vector<double>& v, std::vector<double>& e,
                          const ScatterStats& st, int color, const char* unit) {
        c->cd(pad);
        gPad->SetGrid(); gPad->SetLeftMargin(0.09); gPad->SetRightMargin(0.30);
        gPad->SetTopMargin(0.11); gPad->SetBottomMargin(0.14);

        TGraphErrors* gr = new TGraphErrors(v.size(), &runX[0], &v[0], 0, &e[0]);
        gr->SetTitle(Form("%s;Run number;%s", title, yTitle));
        gr->SetMarkerStyle(20); gr->SetMarkerSize(1.0);
        gr->SetMarkerColor(color); gr->SetLineColor(color); gr->SetLineWidth(1);
        gr->Draw("AP");
        gr->GetXaxis()->SetTitleSize(0.050); gr->GetXaxis()->SetLabelSize(0.043);
        gr->GetYaxis()->SetTitleSize(0.050); gr->GetYaxis()->SetLabelSize(0.043);
        gr->GetYaxis()->SetTitleOffset(0.85);

        double xLo = runX.front() - 1, xHi = runX.back() + 1;
        TLine* lMean = new TLine(xLo, st.mean, xHi, st.mean);
        lMean->SetLineColor(kBlack); lMean->SetLineWidth(3); lMean->Draw();
        // +-1 expected (pure statistics) band, so excess scatter is visible by eye
        TLine* lBand = nullptr;
        for (int sgn : {-1, 1}) {
            TLine* l = new TLine(xLo, st.mean + sgn * st.expSigma, xHi, st.mean + sgn * st.expSigma);
            l->SetLineColor(kGray + 2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
            if (sgn == 1) lBand = l;
        }

        TLegend* lg = new TLegend(0.705, 0.755, 0.995, 0.90);
        lg->SetBorderSize(1); lg->SetFillColor(kWhite); lg->SetTextFont(132); lg->SetTextSize(0.038);
        lg->AddEntry(gr, "run value #pm stat. error", "lp");
        lg->AddEntry(lMean, "mean", "l");
        lg->AddEntry(lBand, "#pm#sigma_{exp} (stat. only)", "l");
        lg->Draw();

        // Results box -- spells out where sigma_exp comes from, since that is the
        // least obvious quantity here: it is just the RMS of the error bars that
        // are already drawn on this very panel, i.e. "how big is a typical bar".
        TPaveText* box = new TPaveText(0.705, 0.16, 0.995, 0.74, "NDC");
        box->SetBorderSize(2); box->SetFillColor(kWhite); box->SetFillStyle(1001);
        box->SetTextFont(132); box->SetTextSize(0.037); box->SetTextAlign(12);
        box->AddText(Form("N runs = %d", st.n));
        box->AddText(Form("mean = %.4f %s", st.mean, unit));
        box->AddText("");
        box->AddText(Form("#sigma_{obs} = %.4f", st.obsSigma));
        box->AddText("   scatter of points");
        box->AddText(Form("#sigma_{exp} = %.4f", st.expSigma));
        box->AddText("   RMS of error bars");
        box->AddText("");
        box->AddText(Form("#sigma_{obs}/#sigma_{exp} = %.2f", st.ratio));
        box->AddText(Form("#sigma_{excess} = %.4f", st.excess));
        if (st.mean > 0) box->AddText(Form("   %.2f %% of mean", st.excess / st.mean * 100));
        box->Draw();
    };

    drawPanel(1, "Monitor PMT SPE Gain", "SPE mean [pC]",
               gainV, gainE, gainS, kBlue + 2, "pC");
    drawPanel(2, "Monitor PMT Counting QE", "Raw QE [%]",
               qeV, qeE, qeS, kRed + 1, "%");

    TString outPath = Form("./Data/image/Uniformity/MonitorStability_%s_%d_%d.pdf",
                            tag.Data(), run_start, run_end);
    c->Print(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
