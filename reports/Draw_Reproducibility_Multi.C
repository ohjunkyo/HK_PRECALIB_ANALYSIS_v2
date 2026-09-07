// Multi-set reproducibility: every scan compared against ONE reference (the
// first tag), instead of the pairwise Draw_Reproducibility.C.
//
// Why a reference rather than all pairs: with N sets there are N(N-1)/2 pairs
// (6 for 4 sets) and no natural ordering to read them in, so a drift over
// time is invisible. Pinning Rep.1 as the reference gives N-1 comparisons
// that share a baseline and can be read as a sequence.
//
// Per (PMT, axis) page:
//   top pad  -- all sets overlaid, one marker shape per set
//   bot pad  -- ratio (Set_n / Ref) per set, x-offset so overlapping error
//               bars stay readable, with the reference band drawn from the
//               reference set's own statistical error
// Then a summary-table page, a combined pull page, and -- when the B-field
// log covers the runs -- a page showing the coil current / |B| each set was
// taken under, since "were these actually taken under the same conditions?"
// is the first question a reproducibility difference raises.
//
// Usage:
//   root -l -b -q 'Draw_Reproducibility_Multi.C({"20260815_0_45","20260816_0_45",
//                   "20260816_100_145","20260817_0_45"},{"EM6400","EL5150"})'
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
#include <TLine.h>
#include <TBox.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TF1.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "repro_common.h"

namespace {

const int kSetColor[6]  = {kBlack, kBlue + 2, kRed + 1, kGreen + 3, kMagenta + 2, kOrange + 8};
// Shape carries the set identity so the plot survives greyscale printing.
const int kSetMarker[6] = {20, 25, 22, 32, 33, 27};

struct SetCurve { std::vector<double> ang, val, err; };

SetCurve LoadCurve(const TString& tag, const TString& serial, const char* axis,
                   const TString& metric) {
    SetCurve c;
    TFile* f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag.Data()));
    if (!f || f->IsZombie()) { if (f) f->Close(); return c; }
    TGraphErrors* g = (TGraphErrors*)f->Get(Form("gr_%s_%s_%s", serial.Data(), axis, metric.Data()));
    if (g) {
        for (int i = 0; i < g->GetN(); ++i) {
            if (g->GetY()[i] <= 0) continue;          // sentinel/invalid guard
            c.ang.push_back(g->GetX()[i]);
            c.val.push_back(g->GetY()[i]);
            c.err.push_back(g->GetEY()[i]);
        }
    }
    f->Close();
    return c;
}

// Match set `b` onto reference `a` by angle, producing the pull/ratio points
// ComputeRepro() expects.
std::vector<MatchedPoint> Match(const SetCurve& a, const SetCurve& b, double tol) {
    std::vector<MatchedPoint> out;
    for (size_t i = 0; i < a.ang.size(); ++i) {
        for (size_t j = 0; j < b.ang.size(); ++j) {
            if (std::abs(a.ang[i] - b.ang[j]) > tol) continue;
            MatchedPoint p;
            p.angle = a.ang[i];
            p.v1 = a.val[i]; p.e1 = a.err[i];
            p.v2 = b.val[j]; p.e2 = b.err[j];
            p.ratio = p.v2 / p.v1;
            p.ratioErr = p.ratio * std::sqrt(std::pow(p.e1 / p.v1, 2) + std::pow(p.e2 / p.v2, 2));
            double den = std::sqrt(p.e1 * p.e1 + p.e2 * p.e2);
            p.pull = den > 0 ? (p.v2 - p.v1) / den : 0;
            out.push_back(p);
            break;
        }
    }
    return out;
}

// ---- B-field conditions per set -------------------------------------------
// Reads the same CSV Draw_Overlay_Uniformity_v7.C uses. Returns false when the
// log doesn't cover the window, so the page can say so instead of drawing
// zeros that look like a real measurement.
struct FieldCond { bool ok = false; double I[4] = {0,0,0,0}; double Isd[4] = {0,0,0,0};
                   double B[3] = {0,0,0}; double Bsd[3] = {0,0,0}; int n = 0; };

double ParseIso(const std::string& s) {
    struct tm tmv = {};
    if (sscanf(s.c_str(), "%d-%d-%d %d:%d:%d", &tmv.tm_year, &tmv.tm_mon, &tmv.tm_mday,
               &tmv.tm_hour, &tmv.tm_min, &tmv.tm_sec) != 6) return -1;
    tmv.tm_year -= 1900; tmv.tm_mon -= 1;
    return (double)mktime(&tmv);
}

double RawMTime(const TString& date, int run) {
    const char* dirs[] = {"/home/precalkor/Data/RAW/Laser",
                          "/home/precalkor/external_HDD_1_4T/Data_Backup/RAW/Laser"};
    for (const char* d : dirs) {
        TString p = Form("%s/precal_raw_kor_run_%s_%03d.root", d, date.Data(), run);
        Long_t id, flags, mt; Long64_t sz;
        if (gSystem->GetPathInfo(p.Data(), &id, &sz, &flags, &mt) == 0) return (double)mt;
    }
    return -1;
}

FieldCond ReadField(const TString& tag) {
    FieldCond fc;
    // tag looks like "<date>_<runA>_<runB>"
    TObjArray* tk = tag.Tokenize("_");
    if (!tk || tk->GetEntries() < 3) { delete tk; return fc; }
    TString date = ((TObjString*)tk->At(0))->GetString();
    int runA = ((TObjString*)tk->At(1))->GetString().Atoi();
    int runB = ((TObjString*)tk->At(2))->GetString().Atoi();
    delete tk;

    double t0 = RawMTime(date, runA), t1 = RawMTime(date, runB);
    if (t0 < 0 || t1 < 0) return fc;

    std::ifstream f("/home/precalkor/ADC/ADC_test/log_bfield_20260627.csv");
    if (!f.is_open()) return fc;
    std::string line; bool header = true;
    std::vector<double> I[4], B[3];
    while (std::getline(f, line)) {
        if (header) { header = false; continue; }
        std::stringstream ss(line); std::string cell; std::vector<std::string> c;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        if (c.size() < 21) continue;
        double t = ParseIso(c[0]);
        if (t < t0 || t > t1) continue;
        for (int k = 0; k < 4; ++k) { try { I[k].push_back(std::stod(c[5 + k])); } catch (...) {} }
        for (int m = 0; m < 3; ++m) {
            try {
                double x = std::stod(c[9 + m*3]), y = std::stod(c[10 + m*3]), z = std::stod(c[11 + m*3]);
                B[m].push_back(std::sqrt(x*x + y*y + z*z) * 1000.0);   // G -> mG
            } catch (...) {}
        }
    }
    // mean AND spread: a set whose coil current merely averaged the same but
    // wandered during the scan is not the same condition as a steady one.
    auto meanSd = [](const std::vector<double>& v, double& m, double& sd) {
        m = sd = 0; if (v.empty()) return;
        for (double x : v) m += x; m /= v.size();
        if (v.size() < 2) return;
        for (double x : v) sd += (x - m) * (x - m);
        sd = std::sqrt(sd / (v.size() - 1));
    };
    for (int k = 0; k < 4; ++k) meanSd(I[k], fc.I[k], fc.Isd[k]);
    for (int m = 0; m < 3; ++m) meanSd(B[m], fc.B[m], fc.Bsd[m]);
    fc.n = (int)I[0].size();
    fc.ok = fc.n > 0;
    return fc;
}

}  // namespace

void Draw_Reproducibility_Multi(std::vector<TString> tags,
                                std::vector<TString> serialList = {"EM6400", "EL5150"},
                                TString metric = "MonNorm_CorrRawQE",
                                double angleTol = 1.0,
                                std::vector<TString> setLabels = {}) {
    if (tags.size() < 2) { std::cout << "[ERROR] need >= 2 tags (first is the reference)." << std::endl; return; }
    const int nSets = (int)tags.size();
    if ((int)setLabels.size() != nSets) {
        setLabels.clear();
        for (int i = 0; i < nSets; ++i) setLabels.push_back(Form("Set %d", i + 1));
    }

    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(22, "");
    gStyle->SetTitleFont(132, "XY"); gStyle->SetLabelFont(132, "XY");
    gStyle->SetGridColor(kGray + 1);

    const char* axes[2] = {"X", "Y"};
    TString out = Form("./Data/image/Uniformity/ReproducibilityMulti_%s_ref.pdf", tags[0].Data());
    TCanvas* c = new TCanvas("cRM", "Reproducibility (multi)", 1600, 950);
    c->Print(out + "[");

    std::vector<double> allPulls;
    // Rows for the summary table: serial, axis, set, stats.
    struct Row { TString serial; const char* axis; TString set; ReproStats st; };
    std::vector<Row> rows;

    for (int a = 0; a < 2; ++a) {
        c->Clear(); c->cd();
        TLatex* pageTitle = new TLatex(0.5, 0.975,
            Form("%s-axis   Reproducibility vs %s   (metric: %s)",
                 axes[a], setLabels[0].Data(), metric.Data()));
        pageTitle->SetNDC(); pageTitle->SetTextAlign(22); pageTitle->SetTextFont(22);
        pageTitle->SetTextSize(0.030); pageTitle->Draw();

        bool drewAny = false;
        const int nSer = (int)serialList.size();
        for (int s = 0; s < nSer; ++s) {
            TString serial = serialList[s];
            SetCurve ref = LoadCurve(tags[0], serial, axes[a], metric);
            if (ref.ang.empty()) continue;
            // Graph points come back in acquisition order, not angle order, so
            // front()/back() are NOT the axis limits -- take real extrema.
            double angLo = *std::min_element(ref.ang.begin(), ref.ang.end());
            double angHi = *std::max_element(ref.ang.begin(), ref.ang.end());

            // Back to the canvas before creating this PMT's pads: the previous
            // iteration left pR current, so without this the next TPad is
            // constructed as a CHILD of that pad and lands as a shrunken inset
            // inside it instead of beside it (2026-08-25).
            c->cd();
            double x0 = 0.02 + s * (0.96 / nSer), x1 = x0 + (0.96 / nSer) - 0.02;
            TPad* pV = new TPad(Form("pV%d%d", a, s), "", x0, 0.42, x1, 0.94);
            TPad* pR = new TPad(Form("pR%d%d", a, s), "", x0, 0.06, x1, 0.42);
            pV->SetBottomMargin(0.02); pV->SetLeftMargin(0.16); pV->SetRightMargin(0.04); pV->SetTopMargin(0.06);
            pR->SetTopMargin(0.02);    pR->SetLeftMargin(0.16); pR->SetRightMargin(0.04); pR->SetBottomMargin(0.30);
            pV->SetGridx(); pV->SetGridy(); pR->SetGridx(); pR->SetGridy();
            pV->Draw(); pR->Draw();
            drewAny = true;

            // ---------- value pad: every set overlaid ----------
            pV->cd();
            double lo = 1e30, hi = -1e30;
            std::vector<SetCurve> curves;
            for (int k = 0; k < nSets; ++k) {
                SetCurve cv = (k == 0) ? ref : LoadCurve(tags[k], serial, axes[a], metric);
                curves.push_back(cv);
                for (size_t i = 0; i < cv.val.size(); ++i) {
                    lo = std::min(lo, cv.val[i] - cv.err[i]);
                    hi = std::max(hi, cv.val[i] + cv.err[i]);
                }
            }
            double span = (hi > lo) ? hi - lo : 1.0;
            // 55% headroom: the legend sits in the top-right of this pad, so
            // the data must be pushed clear of it rather than drawn under it.
            TH1F* fr = pV->DrawFrame(angLo - 5, lo - span * 0.10,
                                     angHi + 5, hi + span * 0.55);
            fr->SetTitle(Form("%s, %s-axis;;Monitor-Normalized QE [a.u.]", serial.Data(), axes[a]));
            fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetLabelSize(0.046);
            fr->GetYaxis()->SetTitleOffset(1.35); fr->GetXaxis()->SetLabelSize(0);
            fr->GetYaxis()->SetNdivisions(508);

            TLegend* leg = new TLegend(0.62, 0.74, 0.955, 0.935);
            leg->SetBorderSize(1); leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
            leg->SetTextFont(132); leg->SetTextSize(0.038);
            leg->SetNColumns(2);   // 4 sets in 2 columns -> half the height

            for (int k = 0; k < nSets; ++k) {
                if (curves[k].ang.empty()) continue;
                TGraphErrors* g = new TGraphErrors(curves[k].ang.size(), &curves[k].ang[0],
                                                   &curves[k].val[0], 0, &curves[k].err[0]);
                g->SetMarkerStyle(kSetMarker[k % 6]); g->SetMarkerColor(kSetColor[k % 6]);
                g->SetLineColor(kSetColor[k % 6]); g->SetMarkerSize(0.8);
                g->Draw("PZ same");
                leg->AddEntry(g, setLabels[k], "lp");
            }
            leg->Draw();

            // ---------- ratio pad: each set / reference ----------
            pR->cd();
            // Ratio range from the data (with a floor), not a hardcoded
            // +-10%: a set that really does sit outside that window was
            // silently clipped off the pad before.
            double rLo = 1.0, rHi = 1.0;
            for (int k = 1; k < nSets; ++k) {
                for (auto& p2 : Match(ref, curves[k], angleTol)) {
                    rLo = std::min(rLo, p2.ratio - p2.ratioErr);
                    rHi = std::max(rHi, p2.ratio + p2.ratioErr);
                }
            }
            double rPad = std::max(0.02, (rHi - rLo) * 0.12);
            TH1F* fr2 = pR->DrawFrame(angLo - 5, rLo - rPad, angHi + 5, rHi + rPad);
            fr2->SetTitle(";Position angle [degree];Set / Ref");
            // pV spans 0.52 of the canvas height, pR only 0.36 -- text size is
            // relative to pad height, so to look the same on paper pR's sizes
            // must be scaled by 0.52/0.36 ~= 1.44 relative to pV's.
            const double kR = 0.52 / 0.36;
            fr2->GetXaxis()->SetTitleSize(0.055 * kR); fr2->GetXaxis()->SetLabelSize(0.046 * kR);
            fr2->GetXaxis()->SetTitleOffset(1.05);
            fr2->GetYaxis()->SetTitleSize(0.055 * kR); fr2->GetYaxis()->SetLabelSize(0.046 * kR);
            fr2->GetYaxis()->SetTitleOffset(0.95); fr2->GetYaxis()->SetNdivisions(505);

            for (int k = 1; k < nSets; ++k) {
                std::vector<MatchedPoint> pts = Match(ref, curves[k], angleTol);
                if (pts.size() < 2) continue;
                ReproStats st = ComputeRepro(pts);
                rows.push_back({serial, axes[a], setLabels[k], st});
                for (auto& p : pts) allPulls.push_back(p.pull);

                // Small sideways nudge so coincident error bars stay readable.
                // Kept well under the ~7 deg spacing between real scan points
                // (was 1.2 deg/step, enough that points looked like they were
                // measured at DIFFERENT angles rather than jittered) -- marker
                // shape already separates the sets, this only unstacks the bars.
                double dx = (k - (nSets - 1) / 2.0) * 0.45;
                std::vector<double> ax, av, ae;
                for (auto& p : pts) { ax.push_back(p.angle + dx); av.push_back(p.ratio); ae.push_back(p.ratioErr); }
                TGraphErrors* gr = new TGraphErrors(ax.size(), &ax[0], &av[0], 0, &ae[0]);
                gr->SetMarkerStyle(kSetMarker[k % 6]); gr->SetMarkerColor(kSetColor[k % 6]);
                gr->SetLineColor(kSetColor[k % 6]); gr->SetMarkerSize(0.7);
                gr->Draw("PZ same");
            }
            TLine* l1 = new TLine(angLo - 5, 1.0, angHi + 5, 1.0);
            l1->SetLineColor(kGray + 2); l1->SetLineWidth(2); l1->Draw();
        }
        if (drewAny) c->Print(out);
    }

    // ---------------- summary table page ----------------
    {
        c->Clear(); c->cd();
        TLatex* t = new TLatex(0.5, 0.95, Form("Reproducibility summary  (reference: %s)", setLabels[0].Data()));
        t->SetNDC(); t->SetTextAlign(22); t->SetTextFont(22); t->SetTextSize(0.030); t->Draw();

        TPaveText* tab = new TPaveText(0.06, 0.08, 0.94, 0.90, "NDC");
        tab->SetBorderSize(1); tab->SetFillColor(kWhite); tab->SetFillStyle(1001);
        tab->SetTextFont(102);   // monospace, so the columns line up
        tab->SetTextSize(0.022); tab->SetTextAlign(12);
        tab->AddText(" PMT      axis  set        N   bias[%]   RMS(pull)      systematic");
        tab->AddText(" ------------------------------------------------------------------------");
        for (auto& r : rows) {
            TString sysTxt;
            if (r.st.sysResolved)          sysTxt = Form("%.2f %%", r.st.sysComp * 100);
            else if (r.st.statOverestimated) sysTxt = Form("< %.2f %% (stat overest.)", r.st.sysUL * 100);
            else                            sysTxt = Form("< %.2f %% (95%% CL)", r.st.sysUL * 100);
            tab->AddText(Form(" %-8s  %-4s  %-8s  %3d  %+6.2f   %.2f+-%.2f   %s",
                              r.serial.Data(), r.axis, r.set.Data(), r.st.n,
                              (r.st.meanRatio - 1) * 100, r.st.pullRMS,
                              r.st.pullRMS / std::sqrt(2.0 * r.st.n), sysTxt.Data()));
        }
        tab->Draw();
        c->Print(out);
    }

    // ---------------- B-field conditions page ----------------
    {
        c->Clear(); c->cd();
        TLatex* t = new TLatex(0.5, 0.95, "Field conditions during each set");
        t->SetNDC(); t->SetTextAlign(22); t->SetTextFont(22); t->SetTextSize(0.030); t->Draw();

        TPaveText* tab = new TPaveText(0.06, 0.30, 0.94, 0.90, "NDC");
        tab->SetBorderSize(1); tab->SetFillColor(kWhite); tab->SetFillStyle(1001);
        tab->SetTextFont(102); tab->SetTextSize(0.024); tab->SetTextAlign(12);
        tab->AddText(" set        tag                  N   I1 [A]        I2 [A]        I3 [A]        |B|1 [mG]     |B|3 [mG]");
        tab->AddText(" ------------------------------------------------------------------------------------------------------");
        bool anyField = false;
        for (int k = 0; k < nSets; ++k) {
            FieldCond fc = ReadField(tags[k]);
            if (!fc.ok) {
                tab->AddText(Form(" %-9s  %-18s  (no B-field log coverage)", setLabels[k].Data(), tags[k].Data()));
                continue;
            }
            anyField = true;
            tab->AddText(Form(" %-9s  %-18s %4d %5.2f+-%-5.2f %5.2f+-%-5.2f %5.2f+-%-5.2f %5.1f+-%-5.1f %5.1f+-%-5.1f",
                              setLabels[k].Data(), tags[k].Data(), fc.n,
                              fc.I[0], fc.Isd[0], fc.I[1], fc.Isd[1], fc.I[2], fc.Isd[2],
                              fc.B[0], fc.Bsd[0], fc.B[2], fc.Bsd[2]));
        }
        tab->Draw();

        TPaveText* note = new TPaveText(0.06, 0.12, 0.94, 0.27, "NDC");
        note->SetBorderSize(0); note->SetFillStyle(0); note->SetTextFont(132);
        note->SetTextSize(0.022); note->SetTextAlign(12);
        if (anyField) {
            note->AddText("Mean over each set's acquisition window (RAW file mtimes).");
            note->AddText("Sets whose currents/|B| differ were NOT taken under the same conditions --");
            note->AddText("a reproducibility difference between them is not necessarily a PMT effect.");
        } else {
            note->AddText("No B-field log rows fell inside any set's window -- conditions unverified.");
        }
        note->Draw();
        c->Print(out);
    }

    // ---------------- combined pull page ----------------
    if (allPulls.size() >= 2) {
        c->Clear(); c->cd();
        gPad->SetGrid(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.13);

        TH1F* h = new TH1F("hPullMulti",
            "Pull distribution, all sets vs reference   (Set - Ref) / #sqrt{#sigma_{ref}^{2}+#sigma_{set}^{2}};Pull;Entries",
            25, -6, 6);
        double sum2 = 0;
        for (double p : allPulls) { h->Fill(p); sum2 += p * p; }
        double rms = std::sqrt(sum2 / allPulls.size());
        h->SetLineColor(kBlue + 2); h->SetLineWidth(2); h->SetFillColorAlpha(kBlue + 2, 0.25);
        h->SetMaximum(h->GetMaximum() * 1.45);
        h->Draw("hist");

        TF1* ref = new TF1("pullRefM", "[0]*TMath::Gaus(x,0,1,kTRUE)", -6, 6);
        ref->SetParameter(0, allPulls.size() * h->GetBinWidth(1));
        ref->SetLineColor(kGray + 2); ref->SetLineStyle(2); ref->SetLineWidth(2);
        ref->Draw("same");

        TLegend* lg = new TLegend(0.16, 0.80, 0.40, 0.88);
        lg->SetBorderSize(1); lg->SetFillColor(kWhite); lg->SetTextFont(132); lg->SetTextSize(0.032);
        lg->AddEntry(h, "Pull", "f"); lg->AddEntry(ref, "Gaussian", "l");
        lg->Draw();

        TPaveText* pb = new TPaveText(0.62, 0.74, 0.93, 0.87, "NDC");
        pb->SetBorderSize(2); pb->SetFillColor(kWhite); pb->SetFillStyle(1001);
        pb->SetTextFont(132); pb->SetTextSize(0.035); pb->SetTextAlign(12);
        pb->AddText(Form("N = %d", (int)allPulls.size()));
        pb->AddText(Form("RMS(pull) = %.2f #pm %.2f", rms, rms / std::sqrt(2.0 * allPulls.size())));
        pb->Draw();
        c->Print(out);

        std::cout << Form("\n[ All sets vs %s ]  N = %d,  RMS(pull) = %.2f",
                          setLabels[0].Data(), (int)allPulls.size(), rms) << std::endl;
    }

    c->Print(out + "]");

    for (auto& r : rows) {
        TString sysTxt;
        if (r.st.sysResolved)            sysTxt = Form("%.2f %%", r.st.sysComp * 100);
        else if (r.st.statOverestimated) sysTxt = Form("< %.2f %% (stat overest.)", r.st.sysUL * 100);
        else                             sysTxt = Form("< %.2f %% (95%% CL)", r.st.sysUL * 100);
        std::cout << Form("  %-8s %-2s %-8s  N=%3d  bias %+6.2f %%  RMS(pull) %.2f  syst %s",
                          r.serial.Data(), r.axis, r.set.Data(), r.st.n,
                          (r.st.meanRatio - 1) * 100, r.st.pullRMS, sysTxt.Data()) << std::endl;
    }
    std::cout << "\n[INFO] Saved: " << out << std::endl;
}
