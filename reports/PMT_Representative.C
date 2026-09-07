// Representative per-PMT numbers from repeated uniformity measurements, plus
// the monitoring conditions each measurement was taken under.
//
// One scan gives one C/S ratio; several scans of the same PMT give a spread.
// The representative value is the mean over repeats, and the spread ACROSS
// repeats is the reproducibility -- which we showed (2026-08) is the dominant
// uncertainty here, larger than the per-point statistical error. So:
//
//   value  = (1/R) * sum_r  x_r
//   sig_repro = sqrt( 1/(R-1) * sum_r (x_r - value)^2 )        [R >= 2]
//   err(value) = sig_repro / sqrt(R)      <- error on the mean
//
// With R = 1 there is no reproducibility estimate at all; that is reported as
// such rather than as "+- 0", which would claim a precision one measurement
// cannot support (same reasoning as the systematic upper limits in
// repro_common.h).
//
// Monitoring conditions (Dark Box T/H, laser TEC temp + drive, coil current,
// |B|) are averaged over each measurement's own acquisition window, so the
// table answers "was this PMT measured under the same conditions as that
// one?" without leaving the page. Previously only a whole-campaign average
// existed, which cannot answer that.
//
// Usage:
//   root -l -b -q 'PMT_Representative.C()'                    // all PMTs
//   root -l -b -q 'PMT_Representative.C("MonNorm_CorrRawQE")' // pick metric
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>

namespace {

struct Meas { TString tag; double cs, sig; int npts; };
struct Key  { TString serial, axis; bool operator<(const Key& o) const {
                  return serial != o.serial ? serial < o.serial : axis < o.axis; } };

struct Rep {
    int    R = 0;
    double cs = 0, csRepro = 0, csErr = 0;
    double sig = 0, sigRepro = 0, sigErr = 0;
    bool   hasRepro = false;
};

Rep Summarize(const std::vector<Meas>& ms) {
    Rep r; r.R = (int)ms.size();
    if (r.R == 0) return r;
    for (auto& m : ms) { r.cs += m.cs; r.sig += m.sig; }
    r.cs /= r.R; r.sig /= r.R;
    if (r.R >= 2) {
        double a = 0, b = 0;
        for (auto& m : ms) { a += (m.cs - r.cs) * (m.cs - r.cs); b += (m.sig - r.sig) * (m.sig - r.sig); }
        r.csRepro  = std::sqrt(a / (r.R - 1));
        r.sigRepro = std::sqrt(b / (r.R - 1));
        r.csErr  = r.csRepro  / std::sqrt((double)r.R);
        r.sigErr = r.sigRepro / std::sqrt((double)r.R);
        r.hasRepro = true;
    }
    return r;
}

// ---- monitoring conditions over one tag's acquisition window --------------
struct Cond {
    bool ok = false;
    double boxT = 0, boxTsd = 0, boxH = 0, boxHsd = 0;
    double ldT  = 0, ldTsd  = 0, ldP = 0;
    double coilI = 0, coilIsd = 0, bmag = 0, bmagsd = 0;
    int nEnv = 0, nLaser = 0, nField = 0;
};

double ParseIso(const std::string& s) {
    struct tm tmv = {}; double sec = 0;
    // Handles both "YYYY-MM-DD HH:MM:SS" and ISO "YYYY-MM-DDTHH:MM:SS.ffffff"
    if (sscanf(s.c_str(), "%d-%d-%d%*c%d:%d:%lf", &tmv.tm_year, &tmv.tm_mon, &tmv.tm_mday,
               &tmv.tm_hour, &tmv.tm_min, &sec) != 6) return -1;
    tmv.tm_year -= 1900; tmv.tm_mon -= 1; tmv.tm_sec = (int)sec;
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

void MeanSd(const std::vector<double>& v, double& m, double& sd) {
    m = sd = 0; if (v.empty()) return;
    for (double x : v) m += x; m /= v.size();
    if (v.size() < 2) return;
    for (double x : v) sd += (x - m) * (x - m);
    sd = std::sqrt(sd / (v.size() - 1));
}

Cond ReadCond(const TString& tag) {
    Cond c;
    TObjArray* tk = tag.Tokenize("_");
    if (!tk || tk->GetEntries() < 3) { delete tk; return c; }
    TString date = ((TObjString*)tk->At(0))->GetString();
    int runA = ((TObjString*)tk->At(1))->GetString().Atoi();
    int runB = ((TObjString*)tk->At(2))->GetString().Atoi();
    delete tk;
    double t0 = RawMTime(date, runA), t1 = RawMTime(date, runB);
    if (t0 < 0 || t1 < 0) return c;

    // Dark Box temperature / humidity, from the HV monitoring sqlite via its
    // CLI (no TSQLiteServer dependency -- same approach Draw_Stability_v1 uses).
    {
        std::vector<double> T, H;
        TString lo = TString(Form("%s", TDatime((UInt_t)t0).AsSQLString())); lo.ReplaceAll(" ", "T");
        TString hi = TString(Form("%s", TDatime((UInt_t)t1).AsSQLString())); hi.ReplaceAll(" ", "T");
        TString cmd = Form("sqlite3 -separator , "
            "/home/precalkor/Integrated_Control_SW/HV_Control_SW/monitoring_log.db "
            "\"SELECT Dark_Box_1_T,Dark_Box_1_H FROM monitoring_data "
            "WHERE timestamp BETWEEN '%s' AND '%s';\"", lo.Data(), hi.Data());
        FILE* p = gSystem->OpenPipe(cmd, "r");
        if (p) {
            char buf[256];
            while (fgets(buf, sizeof(buf), p)) {
                double a, b;
                if (sscanf(buf, "%lf,%lf", &a, &b) == 2) { T.push_back(a); H.push_back(b); }
            }
            gSystem->ClosePipe(p);
        }
        MeanSd(T, c.boxT, c.boxTsd); MeanSd(H, c.boxH, c.boxHsd);
        c.nEnv = (int)T.size();
    }

    // Laser diode TEC temperature + pulse drive, per wavelength file.
    for (int wl : {405, 375, 450, 473}) {
        TString f = Form("/home/precalkor/ADC/ADC_test/LOG/LASER/laser_data_%dnm_%s.csv", wl, date.Data());
        std::ifstream in(f.Data());
        if (!in.is_open()) continue;
        std::string line; bool header = true;
        std::vector<double> T, P;
        while (std::getline(in, line)) {
            if (header) { header = false; continue; }
            std::stringstream ss(line); std::string cell; std::vector<std::string> col;
            while (std::getline(ss, cell, ',')) col.push_back(cell);
            if (col.size() < 3) continue;
            double t = ParseIso(col[0]);
            if (t < t0 || t > t1) continue;
            try { T.push_back(std::stod(col[1])); P.push_back(std::stod(col[2])); } catch (...) {}
        }
        if (!T.empty()) {
            double d;
            MeanSd(T, c.ldT, c.ldTsd); MeanSd(P, c.ldP, d);
            c.nLaser = (int)T.size();
            break;    // first wavelength with coverage is the one in use
        }
    }

    // Coil current + |B|, from the B-field CSV.
    {
        std::ifstream in("/home/precalkor/ADC/ADC_test/log_bfield_20260627.csv");
        std::string line; bool header = true;
        std::vector<double> I, B;
        while (std::getline(in, line)) {
            if (header) { header = false; continue; }
            std::stringstream ss(line); std::string cell; std::vector<std::string> col;
            while (std::getline(ss, cell, ',')) col.push_back(cell);
            if (col.size() < 21) continue;
            double t = ParseIso(col[0]);
            if (t < t0 || t > t1) continue;
            try { I.push_back(std::stod(col[5])); } catch (...) {}
            try {
                double x = std::stod(col[9]), y = std::stod(col[10]), z = std::stod(col[11]);
                B.push_back(std::sqrt(x*x + y*y + z*z) * 1000.0);
            } catch (...) {}
        }
        MeanSd(I, c.coilI, c.coilIsd); MeanSd(B, c.bmag, c.bmagsd);
        c.nField = (int)I.size();
    }

    c.ok = (c.nEnv > 0 || c.nField > 0 || c.nLaser > 0);
    return c;
}

}  // namespace

void PMT_Representative(TString metric = "MonNorm_CorrRawQE",
                        TString summaryFile = "./Data/summary/pmt_uniformity_summary.root") {
    TFile* f = TFile::Open(summaryFile);
    if (!f || f->IsZombie()) { std::cout << "[ERROR] cannot open " << summaryFile << std::endl; return; }
    TTree* t = (TTree*)f->Get("pmt_summary");
    if (!t) { std::cout << "[ERROR] no pmt_summary tree in " << summaryFile << std::endl; return; }

    char b_serial[32], b_tag[128], b_metric[64], b_axis[16];
    int b_npts; double b_cs, b_sig;
    t->SetBranchAddress("serial", b_serial); t->SetBranchAddress("tag", b_tag);
    t->SetBranchAddress("metric", b_metric); t->SetBranchAddress("axis", b_axis);
    t->SetBranchAddress("n_points", &b_npts);
    t->SetBranchAddress("cs_ratio", &b_cs);  t->SetBranchAddress("sigma1000", &b_sig);

    std::map<Key, std::vector<Meas>> byKey;
    std::set<TString> tagsSeen;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (metric != "" && metric != b_metric) continue;
        byKey[{b_serial, b_axis}].push_back({b_tag, b_cs, b_sig, b_npts});
        tagsSeen.insert(b_tag);
    }
    f->Close();
    if (byKey.empty()) { std::cout << "[ERROR] no rows for metric " << metric << std::endl; return; }

    // ---- console + CSV ----
    gSystem->mkdir("./Data/summary", kTRUE);
    std::ofstream out("./Data/summary/pmt_representative.csv");
    out << "serial,axis,metric,n_repeats,cs_ratio,cs_err,cs_reproducibility,"
           "sigma1000,sigma1000_err,sigma1000_reproducibility\n";

    std::cout << "\n=== Representative values (metric: " << metric << ") ===" << std::endl;
    std::cout << Form("%-8s %-7s %3s  %-22s  %-22s", "serial", "axis", "R",
                      "C/S ratio", "sigma x1000") << std::endl;
    std::cout << std::string(74, '-') << std::endl;
    for (auto& kv : byKey) {
        // X and Y only -- the same filter the plot below already applies. The
        // merged "pooled" entry duplicated what the two axis rows already say
        // and only added a third row per PMT to the GUI table (2026-08-26).
        // Filtered on read, not just at the writer, because summary rows
        // recorded before that change still carry axis=="pooled".
        if (kv.first.axis == "pooled") continue;
        Rep r = Summarize(kv.second);
        TString csTxt, sgTxt;
        if (r.hasRepro) {
            csTxt = Form("%.4f +- %.4f (rep %.4f)", r.cs, r.csErr, r.csRepro);
            sgTxt = Form("%.2f +- %.2f (rep %.2f)", r.sig, r.sigErr, r.sigRepro);
        } else {
            csTxt = Form("%.4f (single meas.)", r.cs);
            sgTxt = Form("%.2f (single meas.)", r.sig);
        }
        std::cout << Form("%-8s %-7s %3d  %-22s  %-22s",
                          kv.first.serial.Data(), kv.first.axis.Data(), r.R,
                          csTxt.Data(), sgTxt.Data()) << std::endl;
        out << kv.first.serial << "," << kv.first.axis << "," << metric << "," << r.R << ","
            << Form("%.4f", r.cs) << ",";
        if (r.hasRepro) out << Form("%.4f,%.4f,", r.csErr, r.csRepro); else out << ",,";
        out << Form("%.4f", r.sig) << ",";
        if (r.hasRepro) out << Form("%.4f,%.4f", r.sigErr, r.sigRepro); else out << ",";
        out << "\n";
    }
    out.close();
    std::cout << "\n[INFO] Saved: ./Data/summary/pmt_representative.csv" << std::endl;

    // ---- monitoring conditions per tag ----
    std::cout << "\n=== Monitoring conditions per measurement window ===" << std::endl;
    std::cout << Form("%-20s %-16s %-16s %-14s %-16s", "tag", "DarkBox T [C]",
                      "DarkBox H [%]", "LD TEC T [C]", "|B| [mG]") << std::endl;
    std::cout << std::string(88, '-') << std::endl;
    std::map<TString, Cond> conds;
    for (const TString& tg : tagsSeen) {
        Cond c = ReadCond(tg);
        conds[tg] = c;
        if (!c.ok) { std::cout << Form("%-20s  (no monitoring coverage)", tg.Data()) << std::endl; continue; }
        std::cout << Form("%-20s %6.2f+-%-8.2f %6.2f+-%-8.2f %6.3f+-%-6.3f %6.1f+-%-8.1f",
                          tg.Data(), c.boxT, c.boxTsd, c.boxH, c.boxHsd,
                          c.ldT, c.ldTsd, c.bmag, c.bmagsd) << std::endl;
    }

    // ---- PDF: table + per-PMT comparison plot ----
    gStyle->SetOptStat(0); gStyle->SetTextFont(132);
    gStyle->SetGridColor(kGray + 1);
    gSystem->mkdir("./Data/image/Summary", kTRUE);
    TString pdf = "./Data/image/Summary/PMT_Representative.pdf";
    TCanvas* c = new TCanvas("cRep", "PMT representative", 1500, 900);
    c->Print(pdf + "[");

    // page 1: representative table
    {
        c->Clear(); c->cd();
        TLatex* h = new TLatex(0.5, 0.955, Form("Representative per-PMT values   (metric: %s)", metric.Data()));
        h->SetNDC(); h->SetTextAlign(22); h->SetTextFont(22); h->SetTextSize(0.028); h->Draw();
        TPaveText* p = new TPaveText(0.05, 0.10, 0.95, 0.92, "NDC");
        p->SetBorderSize(1); p->SetFillColor(kWhite); p->SetTextFont(102);
        p->SetTextSize(0.021); p->SetTextAlign(12);
        p->AddText("  serial    axis     R    C/S ratio                    sigma x1000");
        p->AddText("  ---------------------------------------------------------------------------");
        for (auto& kv : byKey) {
            Rep r = Summarize(kv.second);
            if (r.hasRepro)
                p->AddText(Form("  %-8s  %-7s %3d   %.4f +- %.4f (rep %.4f)   %.2f +- %.2f (rep %.2f)",
                                kv.first.serial.Data(), kv.first.axis.Data(), r.R,
                                r.cs, r.csErr, r.csRepro, r.sig, r.sigErr, r.sigRepro));
            else
                p->AddText(Form("  %-8s  %-7s %3d   %.4f  (single meas.)          %.2f  (single meas.)",
                                kv.first.serial.Data(), kv.first.axis.Data(), r.R, r.cs, r.sig));
        }
        p->AddText("");
        p->AddText("  R = number of repeat scans.  'rep' = spread ACROSS repeats (reproducibility);");
        p->AddText("  the quoted +- is the error on the mean, sig_repro/sqrt(R).");
        p->Draw();
        c->Print(pdf);
    }

    // page 2: conditions table
    {
        c->Clear(); c->cd();
        TLatex* h = new TLatex(0.5, 0.955, "Monitoring conditions per measurement window");
        h->SetNDC(); h->SetTextAlign(22); h->SetTextFont(22); h->SetTextSize(0.028); h->Draw();
        TPaveText* p = new TPaveText(0.03, 0.15, 0.97, 0.92, "NDC");
        p->SetBorderSize(1); p->SetFillColor(kWhite); p->SetTextFont(102);
        p->SetTextSize(0.019); p->SetTextAlign(12);
        p->AddText("  tag                  DarkBox T [C]    DarkBox H [%]    LD TEC T [C]    coil I [A]      |B| [mG]");
        p->AddText("  ------------------------------------------------------------------------------------------------------");
        for (auto& kv : conds) {
            const Cond& cd = kv.second;
            if (!cd.ok) { p->AddText(Form("  %-20s (no monitoring coverage)", kv.first.Data())); continue; }
            p->AddText(Form("  %-20s %6.2f+-%-6.2f  %6.2f+-%-6.2f  %6.3f+-%-5.3f  %5.2f+-%-5.2f  %6.1f+-%-6.1f",
                            kv.first.Data(), cd.boxT, cd.boxTsd, cd.boxH, cd.boxHsd,
                            cd.ldT, cd.ldTsd, cd.coilI, cd.coilIsd, cd.bmag, cd.bmagsd));
        }
        p->AddText("");
        p->AddText("  Mean +- stddev over each measurement's own acquisition window (not a campaign-wide average),");
        p->AddText("  so two PMTs can be checked for having been measured under comparable conditions.");
        p->Draw();
        c->Print(pdf);
    }

    // page 3: each measurement vs the PMT's own representative mean
    {
        c->Clear(); c->cd();
        gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.22);

        std::vector<TString> labels;
        std::vector<double> xs, ys, es;
        std::vector<double> mx, my, me;
        int ix = 0;
        for (auto& kv : byKey) {
            if (kv.first.axis == "pooled") continue;    // X and Y only, pooled duplicates them
            Rep r = Summarize(kv.second);
            for (auto& m : kv.second) { xs.push_back(ix); ys.push_back(m.cs); es.push_back(0); }
            mx.push_back(ix); my.push_back(r.cs); me.push_back(r.hasRepro ? r.csRepro : 0);
            labels.push_back(Form("%s %s", kv.first.serial.Data(), kv.first.axis.Data()));
            ++ix;
        }
        if (ix > 0) {
            double lo = *std::min_element(ys.begin(), ys.end());
            double hi = *std::max_element(ys.begin(), ys.end());
            double pad = std::max(0.01, (hi - lo) * 0.35);
            TH1F* fr = gPad->DrawFrame(-0.5, lo - pad, ix - 0.5, hi + pad);
            fr->SetTitle("Individual measurements vs representative mean;;C/S ratio");
            fr->GetYaxis()->SetTitleSize(0.040); fr->GetYaxis()->SetLabelSize(0.034);
            for (int i = 0; i < ix; ++i)
                fr->GetXaxis()->SetBinLabel(fr->GetXaxis()->FindBin(i), labels[i]);
            fr->GetXaxis()->LabelsOption("v");
            fr->GetXaxis()->SetLabelSize(0.030);

            TGraphErrors* gAll = new TGraphErrors(xs.size(), &xs[0], &ys[0], 0, &es[0]);
            gAll->SetMarkerStyle(24); gAll->SetMarkerColor(kGray + 2); gAll->SetMarkerSize(1.0);
            gAll->Draw("P same");
            TGraphErrors* gMean = new TGraphErrors(mx.size(), &mx[0], &my[0], 0, &me[0]);
            gMean->SetMarkerStyle(20); gMean->SetMarkerColor(kRed + 1);
            gMean->SetLineColor(kRed + 1); gMean->SetLineWidth(2); gMean->SetMarkerSize(1.2);
            gMean->Draw("PZ same");

            TLegend* lg = new TLegend(0.62, 0.80, 0.95, 0.90);
            lg->SetBorderSize(1); lg->SetFillColor(kWhite); lg->SetTextFont(132); lg->SetTextSize(0.030);
            lg->AddEntry(gAll, "individual scans", "p");
            lg->AddEntry(gMean, "mean #pm reproducibility", "lp");
            lg->Draw();
        }
        c->Print(pdf);
    }

    c->Print(pdf + "]");
    std::cout << "[INFO] Saved: " << pdf << std::endl;
}
