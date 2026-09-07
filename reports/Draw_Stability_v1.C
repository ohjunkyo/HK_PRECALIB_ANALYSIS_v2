// Usage:
//   root -l -b -q 'Draw_Stability_v1.C("20260813", 800, 899)'
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TDatime.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TString.h>
#include <sys/stat.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

namespace {

struct Point {
    double t;            
    double epoch;        
    double qeRaw, qeRawErr;          // relativeQE_raw        (RawQE)
    double qeCorr, qeCorrErr;        // relativeQE            (CorrRawQE)
    double spe, speErr;              // spe_mean              (RawSPE / Gain)
    double chargeRes, chargeResErr;  // charge_resolution     (ChargeResolution)
    double tts, ttsErr;              // rms_exG * 2           (TTS, ns)
    double darkRate, darkRateErr;    // NoiseCountRate_ch%d   (Dark Rate, Hz) -- no stored err, kept 0
    bool   ok;
};

struct ChannelSeries {
    TString label;
    int     color, marker;
    std::vector<Point> pts;
};

TString EpochToDbTimestamp(double epoch) {
    time_t t = (time_t)epoch;
    struct tm tmv;
    localtime_r(&t, &tmv);
    char buf[32];
    strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tmv);
    return TString(buf);
}

struct EnvPoint { double epoch; double db1t, db1h, db2t, db2h; };

std::vector<EnvPoint> ReadEnvLog(double t0, double t1) {
    std::vector<EnvPoint> out;
    const char* dbPath = "/home/precalkor/Integrated_Control_SW/HV_Control_SW/monitoring_log.db";
    if (gSystem->AccessPathName(dbPath)) {
        std::cout << "[WARNING] monitoring_log.db not found; skipping Temp/Humidity pages." << std::endl;
        return out;
    }
    TString lo = EpochToDbTimestamp(t0), hi = EpochToDbTimestamp(t1);
    TString cmd = Form("sqlite3 -separator , \"%s\" "
                       "\"SELECT timestamp, Dark_Box_1_T, Dark_Box_1_H, Dark_Box_2_T, Dark_Box_2_H "
                       "FROM monitoring_data WHERE timestamp BETWEEN '%s' AND '%s' ORDER BY timestamp;\"",
                       dbPath, lo.Data(), hi.Data());
    FILE* p = gSystem->OpenPipe(cmd, "r");
    if (!p) return out;
    char line[256];
    while (fgets(line, sizeof(line), p)) {
        struct tm tmv = {};
        double t1v, t2v, h1v, h2v;
        char tsbuf[40];
        if (sscanf(line, "%39[^,],%lf,%lf,%lf,%lf", tsbuf, &t1v, &h1v, &t2v, &h2v) != 5) continue;
        if (sscanf(tsbuf, "%d-%d-%dT%d:%d:%d", &tmv.tm_year, &tmv.tm_mon, &tmv.tm_mday,
                   &tmv.tm_hour, &tmv.tm_min, &tmv.tm_sec) != 6) continue;
        tmv.tm_year -= 1900; tmv.tm_mon -= 1; tmv.tm_isdst = -1;
        double epoch = (double)mktime(&tmv);
        out.push_back({epoch, t1v, h1v, t2v, h2v});
    }
    gSystem->ClosePipe(p);
    return out;
}

static const char* kBfieldCSV = "/home/precalkor/ADC/ADC_test/log_bfield_20260627.csv";

static double DataConditionParseTime(const char* timeStr) {
    int year, month, day, hour, minute, second;
    if (sscanf(timeStr, "%d-%d-%d%*c%d:%d:%d", &year, &month, &day, &hour, &minute, &second) != 6) return -1;
    return (double)TDatime(year, month, day, hour, minute, second).Convert();
}

static void DataConditionLoadMonCSV(const TString& file, const std::vector<TString>& want,
        double t0, double t1, std::vector<TGraph*>& out) {
    for (size_t i = 0; i < want.size(); ++i) out.push_back(new TGraph());
    std::ifstream f(file.Data()); if (!f.is_open()) return;
    std::string line; bool header = true; std::vector<int> idx(want.size(), -1);
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string cell; std::vector<std::string> c;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        if (header) {
            for (size_t i = 0; i < c.size(); ++i)
                for (size_t w = 0; w < want.size(); ++w) if (TString(c[i].c_str()) == want[w]) idx[w] = i;
            header = false; continue;
        }
        if (c.empty()) continue;
        double tt = DataConditionParseTime(c[0].c_str());
        if (tt < t0 || tt > t1) continue;
        for (size_t w = 0; w < want.size(); ++w) {
            int i = idx[w]; if (i < 0 || i >= (int)c.size()) continue;
            if (c[i].empty() || c[i].find("N/A") != std::string::npos) continue;
            out[w]->SetPoint(out[w]->GetN(), tt, atof(c[i].c_str()));
        }
    }
}

static void DataConditionLoadBfieldMag(const TString& file, double t0, double t1, TGraph* out[4]) {
    for (int k = 0; k < 4; ++k) out[k] = new TGraph();
    std::ifstream f(file.Data()); if (!f.is_open()) return;
    std::string line; bool header = true;
    int iT = -1, iX[4] = {-1,-1,-1,-1}, iY[4] = {-1,-1,-1,-1}, iZ[4] = {-1,-1,-1,-1};
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string cell; std::vector<std::string> c;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        if (header) {
            for (size_t i = 0; i < c.size(); ++i) {
                if (c[i] == "timestamp") iT = i;
                for (int m = 0; m < 4; ++m) {
                    if (c[i] == std::string(Form("M%d_X", m + 1))) iX[m] = i;
                    else if (c[i] == std::string(Form("M%d_Y", m + 1))) iY[m] = i;
                    else if (c[i] == std::string(Form("M%d_Z", m + 1))) iZ[m] = i;
                }
            }
            header = false; continue;
        }
        if (iT < 0 || (int)c.size() <= iT) continue;
        double tt = DataConditionParseTime(c[iT].c_str());
        if (tt < t0 || tt > t1) continue;
        for (int m = 0; m < 4; ++m) {
            if (iX[m] < 0 || iY[m] < 0 || iZ[m] < 0) continue;
            if ((int)c.size() <= iX[m] || (int)c.size() <= iY[m] || (int)c.size() <= iZ[m]) continue;
            if (c[iX[m]].find("N/A") != std::string::npos || c[iY[m]].find("N/A") != std::string::npos ||
                    c[iZ[m]].find("N/A") != std::string::npos) continue;
            double x = atof(c[iX[m]].c_str()), y = atof(c[iY[m]].c_str()), z = atof(c[iZ[m]].c_str());
            double magMilliGauss = std::sqrt(x*x + y*y + z*z) * 1000.0;   // G -> mG
            out[m]->SetPoint(out[m]->GetN(), tt, magMilliGauss);
        }
    }
}

static void DataConditionMeanStd(TGraph* g, double& mean, double& stddev) {
    mean = 0; stddev = 0;
    int n = g ? g->GetN() : 0;
    if (n == 0) return;
    for (int i = 0; i < n; ++i) mean += g->GetY()[i];
    mean /= n;
    if (n < 2) return;
    for (int i = 0; i < n; ++i) { double dev = g->GetY()[i] - mean; stddev += dev*dev; }
    stddev = std::sqrt(stddev / n);
}

static double DataConditionReserveLegendHeadroom(TMultiGraph* mg, double topFrac = 0.55, double botFrac = 0.05) {
    double ymin = 1e18, ymax = -1e18;
    TList* gl = mg->GetListOfGraphs();
    if (gl) {
        for (TObject* obj : *gl) {
            TGraph* g = (TGraph*)obj;
            for (int i = 0; i < g->GetN(); ++i) {
                double y = g->GetY()[i];
                ymin = std::min(ymin, y); ymax = std::max(ymax, y);
            }
        }
    }
    if (ymin > ymax) return 0.75;
    double range = ymax - ymin; if (range <= 0) range = std::abs(ymax) * 0.1 + 1e-6;
    mg->SetMinimum(ymin - range * botFrac);
    mg->SetMaximum(ymax + range * topFrac);
    double topM = gPad->GetTopMargin(), botM = gPad->GetBottomMargin();
    double frameLo = botM, frameHi = 1.0 - topM;
    double totalSpan = range * (1.0 + topFrac + botFrac);
    double dataTopFracOfAxis = (ymax - (ymin - range * botFrac)) / totalSpan;
    return frameLo + dataTopFracOfAxis * (frameHi - frameLo);
}

static void DataConditionStyleTimeAxis(TMultiGraph* mg) {
    TAxis* ax = mg->GetXaxis();
    ax->SetTimeDisplay(1); ax->SetTimeFormat("%m-%d %H:%M%F1970-01-01 00:00:00");
    ax->SetNdivisions(510); ax->SetLabelSize(0.05); ax->SetTitleSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05); mg->GetYaxis()->SetTitleSize(0.058); mg->GetYaxis()->SetTitleOffset(0.85);
}


static void DataConditionStyleTimeAxisHM(TMultiGraph* mg, int ndiv = 512) {
    TAxis* ax = mg->GetXaxis();
    ax->SetTimeDisplay(1); ax->SetTimeFormat("%H:%M%F1970-01-01 00:00:00");
    ax->SetNdivisions(ndiv); ax->SetLabelSize(0.05); ax->SetTitleSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05); mg->GetYaxis()->SetTitleSize(0.058); mg->GetYaxis()->SetTitleOffset(0.85);
}

// Mean / sample RMS of a metric, ignoring points the fit flagged invalid.
void Stats(const std::vector<double>& v, double& mean, double& rms) {
    mean = rms = 0;
    if (v.empty()) return;
    for (double x : v) mean += x;
    mean /= v.size();
    if (v.size() < 2) return;
    double s2 = 0;
    for (double x : v) s2 += (x - mean) * (x - mean);
    rms = std::sqrt(s2 / (v.size() - 1));
}

std::vector<double> ResampleUniform(const std::vector<double>& t, const std::vector<double>& y,
                                    const std::vector<double>& grid) {
    TGraph g(t.size(), t.data(), y.data());
    std::vector<double> out(grid.size());
    for (size_t i = 0; i < grid.size(); ++i) out[i] = g.Eval(grid[i]);
    return out;
}

// Centered moving-average
std::vector<double> MovingAvgTrend(const std::vector<double>& y, int windowSteps) {
    int n = (int)y.size(), half = windowSteps / 2;
    std::vector<double> trend(n);
    for (int i = 0; i < n; ++i) {
        int lo = std::max(0, i - half), hi = std::min(n, i + half + 1);
        double s = 0; for (int k = lo; k < hi; ++k) s += y[k];
        trend[i] = s / (hi - lo);
    }
    return trend;
}

void Acf(const std::vector<double>& ym, double dt, int maxLagSteps,
        std::vector<double>& lags, std::vector<double>& rk) {
    int n = (int)ym.size();
    double denom = 0; for (double v : ym) denom += v * v;
    lags.push_back(0); rk.push_back(1.0);
    for (int k = 1; k <= maxLagSteps && k < n; ++k) {
        double num = 0; for (int i = 0; i < n - k; ++i) num += ym[i] * ym[i + k];
        lags.push_back(k * dt); rk.push_back(num / denom);
    }
}

}  // namespace

void Draw_Stability_v1(TString tag = "20260813", int run_start = 800, int run_end = 899,
                       TString resultDir = "./Data/FinalResult", double candidatePeriodMin = 40.2) {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132);
    gStyle->SetLabelFont(132, "XYZ");
    gStyle->SetTitleFont(132, "XYZ");
    gStyle->SetLegendFont(132);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetTimeOffset(0);
    gStyle->SetGridColor(kGray+1);

    // CH0 = monitor, CH1/CH2 = PMTs under test -- same channel convention as
    // read_ntp_v7 / Draw_Uniformity_Norm_v7.
    ChannelSeries ch[3] = {
        {"Mon.",  kBlack,   20, {}},
        {"Rot#1", kAzure+2, 21, {}},
        {"Rot#2", kRed+1,   22, {}},
    };

    double t0 = -1, tLast = -1;
    int nFiles = 0;

    for (int run = run_start; run <= run_end; ++run) {
        TString path = Form("%s/precal_result_kor_run_%s_%03d.root",
                            resultDir.Data(), tag.Data(), run);
        if (gSystem->AccessPathName(path)) continue;

   
        static const char* kRawDirsStab[] = {
            "./Data/RAW/Laser",
            "/home/precalkor/external_HDD_1_4T/Data_Backup/RAW/Laser",
        };
        double mtime = 0;
        struct stat st;
        for (const char* dir : kRawDirsStab) {
            TString rawPath = Form("%s/precal_raw_kor_run_%s_%03d.root", dir, tag.Data(), run);
            if (stat(rawPath.Data(), &st) == 0) { mtime = (double)st.st_mtime; break; }
        }
        if (mtime == 0) {
            if (stat(path.Data(), &st) == 0) mtime = (double)st.st_mtime;
            else continue;
        }

        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }

        if (t0 < 0) t0 = mtime;
        tLast = std::max(tLast, mtime);
        double tMin = (mtime - t0) / 60.0;

        bool any = false;
        for (int c = 0; c < 3; ++c) {
            TTree* tr = (TTree*)f->Get(Form("tree_ch%d", c));
            if (!tr) continue;

            double qeCorr = 0, qeCorrErr = 0, qeRaw = 0, qeRawErr = 0;
            double spe = 0, speErr = 0, cRes = 0, cResErr = 0;
            double tts = 0, ttsRelErr = 0;

            tr->SetBranchAddress("relativeQE", &qeCorr);
            tr->SetBranchAddress("relativeQE_err", &qeCorrErr);
            tr->SetBranchAddress("spe_mean", &spe);
            tr->SetBranchAddress("spe_mean_error", &speErr);
            bool hasRaw = tr->GetBranch("relativeQE_raw") != nullptr;
            if (hasRaw) {
                tr->SetBranchAddress("relativeQE_raw", &qeRaw);
                tr->SetBranchAddress("relativeQE_raw_err", &qeRawErr);
            }
            if (tr->GetBranch("charge_resolution")) {
                tr->SetBranchAddress("charge_resolution", &cRes);
                tr->SetBranchAddress("charge_resolution_err", &cResErr);
            }
            if (tr->GetBranch("rms_exG")) {
                tr->SetBranchAddress("rms_exG", &tts);
                if (tr->GetBranch("rms_exG_err")) tr->SetBranchAddress("rms_exG_err", &ttsRelErr);
            }
            tr->GetEntry(0);

            if (!hasRaw) { qeRaw = qeCorr; qeRawErr = qeCorrErr; }

            Point p;
            p.t = tMin;
            p.epoch = mtime;
            p.qeRaw = qeRaw;   p.qeRawErr = qeRawErr;
            p.qeCorr = qeCorr; p.qeCorrErr = qeCorrErr;
            p.spe = spe;       p.speErr = speErr;
            p.chargeRes = cRes; p.chargeResErr = cResErr;

            p.tts = tts * 2.0;
            p.ttsErr = (tts > 0 && ttsRelErr >= 0) ? ttsRelErr * tts * 2.0 : 0.0;
            TParameter<double>* darkPar = (TParameter<double>*)f->Get(Form("NoiseCountRate_ch%d", c));
            p.darkRate = darkPar ? darkPar->GetVal() : 0.0;
            p.darkRateErr = 0.0;
            p.ok = true;
            ch[c].pts.push_back(p);
            any = true;
        }
        if (any) nFiles++;
        f->Close();
    }

    std::cout << "[INFO] " << nFiles << " runs loaded from " << resultDir << std::endl;
    if (nFiles == 0) { std::cout << "[ERROR] nothing to plot." << std::endl; return; }

    struct MetricDef {
        TString key, title;
        double Point::*val;
        double Point::*err;
    };
    std::vector<MetricDef> metrics = {
        {"CorrRawQE",        "Corrected Raw QE [%]",             &Point::qeCorr,    &Point::qeCorrErr},
        {"RawQE",            "Raw QE [%]",                       &Point::qeRaw,     &Point::qeRawErr},
        {"RawSPE",           "Gain / SPE Charge [pC]",           &Point::spe,       &Point::speErr},
        {"ChargeResolution", "Charge Resolution #sigma/#mu [%]", &Point::chargeRes, &Point::chargeResErr},
        {"TTS",              "TTS [ns]",                         &Point::tts,       &Point::ttsErr},
        {"DarkRate",         "Dark Rate [Hz]",                   &Point::darkRate,  &Point::darkRateErr},
    };

    gSystem->mkdir("./Data/image/Stability", kTRUE);
    TString pdf = Form("./Data/image/Stability/Stability_Report_%s_%d_%d.pdf",
                       tag.Data(), run_start, run_end);

    TCanvas* c = new TCanvas("cStab", "Stability", 950, 1300);
    c->Print(pdf + "[");

    for (size_t m = 0; m < metrics.size(); ++m) {
        const MetricDef& md = metrics[m];
        c->Clear();
        c->Divide(1, 3, 0.002, 0.002);

        std::cout << "\n===== " << md.title << " =====" << std::endl;

        bool anyPad = false;
        for (int k = 0; k < 3; ++k) {
            std::vector<double> xs, ys, eys, valid;
            for (const Point& p : ch[k].pts) {
                double v = p.*(md.val);
    
                if (v <= 0) continue;
                xs.push_back(p.epoch);
                ys.push_back(v);
                eys.push_back(p.*(md.err));
                valid.push_back(v);
            }

            c->cd(k + 1);
            gPad->SetGrid();
            gPad->SetTopMargin(0.14);
            gPad->SetBottomMargin(0.20);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.04);
            if (xs.empty()) continue;
            anyPad = true;

            double mean = 0, rms = 0, statErr = 0, dummy = 0;
            Stats(valid, mean, rms);
            Stats(eys, statErr, dummy);   // statErr = mean per-point error
            double spread = (*std::max_element(valid.begin(), valid.end()))
                          - (*std::min_element(valid.begin(), valid.end()));

            TGraphErrors* g = new TGraphErrors(xs.size(), &xs[0], &ys[0], nullptr, &eys[0]);
            g->SetMarkerStyle(ch[k].marker);
            g->SetMarkerColor(ch[k].color);
            g->SetLineColor(ch[k].color);
            g->SetMarkerSize(0.9);
            g->SetFillColorAlpha(ch[k].color, 0.25);
            g->SetTitle(Form(";Date [MM/DD HH:MM];%s", md.title.Data()));
            g->Draw("A3");       
            g->Draw("PZ same");  
            g->GetXaxis()->SetTimeDisplay(1);
            g->GetXaxis()->SetTimeFormat("%m/%d %H:%M");
            g->GetXaxis()->SetTitleSize(0.062);
            g->GetYaxis()->SetTitleSize(0.062);
            g->GetXaxis()->SetLabelSize(0.052);
            g->GetYaxis()->SetLabelSize(0.052);
            g->GetXaxis()->SetLabelOffset(0.022);
            g->GetXaxis()->SetTitleOffset(1.35);
            g->GetYaxis()->SetTitleOffset(1.12);
            g->GetXaxis()->SetNdivisions(508);   

            double xlo = xs.front(), xhi = xs.back();
            TLine* lMean = new TLine(xlo, mean, xhi, mean);
            lMean->SetLineColor(kBlack); lMean->SetLineStyle(1); lMean->SetLineWidth(2);
            lMean->Draw();
            TLine* lHi = new TLine(xlo, mean + rms, xhi, mean + rms);
            TLine* lLo = new TLine(xlo, mean - rms, xhi, mean - rms);
            for (TLine* l : {lHi, lLo}) {
                l->SetLineColor(kBlack); l->SetLineStyle(2); l->SetLineWidth(1);
                l->Draw();
            }

            TLatex lab; lab.SetNDC(); lab.SetTextFont(132); lab.SetTextSize(0.055);
            lab.SetTextColor(ch[k].color);
            lab.DrawLatex(0.20, 0.82, Form("%s  mean=%.4g  (RMS %.2f%%)", ch[k].label.Data(),
                                           mean, mean != 0 ? rms / mean * 100 : 0));

            printf("  %-6s N=%3zu  mean=%8.4f  RMS=%7.4f (%5.2f%%)  spread=%7.4f (%5.2f%%)"
                   "  stat.err=%7.4f  RMS/stat=%.1fx\n",
                   ch[k].label.Data(), valid.size(), mean, rms,
                   mean != 0 ? rms / mean * 100 : 0, spread,
                   mean != 0 ? spread / mean * 100 : 0, statErr,
                   statErr > 0 ? rms / statErr : 0);
        }
        if (!anyPad) continue;

        c->cd(0);
        TLatex hdr; hdr.SetNDC(); hdr.SetTextFont(132); hdr.SetTextSize(0.028);
        hdr.DrawLatex(0.09, 0.965,
                      Form("Stability: %s  |  %s runs %d-%d  (%d points)",
                           md.title.Data(), tag.Data(), run_start, run_end, nFiles));

        c->Print(pdf);
    }


    {
        c->Clear();
        c->Divide(1, 3, 0.002, 0.015);
        bool anyAcf = false;
        for (int k = 0; k < 3; ++k) {
            std::vector<double> t, y;
            for (const Point& p : ch[k].pts) {
                if (p.qeRaw <= 0) continue;
                t.push_back((p.epoch - t0) / 60.0);
                y.push_back(p.qeRaw);
            }
            if ((int)t.size() < 10) continue;

            std::vector<double> dts;
            for (size_t i = 1; i < t.size(); ++i) dts.push_back(t[i] - t[i-1]);
            std::sort(dts.begin(), dts.end());
            double dtMed = dts[dts.size()/2];

            if (dtMed <= 1e-6) continue;
            anyAcf = true;
            std::vector<double> grid;
            for (double g = 0; g < t.back(); g += dtMed) grid.push_back(g);
            int N = (int)grid.size();

            double windowMin = std::min(180.0, t.back() / 3.0);
            int windowSteps = std::max(3, (int)std::round(windowMin / dtMed));

            std::vector<double> yi = ResampleUniform(t, y, grid);
            std::vector<double> trend = MovingAvgTrend(yi, windowSteps);

            int col = ch[k].color;
            c->cd(k + 1);
            gPad->SetGrid(); gPad->SetLeftMargin(0.09); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.14);
            TGraph* gRaw = new TGraph(N, grid.data(), yi.data());
            gRaw->SetMarkerStyle(20); gRaw->SetMarkerSize(0.6); gRaw->SetMarkerColor(col);
            gRaw->SetLineColor(col); gRaw->SetLineWidth(1);
            gRaw->SetTitle(Form("%s: Raw QE + trend;Elapsed [min];QE [%%]", ch[k].label.Data()));
            gRaw->GetXaxis()->SetLabelSize(0.045); gRaw->GetXaxis()->SetTitleSize(0.05);
            gRaw->GetYaxis()->SetLabelSize(0.045); gRaw->GetYaxis()->SetTitleSize(0.05); gRaw->GetYaxis()->SetTitleOffset(0.8);
            gRaw->Draw("APL");
            TGraph* gTrend = new TGraph(N, grid.data(), trend.data());
            gTrend->SetLineColor(kRed+1); gTrend->SetLineWidth(3);
            gTrend->Draw("L same");
        }
        if (anyAcf) c->Print(pdf);
    }

    // ── Temp / Humidity page: Dark Box 1 & 2 
    std::vector<EnvPoint> env = ReadEnvLog(t0, tLast);
    std::cout << "\n[INFO] " << env.size() << " monitoring-log samples in this window." << std::endl;
    if (!env.empty()) {
        c->Clear();
        c->Divide(1, 2, 0.002, 0.002);


        std::vector<double> envEpoch(env.size());
        for (size_t i = 0; i < env.size(); ++i) envEpoch[i] = env[i].epoch;
        double envSpanMin = env.empty() ? 0 : (env.back().epoch - env.front().epoch) / 60.0;
        std::vector<double> envDts;
        for (size_t i = 1; i < env.size(); ++i) envDts.push_back((env[i].epoch - env[i-1].epoch) / 60.0);
        double envDtMed = 1.0;
        if (!envDts.empty()) {
            std::vector<double> sorted = envDts;
            std::sort(sorted.begin(), sorted.end());
            envDtMed = sorted[sorted.size()/2];
            if (envDtMed <= 1e-6) envDtMed = 1.0;
        }
        double envWindowMin = std::min(180.0, envSpanMin / 3.0);
        int envWindowSteps = std::max(3, (int)std::round(envWindowMin / envDtMed));

        auto drawEnv = [&](int pad, const char* yTitle, double EnvPoint::*v1, double EnvPoint::*v2) {
            c->cd(pad);
            gPad->SetGrid();
            gPad->SetTopMargin(0.14);
            gPad->SetBottomMargin(0.20);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.04);

            TGraph* g1 = new TGraph(), *g2 = new TGraph();
            std::vector<double> y1(env.size()), y2(env.size());
            for (size_t i = 0; i < env.size(); ++i) {
                g1->SetPoint(i, env[i].epoch, env[i].*v1);
                g2->SetPoint(i, env[i].epoch, env[i].*v2);
                y1[i] = env[i].*v1; y2[i] = env[i].*v2;
            }
            g1->SetLineColor(kOrange+7); g1->SetLineWidth(2);
            g2->SetLineColor(kGreen+2);  g2->SetLineWidth(2);

            TMultiGraph* mg = new TMultiGraph();
            mg->Add(g1, "L"); mg->Add(g2, "L");
            mg->SetTitle(Form(";Date [MM/DD HH:MM];%s", yTitle));
            mg->Draw("A");
            mg->GetXaxis()->SetTimeDisplay(1);
            mg->GetXaxis()->SetTimeFormat("%m/%d %H:%M");
            mg->GetXaxis()->SetTitleSize(0.062);
            mg->GetYaxis()->SetTitleSize(0.062);
            mg->GetXaxis()->SetLabelSize(0.052);
            mg->GetYaxis()->SetLabelSize(0.052);
            mg->GetXaxis()->SetLabelOffset(0.022);
            mg->GetXaxis()->SetTitleOffset(1.35);
            mg->GetYaxis()->SetTitleOffset(1.12);
            mg->GetXaxis()->SetNdivisions(508);

            TGraph* t1 = nullptr, *t2 = nullptr;
            if ((int)env.size() >= 3) {
                std::vector<double> trend1 = MovingAvgTrend(y1, envWindowSteps);
                std::vector<double> trend2 = MovingAvgTrend(y2, envWindowSteps);
                t1 = new TGraph((int)env.size(), envEpoch.data(), trend1.data());
                t2 = new TGraph((int)env.size(), envEpoch.data(), trend2.data());
                t1->SetLineColor(kRed+1);  t1->SetLineWidth(3);
                t2->SetLineColor(kBlue+2); t2->SetLineWidth(3);
                t1->Draw("L same");
                t2->Draw("L same");
            }

            TLegend* leg = new TLegend(0.66, 0.72, 0.96, 0.90);
            leg->SetBorderSize(1); leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
            leg->SetTextFont(132); leg->SetTextSize(0.036);
            leg->AddEntry(g1, "Dark Box 1", "l");
            leg->AddEntry(g2, "Dark Box 2", "l");
            if (t1) leg->AddEntry(t1, "Dark Box 1 trend", "l");
            if (t2) leg->AddEntry(t2, "Dark Box 2 trend", "l");
            leg->Draw();
        };
        drawEnv(1, "Temperature [#circC]", &EnvPoint::db1t, &EnvPoint::db2t);
        drawEnv(2, "Humidity [%]",         &EnvPoint::db1h, &EnvPoint::db2h);

        c->cd(0);
        TLatex hdr; hdr.SetNDC(); hdr.SetTextFont(132); hdr.SetTextSize(0.028);
        hdr.DrawLatex(0.09, 0.965,
                      Form("Stability: Temp / Humidity  |  %s runs %d-%d",
                           tag.Data(), run_start, run_end));
        c->Print(pdf);
    }


    auto setG = [](TGraph* gr, int color) {
        gr->SetLineColor(color); gr->SetMarkerColor(color);
        gr->SetMarkerStyle(20); gr->SetMarkerSize(0.5); gr->SetLineWidth(2);
    };
    std::vector<TGraph*> gI;
    DataConditionLoadMonCSV(kBfieldCSV, {"I1", "I2", "I3", "I4"}, t0, tLast, gI);
    TGraph* gMag[4]; DataConditionLoadBfieldMag(kBfieldCSV, t0, tLast, gMag);
    std::vector<TGraph*> gBT;
    DataConditionLoadMonCSV(kBfieldCSV, {"T1", "T2", "T3", "T4", "T5", "T6"}, t0, tLast, gBT);
    int nBfieldPts = 0;
    for (auto* gr : gI) nBfieldPts += gr->GetN();
    for (int k = 0; k < 4; ++k) nBfieldPts += gMag[k]->GetN();
    for (auto* gr : gBT) nBfieldPts += gr->GetN();
    std::cout << "\n[INFO] " << nBfieldPts << " B-field log points in this window." << std::endl;

    if (nBfieldPts > 0) {
        c->Clear();
        c->Divide(1, 3, 0.001, 0.001);
        int colorPalette[6] = {kRed+1, kBlue+1, kGreen+2, kOrange+1, kMagenta+1, kCyan+2};

        c->cd(1); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.16); gPad->SetBottomMargin(0.13);
        TMultiGraph* mgI = new TMultiGraph();
        std::vector<TString> lblI;
        for (int k = 0; k < 4; ++k) {
            setG(gI[k], colorPalette[k]);
            if (gI[k]->GetN()) {
                double m, s; DataConditionMeanStd(gI[k], m, s);
                mgI->Add(gI[k], "L");
                lblI.push_back(Form("Coil #%d (%.2f#pm%.2f A)", k + 1, m, s));
            }
        }
        mgI->SetTitle(Form("Stability: B-field Log  |  %s runs %d-%d;;Coil Current [A]",
                           tag.Data(), run_start, run_end));
        double legLoI = DataConditionReserveLegendHeadroom(mgI);
        mgI->Draw("A"); DataConditionStyleTimeAxis(mgI);
        TLegend* lI = new TLegend(0.60, legLoI + 0.02, 0.97, 1.0 - gPad->GetTopMargin() - 0.02);
        lI->SetNColumns(2); lI->SetTextSize(0.038); lI->SetBorderSize(0); lI->SetFillStyle(0);
        for (size_t k = 0; k < lblI.size(); ++k) lI->AddEntry(mgI->GetListOfGraphs()->At(k), lblI[k], "l");
        lI->Draw();

        c->cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.06); gPad->SetBottomMargin(0.13);
        TMultiGraph* mgM = new TMultiGraph();
        std::vector<TString> lblM;
        for (int k = 0; k < 4; ++k) {
            setG(gMag[k], colorPalette[k]);
            if (gMag[k]->GetN()) {
                double m, s; DataConditionMeanStd(gMag[k], m, s);
                mgM->Add(gMag[k], "L");
                lblM.push_back(Form("Mag #%d (%.1f#pm%.1f mG)", k + 1, m, s));
            }
        }
        mgM->SetTitle(";;|B| [mG]");
        double legLoM = DataConditionReserveLegendHeadroom(mgM);
        mgM->Draw("A"); DataConditionStyleTimeAxis(mgM);
        TLegend* lM = new TLegend(0.60, legLoM + 0.02, 0.97, 1.0 - gPad->GetTopMargin() - 0.02);
        lM->SetNColumns(2); lM->SetTextSize(0.038); lM->SetBorderSize(0); lM->SetFillStyle(0);
        for (size_t k = 0; k < lblM.size(); ++k) lM->AddEntry(mgM->GetListOfGraphs()->At(k), lblM[k], "l");
        lM->Draw();

        c->cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.06); gPad->SetBottomMargin(0.18);
        TMultiGraph* mgBT = new TMultiGraph();
        std::vector<TString> lblBT;
        for (int k = 0; k < 6; ++k) {
            setG(gBT[k], colorPalette[k]);
            if (gBT[k]->GetN()) {
                double m, s; DataConditionMeanStd(gBT[k], m, s);
                mgBT->Add(gBT[k], "L");
                lblBT.push_back(Form("T%d (%.1f#pm%.1f #circC)", k + 1, m, s));
            }
        }
        mgBT->SetTitle(";Time [MM-DD HH:MM];Coil Temp [#circC]");
        double legLoBT = DataConditionReserveLegendHeadroom(mgBT);
        mgBT->Draw("A"); DataConditionStyleTimeAxis(mgBT);
        TLegend* lBT = new TLegend(0.60, legLoBT + 0.02, 0.97, 1.0 - gPad->GetTopMargin() - 0.02);
        lBT->SetNColumns(2); lBT->SetTextSize(0.036); lBT->SetBorderSize(0); lBT->SetFillStyle(0);
        for (size_t k = 0; k < lblBT.size(); ++k) lBT->AddEntry(mgBT->GetListOfGraphs()->At(k), lblBT[k], "l");
        lBT->Draw();

        c->Print(pdf);
    }

    if (!env.empty() || nBfieldPts > 0) {
        c->Clear();
        c->Divide(1, 4, 0.001, 0.001);

            auto inPlotLegend = [](double frameTop, int nCols) {
            //TLegend* l = new TLegend(0.46, frameTop - 0.07, 0.60, frameTop - 0.005);
            TLegend* l = new TLegend(0.48, 0.80, 0.78, 0.90);
            l->SetNColumns(nCols); l->SetTextSize(0.065);
            l->SetBorderSize(0); l->SetFillStyle(0);
            return l;
        };

        c->cd(1); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.06);
        TMultiGraph* mgQE = new TMultiGraph();
        for (int k = 0; k < 3; ++k) {
            TGraph* g = new TGraph();
            for (const Point& p : ch[k].pts) {
                if (!p.ok || p.qeRaw <= 0) continue;
                g->SetPoint(g->GetN(), p.epoch, p.qeRaw);
            }
            if (!g->GetN()) continue;
            setG(g, ch[k].color); g->SetMarkerStyle(ch[k].marker); g->SetMarkerSize(0.6);
            mgQE->Add(g, "LP");
        }
        mgQE->SetTitle(";;Raw QE [%]");
        double frameTopQE = DataConditionReserveLegendHeadroom(mgQE, 0.45, 0.05);
        mgQE->Draw("A"); DataConditionStyleTimeAxisHM(mgQE);
        mgQE->GetXaxis()->SetRangeUser(t0, tLast);
        TLegend* lQE = inPlotLegend(frameTopQE, 3);
        for (int k = 0; k < 3 && mgQE->GetListOfGraphs() && k < mgQE->GetListOfGraphs()->GetSize(); ++k)
            lQE->AddEntry(mgQE->GetListOfGraphs()->At(k), ch[k].label.Data(), "lp");
        lQE->Draw();

        c->cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.06);
        TGraph* gT1 = new TGraph(), *gT2 = new TGraph();
        for (const EnvPoint& e : env) { gT1->SetPoint(gT1->GetN(), e.epoch, e.db1t); gT2->SetPoint(gT2->GetN(), e.epoch, e.db2t); }
        setG(gT1, kOrange+7); setG(gT2, kGreen+2);
        TMultiGraph* mgT = new TMultiGraph();
        mgT->Add(gT1, "PL"); mgT->Add(gT2, "PL");
        mgT->SetTitle(";;Temperature [#circC]");
        double frameTopT = DataConditionReserveLegendHeadroom(mgT, 0.45, 0.05);
        mgT->Draw("A"); DataConditionStyleTimeAxisHM(mgT);
        mgT->GetXaxis()->SetRangeUser(t0, tLast);
        TLegend* lT = inPlotLegend(frameTopT, 2);
        lT->AddEntry(gT1, "Dark Box 1", "l"); lT->AddEntry(gT2, "Dark Box 2", "l");
        lT->Draw();

        c->cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.06);
        TGraph* gH1 = new TGraph(), *gH2 = new TGraph();
        for (const EnvPoint& e : env) { gH1->SetPoint(gH1->GetN(), e.epoch, e.db1h); gH2->SetPoint(gH2->GetN(), e.epoch, e.db2h); }
        setG(gH1, kOrange+7); setG(gH2, kGreen+2);
        TMultiGraph* mgH = new TMultiGraph();
        mgH->Add(gH1, "PL"); mgH->Add(gH2, "PL");
        mgH->SetTitle(";;Humidity [%]");
        double frameTopH = DataConditionReserveLegendHeadroom(mgH, 0.45, 0.05);
        mgH->Draw("A"); DataConditionStyleTimeAxisHM(mgH);
        mgH->GetXaxis()->SetRangeUser(t0, tLast);
        TLegend* lH = inPlotLegend(frameTopH, 2);
        lH->AddEntry(gH1, "Dark Box 1", "l"); lH->AddEntry(gH2, "Dark Box 2", "l");
        lH->Draw();

        c->cd(4); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.20);
        TMultiGraph* mgB = new TMultiGraph();
        std::vector<TString> lblB;
        for (int k = 0; k < 2; ++k) {
            if (!gMag[k]->GetN()) continue;
            setG(gMag[k], k == 0 ? kMagenta+1 : kCyan+2);
            mgB->Add(gMag[k], "L");
            lblB.push_back(Form("|B| M%d", k + 1));
        }
        mgB->SetTitle(";Time [HH:MM];|B| [mG]");
        double frameTopB = DataConditionReserveLegendHeadroom(mgB, 0.45, 0.05);
        mgB->Draw("A"); DataConditionStyleTimeAxisHM(mgB);
        mgB->GetXaxis()->SetRangeUser(t0, tLast);
        TLegend* lB = inPlotLegend(frameTopB, 2);
        for (size_t k = 0; k < lblB.size(); ++k) lB->AddEntry(mgB->GetListOfGraphs()->At(k), lblB[k], "l");
        lB->Draw();

        c->Print(pdf);
    }


    if (!env.empty() && !ch[0].pts.empty()) {
        TCanvas* cZoom = new TCanvas("cStabZoom", "Stability Zoom", 1600, 1500);
        cZoom->Divide(1, 3, 0.001, 0.001);

        double tZoomEnd = std::min(tLast, t0 + 3.0 * 3600.0);


        TGraph* gTAll = new TGraph();
        for (const EnvPoint& e : env) {
            if (e.epoch > tZoomEnd) continue;
            gTAll->SetPoint(gTAll->GetN(), e.epoch, e.db1t);
        }
        double tMin = TMath::MinElement(gTAll->GetN(), gTAll->GetY());
        double tMax = TMath::MaxElement(gTAll->GetN(), gTAll->GetY());
        double tSpan = (tMax > tMin) ? (tMax - tMin) : 0.1;
        double tLo = tMin - tSpan * 0.15, tHi = tMax + tSpan * 0.15;

        for (int k = 0; k < 3; ++k) {
            cZoom->cd(k + 1);
            gPad->SetGrid(); gPad->SetLeftMargin(0.08); gPad->SetRightMargin(0.09);
            gPad->SetTopMargin(0.14); gPad->SetBottomMargin(0.14);

            TGraph* gQE = new TGraph();
            for (const Point& p : ch[k].pts) {
                if (!p.ok || p.qeRaw <= 0 || p.epoch > tZoomEnd) continue;
                gQE->SetPoint(gQE->GetN(), p.epoch, p.qeRaw);
            }
            if (!gQE->GetN()) continue;
            setG(gQE, kBlack); gQE->SetMarkerStyle(20); gQE->SetMarkerSize(0.9); gQE->SetLineWidth(2);

            double qeMin = TMath::MinElement(gQE->GetN(), gQE->GetY());
            double qeMax = TMath::MaxElement(gQE->GetN(), gQE->GetY());
            double qeSpan = (qeMax > qeMin) ? (qeMax - qeMin) : 1.0;
            double qeLo = qeMin - qeSpan * 0.15, qeHi = qeMax + qeSpan * 0.15;

            TH1F* frame = gPad->DrawFrame(t0, qeLo, tZoomEnd, qeHi);
            frame->SetTitle(Form(";Time [HH:MM];%s Raw QE [%%]", ch[k].label.Data()));
            frame->GetXaxis()->SetTimeDisplay(1);
            frame->GetXaxis()->SetTimeFormat("%H:%M%F1970-01-01 00:00:00");
            frame->GetXaxis()->SetNdivisions(510);
            frame->GetXaxis()->SetLabelFont(132); frame->GetXaxis()->SetTitleFont(132);
            frame->GetYaxis()->SetLabelFont(132); frame->GetYaxis()->SetTitleFont(132);
            frame->GetXaxis()->SetLabelSize(0.045); frame->GetXaxis()->SetTitleSize(0.05);
            frame->GetYaxis()->SetLabelSize(0.045); frame->GetYaxis()->SetTitleSize(0.05);
            frame->GetYaxis()->SetLabelColor(kBlack); frame->GetYaxis()->SetTitleColor(kBlack);
            gQE->Draw("LP same");

            std::vector<double> tScaled(gTAll->GetN());
            for (int i = 0; i < gTAll->GetN(); ++i)
                tScaled[i] = qeLo + (gTAll->GetY()[i] - tLo) / (tHi - tLo) * (qeHi - qeLo);
            TGraph* gTScaled = new TGraph(gTAll->GetN(), gTAll->GetX(), tScaled.data());
            setG(gTScaled, kRed); gTScaled->SetMarkerStyle(21); gTScaled->SetMarkerSize(0.7); gTScaled->SetLineWidth(2);
            gTScaled->Draw("LP same");

            TGaxis* axT = new TGaxis(tZoomEnd, qeLo, tZoomEnd, qeHi, tLo, tHi, 508, "+L");
            axT->SetLineColor(kRed); axT->SetLabelColor(kRed); axT->SetTitleColor(kRed);
            axT->SetLabelFont(132); axT->SetTitleFont(132);
            axT->SetLabelSize(0.045); axT->SetTitleSize(0.05); axT->SetTitleOffset(1.1);
            axT->SetTitle("Dark Box 1 Temp [#circC]");
            axT->Draw();

            TLegend* lZoom = new TLegend(0.10, 0.86, 0.60, 0.935);
            lZoom->SetNColumns(2);
            lZoom->SetBorderSize(0); lZoom->SetFillStyle(0);
            lZoom->SetTextFont(132); lZoom->SetTextSize(0.045);
            lZoom->AddEntry(gQE, Form("%s Raw QE", ch[k].label.Data()), "lp");
            lZoom->AddEntry(gTScaled, "Dark Box 1 Temp", "lp");
            lZoom->Draw();
        }

        cZoom->cd(0);
        TLatex hdrZoom; hdrZoom.SetNDC(); hdrZoom.SetTextFont(132); hdrZoom.SetTextSize(0.022);
        hdrZoom.DrawLatex(0.08, 0.985, "QE vs Temp");

        cZoom->Print(pdf);
        delete cZoom;
    }

    c->Print(pdf + "]");
    std::cout << "\n[INFO] Saved: " << pdf << std::endl;
}
