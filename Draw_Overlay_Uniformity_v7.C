#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TKey.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TBox.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TDatime.h>
#include <TAxis.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TNamed.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include <array>
#include <TArrow.h>
#include <TMarker.h>

// ==========================================
const bool drawRaw = true;   // "RawQE"/"RelQE" = uncorrected (no PHC/TimingCut) QE, alongside Corrected/Poisson
const bool drawCorrected = false;   // "Count Corr" = PHC + TimingCut (dark-subtracted) QE, shown alongside Poisson
const bool drawPoissonRaw = false;   // also gates MonNorm_RawPoissonQE below
const bool drawPoisson = false;   // also gates MonNorm_PoissonQE below
const bool drawMonNorm = true;   // Monitor(CH0)-normalized QE/Gain (Test_absolute / Monitor_absolute)
const bool drawRelErr = false;   // extra pages: Error/Value x100% for the metrics below, instead of the value itself
const bool drawCenterNorm = false;    // extra pages: re-normalize each tag to its own Center=1 for cross-tag shape overlay

const bool savePNG = false;

const TString monitorSerial = "EM2740";
// ==========================================

void ApplyGlobalStyle() {/*{{{*/
    gStyle->SetTextFont(132);
    gStyle->SetLabelFont(132, "xyz");
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetTitleFont(22, ""); 
    gStyle->SetLegendFont(132);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetGridWidth(1);
    gStyle->SetGridColor(kGray);
    gStyle->SetOptStat(0);
}/*}}}*/

struct MetricDef {/*{{{*/
    TString name;
    TString yTitle;
};/*}}}*/

TGraphErrors* ToRelErrGraph(TGraphErrors* g) {/*{{{*/
    if (!g) return nullptr;
    TGraphErrors* out = new TGraphErrors(g->GetN());
    for (int i = 0; i < g->GetN(); ++i) {
        double x, y; g->GetPoint(i, x, y);
        double ey = g->GetErrorY(i);
        double relPct = (y != 0) ? 100.0 * std::abs(ey / y) : 0.0;
        out->SetPoint(i, x, relPct);
        out->SetPointError(i, 0, 0);
    }
    return out;
}/*}}}*/

// Re-normalizes a graph so its own Center (mean value at |angle| < 20deg,
// same window Draw_QE_UniformityShape.C uses) equals 1 -- errors are scaled
// by the same factor, not zeroed, since the rescaled value still carries
// meaningful uncertainty
TGraphErrors* ToCenterNormGraph(TGraphErrors* g) {/*{{{*/
    if (!g) return nullptr;
    const double centerMax = 20.0;
    double sum = 0; int n = 0;
    for (int i = 0; i < g->GetN(); ++i) {
        double x, y; g->GetPoint(i, x, y);
        if (std::abs(x) < centerMax && y > 0) { sum += y; n++; }
    }
    if (n == 0) return nullptr;
    double center = sum / n;

    TGraphErrors* out = new TGraphErrors(g->GetN());
    for (int i = 0; i < g->GetN(); ++i) {
        double x, y; g->GetPoint(i, x, y);
        double ey = g->GetErrorY(i);
        out->SetPoint(i, x, y / center);
        out->SetPointError(i, 0, ey / center);
    }
    return out;
}/*}}}*/

static const char* kRawDirs[] = {/*{{{*/
    "/home/precalkor/Data/RAW/Laser",
    "/home/precalkor/external_HDD_1_4T/Data_Backup/RAW/Laser",
};
static const char* kHVdb   = "/home/precalkor/Integrated_Control_SW/HV_Control_SW/monitoring_log.db";
static const char* kLaserDir = "/home/precalkor/ADC/ADC_test/LOG/LASER";
static const char* kBfieldCSV = "/home/precalkor/ADC/ADC_test/log_bfield_20260627.csv";/*}}}*/

static double NiceBinWidth(double range, int nEntries) {/*{{{*/
    if (range <= 0) return 1.0;
    int targetBins = std::min(12, std::max(5, (int)std::round(std::sqrt((double)std::max(nEntries, 1)) * 2)));
    double raw = range / targetBins;
    double mag = std::pow(10.0, std::floor(std::log10(raw)));
    double norm = raw / mag;   // in [1,10)
    double nice = (norm < 1.5) ? 1.0 : (norm < 2.25) ? 2.0 : (norm < 3.5) ? 2.5 : (norm < 7.5) ? 5.0 : 10.0;
    return nice * mag;
}/*}}}*/

static double DataConditionParseTime(const char* timeStr) {/*{{{*/
    int year, month, day, hour, minute, second;
    if (sscanf(timeStr, "%d-%d-%d%*c%d:%d:%d", &year, &month, &day, &hour, &minute, &second) != 6) return -1;
    return (double)TDatime(year, month, day, hour, minute, second).Convert();
}/*}}}*/

static TString DataConditionFmtTime(double unixsec, const char* fmt) {/*{{{*/
    time_t tt = (time_t)unixsec;
    struct tm* lt = localtime(&tt);
    char buf[64] = {0};
    if (lt) strftime(buf, sizeof(buf), fmt, lt);
    return TString(buf);
}/*}}}*/

static double DataConditionTakingLogTime(const TString& date, int run) {/*{{{*/
    TString p = Form("/home/precalkor/ADC/ADC_test/LOG/DAQ/TakingLog_%s.txt", date.Data());
    std::ifstream f(p.Data());
    if (!f.is_open()) return -1;
    TString needle = Form("precal_raw_kor_run_%s_%03d.root", date.Data(), run);
    TString line, lastTs;
    while (line.ReadLine(f)) {
        if (!line.Contains("Generated File") || !line.Contains(needle)) continue;
        TObjArray* tk = line.Tokenize("|");
        if (tk->GetEntries() > 0) lastTs = ((TObjString*)tk->At(0))->GetString().Strip(TString::kBoth);
        delete tk;
    }
    if (lastTs.IsNull()) return -1;
    return DataConditionParseTime(Form("%s:00", lastTs.Data()));
}/*}}}*/

static double DataConditionRawMTime(const TString& date, int run) {/*{{{*/
    for (const char* dir : kRawDirs) {
        TString p = Form("%s/precal_raw_kor_run_%s_%03d.root", dir, date.Data(), run);
        Long_t id, flags, mt; Long64_t sz;
        if (gSystem->GetPathInfo(p.Data(), &id, &sz, &flags, &mt) == 0) return (double)mt;
    }
    return DataConditionTakingLogTime(date, run);
}/*}}}*/

struct DataConditionMeta {/*{{{*/
    int wl = 0, laser = 0, hv1 = 0, hv2 = 0, hv3 = 0;
    TString sn1, sn2, sn3, note;
    TString shifter, expert;
    double pulse = -1, bias = -1;   // from laser CSV, median over the run window
    TString bfield = "N/A";
    double t0 = -1, t1 = -1;        // acquisition span (RAW mtimes)
    int runA = 0, runB = 0;
    TString date;
    TString tag;   // original "<date>_<runA>_<runB>" string, for tagLabels lookup
    bool ok = false;
};/*}}}*/

static DataConditionMeta DataConditionReadRunInfo(const TString& date, int run) {/*{{{*/
    DataConditionMeta m; m.date = date; m.runA = run;
    TString p = Form("./Data/FinalResult/precal_result_kor_run_%s_%03d.root", date.Data(), run);
    TFile* f = TFile::Open(p);
    if (!f || f->IsZombie()) { if (f) f->Close(); return m; }
    TTree* t = (TTree*)f->Get("RunInfo");
    if (!t) { f->Close(); return m; }
    char sn1[128] = {0}, sn2[128] = {0}, sn3[128] = {0}, note[256] = {0};
    char shifter[128] = {0}, expert[128] = {0};
    int wl = 0, la = 0, hv1 = 0, hv2 = 0, hv3 = 0;
    t->SetBranchAddress("Wavelength", &wl); t->SetBranchAddress("Laser_mA", &la);
    t->SetBranchAddress("HV1", &hv1); t->SetBranchAddress("HV2", &hv2); t->SetBranchAddress("HV3", &hv3);
    t->SetBranchAddress("SN1", sn1); t->SetBranchAddress("SN2", sn2); t->SetBranchAddress("SN3", sn3);
    t->SetBranchAddress("NOTE", note);
    if (t->GetBranch("Shifter")) t->SetBranchAddress("Shifter", shifter);
    if (t->GetBranch("Expert"))  t->SetBranchAddress("Expert", expert);
    t->GetEntry(0);
    m.wl = wl; m.laser = la; m.hv1 = hv1; m.hv2 = hv2; m.hv3 = hv3;
    m.sn1 = sn1; m.sn2 = sn2; m.sn3 = sn3; m.note = note; m.ok = true;
    m.shifter = (shifter[0] != '\0') ? shifter : "Unknown";
    m.expert  = (expert[0]  != '\0') ? expert  : "Unknown";
    f->Close();
    return m;
}/*}}}*/

static void DataConditionLaserPulseBias(int wl, const TString& date, double t0, double t1,/*{{{*/
        double& pulse, double& bias) {
    pulse = -1; bias = -1;
    TString file = Form("%s/laser_data_%dnm_%s.csv", kLaserDir, wl, date.Data());
    std::ifstream f(file.Data());
    if (!f.is_open()) return;
    std::string line; bool header = true;
    int iT = -1, iP = -1, iB = -1, iOn = -1, color = 0;
    std::vector<double> ps, bs;
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string cell; std::vector<std::string> c;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        if (header) {
            for (size_t i = 0; i < c.size(); ++i) {
                if (c[i] == "timestamp") iT = i;
                else if (c[i] == "pulse_ma") iP = i;
                else if (c[i] == "bias_ma") iB = i;
                else if (c[i] == "ld_on") iOn = i;
            }
            header = false; continue;
        }
        if (iT < 0 || iP < 0 || (int)c.size() <= iP) continue;
        double tt = DataConditionParseTime(c[iT].c_str());
        if (tt < t0 || tt > t1) continue;
        if (iOn >= 0 && (int)c.size() > iOn && atof(c[iOn].c_str()) < 0.5) continue;
        ps.push_back(atof(c[iP].c_str()));
        if (iB >= 0 && (int)c.size() > iB) bs.push_back(atof(c[iB].c_str()));
        (void)color;
    }
    auto median = [](std::vector<double>& v) -> double {
        if (v.empty()) return -1;
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };
    pulse = median(ps); bias = median(bs);
}/*}}}*/

static TString DataConditionBfieldStatus(double t0, double t1) {/*{{{*/
    std::ifstream f(kBfieldCSV); if (!f.is_open()) return "N/A";
    std::string line; bool header = true;
    int iT = -1, iI[4] = {-1, -1, -1, -1};
    double sumI = 0; int n = 0;
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string cell; std::vector<std::string> c;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        if (header) {
            for (size_t i = 0; i < c.size(); ++i) {
                if (c[i] == "timestamp") iT = i;
                else if (c[i] == "I1") iI[0] = i; else if (c[i] == "I2") iI[1] = i;
                else if (c[i] == "I3") iI[2] = i; else if (c[i] == "I4") iI[3] = i;
            }
            header = false; continue;
        }
        if (iT < 0 || (int)c.size() <= iT) continue;
        double tt = DataConditionParseTime(c[iT].c_str());
        if (tt < t0 || tt > t1) continue;
        double mag = 0;
        for (int k = 0; k < 4; ++k) if (iI[k] >= 0 && (int)c.size() > iI[k]) mag += fabs(atof(c[iI[k]].c_str()));
        sumI += mag; n++;
    }
    if (n == 0) return "N/A";               
    return (sumI / n > 0.4) ? "ON" : "OFF"; 
}/*}}}*/

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
            double magMilliGauss = std::sqrt(x*x + y*y + z*z) * 1000.0; // G -> mG
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

static void DataConditionMeanStdRange(TGraph* g, double t0, double t1, double& mean, double& stddev, int& n) {
    mean = 0; stddev = 0; n = 0;
    if (!g) return;
    for (int i = 0; i < g->GetN(); ++i) {
        double x = g->GetX()[i];
        if (x < t0 || x > t1) continue;
        mean += g->GetY()[i]; ++n;
    }
    if (n == 0) return;
    mean /= n;
    if (n < 2) return;
    for (int i = 0; i < g->GetN(); ++i) {
        double x = g->GetX()[i];
        if (x < t0 || x > t1) continue;
        double dev = g->GetY()[i] - mean; stddev += dev*dev;
    }
    stddev = std::sqrt(stddev / n);
}

static void DataConditionDrawTagBands(const std::vector<DataConditionMeta>& metas, double ymin, double ymax) {
    for (const auto& m : metas) {
        if (m.t0 <= 0 || m.t1 <= 0) continue;
        TBox* bx = new TBox(m.t0, ymin, m.t1, ymax);
        bx->SetFillColorAlpha(kBlack, 0.20);
        bx->SetLineWidth(0);
        bx->Draw();
    }
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
    double dataTopFracOfAxis = (ymax - (ymin - range * botFrac)) / totalSpan; // where the real data ends, within [0,1] of the axis
    return frameLo + dataTopFracOfAxis * (frameHi - frameLo);
}

static void DataConditionStyleTimeAxis(TMultiGraph* mg) {
    TAxis* ax = mg->GetXaxis();
    ax->SetTimeDisplay(1); ax->SetTimeFormat("%m-%d %H:%M%F1970-01-01 00:00:00");
    ax->SetNdivisions(510); ax->SetLabelSize(0.05); ax->SetTitleSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05); mg->GetYaxis()->SetTitleSize(0.058); mg->GetYaxis()->SetTitleOffset(0.85);
}

static void DrawDataConditionPages(const std::vector<TString>& tags,
        const std::map<TString, TString>& tagLabels,
        const TString& pdfPath, const TString& pngDir) {
    std::vector<DataConditionMeta> metas;
    double allT0 = 1e18, allT1 = -1e18;
    for (const auto& tg : tags) {
        TObjArray* tk = tg.Tokenize("_");
        if (tk->GetEntries() < 3) { delete tk; continue; }
        TString date = ((TObjString*)tk->At(0))->GetString();
        int runA = TString(((TObjString*)tk->At(1))->GetString()).Atoi();
        int runB = TString(((TObjString*)tk->At(2))->GetString()).Atoi();
        delete tk;
        DataConditionMeta m = DataConditionReadRunInfo(date, runA);
        m.runA = runA; m.runB = runB; m.tag = tg;
        m.t0 = DataConditionRawMTime(date, runA);
        m.t1 = DataConditionRawMTime(date, runB);
        if (m.t0 > 0) allT0 = std::min(allT0, m.t0);
        if (m.t1 > 0) allT1 = std::max(allT1, m.t1);
        if (m.ok && m.t0 > 0 && m.t1 > 0) {
            DataConditionLaserPulseBias(m.wl, date, m.t0, m.t1, m.pulse, m.bias);
            m.bfield = DataConditionBfieldStatus(m.t0, m.t1);
        }
        metas.push_back(m);
    }

    // ================== PAGE 1: run-metadata table ==================
    TCanvas* c1 = new TCanvas("c_datacond_meta", "Data Condition - Metadata", 1600, 2400);
    c1->cd();
    auto Tx = [&](double x, double y, const char* s, int font = 42, double sz = 0.022) {
        TLatex* t = new TLatex(x, y, s); t->SetNDC(); t->SetTextFont(font);
        t->SetTextSize(sz); t->SetTextAlign(12); t->Draw();
    };
    Tx(0.055, 0.955, "A. Data Condition", 62, 0.040);
    Tx(0.055, 0.912, "Run metadata (laser, HV, B-field, run range and acquisition time)", 42, 0.022);

    if (allT1 > 0) {
        Tx(0.055, 0.872, Form("Data period : %s  #rightarrow  %s",
                    DataConditionFmtTime(allT0, "%Y-%m-%d %H:%M").Data(),
                    DataConditionFmtTime(allT1, "%Y-%m-%d %H:%M").Data()), 62, 0.023);
    }
    if (!metas.empty() && metas[0].ok) {
        bool sameSN = true, samePeople = true;
        for (auto& m : metas) {
            if (!m.ok) continue;
            if (m.sn1 != metas[0].sn1 || m.sn2 != metas[0].sn2 || m.sn3 != metas[0].sn3) sameSN = false;
            if (m.shifter != metas[0].shifter || m.expert != metas[0].expert) samePeople = false;
        }
        if (sameSN) {
            Tx(0.055, 0.840, Form("PMTs :  Mon = %s   |   Rot#1 = %s   |   Rot#2 = %s",
                        metas[0].sn1.Data(), metas[0].sn2.Data(), metas[0].sn3.Data()), 42, 0.021);
        } else {
            Tx(0.055, 0.840, "PMTs :  #color[2]{differ across selected blocks -- see RunInfo per block}", 42, 0.021);
        }
        if (samePeople) {
            Tx(0.055, 0.808, Form("Shifter :  %s   |   Expert :  %s",
                        metas[0].shifter.Data(), metas[0].expert.Data()), 42, 0.021);
        } else {
            Tx(0.055, 0.808, "Shifter/Expert :  #color[2]{differ across selected blocks -- see RunInfo per block}", 42, 0.021);
        }
    }

    // table header
    std::vector<double> cx = {0.055, 0.170, 0.315, 0.520, 0.610, 0.700, 0.815};
    double y = 0.748;
    Tx(cx[0], y, "Wavelength", 62, 0.020);
    Tx(cx[1], y, "Laser [mA]", 62, 0.020);
    Tx(cx[2], y, "HV Mon/Rot1/Rot2 [V]", 62, 0.020);
    Tx(cx[3], y, "B-field", 62, 0.020);
    Tx(cx[4], y, "Runs", 62, 0.020);
    Tx(cx[5], y, "Acquisition span", 62, 0.020);
    y -= 0.010;
    { TLine* ln = new TLine(0.05, y, 0.95, y); ln->SetNDC(); ln->SetLineColor(kGray + 2); ln->Draw(); }
    y -= 0.040;

    for (auto& m : metas) {
        TString wlLab = Form("%d nm", m.wl);
        // laser column: pulse + bias if available, else the recorded set current
        TString laser;
        if (m.pulse > 0)
            laser = (m.bias >= 0) ? Form("%.0f + %.0f = %.0f", m.pulse, m.bias, m.pulse + m.bias)
                : Form("%.0f", m.pulse);
        else
            laser = Form("%d (set)", m.laser);
        TString hv = Form("%d / %d / %d", m.hv1, m.hv2, m.hv3);
        TString runs = Form("%03d-%03d", m.runA, m.runB);
        TString span = "n/a";
        if (m.t0 > 0 && m.t1 > 0)
            span = Form("%s #rightarrow %s", DataConditionFmtTime(m.t0, "%m-%d %H:%M").Data(),
                    DataConditionFmtTime(m.t1, "%H:%M").Data());
        Tx(cx[0], y, wlLab.Data(), 42, 0.020);
        Tx(cx[1], y, laser.Data(), 42, 0.020);
        Tx(cx[2], y, hv.Data(), 42, 0.020);
        Tx(cx[3], y, m.bfield.Data(), 42, 0.020);
        Tx(cx[4], y, runs.Data(), 42, 0.020);
        Tx(cx[5], y, span.Data(), 42, 0.020);
        y -= 0.038;
    }
    y -= 0.020;
    Tx(0.055, y, "Laser [mA] = pulse + bias (median over the run window, diode on).", 42, 0.016);
    y -= 0.028;
    Tx(0.055, y, "B-field : ON/OFF from coil currents over the window; N/A = no B-field log in that window.", 42, 0.016);
    c1->Print(pdfPath + "(", "pdf Title:A. Data Condition - Metadata");   // opens the PDF
    if (savePNG) c1->SaveAs(pngDir + "/DataCondition_Metadata.png");
    delete c1;

    // ================== PAGE 2: environmental monitoring ==================
    if (allT1 < 0) { allT0 = 0; allT1 = (double)TDatime().Convert(); }
    TString slice = "/tmp/dc_hv_slice.csv";
    gSystem->Exec(Form(
                "sqlite3 -header -csv %s "
                "\"SELECT timestamp,Ch1_V,Ch2_V,Dark_Box_1_T,Dark_Box_2_T,Dark_Box_3_T,"
                "Dark_Box_1_H,Dark_Box_2_H,Dark_Box_3_H FROM monitoring_data "
                "WHERE timestamp>='%s' AND timestamp<='%s' ORDER BY timestamp;\" > %s",
                kHVdb, DataConditionFmtTime(allT0, "%Y-%m-%dT%H:%M:%S").Data(),
                DataConditionFmtTime(allT1, "%Y-%m-%dT%H:%M:%S").Data(), slice.Data()));

    std::vector<TGraph*> g;
    DataConditionLoadMonCSV(slice, {"Ch1_V", "Ch2_V", "Dark_Box_3_T", "Dark_Box_1_T", "Dark_Box_2_T",
            "Dark_Box_3_H", "Dark_Box_1_H", "Dark_Box_2_H"}, allT0, allT1, g);
    TGraph *gV1 = g[0], *gV2 = g[1];
    TGraph *gTmon = g[2], *gT1 = g[3], *gT2 = g[4];
    TGraph *gHmon = g[5], *gH1 = g[6], *gH2 = g[7];

    auto setG = [](TGraph* gr, int color) { gr->SetLineColor(color); gr->SetMarkerColor(color);
        gr->SetLineWidth(2); gr->SetMarkerStyle(20); gr->SetMarkerSize(0.4); };

    TCanvas* c2 = new TCanvas("c_datacond_mon", "Data Condition - Monitoring", 1600, 2000);
    c2->Divide(1, 3, 0.001, 0.001);

    // Panel 1: HV
    c2->cd(1); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.16); gPad->SetBottomMargin(0.13);
    setG(gV1, kRed + 1); setG(gV2, kBlue + 1);
    TMultiGraph* mgHV = new TMultiGraph();
    if (gV1->GetN()) mgHV->Add(gV1, "L"); if (gV2->GetN()) mgHV->Add(gV2, "L");
    mgHV->SetTitle("A. Data Condition #minus Monitoring over data period;;Voltage [V]");
    mgHV->Draw("A"); DataConditionStyleTimeAxis(mgHV);
    DataConditionDrawTagBands(metas, mgHV->GetYaxis()->GetXmin(), mgHV->GetYaxis()->GetXmax());
    TLegend* lHV = new TLegend(0.60, 0.40, 0.97, 0.58); lHV->SetNColumns(2);
    lHV->SetTextSize(0.045); lHV->SetBorderSize(0); lHV->SetFillStyle(0);
    lHV->AddEntry(gV1, "Rot#1 EM5370", "l"); lHV->AddEntry(gV2, "Rot#2 EL9590", "l"); lHV->Draw();

    // Panel 2: Temperature
    c2->cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.06); gPad->SetBottomMargin(0.13);
    setG(gTmon, kGreen + 2); setG(gT1, kRed + 1); setG(gT2, kBlue + 1);
    TMultiGraph* mgT = new TMultiGraph();
    if (gTmon->GetN()) mgT->Add(gTmon, "L"); if (gT1->GetN()) mgT->Add(gT1, "L"); if (gT2->GetN()) mgT->Add(gT2, "L");
    mgT->SetTitle(";;Temperature [#circC]"); mgT->Draw("A"); DataConditionStyleTimeAxis(mgT);
    DataConditionDrawTagBands(metas, mgT->GetYaxis()->GetXmin(), mgT->GetYaxis()->GetXmax());
    TLegend* lT = new TLegend(0.60, 0.75, 0.97, 0.93); lT->SetNColumns(3);
    lT->SetTextSize(0.045); lT->SetBorderSize(0); lT->SetFillStyle(0);
    lT->AddEntry(gTmon, "Mon box", "l"); lT->AddEntry(gT1, "Rot#1 box", "l"); lT->AddEntry(gT2, "Rot#2 box", "l"); lT->Draw();

    // Panel 3: Humidity
    c2->cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.06); gPad->SetBottomMargin(0.18);
    setG(gHmon, kGreen + 2); setG(gH1, kRed + 1); setG(gH2, kBlue + 1);
    TMultiGraph* mgH = new TMultiGraph();
    if (gHmon->GetN()) mgH->Add(gHmon, "L"); if (gH1->GetN()) mgH->Add(gH1, "L"); if (gH2->GetN()) mgH->Add(gH2, "L");
    mgH->SetTitle(";Time [MM-DD HH:MM];Humidity [%]"); mgH->Draw("A"); DataConditionStyleTimeAxis(mgH);
    DataConditionDrawTagBands(metas, mgH->GetYaxis()->GetXmin(), mgH->GetYaxis()->GetXmax());
    TLegend* lH = new TLegend(0.60, 0.78, 0.97, 0.96); lH->SetNColumns(3);
    lH->SetTextSize(0.045); lH->SetBorderSize(0); lH->SetFillStyle(0);
    lH->AddEntry(gHmon, "Mon box", "l"); lH->AddEntry(gH1, "Rot#1 box", "l"); lH->AddEntry(gH2, "Rot#2 box", "l"); lH->Draw();

    c2->Print(pdfPath, "pdf Title:A. Data Condition - Monitoring");
    if (savePNG) c2->SaveAs(pngDir + "/DataCondition_Monitoring.png");
    std::cout << "[INFO] Data Condition pages: HV=" << gV1->GetN()
        << " T=" << gT1->GetN() << " H=" << gH1->GetN() << " rows sliced." << std::endl;

    // ================== PAGE 2b: per-dataset monitoring averages ==================
    {
        TCanvas* c2b = new TCanvas("c_datacond_mon_avg", "Data Condition - Monitoring Averages", 1600, 2400);
        c2b->cd();
        auto Tx2 = [&](double x, double y, const char* s, int font = 42, double sz = 0.020) {
            TLatex* t = new TLatex(x, y, s); t->SetNDC(); t->SetTextFont(font);
            t->SetTextSize(sz); t->SetTextAlign(12); t->Draw();
        };
        Tx2(0.055, 0.965, "A. Data Condition #minus Per-dataset Monitoring Averages", 62, 0.030);
        Tx2(0.055, 0.935, "Mean over each dataset's own acquisition window [t0,t1]", 42, 0.018);
        double cx2[5] = {0.05, 0.24, 0.36, 0.48, 0.72};
        const char* hdr2[5] = {"Dataset", "Rot#1 [V]", "Rot#2 [V]", "Mon/Rot1/Rot2 T [#circC]", "Mon/Rot1/Rot2 H [%]"};
        double y2 = 0.90;
        for (int c = 0; c < 5; ++c) Tx2(cx2[c], y2, hdr2[c], 62, 0.017);
        y2 -= 0.020;
        TLine* div2 = new TLine(0.05, y2, 0.95, y2); div2->SetNDC(); div2->Draw();
        y2 -= 0.024;
        for (const auto& m : metas) {
            if (m.t0 <= 0 || m.t1 <= 0) continue;
            TString lbl = tagLabels.count(m.tag) ? tagLabels.at(m.tag) : m.date;
            double mv1, sv1, mv2, sv2; int idxRef, idxTgt;
            DataConditionMeanStdRange(gV1, m.t0, m.t1, mv1, sv1, idxRef); DataConditionMeanStdRange(gV2, m.t0, m.t1, mv2, sv2, idxTgt);
            double mtm, stm, mt1, st1, mt2, st2; int ntm, nt1, nt2;
            DataConditionMeanStdRange(gTmon, m.t0, m.t1, mtm, stm, ntm);
            DataConditionMeanStdRange(gT1, m.t0, m.t1, mt1, st1, nt1);
            DataConditionMeanStdRange(gT2, m.t0, m.t1, mt2, st2, nt2);
            double mhm, shm, mh1, sh1, mh2, sh2; int nhm, nh1, nh2;
            DataConditionMeanStdRange(gHmon, m.t0, m.t1, mhm, shm, nhm);
            DataConditionMeanStdRange(gH1, m.t0, m.t1, mh1, sh1, nh1);
            DataConditionMeanStdRange(gH2, m.t0, m.t1, mh2, sh2, nh2);
            Tx2(cx2[0], y2, Form("%s (%s)", lbl.Data(), m.date.Data()), 42, 0.017);
            Tx2(cx2[1], y2, idxRef ? Form("%.1f#pm%.1f", mv1, sv1) : "N/A", 42, 0.017);
            Tx2(cx2[2], y2, idxTgt ? Form("%.1f#pm%.1f", mv2, sv2) : "N/A", 42, 0.017);
            Tx2(cx2[3], y2, (ntm || nt1 || nt2) ? Form("%.2f / %.2f / %.2f", mtm, mt1, mt2) : "N/A", 42, 0.017);
            Tx2(cx2[4], y2, (nhm || nh1 || nh2) ? Form("%.1f / %.1f / %.1f", mhm, mh1, mh2) : "N/A", 42, 0.017);
            y2 -= 0.022;
        }
        c2b->Print(pdfPath, "pdf Title:A. Data Condition - Monitoring Averages");
        if (savePNG) c2b->SaveAs(pngDir + "/DataCondition_MonitoringAverages.png");
        delete c2b;
    }
    delete c2;   

    // ================== PAGE 3: B-field coil log ==================
    double bfieldT1 = std::max(allT1, (double)TDatime().Convert());
    std::vector<TGraph*> gI;
    DataConditionLoadMonCSV(kBfieldCSV, {"I1", "I2", "I3", "I4"}, allT0, bfieldT1, gI);
    TGraph* gMag[4]; DataConditionLoadBfieldMag(kBfieldCSV, allT0, bfieldT1, gMag);
    std::vector<TGraph*> gBT;
    DataConditionLoadMonCSV(kBfieldCSV, {"T1", "T2", "T3", "T4", "T5", "T6"}, allT0, bfieldT1, gBT);
    int nBfieldPts = 0;
    for (auto* gr : gI) nBfieldPts += gr->GetN();
    for (int k = 0; k < 4; ++k) nBfieldPts += gMag[k]->GetN();
    for (auto* gr : gBT) nBfieldPts += gr->GetN();

    TCanvas* c3 = new TCanvas("c_datacond_bfield", "Data Condition - B-field Log", 1600, 2000);
    if (nBfieldPts > 0) {
        c3->Divide(1, 3, 0.001, 0.001);
        int colorPalette[6] = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 1, kMagenta + 1, kCyan + 2};

        c3->cd(1); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
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
        mgI->SetTitle("A. Data Condition #minus B-field Log over data period;;Coil Current [A]");
        double legLoI = DataConditionReserveLegendHeadroom(mgI);
        mgI->Draw("A"); DataConditionStyleTimeAxis(mgI);
        DataConditionDrawTagBands(metas, mgI->GetYaxis()->GetXmin(), mgI->GetYaxis()->GetXmax());
        TLegend* lI = new TLegend(0.60, legLoI + 0.02, 0.97, 1.0 - gPad->GetTopMargin() - 0.02);
        lI->SetNColumns(2); lI->SetTextSize(0.038); lI->SetBorderSize(0); lI->SetFillStyle(0);
        for (size_t k = 0; k < lblI.size(); ++k) lI->AddEntry(mgI->GetListOfGraphs()->At(k), lblI[k], "l");
        lI->Draw();

        c3->cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
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
        DataConditionDrawTagBands(metas, mgM->GetYaxis()->GetXmin(), mgM->GetYaxis()->GetXmax());
        TLegend* lM = new TLegend(0.60, legLoM + 0.02, 0.97, 1.0 - gPad->GetTopMargin() - 0.02);
        lM->SetNColumns(2); lM->SetTextSize(0.038); lM->SetBorderSize(0); lM->SetFillStyle(0);
        for (size_t k = 0; k < lblM.size(); ++k) lM->AddEntry(mgM->GetListOfGraphs()->At(k), lblM[k], "l");
        lM->Draw();

        c3->cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03);
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
        DataConditionDrawTagBands(metas, mgBT->GetYaxis()->GetXmin(), mgBT->GetYaxis()->GetXmax());
        TLegend* lBT = new TLegend(0.60, legLoBT + 0.02, 0.97, 1.0 - gPad->GetTopMargin() - 0.02);
        lBT->SetNColumns(2); lBT->SetTextSize(0.036); lBT->SetBorderSize(0); lBT->SetFillStyle(0);
        for (size_t k = 0; k < lblBT.size(); ++k) lBT->AddEntry(mgBT->GetListOfGraphs()->At(k), lblBT[k], "l");
        lBT->Draw();
    } else {
        c3->cd();
        TLatex* t = new TLatex(0.5, 0.5, "#bf{N/A}");
        t->SetNDC(); t->SetTextAlign(22); t->SetTextSize(0.08); t->SetTextColor(kGray + 2); t->Draw();
        TLatex* t2 = new TLatex(0.5, 0.42, "(B-field log has no coverage for this data period)");
        t2->SetNDC(); t2->SetTextAlign(22); t2->SetTextSize(0.025); t2->SetTextColor(kGray + 2); t2->Draw();
    }
    c3->Print(pdfPath, "pdf Title:A. Data Condition - B-field Log");
    if (savePNG) c3->SaveAs(pngDir + "/DataCondition_Bfield.png");
    std::cout << "[INFO] B-field log page: " << nBfieldPts << " points in window." << std::endl;

    // ================== PAGE 3b: per-dataset B-field averages ==================
    if (nBfieldPts > 0) {
        TCanvas* c3b = new TCanvas("c_datacond_bfield_avg", "Data Condition - B-field Averages", 1600, 2400);
        c3b->cd();
        auto Tx3 = [&](double x, double y, const char* s, int font = 42, double sz = 0.020) {
            TLatex* t = new TLatex(x, y, s); t->SetNDC(); t->SetTextFont(font);
            t->SetTextSize(sz); t->SetTextAlign(12); t->Draw();
        };
        Tx3(0.055, 0.965, "A. Data Condition #minus Per-dataset B-field Averages", 62, 0.030);
        Tx3(0.055, 0.935, "Mean over each dataset's own acquisition window [t0,t1]", 42, 0.018);
        double cx3[4] = {0.05, 0.26, 0.48, 0.68};
        const char* hdr3[4] = {"Dataset", "Coil #1-4 [A]", "|B| #1-3 [mG]", "Coil T1/T2/T3/T5/T6 [#circC]"};
        double y3 = 0.90;
        for (int c = 0; c < 4; ++c) Tx3(cx3[c], y3, hdr3[c], 62, 0.017);
        y3 -= 0.020;
        TLine* div3 = new TLine(0.05, y3, 0.95, y3); div3->SetNDC(); div3->Draw();
        y3 -= 0.024;
        for (const auto& m : metas) {
            if (m.t0 <= 0 || m.t1 <= 0) continue;
            TString lbl = tagLabels.count(m.tag) ? tagLabels.at(m.tag) : m.date;
            TString iStr = "N/A"; bool anyI = false;
            for (int k = 0; k < 4; ++k) {
                double mv, sv; int n; DataConditionMeanStdRange(gI[k], m.t0, m.t1, mv, sv, n);
                if (n) { iStr = anyI ? iStr + Form(" / %.2f", mv) : Form("%.2f", mv); anyI = true; }
            }
            TString mgStr = "N/A"; bool anyMg = false;
            for (int k = 0; k < 4; ++k) {
                double mv, sv; int n; DataConditionMeanStdRange(gMag[k], m.t0, m.t1, mv, sv, n);
                if (n) { mgStr = anyMg ? mgStr + Form(" / %.0f", mv) : Form("%.0f", mv); anyMg = true; }
            }
            TString tStr = "N/A"; bool anyT = false;
            for (int k = 0; k < 6; ++k) {
                double mv, sv; int n; DataConditionMeanStdRange(gBT[k], m.t0, m.t1, mv, sv, n);
                if (n) { tStr = anyT ? tStr + Form(" / %.1f", mv) : Form("%.1f", mv); anyT = true; }
            }
            Tx3(cx3[0], y3, Form("%s (%s)", lbl.Data(), m.date.Data()), 42, 0.017);
            Tx3(cx3[1], y3, iStr.Data(), 42, 0.017);
            Tx3(cx3[2], y3, mgStr.Data(), 42, 0.017);
            Tx3(cx3[3], y3, tStr.Data(), 42, 0.017);
            y3 -= 0.022;
        }
        c3b->Print(pdfPath, "pdf Title:A. Data Condition - B-field Averages");
        if (savePNG) c3b->SaveAs(pngDir + "/DataCondition_BfieldAverages.png");
        delete c3b;
    }

    // ================== PAGE 3c: combined per-dataset averages (Monitoring + B-field) ==================
    {
        TCanvas* c3c = new TCanvas("c_datacond_combined_avg", "Data Condition - Combined Averages", 2000, 1300);
        c3c->cd();
        auto Tx4 = [&](double x, double y, const char* s, int font = 42, double sz = 0.020) {
            TLatex* t = new TLatex(x, y, s); t->SetNDC(); t->SetTextFont(font);
            t->SetTextSize(sz); t->SetTextAlign(12); t->Draw();
        };
        Tx4(0.03, 0.965, "A. Data Condition #minus Combined Per-dataset Averages", 62, 0.026);
        Tx4(0.03, 0.938, "Mean over each dataset's own acquisition window [t0,t1] -- Monitoring + B-field side by side", 42, 0.016);
        double cx4[8] = {0.03, 0.145, 0.225, 0.305, 0.435, 0.565, 0.685, 0.795};
        const char* hdr4[8] = {"Dataset", "Rot#1 [V]", "Rot#2 [V]", "Mon/Rot1/Rot2 T [#circC]",
            "Mon/Rot1/Rot2 H [%]", "Coil #1-4 [A]", "|B| #1-3 [mG]", "Coil T1-T6 [#circC]"};
        double y4 = 0.90;
        for (int c = 0; c < 8; ++c) Tx4(cx4[c], y4, hdr4[c], 62, 0.014);
        y4 -= 0.020;
        TLine* div4 = new TLine(0.03, y4, 0.97, y4); div4->SetNDC(); div4->Draw();
        y4 -= 0.026;
        for (const auto& m : metas) {
            if (m.t0 <= 0 || m.t1 <= 0) continue;
            TString lbl = tagLabels.count(m.tag) ? tagLabels.at(m.tag) : m.date;

            double mv1, sv1, mv2, sv2; int idxRef, idxTgt;
            DataConditionMeanStdRange(gV1, m.t0, m.t1, mv1, sv1, idxRef);
            DataConditionMeanStdRange(gV2, m.t0, m.t1, mv2, sv2, idxTgt);
            double mtm, stm, mt1, st1, mt2, st2; int ntm, nt1, nt2;
            DataConditionMeanStdRange(gTmon, m.t0, m.t1, mtm, stm, ntm);
            DataConditionMeanStdRange(gT1, m.t0, m.t1, mt1, st1, nt1);
            DataConditionMeanStdRange(gT2, m.t0, m.t1, mt2, st2, nt2);
            double mhm, shm, mh1, sh1, mh2, sh2; int nhm, nh1, nh2;
            DataConditionMeanStdRange(gHmon, m.t0, m.t1, mhm, shm, nhm);
            DataConditionMeanStdRange(gH1, m.t0, m.t1, mh1, sh1, nh1);
            DataConditionMeanStdRange(gH2, m.t0, m.t1, mh2, sh2, nh2);

            TString iStr = "N/A"; bool anyI = false;
            if (nBfieldPts > 0) for (int k = 0; k < 4; ++k) {
                double mv, sv; int n; DataConditionMeanStdRange(gI[k], m.t0, m.t1, mv, sv, n);
                if (n) { iStr = anyI ? iStr + Form(" / %.2f", mv) : Form("%.2f", mv); anyI = true; }
            }
            TString mgStr = "N/A"; bool anyMg = false;
            if (nBfieldPts > 0) for (int k = 0; k < 4; ++k) {
                double mv, sv; int n; DataConditionMeanStdRange(gMag[k], m.t0, m.t1, mv, sv, n);
                if (n) { mgStr = anyMg ? mgStr + Form(" / %.0f", mv) : Form("%.0f", mv); anyMg = true; }
            }
            TString tbStr = "N/A"; bool anyTb = false;
            if (nBfieldPts > 0) for (int k = 0; k < 6; ++k) {
                double mv, sv; int n; DataConditionMeanStdRange(gBT[k], m.t0, m.t1, mv, sv, n);
                if (n) { tbStr = anyTb ? tbStr + Form(" / %.1f", mv) : Form("%.1f", mv); anyTb = true; }
            }

            Tx4(cx4[0], y4, Form("%s (%s)", lbl.Data(), m.date.Data()), 42, 0.014);
            Tx4(cx4[1], y4, idxRef ? Form("%.1f#pm%.1f", mv1, sv1) : "N/A", 42, 0.014);
            Tx4(cx4[2], y4, idxTgt ? Form("%.1f#pm%.1f", mv2, sv2) : "N/A", 42, 0.014);
            Tx4(cx4[3], y4, (ntm || nt1 || nt2) ? Form("%.2f / %.2f / %.2f", mtm, mt1, mt2) : "N/A", 42, 0.014);
            Tx4(cx4[4], y4, (nhm || nh1 || nh2) ? Form("%.1f / %.1f / %.1f", mhm, mh1, mh2) : "N/A", 42, 0.014);
            Tx4(cx4[5], y4, iStr.Data(), 42, 0.014);
            Tx4(cx4[6], y4, mgStr.Data(), 42, 0.014);
            Tx4(cx4[7], y4, tbStr.Data(), 42, 0.014);
            y4 -= 0.026;
        }
        c3c->Print(pdfPath, "pdf Title:A. Data Condition - Combined Averages");
        if (savePNG) c3c->SaveAs(pngDir + "/DataCondition_CombinedAverages.png");
        delete c3c;
    }

    delete c3;
}
// ===========================================================================
// ======================= MAIN FUNCTION /====================================//
static bool kOverlayRawStageAxis = false;
static const char* OverlayAngleAxisTitle() {
    return kOverlayRawStageAxis ? "Raw Stage angle [degree]" : "Position angle [degree]";
}

void Draw_Overlay_Uniformity_v7() {
    ApplyGlobalStyle();

    std::vector<TString> tags = {
        "20260331_100_145", // Reference
                            //"20260331_0_45",
        "20260401_0_45",
        "20260401_100_145",
        //  "20260407_0_45",
        // "20260426_0_45"
    };

    std::map<TString, TString> tagLabels;
    tagLabels["20260226"] = "Feb 26 (Ref)";
    tagLabels["20260319"] = "Mar 19";
    tagLabels["20260331_0_45"] = "Mar 31 (1)";
    tagLabels["20260331_100_145"] = "Mar 31 (Ref)";
    tagLabels["20260401_0_45"] = "Apr 1";
    tagLabels["20260401_100_145"] = "Apr 1";
    tagLabels["20260407_0_45"] = "Apr 7";
    tagLabels["20260426_0_45"] = "Apr 26";
    tagLabels["20260713_0_45"] = "405 nm";
    tagLabels["20260713_100_145"] = "375 nm";
    tagLabels["20260713_200_245"] = "450 nm";
    tagLabels["20260713_300_345"] = "473 nm";
    tagLabels["20260715_0_45"] = "405 nm";
    tagLabels["20260715_100_145"] = "375 nm";
    tagLabels["20260715_200_245"] = "450 nm";
    tagLabels["20260715_300_345"] = "473 nm";

    int refIdx = 0;
    {
        std::ifstream cfg("./Data/UNIFORMITY/overlay_tags.txt");
        if (cfg.is_open()) {
            std::vector<TString> fileTags;
            int fileRefIdx = -1;
            std::string line;
            while (std::getline(cfg, line)) {
                while (!line.empty() && (line.back()=='\r' || line.back()==' ' || line.back()=='\t'))
                    line.pop_back();
                size_t b = line.find_first_not_of(" \t");
                if (b == std::string::npos) continue;
                line = line.substr(b);
                if (line.empty() || line[0] == '#') continue;
                std::string tg = line, lab;
                size_t comma = line.find(',');
                if (comma != std::string::npos) { tg = line.substr(0, comma); lab = line.substr(comma + 1); }
                size_t comma2 = lab.find(',');
                if (comma2 != std::string::npos) {
                    std::string flag = lab.substr(comma2 + 1);
                    lab = lab.substr(0, comma2);
                    while (!flag.empty() && (flag.front()==' ' || flag.front()=='\t')) flag.erase(flag.begin());
                    while (!flag.empty() && (flag.back()==' '  || flag.back()=='\t'))  flag.pop_back();
                    if (flag == "ref" || flag == "REF" || flag == "Ref") fileRefIdx = (int)fileTags.size();
                }
                TString tgs(tg.c_str());
                fileTags.push_back(tgs);
                if (!lab.empty()) tagLabels[tgs] = TString(lab.c_str());
            }
            if (!fileTags.empty()) {
                tags = fileTags;
                if (fileRefIdx >= 0) refIdx = fileRefIdx;
                std::cout << "[INFO] Overlay tag list loaded from overlay_tags.txt ("
                    << tags.size() << " entries), reference = "
                    << (tagLabels.count(tags[refIdx]) ? tagLabels[tags[refIdx]] : tags[refIdx]) << std::endl;
            }
        }
    }
    // =========================================================================

    std::set<int> selCh;
    {
        std::ifstream chf("./Data/UNIFORMITY/overlay_channels.txt");
        if (chf.is_open()) {
            std::string line;
            while (std::getline(chf, line)) {
                for (char& c : line) if (c == ',') c = ' ';
                std::istringstream iss(line);
                int v;
                while (iss >> v) if (v >= 0 && v < 3) selCh.insert(v);
            }
        }
    }

    std::set<TString> uniqueSerials;

    for (auto& tg : tags) {
        TString fileName = Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tg.Data());
        TFile *f = TFile::Open(fileName);
        if (!f || f->IsZombie()) {
            std::cout << "[WARNING] Cannot open file: " << fileName << std::endl;
            continue;
        }
        TIter next(f->GetListOfKeys()); TKey *key;
        while ((key = (TKey*)next())) {
            TString name = key->GetName();
            if (!name.Contains("gr_", TString::kIgnoreCase) || name.Contains(monitorSerial, TString::kIgnoreCase)) continue;
            TObjArray *tokens = name.Tokenize("_");
            if (tokens->GetEntries() > 1) {
                TString serial = ((TObjString*)tokens->At(1))->GetString();
                serial.ToUpper(); 
                uniqueSerials.insert(serial);
            }
            delete tokens;
        }
        f->Close();
    }

    std::vector<TString> testSerials(uniqueSerials.begin(), uniqueSerials.end());
    if (testSerials.empty()) {
        std::cout << "[ERROR] No valid serials found. Check ROOT files." << std::endl;
        return;
    }


    {
        std::set<TString> axesSeen;
        for (auto& tg : tags) {
            TFile* f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tg.Data()));
            if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
            TNamed* ax = (TNamed*)f->Get("AngleAxis");
            axesSeen.insert(ax ? TString(ax->GetTitle()) : TString("hamamatsu"));
            f->Close();
        }
        if (axesSeen.size() > 1) {
            std::cout << "[ERROR] The selected tags do not share one angle convention "
                "(found:";
            for (auto& a : axesSeen) std::cout << " " << a;
            std::cout << "). Re-run Draw_Uniformity_Norm_v7 so every tag uses the same "
                "xaxis before overlaying." << std::endl;
            return;
        }
        kOverlayRawStageAxis = (!axesSeen.empty() && *axesSeen.begin() == "rawstage");
        std::cout << "[INFO] X axis: " << OverlayAngleAxisTitle() << std::endl;
    }

    if (!selCh.empty()) {
        std::set<TString> allowedSerials;
        for (auto& tg : tags) {
            TFile* f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tg.Data()));
            if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
            TNamed* cm = (TNamed*)f->Get("ChannelMap");
            if (cm) {
                TString map = cm->GetTitle();           // "0:SN0;1:SN1;2:SN2;"
                TObjArray* pairs = map.Tokenize(";");
                for (int i = 0; i < pairs->GetEntries(); ++i) {
                    TString pr = ((TObjString*)pairs->At(i))->GetString();
                    Ssiz_t colon = pr.First(':');
                    if (colon == kNPOS) continue;
                    int ch = TString(pr(0, colon)).Atoi();
                    TString sn = pr(colon + 1, pr.Length()); sn.ToUpper();
                    if (selCh.count(ch)) allowedSerials.insert(sn);
                }
                delete pairs;
            }
            f->Close();
            if (!allowedSerials.empty()) break;   // one valid map is enough
        }
        if (!allowedSerials.empty()) {
            std::vector<TString> filtered;
            for (auto& s : testSerials) if (allowedSerials.count(s)) filtered.push_back(s);
            if (!filtered.empty()) {
                testSerials = filtered;
                std::cout << "[INFO] Channel filter applied -> " << testSerials.size()
                    << " serial(s) kept." << std::endl;
            } else {
                std::cout << "[WARNING] Channel filter matched no serials; drawing all." << std::endl;
            }
        }
    }

    {
        std::map<TString, int> serialToChannel;
        for (auto& tg : tags) {
            TFile* f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tg.Data()));
            if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
            TNamed* cm = (TNamed*)f->Get("ChannelMap");
            if (cm) {
                TString map = cm->GetTitle();
                TObjArray* pairs = map.Tokenize(";");
                for (int i = 0; i < pairs->GetEntries(); ++i) {
                    TString pr = ((TObjString*)pairs->At(i))->GetString();
                    Ssiz_t colon = pr.First(':');
                    if (colon == kNPOS) continue;
                    int ch = TString(pr(0, colon)).Atoi();
                    TString sn = pr(colon + 1, pr.Length()); sn.ToUpper();
                    serialToChannel[sn] = ch;
                }
                delete pairs;
            }
            f->Close();
            if (!serialToChannel.empty()) break;
        }
        if (!serialToChannel.empty()) {
            std::stable_sort(testSerials.begin(), testSerials.end(),
                    [&](const TString& a, const TString& b) {
                    int ca = serialToChannel.count(a) ? serialToChannel[a] : 999;
                    int cb = serialToChannel.count(b) ? serialToChannel[b] : 999;
                    return ca < cb;
                    });
        }
    }


    bool includeMonitor = selCh.empty() || selCh.count(0);
    std::vector<TString> allSerials;
    if (includeMonitor) allSerials.push_back(monitorSerial);
    for (auto s : testSerials) allSerials.push_back(s);

    TString sList = ""; for(auto s : testSerials) sList += "_" + s;
    TString dList = ""; for(auto tg : tags) dList += "_" + tg;
    TString pdfPath = Form("./Data/image/Uniformity/Overlay_Report_v7%s%s.pdf", sList.Data(), dList.Data());

    int dateColors[] = {kBlack, kRed, kBlue, kMagenta+2, kGreen+2, kOrange+1,
        kCyan+2, kViolet+1, kAzure+2, kSpring+4, kPink+7};
    int markerStyles[] = {20, 21, 22, 34, 29, 23, 33, 24, 25, 26, 27};
    int nStyles = sizeof(markerStyles) / sizeof(int);
    int nDateColors = sizeof(dateColors) / sizeof(int);


    auto PaletteSlot  = [&](int i) { return (i == refIdx) ? 0 : (i == 0 ? refIdx : i); };
    auto ColorForTag  = [&](int i) { return dateColors[PaletteSlot(i) % nDateColors]; };
    auto MarkerForTag = [&](int i) { return markerStyles[PaletteSlot(i) % nStyles]; };
    // series[] on the deviation/pull pages skips the reference, so map its
    // running index back onto the real tag index before picking a colour.
    auto ColorForSeries = [&](int s) { return ColorForTag(s < refIdx ? s : s + 1); };

    std::vector<MetricDef> metricsToDraw;
    if (drawRaw) {
        metricsToDraw.push_back({"RawQE", "Raw QE [%]"});
        // "RelQE" is RawQE divided by its own average over |raw tilt| < 48 deg
        metricsToDraw.push_back({"RelQE", "Normalized QE [a.u.]"});
    }
    if (drawCorrected) {
        metricsToDraw.push_back({"CorrRawQE", "Corrected Raw QE [%]"});
        metricsToDraw.push_back({"CorrRelQE", "Normalized Corr QE [a.u.]"});
    }
    if (drawPoissonRaw) {
        metricsToDraw.push_back({"RawPoissonQE", "Poisson Raw QE [%]"});
    }
    if (drawPoisson) {
        metricsToDraw.push_back({"PoissonQE", "Normalized Poisson QE [a.u.]"});
    }
    if (drawMonNorm) {
        // Monitor(CH0)-normalized: Test_absolute / Monitor_absolute at the same
        metricsToDraw.push_back({"MonNorm_RawQE", "Mon-Norm Raw QE [a.u.]"});
        metricsToDraw.push_back({"MonNorm_CorrRawQE", "QE_{corr}^{Test} / QE_{corr}^{Mon}  [a.u.]"});
        if (drawCenterNorm) {
            metricsToDraw.push_back({"MonNorm_CorrRawQE_CenterNorm", "QE (Mon-Norm, Center=1) [a.u.]"});
        }
        if (drawPoissonRaw) metricsToDraw.push_back({"MonNorm_RawPoissonQE", "Mon-Norm Poisson Raw QE [a.u.]"});
        if (drawPoisson) metricsToDraw.push_back({"MonNorm_PoissonQE", "Mon-Norm Poisson QE [a.u.]"});
        metricsToDraw.push_back({"MonNorm_RelGain", "Mon-Norm Gain [a.u.]"});
    }

    metricsToDraw.push_back({"RawSPE", "Charge [pC]"});

    metricsToDraw.push_back({"ChargeResolution", "Charge Resolution #sigma/#mu [%]"});
    metricsToDraw.push_back({"RelGain", "Relative Gain [a.u.]"});
    metricsToDraw.push_back({"DarkCount", "Dark Count [entries]"});
    // Rate = Count / (this dataset's own livetime), already computed correctly
    metricsToDraw.push_back({"DarkRate", "Dark Rate [Hz]"});
    metricsToDraw.push_back({"TTS", "TTS [ns]"});
    metricsToDraw.push_back({"FWHM", "FWHM [ns]"});
    metricsToDraw.push_back({"Sigma", "Sigma (exGaus fit core) [ns]"});
    metricsToDraw.push_back({"PedStability", "Pedestal Mean #pm #sigma [mV]"});

    if (drawRelErr) {

        if (drawMonNorm) {
            metricsToDraw.push_back({"MonNorm_CorrRawQE_RelErr", "QE (Mon-Norm Corr) Rel.Err [%]"});
        }
        metricsToDraw.push_back({"TTS_RelErr", "TTS Rel.Err [%]"});
    }

    int totalCanvases = metricsToDraw.size() * 2;
    int canvasCount = 0;

    DrawDataConditionPages(tags, tagLabels, pdfPath, "./Data/image/Uniformity");

    // ======================= Page 3: Table of Contents =========================

    {
        TCanvas* ctoc = new TCanvas("c_toc", "Table of Contents", 1600, 2000);
        ctoc->cd();
        auto Tx = [&](double x, double y, const char* s, int font = 42, double sz = 0.021) {
            TLatex* t = new TLatex(x, y, s); t->SetNDC(); t->SetTextFont(font);
            t->SetTextSize(sz); t->SetTextAlign(12); t->Draw();
        };
        const int pageDataCond = 1;   
        const int pageTOC = 4;
        const int pageMetricsStart = 5;
        int nMetrics = (int)metricsToDraw.size();
        int pageRatioSummary = pageMetricsStart + 2 * nMetrics;

        Tx(0.055, 0.955, "Table of Contents", 62, 0.038);
        double y = 0.900;
        Tx(0.06, y, Form("A. Data Condition ................................................  p.%d - %d", pageDataCond, pageDataCond + 2), 42, 0.023);
        y -= 0.040;
        Tx(0.06, y, "B. Overlay Plots (per metric: X-axis then Y-axis)", 62, 0.023);
        y -= 0.038;
        for (int k = 0; k < nMetrics; ++k) {
            int pX = pageMetricsStart + 2 * k, pY = pX + 1;
            Tx(0.10, y, Form("%-24s  %s", metricsToDraw[k].name.Data(), metricsToDraw[k].yTitle.Data()),
                    42, 0.0185);
            Tx(0.72, y, Form("p.%d - %d", pX, pY), 42, 0.0185);
            y -= 0.0265;
        }
        y -= 0.014;
        Tx(0.06, y, Form("C. Ratio Summary  (mean/RMS vs %s reference) ..........  p.%d",
                    (tagLabels.count(tags[refIdx]) ? tagLabels.at(tags[refIdx]) : tags[refIdx]).Data(), pageRatioSummary), 62, 0.023);
        y -= 0.045;
        Tx(0.06, y, Form("Wavelength blocks compared: %d   |   PMT/SN mapping: see p.1", (int)tags.size()), 42, 0.018);

        ctoc->Print(pdfPath, "pdf Title:Table of Contents");
        delete ctoc;
    }

    for (const auto& metricDef : metricsToDraw) {
        TString metric = metricDef.name;

        bool isRelErr = metric.EndsWith("_RelErr");
        bool isCenterNorm = metric.EndsWith("_CenterNorm");
        TString baseMetric = isRelErr ? TString(metric(0, metric.Length() - 7)) :
            isCenterNorm ? TString(metric(0, metric.Length() - 11)) : metric;
        TString yTitle = metricDef.yTitle;
        TString metricLabel = yTitle;
        if (metricLabel.Contains(" [")) metricLabel = metricLabel(0, metricLabel.Index(" ["));
        for (TString axis : {"X", "Y"}) {
            canvasCount++;

            auto drawRatioHistPage = [&](const std::vector<std::vector<std::pair<TString, std::vector<double>>>>& series,
                    int numItems, int nCols, int nRows,
                    const std::vector<TString>& serialsForTitle) {
                bool anyData = false;
                for (auto& v : series) if (!v.empty()) { anyData = true; break; }
                if (!anyData) return;

                TCanvas* cHist = new TCanvas(Form("c_%d_hist", canvasCount), metric + " " + axis + " RatioHist",
                        900 * nCols, 850 * nRows);
                cHist->Divide(nCols, nRows);
                for (int i = 0; i < numItems; ++i) {
                    cHist->cd(i + 1);
                    gPad->SetGrid();
                    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.05);
                    gPad->SetTopMargin(0.14); gPad->SetBottomMargin(0.14);
                    if (series[i].empty()) {
                        TLatex* t = new TLatex(0.5, 0.5, "No Data"); t->SetNDC(); t->SetTextAlign(22);
                        t->SetTextFont(132); t->SetTextSize(0.05); t->Draw();
                        continue;
                    }

                    int nEntriesR = 0;
                    double lo = 1e9, hi = -1e9;
                    for (auto& s : series[i]) for (double v : s.second) {
                        double pct = (v - 1.0) * 100.0;
                        if (pct < lo) lo = pct; if (pct > hi) hi = pct;
                        nEntriesR++;
                    }
                    const double binWidthR = NiceBinWidth(hi - lo, nEntriesR);
                    double loR = std::floor((lo - binWidthR / 2) / binWidthR) * binWidthR;
                    double hiR = std::ceil((hi + binWidthR / 2) / binWidthR) * binWidthR;
                    loR = std::min(loR, -binWidthR); hiR = std::max(hiR, binWidthR);
                    const int nBinsR = (int)std::lround((hiR - loR) / binWidthR);

                    TLegend* legH = new TLegend(0.66, 0.58, 0.94, 0.76);
                    legH->SetBorderSize(0); legH->SetFillStyle(0); legH->SetTextFont(132); legH->SetTextSize(0.036);
                    TPaveText* ptH = new TPaveText(0.20, 0.58, 0.64, 0.76, "NDC");
                    ptH->SetFillColor(kWhite); ptH->SetFillStyle(1001); ptH->SetBorderSize(1);
                    ptH->SetTextFont(132); ptH->SetTextAlign(12); ptH->SetTextSize(0.024);

                    bool first = true;
                    std::vector<double> perTagMeans;   // each tag's own Mean(%), for mu/sigma below
                    double maxCount = 0;
                    for (size_t s = 0; s < series[i].size(); ++s) {
                        const TString& label = series[i][s].first;
                        const std::vector<double>& vals = series[i][s].second;
                        if (vals.empty()) continue;
                        int color = ColorForSeries((int)s);
                        TH1F* h = new TH1F(Form("h_%d_%d_%zu", canvasCount, i, s), "", nBinsR, loR, hiR);
                        for (double v : vals) h->Fill((v - 1.0) * 100.0);
                        h->SetLineColor(color); h->SetLineWidth(3);
                        h->SetFillColorAlpha(color, 0.22); h->SetFillStyle(1001);
                        h->SetTitle(Form(";Deviation from Ref [%%];Entries / %.0f%%", binWidthR));
                        h->GetXaxis()->SetLabelSize(0.050); h->GetXaxis()->SetTitleSize(0.056);
                        h->GetXaxis()->SetTitleOffset(1.05); h->GetXaxis()->SetLabelFont(132); h->GetXaxis()->SetTitleFont(132);
                        h->GetXaxis()->SetNdivisions(505);
                        h->GetYaxis()->SetLabelSize(0.050); h->GetYaxis()->SetTitleSize(0.056);
                        h->GetYaxis()->SetTitleOffset(1.25); h->GetYaxis()->SetLabelFont(132); h->GetYaxis()->SetTitleFont(132);
                        h->GetYaxis()->SetNdivisions(505);
                        if (h->GetMaximum() > maxCount) maxCount = h->GetMaximum();
                        h->Draw(first ? "HIST" : "HIST SAME");
                        first = false;
                        legH->AddEntry(h, label, "f");
                        ptH->AddText(Form("%s: M=%+.1f%%, R=%.2f%%", label.Data(), h->GetMean(), h->GetRMS()));
                        perTagMeans.push_back(h->GetMean());
                    }

                    if (perTagMeans.size() > 1) {
                        double mu = 0; for (double v : perTagMeans) mu += v; mu /= perTagMeans.size();
                        double sigma = 0; for (double v : perTagMeans) sigma += (v - mu) * (v - mu);
                        sigma = std::sqrt(sigma / perTagMeans.size());
                        ptH->AddText(Form("#mu=%+.1f%%, #sigma=%.1f%%", mu, sigma));
                    }
                    gPad->cd();
                    TH1F* frameH = (TH1F*)gPad->GetListOfPrimitives()->FindObject(Form("h_%d_%d_%zu", canvasCount, i, (size_t)0));
                    if (frameH) frameH->SetMaximum(maxCount * 1.8);
                    if (loR < 0.0 && hiR > 0.0) {
                        TLine* refLine = new TLine(0.0, 0, 0.0, maxCount * 1.8);
                        refLine->SetLineStyle(2); refLine->SetLineColor(kGray + 2); refLine->SetLineWidth(2);
                        refLine->Draw();
                    }
                    legH->Draw();
                    ptH->Draw();
                    TLatex* tt = new TLatex(0.5, 0.965,
                            Form("%s: %s (%s-axis)  #minus  Ratio Distribution",
                                metricLabel.Data(), serialsForTitle[i].Data(), axis.Data()));
                    tt->SetNDC(); tt->SetTextAlign(22); tt->SetTextFont(132); tt->SetTextSize(0.040); tt->Draw();
                }
                cHist->Print(pdfPath, Form("pdf Title:%s (%s-axis) Ratio Distribution", metricLabel.Data(), axis.Data()));
                if (savePNG) cHist->SaveAs(Form("./Data/image/Uniformity/Overlay_%s_%sRatioHist%s%s.png",
                            metric.Data(), axis.Data(), sList.Data(), dList.Data()));
                delete cHist;
            };


            auto drawPullHistPage = [&](const std::vector<std::vector<std::pair<TString, std::vector<double>>>>& series,
                    int numItems, int nCols, int nRows,
                    const std::vector<TString>& serialsForTitle) {
                bool anyData = false;
                for (auto& v : series) if (!v.empty()) { anyData = true; break; }
                if (!anyData) return;

                TCanvas* cPull = new TCanvas(Form("c_%d_pull", canvasCount), metric + " " + axis + " PullHist",
                        900 * nCols, 850 * nRows);
                cPull->Divide(nCols, nRows);
                for (int i = 0; i < numItems; ++i) {
                    cPull->cd(i + 1);
                    gPad->SetGrid();
                    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.05);
                    gPad->SetTopMargin(0.14); gPad->SetBottomMargin(0.14);
                    if (series[i].empty()) {
                        TLatex* t = new TLatex(0.5, 0.5, "No Data"); t->SetNDC(); t->SetTextAlign(22);
                        t->SetTextFont(132); t->SetTextSize(0.05); t->Draw();
                        continue;
                    }
                    int nEntriesP = 0;
                    double lo = 1e9, hi = -1e9;
                    for (auto& s : series[i]) for (double v : s.second) { if (v < lo) lo = v; if (v > hi) hi = v; nEntriesP++; }

                    const double binWidthP = std::min(1.0, NiceBinWidth(8.0, nEntriesP));
                    double loP = std::floor((lo - binWidthP / 2) / binWidthP) * binWidthP;
                    double hiP = std::ceil((hi + binWidthP / 2) / binWidthP) * binWidthP;
                    loP = std::min(loP, -4.0); hiP = std::max(hiP, 4.0);   // always show the +-4sigma frame
                    const int nBinsP = (int)std::lround((hiP - loP) / binWidthP);

                                        const int legPCols = 2;
                    int legPRows = (int)std::ceil((double)(series[i].size() + 1) / legPCols);
                    double legPRowH = 0.045;
                    double legPTop = 0.90;
                    TLegend* legP = new TLegend(0.60, legPTop - legPRows * legPRowH, 0.94, legPTop);
                    legP->SetNColumns(std::min((int)series[i].size() + 1, legPCols));
                    legP->SetBorderSize(0); legP->SetFillStyle(0); legP->SetTextFont(132); legP->SetTextSize(0.030);

                    bool first = true;
                    double maxCount = 0;
                    int totalN = 0;
                    for (size_t s = 0; s < series[i].size(); ++s) {
                        const TString& label = series[i][s].first;
                        const std::vector<double>& vals = series[i][s].second;
                        if (vals.empty()) continue;
                        int color = ColorForSeries((int)s);
                        TH1F* h = new TH1F(Form("hp_%d_%d_%zu", canvasCount, i, s), "", nBinsP, loP, hiP);
                        for (double v : vals) { h->Fill(v); totalN++; }
                        h->SetLineColor(color); h->SetLineWidth(3);
                        h->SetFillColorAlpha(color, 0.22); h->SetFillStyle(1001);
                        h->SetTitle(";Pull;Entries");
                        h->GetXaxis()->SetLabelSize(0.050); h->GetXaxis()->SetTitleSize(0.056);
                        h->GetXaxis()->SetTitleOffset(1.05); h->GetXaxis()->SetLabelFont(132); h->GetXaxis()->SetTitleFont(132);
                        h->GetXaxis()->SetNdivisions(505);
                        h->GetYaxis()->SetLabelSize(0.050); h->GetYaxis()->SetTitleSize(0.056);
                        h->GetYaxis()->SetTitleOffset(1.25); h->GetYaxis()->SetLabelFont(132); h->GetYaxis()->SetTitleFont(132);
                        h->GetYaxis()->SetNdivisions(505);
                        if (h->GetMaximum() > maxCount) maxCount = h->GetMaximum();
                        h->Draw(first ? "HIST" : "HIST SAME");
                        first = false;
                        legP->AddEntry(h, label, "f");
                    }

                    double gausPeak = 0;
                    if (totalN > 0) {
                        gausPeak = totalN * binWidthP / std::sqrt(2 * TMath::Pi());
                        TF1* fGaus = new TF1(Form("fGausPull_%d_%d", canvasCount, i),
                                Form("%d*%f/sqrt(2*TMath::Pi())*exp(-0.5*x*x)", totalN, binWidthP), loP, hiP);
                        fGaus->SetLineColor(kBlack); fGaus->SetLineStyle(1); fGaus->SetLineWidth(2);
                        fGaus->Draw("same");
                        legP->AddEntry(fGaus, "Gaussian", "l");
                    }

                    gPad->cd();
                    TH1F* frameP = (TH1F*)gPad->GetListOfPrimitives()->FindObject(Form("hp_%d_%d_%zu", canvasCount, i, (size_t)0));
                    if (frameP) frameP->SetMaximum(std::max(maxCount, gausPeak) * 1.2);
                    TLine* zeroLine = new TLine(0, 0, 0, std::max(maxCount, gausPeak) * 1.2);
                    zeroLine->SetLineStyle(2); zeroLine->SetLineColor(kGray + 2); zeroLine->SetLineWidth(2);
                    zeroLine->Draw();
                    legP->Draw();
                    TLatex* tt = new TLatex(0.5, 0.965,
                            Form("%s: %s (%s-axis)  #minus  Pull Distribution",
                                metricLabel.Data(), serialsForTitle[i].Data(), axis.Data()));
                    tt->SetNDC(); tt->SetTextAlign(22); tt->SetTextFont(132); tt->SetTextSize(0.040); tt->Draw();
                }
                cPull->Print(pdfPath, Form("pdf Title:%s (%s-axis) Pull Distribution", metricLabel.Data(), axis.Data()));
                if (savePNG) cPull->SaveAs(Form("./Data/image/Uniformity/Overlay_%s_%sPullHist%s%s.png",
                            metric.Data(), axis.Data(), sList.Data(), dList.Data()));
                delete cPull;
            };

            if (metric == "PedStability") {
                int numItems = allSerials.size();
                int nTags = (int)tags.size();
                int nCols = numItems, nRows = 1;   // one row (Monitor | PMT1 | PMT2)
                TCanvas *c = new TCanvas(Form("c_%d", canvasCount), metric + " " + axis, 800 * nCols, 850 * nRows);

                c->cd();
                TLatex* pageTitle = new TLatex(0.5, 0.985, Form("%s vs Angle (%s-axis)", yTitle.Data(), axis.Data()));
                pageTitle->SetNDC(); pageTitle->SetTextAlign(23); pageTitle->SetTextFont(132); pageTitle->SetTextSize(0.022);
                pageTitle->Draw();

                TPad* gridPad = new TPad(Form("gridPad_%d", canvasCount), "grid", 0.0, 0.0, 1.0, 0.96);
                gridPad->Draw();
                gridPad->Divide(nCols, nRows);

                for (int i = 0; i < numItems; ++i) {
                    TString currentSerial = allSerials[i];

                    double yMax = -1e9, yMin = 1e9;
                    std::vector<TGraphErrors*> wgraphs(nTags, nullptr);
                    for (int tagIdx = 0; tagIdx < nTags; ++tagIdx) {
                        TFile *f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tags[tagIdx].Data()));
                        if (!f) continue;
                        TGraphErrors *gr = (TGraphErrors*)f->Get(Form("gr_%s_%s_%s", currentSerial.Data(), axis.Data(), metric.Data()));
                        if (gr) {
                            wgraphs[tagIdx] = (TGraphErrors*)gr->Clone();
                            for (int k = 0; k < wgraphs[tagIdx]->GetN(); ++k) {
                                double y = wgraphs[tagIdx]->GetY()[k], ey = wgraphs[tagIdx]->GetEY()[k];
                                if (y + ey > yMax) yMax = y + ey;
                                if (y - ey < yMin) yMin = y - ey;
                            }
                        }
                        f->Close();
                    }

                    TVirtualPad* pmtPad = gridPad->cd(i + 1);
                    pmtPad->SetLeftMargin(0.15);
                    pmtPad->Divide(1, nTags, 0.0, 0.006);

                    for (int tagIdx = 0; tagIdx < nTags; ++tagIdx) {
                        bool isLast = (tagIdx == nTags - 1);
                        pmtPad->cd(tagIdx + 1);
                        gPad->SetFrameLineWidth(1);

                        gPad->SetLeftMargin(0.15);
                        gPad->SetTopMargin(0.01);
                        gPad->SetBottomMargin(isLast ? 0.28 : 0.01);
                        gPad->SetRightMargin(0.04);

                        TString label = tags[tagIdx];
                        if (tagLabels.count(tags[tagIdx])) label = tagLabels[tags[tagIdx]];

                        if (!wgraphs[tagIdx] || wgraphs[tagIdx]->GetN() == 0) {
                            TLatex* t = new TLatex(0.5, 0.5, "No Data"); t->SetNDC(); t->SetTextAlign(22); t->Draw();
                            continue;
                        }

                        wgraphs[tagIdx]->SetMarkerStyle(MarkerForTag(tagIdx));
                        wgraphs[tagIdx]->SetMarkerColor(ColorForTag(tagIdx));
                        wgraphs[tagIdx]->SetLineColor(ColorForTag(tagIdx));
                        wgraphs[tagIdx]->SetMarkerSize(1.0);
                        wgraphs[tagIdx]->SetTitle(isLast ? Form(";%s;%s", OverlayAngleAxisTitle(), yTitle.Data())
                                : Form(";;%s", yTitle.Data()));

                        if (yMax > yMin) {
                            double diff = yMax - yMin;
                            wgraphs[tagIdx]->SetMinimum(yMin - diff * 0.15);
                            wgraphs[tagIdx]->SetMaximum(yMax + diff * 0.15);
                        }
                        wgraphs[tagIdx]->GetXaxis()->SetTitleSize(isLast ? 0.11 : 0.0);
                        wgraphs[tagIdx]->GetXaxis()->SetLabelSize(isLast ? 0.10 : 0.0);
                        wgraphs[tagIdx]->GetXaxis()->SetTitleOffset(1.1);
                        wgraphs[tagIdx]->GetYaxis()->SetTitleSize(0.075); wgraphs[tagIdx]->GetYaxis()->SetLabelSize(0.070);
                        wgraphs[tagIdx]->GetYaxis()->SetTitleOffset(0.95);
                        wgraphs[tagIdx]->Draw("APZ");

                        TLatex* wl = new TLatex(0.20, 0.85, Form("%s (%s)", label.Data(), currentSerial.Data()));
                        wl->SetNDC(); wl->SetTextFont(132); wl->SetTextSize(0.10); wl->SetTextColor(ColorForTag(tagIdx));
                        wl->Draw();
                    }
                }

                c->Print(pdfPath, Form("pdf Title:%s (%s-axis)", metric.Data(), axis.Data()));
                if (savePNG) c->SaveAs(Form("./Data/image/Uniformity/Overlay_%s_%s%s%s.png", metric.Data(), axis.Data(), sList.Data(), dList.Data()));
                delete c;
                continue;
            }

            if ((baseMetric == "RelQE" || baseMetric == "RawQE" || baseMetric == "CorrRelQE" || baseMetric == "CorrRawQE" ||
                        baseMetric == "PoissonQE" || baseMetric == "RawPoissonQE" ||
                        baseMetric == "RelGain" || baseMetric == "RawSPE" || baseMetric == "TTS" ||
                        baseMetric == "FWHM" || baseMetric == "Sigma" ||
                        baseMetric == "ChargeResolution" ||
                        baseMetric == "MonNorm_RawQE" || baseMetric == "MonNorm_CorrRawQE" || baseMetric == "MonNorm_RawPoissonQE" ||
                        baseMetric == "MonNorm_PoissonQE" || baseMetric == "MonNorm_RelGain")) {
                const std::vector<TString>& displaySerials = baseMetric.BeginsWith("MonNorm_") ? testSerials : allSerials;
                int numItems = displaySerials.size();
                int nTags = (int)tags.size();
                int nCols = numItems, nRows = 1;   // one row (Monitor | PMT1 | PMT2), not a 2x2-style grid
                bool single = (nTags == 1);

                double fscale = 2.0 / nCols;

                double tscale = std::min(1.0, std::max(0.6, nTags / 4.0));

                auto loadGr = [&](const TString& serial, int tagIdx) -> TGraphErrors* {
                    TFile* f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tags[tagIdx].Data()));
                    if (!f) return nullptr;
                    TGraphErrors* gr = (TGraphErrors*)f->Get(Form("gr_%s_%s_%s", serial.Data(), axis.Data(), baseMetric.Data()));
                    TGraphErrors* out = gr ? (TGraphErrors*)gr->Clone() : nullptr;
                    if (out && isRelErr) { TGraphErrors* re = ToRelErrGraph(out); delete out; out = re; }
                    if (out && isCenterNorm) { TGraphErrors* cn = ToCenterNormGraph(out); delete out; out = cn; }
                    if (out) out->Sort();   // acquisition order, not angle order -- band fill needs monotonic X
                    f->Close();
                    return out;
                };

                TCanvas* c = new TCanvas(Form("c_%d", canvasCount), metric + " " + axis, 800 * nCols, 1200 * nRows);
                c->Divide(nCols, nRows);

                double sharedYmin = 1e9, sharedYmax = -1e9;
                for (int i2 = 0; i2 < numItems; ++i2) {
                    if (displaySerials[i2] == monitorSerial) continue;
                    for (int tagIdx2 = 0; tagIdx2 < nTags; ++tagIdx2) {
                        TGraphErrors* g2 = loadGr(displaySerials[i2], tagIdx2);
                        if (!g2) continue;
                        for (int k2 = 0; k2 < g2->GetN(); ++k2) {
                            double y2 = g2->GetY()[k2], ey2 = g2->GetEY()[k2];
                            if (y2 <= 0) continue;
                            if (y2 + ey2 > sharedYmax) sharedYmax = y2 + ey2;
                            if (y2 - ey2 < sharedYmin) sharedYmin = y2 - ey2;
                        }
                        delete g2;
                    }
                }

                std::vector<std::vector<std::pair<TString, std::vector<double>>>> itemRatioSeries(numItems);
                std::vector<std::vector<std::pair<TString, std::vector<double>>>> itemPullSeries(numItems);

                for (int i = 0; i < numItems; ++i) {
                    TString serial = displaySerials[i];
                    TVirtualPad* pmtPad = c->cd(i + 1);

                    TString metricName = yTitle; 
                    if (metricName.Contains(" [")) metricName = metricName(0, metricName.Index(" ["));
                    TLatex* hdr = new TLatex(0.5, 0.985, Form("%s;  %s, %s-axis", metricName.Data(), serial.Data(), axis.Data()));
                    hdr->SetNDC(); hdr->SetTextAlign(23); hdr->SetTextFont(62); hdr->SetTextSize((single ? 0.045 : 0.040) * fscale); hdr->Draw();

                    std::vector<TGraphErrors*> grs(nTags, nullptr);
                    for (int tagIdx = 0; tagIdx < nTags; ++tagIdx) grs[tagIdx] = loadGr(serial, tagIdx);
                    TGraphErrors* grRef = grs[0];

                    auto rangeOf = [&](TGraphErrors* g, double& lo, double& hi) {
                        lo = 1e9; hi = -1e9;
                        if (!g) return;
                        for (int k = 0; k < g->GetN(); ++k) {
                            double y = g->GetY()[k], ey = g->GetEY()[k];
                            if (y <= 0) continue;
                            if (y + ey > hi) hi = y + ey; if (y - ey < lo) lo = y - ey;
                        }
                    };

                    double gxLo = 1e9, gxHi = -1e9;
                    for (int tagIdx = 0; tagIdx < nTags; ++tagIdx) {
                        if (!grs[tagIdx]) continue;
                        for (int k = 0; k < grs[tagIdx]->GetN(); ++k) {
                            double x = grs[tagIdx]->GetX()[k];
                            if (x < gxLo) gxLo = x; if (x > gxHi) gxHi = x;
                        }
                    }
                    if (gxHi > gxLo) { double mx = (gxHi - gxLo) * 0.05; gxLo -= mx; gxHi += mx; }
                    else if (gxHi == -1e9) { gxLo = 0; gxHi = 1; }   
                    else { gxLo -= 1; gxHi += 1; }

                    double headTop  = 0.955;                    
                    double ratioH   = single ? 0.0 : 0.20;      
                    double valZone  = headTop - ratioH;         
                    double vH       = valZone / nTags;

                    for (int tagIdx = 0; tagIdx < nTags; ++tagIdx) {
                        double vTop = headTop - tagIdx * vH;
                        double vBot = headTop - (tagIdx + 1) * vH;
                        TString wl = tags[tagIdx]; if (tagLabels.count(tags[tagIdx])) wl = tagLabels[tags[tagIdx]];
                        int color = ColorForTag(tagIdx);

                        pmtPad->cd();
                        TPad* vPad = new TPad(Form("v_%d_%d_%d", canvasCount, i, tagIdx), "", 0.0, vBot, 1.0, vTop);
                        vPad->SetLeftMargin(0.22); vPad->SetRightMargin(0.05);
                        vPad->SetTopMargin(0.02); vPad->SetBottomMargin(single ? 0.16 : 0.03);
                        vPad->SetGridx(); vPad->SetGridy();
                        vPad->Draw(); vPad->cd();

                        if (!grs[tagIdx] || grs[tagIdx]->GetN() == 0) {
                            TLatex* t = new TLatex(0.5, 0.5, "No Data"); t->SetNDC(); t->SetTextAlign(22); t->Draw();
                            continue;
                        }

                        grs[tagIdx]->SetMarkerStyle(MarkerForTag(tagIdx)); grs[tagIdx]->SetMarkerColor(color); grs[tagIdx]->SetLineColor(color); grs[tagIdx]->SetMarkerSize(0.7);
                        double ymin, ymax;
                        if (serial == monitorSerial || sharedYmax <= sharedYmin) {
                            rangeOf(grs[tagIdx], ymin, ymax);
                        } else {
                            ymin = sharedYmin; ymax = sharedYmax;
                        }
                        double dy = (ymax > ymin) ? ymax - ymin : std::abs(ymax) * 0.2 + 1e-6;
                        double vMin = ymin - dy * 0.10, vMax = ymax + dy * 0.55;
                        grs[tagIdx]->SetMinimum(vMin); grs[tagIdx]->SetMaximum(vMax);
                        grs[tagIdx]->SetTitle(single ? Form(";%s;%s", OverlayAngleAxisTitle(), yTitle.Data()) : ";;");
                        grs[tagIdx]->GetXaxis()->SetLabelSize((single ? 0.05 : 0.0) * fscale);
                        grs[tagIdx]->GetXaxis()->SetTitleSize((single ? 0.06 : 0.0) * fscale);
                        grs[tagIdx]->GetYaxis()->SetTitleSize((single ? 0.050 : 0.0) * fscale);

                        grs[tagIdx]->GetYaxis()->SetLabelSize(single ? 0.060 * fscale : 0.20 * tscale);
                        grs[tagIdx]->GetYaxis()->SetTitleOffset(single ? 1.3 : 0.95);
                        grs[tagIdx]->GetYaxis()->SetNdivisions(505);
                        grs[tagIdx]->SetFillColorAlpha(color, 0.25);
                        grs[tagIdx]->Draw("A3");     // axis + shaded +-1sigma band, no points yet
                        grs[tagIdx]->Draw("PZ same"); // points+bars on top of the band
                        grs[tagIdx]->GetXaxis()->SetLimits(gxLo, gxHi);   // shared angular window
                        TLatex* lab = new TLatex(0.28, 0.80, wl.Data());
                        lab->SetNDC(); lab->SetTextFont(22); lab->SetTextSize(0.16 * tscale); lab->SetTextColor(color); lab->Draw();
                    }

                    if (!single) {
                        pmtPad->cd();
                        TLatex* sharedYTitle = new TLatex(0.07, headTop - valZone / 2.0, yTitle.Data());
                        sharedYTitle->SetNDC(); sharedYTitle->SetTextAlign(22); sharedYTitle->SetTextAngle(90);

                        sharedYTitle->SetTextFont(132); sharedYTitle->SetTextSize(0.065 * std::sqrt(fscale)); sharedYTitle->Draw();
                    }

                    if (single) continue;

                    pmtPad->cd();
                    TPad* rPad = new TPad(Form("r_%d_%d", canvasCount, i), "", 0.0, 0.0, 1.0, ratioH);
                    rPad->SetLeftMargin(0.22); rPad->SetRightMargin(0.05);
                    rPad->SetTopMargin(0.03); rPad->SetBottomMargin(0.40);
                    rPad->SetGridx(); rPad->SetGridy();
                    rPad->Draw(); rPad->cd();

                    TString refLab = tags[refIdx]; if (tagLabels.count(tags[refIdx])) refLab = tagLabels[tags[refIdx]];
                    TMultiGraph* mgR = new TMultiGraph();
                    double rmin = 1e9, rmax = -1e9;

                    std::vector<std::array<double,2>> allDevXY;
                    std::vector<int> allDevColor;
                    for (int tagIdx = 1; tagIdx < nTags; ++tagIdx) {   // skip reference
                        if (!grRef || !grs[tagIdx]) continue;
                        TGraphErrors* grRatio = new TGraphErrors();
                        int color = ColorForTag(tagIdx);
                        grRatio->SetMarkerStyle(MarkerForTag(tagIdx)); grRatio->SetMarkerColor(color); grRatio->SetLineColor(color); grRatio->SetMarkerSize(0.7);
                        int pointCount = 0;
                        std::vector<double> ratioVals;   // raw values -- projected onto 1D below
                        std::vector<double> pullVals;    // (target-ref)/sqrt(errT^2+errR^2), absolute-error version of ratioVals
                        for (int idxRef = 0; idxRef < grRef->GetN(); ++idxRef) {
                            double angleRef = grRef->GetX()[idxRef], valueRef = grRef->GetY()[idxRef];
                            if (valueRef <= 0) continue;   // skip invalid/sentinel reference points (TTS -1/-2)
                            for (int idxTgt = 0; idxTgt < grs[tagIdx]->GetN(); ++idxTgt) {
                                if (std::abs(grs[tagIdx]->GetX()[idxTgt] - angleRef) < 0.1) {
                                    double valueTgt = grs[tagIdx]->GetY()[idxTgt];
                                    if (valueTgt <= 0) break;   // skip invalid target point
                                    double ratioVal = valueTgt / valueRef;
                                    // Plotted as % deviation from Ref, zero-centered -- see the
                                    // equivalent conversion/comment at the other ratio-pad block.
                                    grRatio->SetPoint(pointCount, angleRef, (ratioVal - 1.0) * 100.0);
                                    // Propagated stat error of the ratio, treating ref and
                                    // target as independent: sigma_r/r = sqrt((sigma_t/t)^2+(sigma_r/r)^2).
                                    // Ref is SHARED across every non-ref tag on this page, so its
                                    // fluctuation is actually common to all points here -- this
                                    // slightly overestimates the point-to-point (relative) scatter,
                                    // but is the right size for judging "is this point consistent
                                    // with 1 given its own statistics".
                                    double errRef = grRef->GetErrorY(idxRef), errTgt = grs[tagIdx]->GetErrorY(idxTgt);
                                    double relErrRef = (valueRef != 0) ? errRef / valueRef : 0;
                                    double relErrTgt = (valueTgt != 0) ? errTgt / valueTgt : 0;
                                    grRatio->SetPointError(pointCount, 0, ratioVal * 100.0 * std::sqrt(relErrRef*relErrRef + relErrTgt*relErrTgt));
                                    pointCount++;
                                    ratioVals.push_back(ratioVal);   // all angles, not restricted to inner +-60deg
                                    double denom = std::sqrt(errRef*errRef + errTgt*errTgt);
                                    if (denom > 0) pullVals.push_back((valueTgt - valueRef) / denom);
                                    double ratioPct = (ratioVal - 1.0) * 100.0;
                                    if (ratioPct > rmax) rmax = ratioPct; if (ratioPct < rmin) rmin = ratioPct;
                                    allDevXY.push_back({angleRef, ratioPct});
                                    allDevColor.push_back(color);
                                    break;
                                }
                            }
                        }
                        if (pointCount > 0) {
                            mgR->Add(grRatio, "pz");
                            TString wlLabel = tags[tagIdx]; if (tagLabels.count(tags[tagIdx])) wlLabel = tagLabels[tags[tagIdx]];
                            itemRatioSeries[i].push_back({wlLabel, ratioVals});
                            itemPullSeries[i].push_back({wlLabel, pullVals});
                        }
                    }
                    if (mgR->GetListOfGraphs()) {
                        mgR->SetTitle(Form(";%s;Dev [%%]", OverlayAngleAxisTitle()));
                        mgR->Draw("AP");
                        mgR->GetXaxis()->SetLimits(gxLo, gxHi);

                        double coreLo = rmin, coreHi = rmax;
                        if (allDevXY.size() >= 3) {
                            std::vector<double> devs; devs.reserve(allDevXY.size());
                            for (auto& xy : allDevXY) devs.push_back(xy[1]);
                            std::vector<double> sorted = devs;
                            std::sort(sorted.begin(), sorted.end());
                            double median = sorted[sorted.size() / 2];
                            std::vector<double> absDev; absDev.reserve(devs.size());
                            for (double v : devs) absDev.push_back(std::abs(v - median));
                            std::sort(absDev.begin(), absDev.end());
                            double mad = absDev[absDev.size() / 2];
                            if (mad > 1e-9) {
                                coreLo = median - 4.0 * mad;
                                coreHi = median + 4.0 * mad;
                                // Never crop tighter than the raw range needs, and
                                // never wider either -- clamp to what's actually there.
                                coreLo = std::max(coreLo, rmin);
                                coreHi = std::min(coreHi, rmax);
                            }
                        }
                        double padLo = coreLo, padHi = coreHi;
                        if (coreHi > coreLo) {
                            double dr = coreHi - coreLo;
                            padLo = coreLo - dr * 0.15; padHi = coreHi + dr * 0.15;
                        } else if (rmax > rmin) {
                            double dr = rmax - rmin;
                            padLo = rmin - dr * 0.15; padHi = rmax + dr * 0.15;
                        }
                        mgR->SetMinimum(padLo);
                        mgR->SetMaximum(padHi);

                        mgR->GetXaxis()->SetLabelSize(0.10); mgR->GetXaxis()->SetTitleSize(0.12); mgR->GetXaxis()->SetTitleOffset(1.1);
                        mgR->GetYaxis()->SetLabelSize(0.08); mgR->GetYaxis()->SetTitleSize(0.10); mgR->GetYaxis()->SetTitleOffset(0.6);
                        mgR->GetYaxis()->SetNdivisions(505);
                        gPad->Modified(); gPad->Update();
                        TLine* one = new TLine(mgR->GetXaxis()->GetXmin(), 0.0, mgR->GetXaxis()->GetXmax(), 0.0);
                        one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->SetLineWidth(2); one->Draw();


                        double yLo = padLo, yHi = padHi;
                        for (size_t k = 0; k < allDevXY.size(); ++k) {
                            double ax = allDevXY[k][0], ay = allDevXY[k][1];
                            if (ay >= yLo && ay <= yHi) continue;   // inside the visible range already
                            bool goesUp = (ay > yHi);
                            double edgeY = goesUp ? yHi : yLo;
                            TMarker* tri = new TMarker(ax, edgeY, goesUp ? kFullTriangleUp : kFullTriangleDown);
                            tri->SetMarkerColor(allDevColor[k]); tri->SetMarkerSize(1.1);
                            tri->Draw();
                            double labY = goesUp ? (yHi - (yHi - yLo) * 0.09) : (yLo + (yHi - yLo) * 0.09);
                            TLatex* lab = new TLatex(ax, labY, Form("%.0f%%", ay));
                            lab->SetTextFont(132); lab->SetTextSize(0.06); lab->SetTextAlign(22);
                            lab->SetTextColor(allDevColor[k]); lab->Draw();
                        }
                    }
                }

                c->Print(pdfPath, Form("pdf Title:%s (%s-axis)", metric.Data(), axis.Data()));
                if (savePNG) c->SaveAs(Form("./Data/image/Uniformity/Overlay_%s_%s%s%s.png", metric.Data(), axis.Data(), sList.Data(), dList.Data()));
                delete c;


            }


            canvasCount++;


            const std::vector<TString>& generalDisplaySerials =
                baseMetric.BeginsWith("MonNorm_") ? testSerials : allSerials;
            int numItems = generalDisplaySerials.size();
            int nCols = numItems, nRows = 1;   // one row (Monitor | PMT1 | PMT2), not a 2x2-style grid

            std::vector<double> itemTopMax(numItems, -1e9), itemTopMin(numItems, 1e9);
            std::vector<double> itemBotMax(numItems, -1e9), itemBotMin(numItems, 1e9);

            for (int i = 0; i < numItems; ++i) {
                TString currentSerial = generalDisplaySerials[i];
                std::vector<TGraphErrors*> tempGraphs(tags.size(), nullptr);

                for (size_t tagIdx = 0; tagIdx < tags.size(); ++tagIdx) {
                    TFile *fUniformity = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tags[tagIdx].Data()));
                    if (!fUniformity) continue;

                    TGraphErrors *gr = (TGraphErrors*)fUniformity->Get(Form("gr_%s_%s_%s", currentSerial.Data(), axis.Data(), baseMetric.Data()));
                    if (gr) {
                        tempGraphs[tagIdx] = isRelErr ? ToRelErrGraph(gr) :
                            isCenterNorm ? ToCenterNormGraph(gr) : (TGraphErrors*)gr->Clone();
                        if (!tempGraphs[tagIdx]) continue;   // ToCenterNormGraph found no valid Center points
                        for (int k = 0; k < tempGraphs[tagIdx]->GetN(); ++k) {
                            double y = tempGraphs[tagIdx]->GetY()[k];
                            double ey = tempGraphs[tagIdx]->GetEY()[k];
                            if (y + ey > itemTopMax[i]) itemTopMax[i] = y + ey;
                            if (y - ey < itemTopMin[i]) itemTopMin[i] = y - ey;
                        }
                    }
                    fUniformity->Close();
                }

                TGraphErrors *grRef = tempGraphs[0];
                if (grRef) {
                    for (size_t tagIdx = 1; tagIdx < tags.size(); ++tagIdx) {
                        TGraphErrors *grTgt = tempGraphs[tagIdx];
                        if (!grTgt) continue;

                        for (int idxRef = 0; idxRef < grRef->GetN(); ++idxRef) {
                            for (int idxTgt = 0; idxTgt < grTgt->GetN(); ++idxTgt) {
                                if (std::abs(grRef->GetX()[idxRef] - grTgt->GetX()[idxTgt]) < 0.1) {
                                    double valRef = grRef->GetY()[idxRef];
                                    double valTgt = grTgt->GetY()[idxTgt];
                                    if (valRef != 0) {
                                        double ratio = (valTgt / valRef - 1.0) * 100.0;
                                        if (ratio > itemBotMax[i]) itemBotMax[i] = ratio;
                                        if (ratio < itemBotMin[i]) itemBotMin[i] = ratio;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
                for (auto g : tempGraphs) { if (g) delete g; }
            }


            {
                double sharedTopMin = 1e9, sharedTopMax = -1e9;
                for (int i2 = 0; i2 < numItems; ++i2) {
                    if (generalDisplaySerials[i2] == monitorSerial) continue;
                    if (itemTopMax[i2] > itemTopMin[i2]) {
                        sharedTopMin = std::min(sharedTopMin, itemTopMin[i2]);
                        sharedTopMax = std::max(sharedTopMax, itemTopMax[i2]);
                    }
                }
                if (sharedTopMax > sharedTopMin) {
                    for (int i2 = 0; i2 < numItems; ++i2) {
                        if (generalDisplaySerials[i2] == monitorSerial) continue;
                        itemTopMin[i2] = sharedTopMin; itemTopMax[i2] = sharedTopMax;
                    }
                }
            }
            // =========================================================================
            TCanvas *c = new TCanvas(Form("c_%d", canvasCount), metric + " " + axis, 800 * nCols, 800 * nRows);
            c->Divide(nCols, nRows);

            std::vector<std::vector<std::pair<TString, std::vector<double>>>> itemRatioSeries(numItems);
            std::vector<std::vector<std::pair<TString, std::vector<double>>>> itemPullSeries(numItems);

            for (int i = 0; i < numItems; ++i) {
                TString currentSerial = generalDisplaySerials[i];
                c->cd(i + 1);

                TPad *padTop = new TPad(Form("pad_top_%d", i), "top", 0.0, 0.35, 1.0, 1.0);
                TPad *padBot = new TPad(Form("pad_bot_%d", i), "bot", 0.0, 0.0, 1.0, 0.35);

                double rMargin = (metric == "DarkCount") ? 0.12 : 0.05;
                padTop->SetBottomMargin(0.02);
                padTop->SetLeftMargin(0.15);
                padTop->SetRightMargin(rMargin);
                padTop->SetTopMargin(0.10);
                padBot->SetTopMargin(0.02);
                padBot->SetBottomMargin(0.3);
                padBot->SetLeftMargin(0.15);
                padBot->SetRightMargin(rMargin);

                padTop->Draw();
                padBot->Draw();

                // ----- TOP PAD -----
                padTop->cd();
                gPad->SetGrid();

                TMultiGraph *mgTop = new TMultiGraph();

                const int legCols = 4;
                int legRows = (int)std::ceil((double)tags.size() / legCols);

                double legRowH = 0.055;
                double legTopY = 1.0 - padTop->GetTopMargin() - 0.02;
                double legRightX = 1.0 - rMargin - 0.01;
                TLegend *leg = new TLegend(0.15, legTopY - legRows * legRowH, legRightX, legTopY);
                leg->SetNColumns(std::min((int)tags.size(), legCols));
                leg->SetBorderSize(1);
                leg->SetTextFont(132);
                leg->SetTextSize(0.040);

                double statsX2 = (metric == "DarkCount") ? 0.74 : 0.93;
                TPaveText *ptStats = new TPaveText(0.35, 0.65, statsX2, 0.88, "NDC");
                ptStats->SetFillColor(kWhite);
                ptStats->SetBorderSize(0);
                ptStats->SetTextAlign(12);
                ptStats->SetTextSize(0.035);

                TString labelRef = tags[refIdx];
                if (tagLabels.count(tags[refIdx])) labelRef = tagLabels[tags[refIdx]];
                ptStats->AddText(Form("#bf{Serial: %s | Ref: %s}", currentSerial.Data(), labelRef.Data()));

                std::vector<TGraphErrors*> graphs(tags.size(), nullptr);

                for (size_t tagIdx = 0; tagIdx < tags.size(); ++tagIdx) {
                    TFile *f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tags[tagIdx].Data()));
                    if (!f) continue;

                    TString label = tags[tagIdx];
                    if (tagLabels.count(tags[tagIdx])) label = tagLabels[tags[tagIdx]];

                    TGraphErrors *gr = (TGraphErrors*)f->Get(Form("gr_%s_%s_%s", currentSerial.Data(), axis.Data(), baseMetric.Data()));
                    if (gr) {
                        graphs[tagIdx] = isRelErr ? ToRelErrGraph(gr) :
                            isCenterNorm ? ToCenterNormGraph(gr) : (TGraphErrors*)gr->Clone();
                        if (!graphs[tagIdx]) { f->Close(); continue; }   // ToCenterNormGraph found no valid Center points
                        graphs[tagIdx]->Sort();   // points are stored in acquisition order, not angle order --
                                                  // the "3" band fill needs monotonic X or it self-crosses
                        graphs[tagIdx]->SetLineColor(ColorForTag(tagIdx));
                        graphs[tagIdx]->SetMarkerColor(ColorForTag(tagIdx));
                        graphs[tagIdx]->SetMarkerStyle(MarkerForTag(tagIdx));
                        graphs[tagIdx]->SetMarkerSize(1.0);

                        mgTop->Add(graphs[tagIdx], (metric == "PedStability") ? "pz" : "p");
                        leg->AddEntry(graphs[tagIdx], label, "p");
                    }
                    f->Close();
                }

                if (mgTop->GetListOfGraphs()) {
                    mgTop->SetTitle(Form("%s (%s-axis);;%s", currentSerial.Data(), axis.Data(), yTitle.Data()));

                    if (itemTopMax[i] > itemTopMin[i]) {
                        double diff = itemTopMax[i] - itemTopMin[i];
                        mgTop->SetMinimum(itemTopMin[i] - diff * 0.45);
                        mgTop->SetMaximum(itemTopMax[i] + diff * 0.90);
                    } else {
                        mgTop->SetMinimum(itemTopMin[i] * 0.5);
                        mgTop->SetMaximum(itemTopMax[i] * 1.5);
                    }

                    mgTop->Draw("AP");
                    mgTop->GetXaxis()->SetLabelSize(0);
                    mgTop->GetYaxis()->SetTitleSize(0.05);
                    mgTop->GetYaxis()->SetLabelSize(0.05);

                    // +-1sigma shaded band per tag, drawn under the points.
                    for (size_t tagIdx = 0; tagIdx < tags.size(); ++tagIdx) {
                        if (!graphs[tagIdx]) continue;
                        graphs[tagIdx]->SetFillColorAlpha(ColorForTag(tagIdx), 0.18);
                        graphs[tagIdx]->Draw("3 same");
                    }
                    mgTop->Draw("P same");   // points back on top of the bands

                    if (metric == "DarkCount") {
                        double refFactor = 0.0;
                        bool factorsConsistent = true;
                        int nFactors = 0;
                        for (size_t tagIdx = 0; tagIdx < tags.size(); ++tagIdx) {
                            TFile* fr = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tags[tagIdx].Data()));
                            if (!fr || fr->IsZombie()) { if (fr) fr->Close(); continue; }
                            TGraphErrors* gcnt = (TGraphErrors*)fr->Get(Form("gr_%s_%s_DarkCount", currentSerial.Data(), axis.Data()));
                            TGraphErrors* grat = (TGraphErrors*)fr->Get(Form("gr_%s_%s_DarkRate",  currentSerial.Data(), axis.Data()));
                            if (gcnt && grat && gcnt->GetN() > 0 && grat->GetN() > 0 && gcnt->GetY()[0] != 0) {
                                double f = grat->GetY()[0] / gcnt->GetY()[0];
                                if (nFactors == 0) refFactor = f;
                                else if (std::abs(f - refFactor) > 0.02 * refFactor) factorsConsistent = false;   // >2% off = a real livetime difference, not float noise
                                ++nFactors;
                            }
                            fr->Close();
                        }
                        if (nFactors > 0 && factorsConsistent) {
                            gPad->Update();
                            double angleRef = gPad->GetUxmax();
                            double ylo = gPad->GetUymin(), yhi = gPad->GetUymax();
                            TGaxis* axHz = new TGaxis(angleRef, ylo, angleRef, yhi, ylo * refFactor, yhi * refFactor, 510, "+L");
                            axHz->SetTitle("Dark Rate [Hz]");
                            axHz->SetLabelFont(132); axHz->SetTitleFont(132);
                            axHz->SetLabelSize(0.045); axHz->SetTitleSize(0.05); axHz->SetTitleOffset(1.15);
                            axHz->SetLineColor(kGray + 3); axHz->SetLabelColor(kGray + 3); axHz->SetTitleColor(kGray + 3);
                            axHz->Draw();
                        } else if (nFactors > 1 && !factorsConsistent) {
                            std::cout << Form("[WARNING] %s (%s-axis): DarkCount->Rate factor differs across tags "
                                    "(livetime mismatch) -- Dark Rate axis omitted, left as Counts only.",
                                    currentSerial.Data(), axis.Data()) << std::endl;
                            TLatex* noteHz = new TLatex(0.985, 0.5, "Dark Rate axis omitted: livetime differs across datasets");
                            noteHz->SetNDC(); noteHz->SetTextAngle(90); noteHz->SetTextAlign(21);
                            noteHz->SetTextFont(132); noteHz->SetTextSize(0.03); noteHz->SetTextColor(kGray + 2);
                            noteHz->Draw();
                        }
                    }
                }
                leg->Draw();


                // ----- BOTTOM PAD -----
                padBot->cd();
                gPad->SetGrid();
                TMultiGraph *mgBot = new TMultiGraph();

                TGraphErrors *grRef = graphs[0]; 

                for (size_t tagIdx = 1; tagIdx < tags.size(); ++tagIdx) {
                    TGraphErrors *grTgt = graphs[tagIdx];

                    TString labelTgt = tags[tagIdx];
                    if (tagLabels.count(tags[tagIdx])) labelTgt = tagLabels[tags[tagIdx]];

                    if (!grRef || !grTgt) {
                        ptStats->AddText(Form("vs %s: Missing Data", labelTgt.Data()));
                        continue;
                    }

                    TGraphErrors *grRatio = new TGraphErrors();
                    grRatio->SetLineColor(ColorForTag(tagIdx));
                    grRatio->SetMarkerColor(ColorForTag(tagIdx));
                    grRatio->SetMarkerStyle(MarkerForTag(tagIdx)); 
                    grRatio->SetMarkerSize(1.0);

                    double sumRatio = 0;
                    double maxDiff = 0;
                    int pCount = 0;
                    std::vector<double> ratioVals;   // raw values -- projected onto 1D below
                    std::vector<double> pullVals;    // (target-ref)/sqrt(errT^2+errR^2), absolute-error version of ratioVals

                    for (int idxRef = 0; idxRef < grRef->GetN(); ++idxRef) {
                        for (int idxTgt = 0; idxTgt < grTgt->GetN(); ++idxTgt) {
                            if (std::abs(grRef->GetX()[idxRef] - grTgt->GetX()[idxTgt]) < 0.1) {
                                double valRef = grRef->GetY()[idxRef];
                                double valTgt = grTgt->GetY()[idxTgt];

                                if (valRef != 0) {
                                    double ratio = valTgt / valRef;
                                    double diff = std::abs(valTgt - valRef);

                                    // Plotted as % deviation from Ref ((ratio-1)*100), zero-centered
                                    // so "how many % off" reads directly off the y-axis. ratioVals
                                    // below stays the raw ratio (used for the RMS/text stats and the
                                    // ratio-histogram page).
                                    grRatio->SetPoint(pCount, grRef->GetX()[idxRef], (ratio - 1.0) * 100.0);
                                    // Propagated stat error of the ratio (independent-ref/target
                                    // approximation -- see the equivalent comment at the top-pad
                                    // ratio construction above). The error's magnitude is the same
                                    // whether centered on 1 or 0, only the point itself shifts.
                                    double errRef = grRef->GetErrorY(idxRef), errTgt = grTgt->GetErrorY(idxTgt);
                                    double relErrRef = (valRef != 0) ? errRef / valRef : 0;
                                    double relErrTgt = (valTgt != 0) ? errTgt / valTgt : 0;
                                    grRatio->SetPointError(pCount, 0, ratio * 100.0 * std::sqrt(relErrRef*relErrRef + relErrTgt*relErrTgt));
                                    ratioVals.push_back(ratio);
                                    double denom = std::sqrt(errRef*errRef + errTgt*errTgt);
                                    if (denom > 0) pullVals.push_back((valTgt - valRef) / denom);

                                    sumRatio += ratio;
                                    if (diff > maxDiff) maxDiff = diff;
                                    pCount++;
                                }
                                break;
                            }
                        }
                    }

                    if (pCount > 0) {
                        mgBot->Add(grRatio, "pz");
                        double avgRatio = sumRatio / pCount;
                        // RMS of the ratio itself (point-to-point scatter)
                        double rmsRatio = 0;
                        for (double v : ratioVals) rmsRatio += (v - avgRatio) * (v - avgRatio);
                        rmsRatio = std::sqrt(rmsRatio / pCount);
                        ptStats->AddText(Form("vs %s : Avg Ratio = %.3f, RMS = %.2f%%, Max Diff = %.3f",
                                    labelTgt.Data(), avgRatio, rmsRatio * 100, maxDiff));
                        itemRatioSeries[i].push_back({labelTgt, ratioVals});
                        itemPullSeries[i].push_back({labelTgt, pullVals});
                    } else {
                        ptStats->AddText(Form("vs %s : X-axis mismatched", labelTgt.Data()));
                    }
                }

                padBot->cd();
                if (mgBot->GetListOfGraphs()) {
                    mgBot->SetTitle(Form(";%s;Dev [%%]", OverlayAngleAxisTitle()));

                    if (itemBotMax[i] > itemBotMin[i]) {
                        double diff = itemBotMax[i] - itemBotMin[i];
                        mgBot->SetMinimum(itemBotMin[i] - diff * 0.25);
                        mgBot->SetMaximum(itemBotMax[i] + diff * 0.25);
                    } else {
                        mgBot->SetMinimum(-20);
                        mgBot->SetMaximum(20);
                    }

                    mgBot->Draw("AP");

                    TLine *line = new TLine(mgBot->GetXaxis()->GetXmin(), 0.0, mgBot->GetXaxis()->GetXmax(), 0.0);
                    line->SetLineStyle(2); line->SetLineColor(kBlack); line->SetLineWidth(1);
                    line->Draw("same");

                    mgBot->GetXaxis()->SetTitleSize(0.12);
                    mgBot->GetXaxis()->SetLabelSize(0.10);
                    mgBot->GetYaxis()->SetTitleSize(0.10);
                    mgBot->GetYaxis()->SetTitleOffset(0.6);
                    mgBot->GetYaxis()->SetLabelSize(0.08);
                }
            }

            c->Print(pdfPath, Form("pdf Title:%s (%s-axis)", metric.Data(), axis.Data()));
            if (savePNG) c->SaveAs(Form("./Data/image/Uniformity/Overlay_%s_%s%s%s.png", metric.Data(), axis.Data(), sList.Data(), dList.Data()));
            delete c;

            drawPullHistPage(itemPullSeries, numItems, numItems, 1, generalDisplaySerials);
        }

    }

    auto noNegZero = [](double v) { return (std::abs(v) < 0.005) ? 0.0 : v; };

    // ================= Final PDF page: mean / RMS vs reference =================
    {
        auto ratioStat = [&](const TString& sn, const TString& metric, const TString& axis,
                const TString& tag, bool inner, double& mean, double& rms, int& npts)->bool {
            mean = 0; rms = 0; npts = 0;
            TFile* fr = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tags[refIdx].Data()));
            TFile* ft = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag.Data()));
            if (!fr || !ft || fr->IsZombie() || ft->IsZombie()) { if(fr)fr->Close(); if(ft)ft->Close(); return false; }
            TGraphErrors* gr = (TGraphErrors*)fr->Get(Form("gr_%s_%s_%s", sn.Data(), axis.Data(), metric.Data()));
            TGraphErrors* gt = (TGraphErrors*)ft->Get(Form("gr_%s_%s_%s", sn.Data(), axis.Data(), metric.Data()));
            std::vector<double> ratioVals;
            if (gr && gt) {
                for (int i = 0; i < gr->GetN(); ++i) {
                    double ang = gr->GetX()[i], vref = gr->GetY()[i];
                    // <=0 (not just ==0) also excludes the -1/-2 invalid-fit
                    // sentinels TTS uses (e.g. a failed ex-Gaussian fit) -- QE/Gain
                    // are never negative in practice, so this is a no-op for them.
                    if (vref <= 0) continue;
                    if (inner && std::abs(ang) > 60.0) continue;
                    for (int j = 0; j < gt->GetN(); ++j) {
                        if (std::abs(gt->GetX()[j] - ang) < 0.5) {
                            double vtgt = gt->GetY()[j];
                            if (vtgt <= 0) break;
                            ratioVals.push_back(vtgt / vref);
                            break;
                        }
                    }
                }
            }
            fr->Close(); ft->Close();
            if (ratioVals.empty()) return false;
            for (double v : ratioVals) mean += v; mean /= ratioVals.size();
            for (double v : ratioVals) rms += (v - mean) * (v - mean); rms = std::sqrt(rms / ratioVals.size());
            npts = ratioVals.size();
            return true;
        };

        TString refLabel = tags[refIdx]; if (tagLabels.count(tags[refIdx])) refLabel = tagLabels[tags[refIdx]];
        std::vector<std::pair<TString,TString>> sumMetrics = {{"RelQE","QE (PHC)"}, {"PoissonQE","QE (Poisson)"}, {"MonNorm_CorrRawQE","QE (Mon-Norm Corr)"}, {"RelGain","Gain"}, {"TTS","TTS"}};

        TCanvas* csum = new TCanvas("c_ratio_summary", "Ratio Summary", 1600, 2000);
        csum->cd();
        auto Tx = [&](double x, double y, const char* s, int font = 42, double sz = 0.020) {
            TLatex* t = new TLatex(x, y, s); t->SetNDC(); t->SetTextFont(font);
            t->SetTextSize(sz); t->SetTextAlign(12); t->Draw();
        };
        std::vector<double> cx = {0.045, 0.175, 0.280};
        int nComp = (int)tags.size() - 1;
        double span = (nComp > 0) ? std::min(0.30, (0.955 - 0.35) / nComp) : 0.30;
        double dataFont = (nComp >= 3) ? 0.0135 : 0.017;
        for (size_t tagIdx = 1; tagIdx < tags.size(); ++tagIdx) cx.push_back(0.35 + (tagIdx - 1) * span);

        Tx(0.055, 0.955, Form("Ratio vs %s reference    ( mean%% / RMS%% ,  inner #pm60#circ )",
                    refLabel.Data()), 62, 0.026);
        double y = 0.895;
        Tx(cx[0], y, "Metric", 62, 0.021); Tx(cx[1], y, "PMT", 62, 0.021); Tx(cx[2], y, "Axis", 62, 0.021);
        for (size_t tagIdx = 1; tagIdx < tags.size(); ++tagIdx) {
            TString tl = tags[tagIdx]; if (tagLabels.count(tags[tagIdx])) tl = tagLabels[tags[tagIdx]];
            Tx(cx[2 + tagIdx], y, Form("vs %s", tl.Data()), 62, 0.021);
        }
        y -= 0.012;
        { TLine* ln = new TLine(0.05, y, 0.96, y); ln->SetNDC(); ln->SetLineColor(kGray + 2); ln->Draw(); }
        y -= 0.032;
        for (auto& mm : sumMetrics) {
            for (auto& sn : allSerials) {
                for (TString axis : {"X", "Y"}) {
                    Tx(cx[0], y, mm.second.Data()); Tx(cx[1], y, sn.Data()); Tx(cx[2], y, axis.Data());
                    for (size_t tagIdx = 1; tagIdx < tags.size(); ++tagIdx) {
                        double mf, rf, mi, ri; int nf, ni;
                        bool okf = ratioStat(sn, mm.first, axis, tags[tagIdx], false, mf, rf, nf);
                        ratioStat(sn, mm.first, axis, tags[tagIdx], true, mi, ri, ni);

                        if (okf) Tx(cx[2 + tagIdx], y, Form("%+.2f/%.2f  (%+.2f/%.2f)",
                                    noNegZero((mf-1)*100), noNegZero(rf*100),
                                    noNegZero((mi-1)*100), noNegZero(ri*100)), 42, dataFont);
                        else     Tx(cx[2 + tagIdx], y, "n/a", 42, dataFont);
                    }
                    y -= 0.030;
                }
            }
            y -= 0.022;
        }
        y -= 0.010;
        Tx(0.055, y, "mean = average of (target / ref) over angles ;   RMS = point-to-point scatter.", 42, 0.017);
        y -= 0.030;
        Tx(0.055, y, "inner #pm60#circ excludes grazing-angle fit artifacts (dominant in QE / edge Gain).", 42, 0.017);
        csum->Print(pdfPath, "pdf Title:C. Ratio Summary");   // not the last page anymore -- the methodology page below closes the PDF
        if (savePNG) csum->SaveAs(Form("./Data/image/Uniformity/Overlay_RatioSummary%s%s.png", sList.Data(), dList.Data()));
        delete csum;
    }

    // ========== Final PDF page: PHC-cut vs Poisson QE methodology comparison =========={{{
    // The table above compares each wavelength against the reference wavelength
    // WITHIN one cut methodology. This page instead compares the two
    // methodologies AGAINST EACH OTHER, point-by-point at the same wavelength/
    // angle: PHC (threshold-only cut, "RelQE") vs Poisson (no cut, shape-fit of
    // the full charge distribution, "PoissonQE"). PHC is being adopted as the
    // primary QE method, so this quantifies how much that choice would move
    // the reported QE relative to the Poisson-fit alternative.}}}
    {
        auto methodStat = [&](const TString& sn, const TString& axis, const TString& tag,
                bool inner, double& mean, double& rms, int& npts)->bool {
            mean = 0; rms = 0; npts = 0;
            TFile* f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag.Data()));
            if (!f || f->IsZombie()) { if (f) f->Close(); return false; }
            TGraphErrors* gPHC = (TGraphErrors*)f->Get(Form("gr_%s_%s_RelQE", sn.Data(), axis.Data()));
            TGraphErrors* gPoi = (TGraphErrors*)f->Get(Form("gr_%s_%s_PoissonQE", sn.Data(), axis.Data()));
            std::vector<double> ratioVals;
            if (gPHC && gPoi) {
                for (int i = 0; i < gPHC->GetN(); ++i) {
                    double ang = gPHC->GetX()[i], vphc = gPHC->GetY()[i];
                    if (vphc <= 0) continue;
                    if (inner && std::abs(ang) > 60.0) continue;
                    for (int j = 0; j < gPoi->GetN(); ++j) {
                        if (std::abs(gPoi->GetX()[j] - ang) < 0.5) {
                            double vpoi = gPoi->GetY()[j];
                            if (vpoi <= 0) break;
                            ratioVals.push_back(vphc / vpoi);
                            break;
                        }
                    }
                }
            }
            f->Close();
            if (ratioVals.empty()) return false;
            for (double v : ratioVals) mean += v; mean /= ratioVals.size();
            for (double v : ratioVals) rms += (v - mean) * (v - mean); rms = std::sqrt(rms / ratioVals.size());
            npts = ratioVals.size();
            return true;
        };

        TCanvas* cmeth = new TCanvas("c_method_compare", "PHC vs Poisson QE", 1600, 2000);
        cmeth->cd();
        auto Tx2 = [&](double x, double y, const char* s, int font = 42, double sz = 0.020) {
            TLatex* t = new TLatex(x, y, s); t->SetNDC(); t->SetTextFont(font);
            t->SetTextSize(sz); t->SetTextAlign(12); t->Draw();
        };
        std::vector<double> cx = {0.055, 0.170};   
        int nTagsC = (int)tags.size();
        double span = std::min(0.32, (0.955 - 0.30) / std::max(1, nTagsC));
        double dataFont = (nTagsC >= 4) ? 0.0120 : 0.015;
        for (int tagIdx = 0; tagIdx < nTagsC; ++tagIdx) cx.push_back(0.30 + tagIdx * span);

        Tx2(0.055, 0.955, "PHC cut vs Poisson QE  ( PHC / Poisson,  mean% / RMS% ,  inner #pm60#circ )", 62, 0.026);
        Tx2(0.055, 0.912, "PHC = threshold-only cut (RelQE)   |   Poisson = no cut, full charge-shape fit (PoissonQE)", 42, 0.018);
        double y = 0.865;
        Tx2(cx[0], y, "PMT", 62, 0.021); Tx2(cx[1], y, "Axis", 62, 0.021);
        for (int tagIdx = 0; tagIdx < nTagsC; ++tagIdx) {
            TString tl = tags[tagIdx]; if (tagLabels.count(tags[tagIdx])) tl = tagLabels[tags[tagIdx]];
            Tx2(cx[2 + tagIdx], y, tl.Data(), 62, 0.021);
        }
        y -= 0.014;
        { TLine* ln = new TLine(0.05, y, 0.96, y); ln->SetNDC(); ln->SetLineColor(kGray + 2); ln->Draw(); }
        y -= 0.034;
        for (auto& sn : allSerials) {
            for (TString axis : {"X", "Y"}) {
                Tx2(cx[0], y, sn.Data()); Tx2(cx[1], y, axis.Data());
                for (int tagIdx = 0; tagIdx < nTagsC; ++tagIdx) {
                    double mf, rf, mi, ri; int nf, ni;
                    bool okf = methodStat(sn, axis, tags[tagIdx], false, mf, rf, nf);
                    methodStat(sn, axis, tags[tagIdx], true, mi, ri, ni);
                    if (okf) Tx2(cx[2 + tagIdx], y, Form("%+.2f/%.2f  (%+.2f/%.2f)",
                                noNegZero((mf-1)*100), noNegZero(rf*100),
                                noNegZero((mi-1)*100), noNegZero(ri*100)), 42, dataFont);
                    else     Tx2(cx[2 + tagIdx], y, "n/a", 42, dataFont);
                }
                y -= 0.032;
            }
        }
        y -= 0.020;
        Tx2(0.055, y, "mean = average of (PHC / Poisson) over angles ;   RMS = point-to-point scatter.", 42, 0.017);
        y -= 0.030;
        Tx2(0.055, y, "A mean near 0% means the two cut methodologies agree on absolute QE level.", 42, 0.017);
        cmeth->Print(pdfPath + ")", "pdf Title:D. PHC vs Poisson QE");   // final page, closes the PDF
        if (savePNG) cmeth->SaveAs(Form("./Data/image/Uniformity/Overlay_MethodCompare%s%s.png", sList.Data(), dList.Data()));
        delete cmeth;
    }

    std::cout << "[INFO] Analysis Complete. Overlay graphs (Fixed Top & Bottom Overlap) saved." << std::endl;
}
