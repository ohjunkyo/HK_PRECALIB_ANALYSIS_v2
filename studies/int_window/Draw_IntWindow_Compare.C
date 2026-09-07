// Usage: root -l -b -q Draw_IntWindow_Compare.C
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

struct Row {
    TString repTag; int run; int ch; TString pmtSN; TString axis;
    double tilt, rot; int intWindowNs;
    double qeFinal, qeRaw, tts, gain, chargeRes;
};

void Draw_IntWindow_Compare(const char* csvPath = "/home/precalkor/ADC/ADC_test/Data/summary/intwin_qe_tts_gain.csv",
                             const char* outDir = "/home/precalkor/ADC/ADC_test/Data/FinalResult/IntWindowStudy") {
    std::ifstream in(csvPath);
    if (!in.is_open()) { std::cerr << "[ERROR] cannot open " << csvPath << std::endl; return; }
    std::string line;
    std::getline(in, line);   // header
    std::vector<Row> rows;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string tok;
        Row r;
        std::getline(ss, tok, ','); r.repTag = tok;
        std::getline(ss, tok, ','); r.run = atoi(tok.c_str());
        std::getline(ss, tok, ','); r.ch = atoi(tok.c_str());
        std::getline(ss, tok, ','); r.pmtSN = tok;
        std::getline(ss, tok, ','); r.axis = tok;
        std::getline(ss, tok, ','); r.tilt = atof(tok.c_str());
        std::getline(ss, tok, ','); r.rot = atof(tok.c_str());
        std::getline(ss, tok, ','); r.intWindowNs = atoi(tok.c_str());
        std::getline(ss, tok, ','); r.qeFinal = atof(tok.c_str());
        std::getline(ss, tok, ','); r.qeRaw = atof(tok.c_str());
        std::getline(ss, tok, ','); r.tts = atof(tok.c_str());
        std::getline(ss, tok, ','); r.gain = atof(tok.c_str());
        std::getline(ss, tok, ','); r.chargeRes = atof(tok.c_str());
        if (r.ch == 0) continue;   // Monitor is fixed, not angle-scanned
        rows.push_back(r);
    }
    std::cout << "[INFO] " << rows.size() << " rows loaded (Monitor channel excluded)." << std::endl;

    gSystem->mkdir(outDir, kTRUE);

    struct MetricDef { double Row::*field; const char* label; const char* fname; };
    std::vector<MetricDef> metrics = {
        {&Row::qeFinal, "QE (dark-subtracted) [%]", "QE"},
        {&Row::tts,     "TTS [ns]",                  "TTS"},
        {&Row::gain,    "Gain (SPE mean) [pC]",       "Gain"},
    };
    std::vector<int> widths = {100, 200, 300};
    int widthColor[3] = {kBlue + 2, kGreen + 2, kRed + 1};

    TString pdfPath = Form("%s/QE_TTS_Gain_Compare.pdf", outDir);
    TCanvas c("c", "c", 1000, 700);
    bool first = true;

    // repTags found in the CSV itself, not hardcoded Rep1/Rep2 -- any tag
    // string Run_IntWindowStudy.C was called with shows up here too.
    std::vector<TString> repTags;
    for (auto& r : rows)
        if (std::find(repTags.begin(), repTags.end(), r.repTag) == repTags.end()) repTags.push_back(r.repTag);
    std::sort(repTags.begin(), repTags.end());

    for (auto& m : metrics) {

        for (TString repTag : repTags) {
            // Group by (pmtSN) so Rot#1 / Rot#2 each get their own page.
            std::vector<TString> pmts;
            for (auto& r : rows) if (r.repTag == repTag &&
                    std::find(pmts.begin(), pmts.end(), r.pmtSN) == pmts.end()) pmts.push_back(r.pmtSN);
            std::sort(pmts.begin(), pmts.end());

            for (auto& pmt : pmts) {
                c.cd(); c.Clear();
                TLegend leg(0.72, 0.75, 0.93, 0.90);
                leg.SetBorderSize(1); leg.SetFillColor(kWhite);
                bool drewAny = false;
                double yMin = 1e18, yMax = -1e18;

                std::vector<TGraph*> graphs;
                for (size_t wi = 0; wi < widths.size(); ++wi) {
                    std::vector<double> ang, val;
                    for (auto& r : rows) {
                        if (r.repTag != repTag || r.pmtSN != pmt || r.intWindowNs != widths[wi]) continue;
                        ang.push_back(r.tilt);
                        val.push_back(r.*(m.field));
                    }
                    if (ang.empty()) continue;
                    // sort by angle
                    std::vector<int> idx(ang.size()); for (size_t k = 0; k < idx.size(); ++k) idx[k] = k;
                    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return ang[a] < ang[b]; });
                    std::vector<double> angS, valS;
                    for (int k : idx) { angS.push_back(ang[k]); valS.push_back(val[k]); }
                    for (double v : valS) { yMin = std::min(yMin, v); yMax = std::max(yMax, v); }

                    TGraph* g = new TGraph((int)angS.size(), &angS[0], &valS[0]);
                    g->SetLineColor(widthColor[wi]); g->SetMarkerColor(widthColor[wi]);
                    g->SetMarkerStyle(20 + wi); g->SetLineWidth(2); g->SetMarkerSize(0.9);
                    graphs.push_back(g);
                    leg.AddEntry(g, Form("%d ns", widths[wi]), "lp");
                    drewAny = true;
                }
                if (!drewAny) continue;

                double pad = (yMax - yMin) * 0.15 + 1e-6;
                TGraph* frame = graphs[0];
                frame->SetTitle(Form("%s   %s   (%s);Position angle [deg];%s", pmt.Data(), repTag.Data(), m.fname, m.label));
                frame->GetYaxis()->SetRangeUser(yMin - pad, yMax + pad);
                frame->Draw("ALP");
                for (size_t k = 1; k < graphs.size(); ++k) graphs[k]->Draw("LP SAME");
                leg.Draw();
                gPad->SetGrid();

                c.Print(pdfPath + (first ? "(" : ""), "pdf");
                first = false;
            }
        }
    }
    if (!first) {
        c.Print(pdfPath + ")", "pdf");
        std::cout << "  -> " << pdfPath << std::endl;
    }
}
