// Standalone: Hamamatsu-style Cathode-Uniformity shape classification, applied
// to our QE (CorrRelQE) metric. Center/Side/sigma definitions match the
// R12860-22 slide convention:
//   Center = mean(value) for |angle| < 20 deg
//   Side   = mean(value) for 30 < |angle| < 70 deg (both sides pooled)
//   sigma  = population std(value) for |angle| < 70 deg, x1000
//   C/S Ratio = Center / Side
// Our scan only covers |angle| <= 55 deg, so the "Side"/"sigma" windows above
// 55 deg simply have no points to average -- not a bug, just a narrower scan
// range than Hamamatsu's bench (which goes to ~90 deg).
// No shape-classification thresholds are hardcoded yet (C/S Ratio and sigma
// cut values are still TBD) -- this only computes and plots the raw
// (C/S Ratio, sigma) points so those thresholds can be chosen by eye first.

void Draw_QE_UniformityShape(std::vector<TString> tags, std::vector<TString> serialList = {"EM5370", "EL9590"},
                              TString metric = "CorrRelQE") {
    // Classification is a per-PMT-UNIT property (cathode shape of that
    // physical tube), so tags mixed into one call must all use the SAME
    // physical serials in serialList -- comparing EM5370 across dates is
    // "how did this one tube's shape hold up over time", not a cross-unit
    // comparison. To compare different tubes (e.g. EL1635 vs EM5370, both
    // slot-2 units at different times), call this twice with each one's own
    // tags+serial and compare the two printed tables / overlay the two PDFs.
    // EM2740 (Monitor PMT) is excluded by default -- it doesn't rotate, so
    // "uniformity vs angle" isn't a physically meaningful shape for it the
    // way it is for a rotating PMT.
    int nSerials = serialList.size();
    std::vector<const char*> serials; for (auto& s : serialList) serials.push_back(s.Data());
    int colorPalette[6] = {kRed+1, kBlue+1, kGreen+2, kOrange+1, kMagenta+1, kCyan+2};
    int markerPalette[6] = {21, 22, 20, 23, 33, 34};
    const char* axes[2] = {"X", "Y"};

    const double centerMax = 20.0;
    const double sideMin = 30.0, sideMax = 70.0;

    struct Point { TString tag, serial, axis; double center, side, csRatio, sigma; int n; };
    std::vector<Point> results;

    for (auto& tag : tags) {
        TString path = Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag.Data());
        TFile* f = TFile::Open(path);
        if (!f || f->IsZombie()) { std::cout << "[WARN] missing: " << path << std::endl; continue; }

        for (int s = 0; s < nSerials; ++s) {
            for (int a = 0; a < 2; ++a) {
                TGraphErrors* gr = (TGraphErrors*)f->Get(Form("gr_%s_%s_%s", serials[s], axes[a], metric.Data()));
                if (!gr || gr->GetN() == 0) continue;

                std::vector<double> centerVals, sideVals, sigmaVals;
                for (int i = 0; i < gr->GetN(); ++i) {
                    double ang = gr->GetX()[i], val = gr->GetY()[i];
                    if (val <= 0) continue;   // sentinel/invalid guard, same convention as elsewhere
                    double absAng = std::abs(ang);
                    if (absAng < centerMax) centerVals.push_back(val);
                    if (absAng > sideMin && absAng < sideMax) sideVals.push_back(val);
                    if (absAng < sideMax) sigmaVals.push_back(val);
                }
                if (centerVals.empty() || sideVals.empty() || sigmaVals.size() < 2) continue;

                double center = 0; for (double v : centerVals) center += v; center /= centerVals.size();
                double side = 0; for (double v : sideVals) side += v; side /= sideVals.size();
                double sigMean = 0; for (double v : sigmaVals) sigMean += v; sigMean /= sigmaVals.size();
                double sigma = 0; for (double v : sigmaVals) sigma += (v - sigMean) * (v - sigMean);
                sigma = std::sqrt(sigma / sigmaVals.size()) * 1000.0;

                results.push_back({tag, serials[s], axes[a], center, side, center / side, sigma, (int)sigmaVals.size()});
            }
        }
        f->Close();
    }

    if (results.empty()) { std::cout << "[ERROR] No usable points found." << std::endl; return; }

    std::cout << Form("%-22s %-8s %-4s %8s %8s %10s %10s %5s", "Tag", "Serial", "Axis", "Center", "Side", "C/S Ratio", "sigma*1000", "N") << std::endl;
    for (auto& r : results) {
        std::cout << Form("%-22s %-8s %-4s %8.4f %8.4f %10.4f %10.2f %5d",
                           r.tag.Data(), r.serial.Data(), r.axis.Data(), r.center, r.side, r.csRatio, r.sigma, r.n) << std::endl;
    }

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c_qe_shape", "QE Uniformity Shape (C/S Ratio vs sigma)", 1200, 900);
    c->SetGrid(); c->SetLeftMargin(0.13); c->SetBottomMargin(0.12); c->SetTopMargin(0.08); c->SetRightMargin(0.04);

    double csMin = 1e9, csMax = -1e9, sigMax = -1e9;
    for (auto& r : results) {
        csMin = std::min(csMin, r.csRatio); csMax = std::max(csMax, r.csRatio);
        sigMax = std::max(sigMax, r.sigma);
    }
    double csPad = std::max(0.02, (csMax - csMin) * 0.15);

    TH1F* frame = c->DrawFrame(csMin - csPad, 0, csMax + csPad, sigMax * 1.2);
    frame->SetTitle(Form("QE (%s) Uniformity Shape;Center / Side Ratio;#sigma #times 1000", metric.Data()));
    frame->GetXaxis()->SetTitleSize(0.045); frame->GetXaxis()->SetLabelSize(0.035); frame->GetXaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetLabelSize(0.035); frame->GetYaxis()->SetTitleOffset(1.3);

    TLegend* leg = new TLegend(0.70, 0.78, 0.94, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(132); leg->SetTextSize(0.03);
    std::vector<bool> shown(nSerials, false);

    for (auto& r : results) {
        int si = -1;
        for (int s = 0; s < nSerials; ++s) if (r.serial == serials[s]) si = s;
        if (si < 0) continue;
        TGraph* pt = new TGraph(1);
        pt->SetPoint(0, r.csRatio, r.sigma);
        pt->SetMarkerStyle(markerPalette[si % 6]); pt->SetMarkerColor(colorPalette[si % 6]); pt->SetMarkerSize(1.4);
        pt->Draw("P SAME");
        if (!shown[si]) { leg->AddEntry(pt, serials[si], "p"); shown[si] = true; }
    }
    leg->Draw();

    TString tagList = ""; for (auto& t : tags) tagList += "_" + t;
    TString outPath = Form("./Data/image/Uniformity/QE_UniformityShape%s.pdf", tagList.Data());
    c->Print(outPath);
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
