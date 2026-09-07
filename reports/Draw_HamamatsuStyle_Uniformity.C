// Recreates the Hamamatsu "R12860-22 Cathode Uniformity" slide style for our
// own uniformity-scan data: Normalized value [a.u.] vs Position angle, X-axis
// & Y-axis overlaid (with error bars), 0.9/1.1 reference lines, and a
// Center/Side/C-S-Ratio/sigma summary box -- one page per PMT serial, all
// pages in a single PDF.
//
// Default metric is "MonNorm_CorrRawQE" (Monitor(CH0)-normalized, dark-
// subtracted QE: Test_QEcorr / Monitor_QEcorr at the matching raw tilt),
// plotted on its native scale, Center NOT forced to 1.0 (centerNormalize
// defaults to false). Hamamatsu's own slide doesn't force Center=100% either
// (their 3 examples read 108%/99%/95%) -- their Center carries real
// information (this unit vs. a reference), and so does ours (this PMT's QE
// vs. the Monitor PMT's QE). Forcing Center=1 would erase exactly that
// absolute-vs-Monitor information and leave only the angular shape; set
// centerNormalize=true only when doing a deliberate shape-only comparison
// (e.g. across repeats/dates where each side's own absolute level differs
// for unrelated reasons). Either way C/S Ratio is unchanged -- it's a ratio,
// so uniform rescaling cancels; only sigma*1000 and the plotted curve height
// depend on this choice.
//
// Center/Side/sigma windows and the C/S-Ratio definition match the R12860-22
// slide (2026/7/29):
//   Center = mean(value) for |angle| < 20 deg
//   Side   = mean(value) for 30 < |angle| < 70 deg (both sides pooled)
//   sigma  = population std(value) for |angle| < 70 deg, x1000 (fractional scale)
//   C/S Ratio = Center / Side
// "angle" here is the graph's plotted X value, already the raw stage tilt run
// through the KR->Hamamatsu incidence conversion (angle_convert.h,
// ConvertKRtoHamamatsu()), NOT the raw stage angle. Our raw stage scan goes to
// |tilt|=55deg, which maps to |angle| ~= 81.5deg -- close to (slightly short
// of) Hamamatsu's own ~90deg bench range. Points beyond that simply don't
// contribute to Side/sigma -- not a bug, just the residual coverage gap.

// separateAxes: false (default) keeps X and Y overlaid on one page, as the
//   Hamamatsu slide does. true gives each axis its own page instead, which is
//   what you want when the two axes differ enough that overlaying them hides
//   the shape of each (2026-08-25).
// A per-PMT summary row is always appended to
//   ./Data/summary/pmt_uniformity_summary.csv
// so results accumulate across PMTs/tags into one table instead of living
// only inside a PDF -- see the header written below for the columns.
void Draw_HamamatsuStyle_Uniformity(std::vector<TString> tags, std::vector<TString> serialList = {"EM5370", "EL9590"},
                                     TString metric = "MonNorm_CorrRawQE", bool centerNormalize = false,
                                     bool separateAxes = false) {
    int nSerials = serialList.size();
    const char* axes[2] = {"X", "Y"};
    int axisColor[2] = {kBlue+2, kRed+1};
    int axisMarker[2] = {21, 20};

    const double centerMax = 20.0;
    const double sideMin = 30.0, sideMax = 70.0;

    // Human-readable label for the info box. The internal graph names are
    // confusing on purpose-by-accident: "Raw" there means "not angle-normalized"
    // (i.e. still on the absolute scale), NOT "not dark-subtracted" -- the
    // dark-subtraction state is spelled by the presence/absence of "Corr".
    // Split over two lines so it fits the narrow info box.
    TString metricLabel = metric, metricLabel2 = "";
    if (metric == "MonNorm_CorrRawQE")      { metricLabel = "Monitor-Normalized QE"; metricLabel2 = "(dark-corrected)"; }
    else if (metric == "MonNorm_RawQE")     { metricLabel = "Monitor-Normalized QE"; metricLabel2 = "(no dark corr.)"; }
    else if (metric == "CorrRelQE")         { metricLabel = "Angle-Normalized QE";   metricLabel2 = "(dark-corrected)"; }
    else if (metric == "RelQE")             { metricLabel = "Angle-Normalized QE";   metricLabel2 = "(no dark corr.)"; }
    if (centerNormalize) metricLabel2 += " Center=1";

    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(22, "");     // pad title: Times New Roman Bold
    gStyle->SetTitleSize(0.045, "");
    gStyle->SetTitleFont(132, "X");   // axis titles/labels: Times New Roman (regular)
    gStyle->SetTitleFont(132, "Y");
    gStyle->SetLabelFont(132, "X");
    gStyle->SetLabelFont(132, "Y");

    TString tagList = ""; for (auto& t : tags) tagList += "_" + t;
    TString outPath = Form("./Data/image/Uniformity/HamamatsuStyle_Uniformity%s.pdf", tagList.Data());

    TCanvas* c = new TCanvas("c_hama", "Hamamatsu-style Cathode Uniformity", 1400, 850);
    c->Print(outPath + "[");

    for (int s = 0; s < nSerials; ++s) {
        TString serial = serialList[s];

        std::vector<double> ang[2], val[2], err[2];   // [0]=X-axis, [1]=Y-axis, val/err in fractional units (Center~1.0)
        for (auto& tag : tags) {
            TString path = Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tag.Data());
            TFile* f = TFile::Open(path);
            if (!f || f->IsZombie()) { std::cout << "[WARN] missing: " << path << std::endl; continue; }
            for (int a = 0; a < 2; ++a) {
                TGraphErrors* gr = (TGraphErrors*)f->Get(Form("gr_%s_%s_%s", serial.Data(), axes[a], metric.Data()));
                if (!gr) continue;
                for (int i = 0; i < gr->GetN(); ++i) {
                    double an = gr->GetX()[i], v = gr->GetY()[i], e = gr->GetEY()[i];
                    if (v <= 0) continue;   // sentinel/invalid guard
                    ang[a].push_back(an); val[a].push_back(v); err[a].push_back(e);
                }
            }
            f->Close();
        }
        if (ang[0].empty() && ang[1].empty()) { std::cout << "[WARN] no data for " << serial << std::endl; continue; }

        // Re-normalize so this PMT's own Center(|angle|<20deg) average == 1.0,
        // pooling both axes -- matches Hamamatsu's own plotting convention.
        if (centerNormalize) {
            std::vector<double> preCenter;
            for (int a = 0; a < 2; ++a)
                for (size_t i = 0; i < ang[a].size(); ++i)
                    if (std::abs(ang[a][i]) < centerMax) preCenter.push_back(val[a][i]);
            if (!preCenter.empty()) {
                double centerAvg = 0; for (double v : preCenter) centerAvg += v; centerAvg /= preCenter.size();
                if (centerAvg != 0) {
                    for (int a = 0; a < 2; ++a) {
                        for (size_t i = 0; i < val[a].size(); ++i) { val[a][i] /= centerAvg; err[a][i] /= centerAvg; }
                    }
                }
            }
        }

        // Center/Side/sigma -- computed per axis AND pooled. The pooled numbers
        // follow the Hamamatsu convention (their bench sweeps one diameter at a
        // time but quotes a single figure per tube); the per-axis split is shown
        // alongside because X and Y can genuinely differ (e.g. under a residual
        // field the two scan planes see different Lorentz deflection).
        auto shapeStats = [&](const std::vector<int>& useAxes,
                               double& csRatio, double& sigma) -> bool {
            std::vector<double> centerVals, sideVals, sigmaVals;
            for (int a : useAxes) {
                for (size_t i = 0; i < ang[a].size(); ++i) {
                    double absAng = std::abs(ang[a][i]), v = val[a][i];
                    if (absAng < centerMax) centerVals.push_back(v);
                    if (absAng > sideMin && absAng < sideMax) sideVals.push_back(v);
                    if (absAng < sideMax) sigmaVals.push_back(v);
                }
            }
            if (centerVals.empty() || sideVals.empty() || sigmaVals.size() < 2) return false;
            double center = 0, side = 0;
            for (double v : centerVals) center += v; center /= centerVals.size();
            for (double v : sideVals) side += v; side /= sideVals.size();
            csRatio = (side != 0) ? center / side : 0;
            double m = 0; for (double v : sigmaVals) m += v; m /= sigmaVals.size();
            double s2 = 0; for (double v : sigmaVals) s2 += (v - m) * (v - m);
            sigma = std::sqrt(s2 / sigmaVals.size()) * 1000.0;
            return true;
        };

        double csRatio = 0, sigma = 0;
        double csX = 0, sigX = 0, csY = 0, sigY = 0;
        bool haveStats = shapeStats({0, 1}, csRatio, sigma);
        bool haveX = shapeStats({0}, csX, sigX);
        bool haveY = shapeStats({1}, csY, sigY);

        // ---- Draw page(s) ----
        // separateAxes=false -> one page with both axes; true -> one page per
        // axis. `pages` is the list of axis sets to put on each page, so the
        // body below is written once either way.
        std::vector<std::vector<int>> pages;
        if (separateAxes) { pages.push_back({0}); pages.push_back({1}); }
        else              { pages.push_back({0, 1}); }

        for (auto& pageAxes : pages) {
            bool pageHasData = false;
            for (int a : pageAxes) if (!ang[a].empty()) pageHasData = true;
            if (!pageHasData) continue;

            c->cd(); c->Clear();
            c->SetGrid(); c->SetLeftMargin(0.10); c->SetRightMargin(0.28);
            c->SetTopMargin(0.10); c->SetBottomMargin(0.12);

            double xMin = 1e9, xMax = -1e9;
            for (int a : pageAxes)
                for (double an : ang[a]) { xMin = std::min(xMin, an); xMax = std::max(xMax, an); }
            double xPad = std::max(5.0, (xMax - xMin) * 0.08);

            TString axisTitle = "";
            if (pageAxes.size() == 1) axisTitle = Form("  (%s-Axis)", axes[pageAxes[0]]);

            TH1F* frame = c->DrawFrame(xMin - xPad, 0, xMax + xPad, 1.5);
            frame->SetTitle(Form("%s Cathode Uniformity%s;Position angle [degree];Relative QE [a.u.]",
                                 serial.Data(), axisTitle.Data()));
            frame->GetXaxis()->SetTitleSize(0.040); frame->GetXaxis()->SetLabelSize(0.032); frame->GetXaxis()->SetTitleOffset(1.2);
            frame->GetYaxis()->SetTitleSize(0.040); frame->GetYaxis()->SetLabelSize(0.032); frame->GetYaxis()->SetTitleOffset(1.1);
            frame->GetYaxis()->SetNdivisions(515);

            TLine* l110 = new TLine(xMin - xPad, 1.1, xMax + xPad, 1.1);
            l110->SetLineColor(kAzure+7); l110->SetLineWidth(2); l110->Draw();
            TLine* l90 = new TLine(xMin - xPad, 0.9, xMax + xPad, 0.9);
            l90->SetLineColor(kOrange+1); l90->SetLineWidth(2); l90->Draw();

            TLegend* leg = new TLegend(0.735, 0.80, 0.985, 0.90);
            leg->SetBorderSize(1); leg->SetTextFont(132); leg->SetTextSize(0.028);

            for (int a : pageAxes) {
                if (ang[a].empty()) continue;
                std::vector<size_t> idx(ang[a].size());
                for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
                std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j){ return ang[a][i] < ang[a][j]; });

                TGraphErrors* gr = new TGraphErrors(ang[a].size());
                for (size_t k = 0; k < idx.size(); ++k)
                    gr->SetPoint(k, ang[a][idx[k]], val[a][idx[k]]), gr->SetPointError(k, 0, err[a][idx[k]]);
                gr->SetMarkerStyle(axisMarker[a]); gr->SetMarkerColor(axisColor[a]); gr->SetMarkerSize(1.1);
                gr->SetLineColor(axisColor[a]); gr->SetLineWidth(2);
                gr->Draw("PL SAME");
                leg->AddEntry(gr, Form("%s-Axis", axes[a]), "lp");
            }
            leg->Draw();

            // Info box (top-right, outside plot area)
            TPaveText* info = new TPaveText(0.735, 0.55, 0.985, 0.78, "NDC");
            info->SetBorderSize(1); info->SetFillColor(kWhite); info->SetFillStyle(1001); info->SetTextFont(132); info->SetTextSize(0.026); info->SetTextAlign(12);
            info->AddText(Form("Serial No. : %s", serial.Data()));
            info->AddText(Form("Metric : %s", metricLabel.Data()));
            if (metricLabel2 != "") info->AddText(Form("    %s", metricLabel2.Data()));
            TString tagsShown = ""; for (size_t i = 0; i < tags.size(); ++i) tagsShown += (i? ", " : "") + tags[i];
            info->AddText(Form("Tag(s) : %s", tagsShown.Data()));
            info->Draw();

            // Stats box. Collect the lines FIRST, then size the box to them --
            // it used to be a fixed 0.13..0.37 NDC rectangle regardless of
            // content, so a single-axis page (one line) got a box sized for
            // three and that oversized box reached far enough left to cover
            // the curve's descending edge tail (2026-08-25).
            std::vector<TString> statLines;
            if (pageAxes.size() == 1) {
                int a = pageAxes[0];
                bool ok = a ? haveY : haveX;
                double cs = a ? csY : csX, sg = a ? sigY : sigX;
                if (ok) statLines.push_back(Form("%s  %s-Axis :  C/S = %.3f,  #sigma#times1000 = %.2f",
                                                 serial.Data(), axes[a], cs, sg));
                else    statLines.push_back(Form("%s  %s-Axis : insufficient points", serial.Data(), axes[a]));
            } else if (haveStats) {
                statLines.push_back(Form("%s   C/S = %.3f,  #sigma#times1000 = %.2f",
                                         serial.Data(), csRatio, sigma));
                if (haveX) statLines.push_back(Form("   X-Axis :  C/S = %.3f,  #sigma#times1000 = %.2f", csX, sigX));
                if (haveY) statLines.push_back(Form("   Y-Axis :  C/S = %.3f,  #sigma#times1000 = %.2f", csY, sigY));
            } else {
                statLines.push_back(Form("%s", serial.Data()));
                statLines.push_back("Insufficient points for Center/Side stats");
            }

            // Bottom-centre: the uniformity curve is high across the middle and
            // only dips at |angle| > ~60, so the centre of the lower half is the
            // one region that stays clear on every PMT. Height follows the line
            // count; width follows the longest line.
            const double lineH = 0.045;
            size_t maxLen = 0;
            for (auto& s2 : statLines) maxLen = std::max(maxLen, (size_t)s2.Length());
            double boxW = std::min(0.56, 0.055 + maxLen * 0.0092);
            double boxH = statLines.size() * lineH;
            double xC = 0.50;                       // centred within the pad area
            double x1b = std::max(0.12, xC - boxW / 2), x2b = std::min(0.70, xC + boxW / 2);
            double y1b = 0.14, y2b = y1b + boxH;

            TPaveText* stats = new TPaveText(x1b, y1b, x2b, y2b, "NDC");
            stats->SetBorderSize(1); stats->SetFillColor(kWhite); stats->SetFillStyle(1001);
            stats->SetTextFont(132); stats->SetTextSize(0.026); stats->SetTextAlign(12);
            stats->SetMargin(0.03);
            for (auto& s2 : statLines) stats->AddText(s2);
            stats->Draw();

            c->Print(outPath);
        }

        std::cout << Form("[INFO] %-8s  pooled C/S=%.3f sig=%.2f | X: C/S=%.3f sig=%.2f | Y: C/S=%.3f sig=%.2f",
                           serial.Data(), csRatio, sigma, csX, sigX, csY, sigY) << std::endl;

        // ---- Persist this PMT's numbers ----
        // A PDF page is fine for looking at one PMT, but useless for "compare
        // every PMT we have measured". Written in LONG format -- one row per
        // (PMT, tag, metric, AXIS) rather than one wide row with cs_X/cs_Y/
        // cs_pooled columns. Long format is what makes "store together, split
        // X vs Y later" a filter (axis=="X") instead of a schema change, and
        // it is also the shape a TTree wants.
        {
            gSystem->mkdir("./Data/summary", kTRUE);
            TString tagsJoined = "";
            for (size_t i = 0; i < tags.size(); ++i) tagsJoined += (i ? "+" : "") + tags[i];
            TDatime now;

            // X and Y only. A "pooled" row (both axes merged) used to be stored
            // alongside them, but X and Y are always wanted separately -- the
            // pooled number just duplicated information already in the two axis
            // rows and added a third row per PMT to the GUI's PMT Info table
            // that nothing acted on (2026-08-26). The pooled C/S and sigma are
            // still printed in the console [INFO] line above if needed ad hoc.
            struct AxRow { const char* axis; bool ok; double cs, sig; int npts; };
            AxRow axRows[2] = {
                {"X",      haveX,     csX,     sigX,  (int)ang[0].size()},
                {"Y",      haveY,     csY,     sigY,  (int)ang[1].size()},
            };

            TString csvPath = "./Data/summary/pmt_uniformity_summary.csv";
            bool needHeader = gSystem->AccessPathName(csvPath.Data());
            std::ofstream csv(csvPath.Data(), std::ios::app);
            if (csv.is_open()) {
                if (needHeader)
                    csv << "timestamp,serial,tag,metric,axis,center_normalized,"
                           "n_points,cs_ratio,sigma1000\n";
                for (auto& r : axRows) {
                    if (!r.ok) continue;
                    csv << now.AsSQLString() << "," << serial << "," << tagsJoined << ","
                        << metric << "," << r.axis << "," << (centerNormalize ? 1 : 0) << ","
                        << r.npts << "," << Form("%.4f", r.cs) << "," << Form("%.4f", r.sig) << "\n";
                }
                csv.close();
            }

            // Same rows into a TTree. Opened UPDATE so entries accumulate
            // across runs of this macro; the tree is read back by the GUI's
            // PMT Info tab and by any ROOT session wanting
            // t->Draw("cs_ratio","axis==\"X\"").
            TString rootPath = "./Data/summary/pmt_uniformity_summary.root";
            TFile* fs = TFile::Open(rootPath, "UPDATE");
            if (fs && !fs->IsZombie()) {
                TTree* t = (TTree*)fs->Get("pmt_summary");
                char b_serial[32], b_tag[128], b_metric[64], b_axis[16], b_time[32];
                int   b_npts = 0, b_cnorm = 0;
                double b_cs = 0, b_sig = 0;
                if (!t) {
                    t = new TTree("pmt_summary", "Per-PMT cathode uniformity summary");
                    t->Branch("timestamp", b_time,   "timestamp/C");
                    t->Branch("serial",    b_serial, "serial/C");
                    t->Branch("tag",       b_tag,    "tag/C");
                    t->Branch("metric",    b_metric, "metric/C");
                    t->Branch("axis",      b_axis,   "axis/C");
                    t->Branch("center_normalized", &b_cnorm, "center_normalized/I");
                    t->Branch("n_points",  &b_npts,  "n_points/I");
                    t->Branch("cs_ratio",  &b_cs,    "cs_ratio/D");
                    t->Branch("sigma1000", &b_sig,   "sigma1000/D");
                } else {
                    t->SetBranchAddress("timestamp", b_time);
                    t->SetBranchAddress("serial",    b_serial);
                    t->SetBranchAddress("tag",       b_tag);
                    t->SetBranchAddress("metric",    b_metric);
                    t->SetBranchAddress("axis",      b_axis);
                    t->SetBranchAddress("center_normalized", &b_cnorm);
                    t->SetBranchAddress("n_points",  &b_npts);
                    t->SetBranchAddress("cs_ratio",  &b_cs);
                    t->SetBranchAddress("sigma1000", &b_sig);
                }
                for (auto& r : axRows) {
                    if (!r.ok) continue;
                    snprintf(b_time,   sizeof(b_time),   "%s", now.AsSQLString());
                    snprintf(b_serial, sizeof(b_serial), "%s", serial.Data());
                    snprintf(b_tag,    sizeof(b_tag),    "%s", tagsJoined.Data());
                    snprintf(b_metric, sizeof(b_metric), "%s", metric.Data());
                    snprintf(b_axis,   sizeof(b_axis),   "%s", r.axis);
                    b_cnorm = centerNormalize ? 1 : 0;
                    b_npts = r.npts; b_cs = r.cs; b_sig = r.sig;
                    t->Fill();
                }
                fs->cd();
                t->Write("pmt_summary", TObject::kOverwrite);
                fs->Close();
            }
            std::cout << "[INFO]   -> appended summary rows to " << csvPath
                      << " and " << rootPath << std::endl;
        }
    }

    c->Print(outPath + "]");
    std::cout << "[INFO] Saved: " << outPath << std::endl;
}
