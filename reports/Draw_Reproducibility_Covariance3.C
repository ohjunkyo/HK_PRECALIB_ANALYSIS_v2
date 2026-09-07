// 3-repeat covariance / chi2_r analysis (steps 3-4 of the reproducibility
// framework) -- needs >=3 repeats to estimate an empirical point-to-point
// covariance matrix, which Rep.1/Rep.2-only data (2 repeats) cannot give.
// 20260713_100_145 / _200_245 / _300_345 are three full repeat scans of the
// same PMTs/angles on the same day, so covariance IS estimable here.
//
// Method
// ------
// For a given axis/PMT, let Y_r,i be the value at angle i in repeat r
// (r=1..3). The point-to-point covariance matrix over the 3 repeats:
//
//   C_ij = (1/(R-1)) * sum_r (Y_r,i - Ybar_i)(Y_r,j - Ybar_j)
//
// This mixes real per-point stat noise (diagonal-dominated) with any
// common-mode drift shared across nearby angles (off-diagonal, e.g. laser
// intensity drift shifting the whole curve together within a repeat). Then
// for each repeat r, the deviation vector d_r = Y_r - Ybar is tested against
// the full (all-repeat) covariance via
//
//   chi2_r = d_r^T C^-1 d_r,  dof = N_points
//
// chi2_r/dof ~ 1 means that repeat's curve-to-curve deviation is fully
// consistent with the same covariance structure seen across all repeats
// (i.e., no single repeat is an outlier beyond what the ensemble already
// shows). This is a genuinely different question from the point-by-point
// pull test in Draw_Reproducibility.C: it tests SELF-CONSISTENCY of the
// whole curve as a correlated object, not just point-by-point agreement.
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

void Draw_Reproducibility_Covariance3(
        std::vector<TString> tags = {"20260713_100_145", "20260713_200_245", "20260713_300_345"},
        std::vector<TString> serialList = {"EM5370", "EL9590"},
        TString metric = "CorrRelQE", double angleTol = 1.0) {

    const char* axes[2] = {"X", "Y"};
    int R = (int)tags.size();
    if (R < 3) { std::cout << "[ERROR] need >=3 repeat tags for a covariance estimate." << std::endl; return; }

    std::vector<TFile*> files;
    for (auto& t : tags) {
        TFile* f = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", t.Data()));
        if (!f || f->IsZombie()) { std::cout << "[ERROR] cannot open tag " << t << std::endl; return; }
        files.push_back(f);
    }

    std::cout << "\n================  3-repeat covariance / chi2_r  ================" << std::endl;
    std::cout << Form("  Repeats: %s, %s, %s   metric=%s", tags[0].Data(), tags[1].Data(), tags[2].Data(), metric.Data()) << std::endl;

    for (int a_idx = 0; a_idx < 2; ++a_idx) {
        for (auto& serial : serialList) {
            // Load the 3 graphs for this axis/PMT.
            std::vector<TGraphErrors*> g(R);
            bool ok = true;
            for (int r = 0; r < R; ++r) {
                g[r] = (TGraphErrors*)files[r]->Get(Form("gr_%s_%s_%s", serial.Data(), axes[a_idx], metric.Data()));
                if (!g[r]) { ok = false; break; }
            }
            if (!ok) { std::cout << Form("[WARN] missing graph for %s %s-axis", serial.Data(), axes[a_idx]) << std::endl; continue; }

            // Match angles common to ALL 3 repeats (within angleTol), using
            // repeat 0's angles as the reference grid.
            std::vector<double> angles;
            std::vector<std::vector<double>> Y(R);   // Y[r][i]
            for (int i = 0; i < g[0]->GetN(); ++i) {
                double an0 = g[0]->GetX()[i], v0 = g[0]->GetY()[i];
                if (v0 <= 0) continue;
                std::vector<double> vals(R, 0.0);
                vals[0] = v0;
                bool matched = true;
                for (int r = 1; r < R; ++r) {
                    bool found = false;
                    for (int j = 0; j < g[r]->GetN(); ++j) {
                        double anr = g[r]->GetX()[j], vr = g[r]->GetY()[j];
                        if (vr <= 0) continue;
                        if (std::abs(an0 - anr) <= angleTol) { vals[r] = vr; found = true; break; }
                    }
                    if (!found) { matched = false; break; }
                }
                if (!matched) continue;
                angles.push_back(an0);
                for (int r = 0; r < R; ++r) Y[r].push_back(vals[r]);
            }
            int N = (int)angles.size();
            if (N < 4) { std::cout << Form("[WARN] %s %s-axis: only %d common points, skipping", serial.Data(), axes[a_idx], N) << std::endl; continue; }

            // Mean curve and covariance matrix C_ij over the R repeats.
            std::vector<double> Ybar(N, 0.0);
            for (int i = 0; i < N; ++i) { for (int r = 0; r < R; ++r) Ybar[i] += Y[r][i]; Ybar[i] /= R; }

            TMatrixDSym C(N);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    double c = 0;
                    for (int r = 0; r < R; ++r) c += (Y[r][i] - Ybar[i]) * (Y[r][j] - Ybar[j]);
                    C(i, j) = c / (R - 1);
                }

            // Regularize: covariance from only 3 repeats is noisy/near-singular
            // off the diagonal-adjacent band. Add a tiny diagonal floor (1e-6 of
            // mean variance) purely for numerical invertibility.
            double diagMean = 0; for (int i = 0; i < N; ++i) diagMean += C(i, i); diagMean /= N;
            for (int i = 0; i < N; ++i) C(i, i) += diagMean * 1e-6;

            TDecompChol chol(C);
            if (!chol.Decompose()) { std::cout << Form("[WARN] %s %s-axis: covariance not invertible, skipping", serial.Data(), axes[a_idx]) << std::endl; continue; }

            std::cout << Form("\n[ %s , %s-axis ]  (N=%d common points, R=%d repeats)", serial.Data(), axes[a_idx], N, R) << std::endl;
            double meanVal = 0; for (double v : Ybar) meanVal += v; meanVal /= N;
            double meanSigma = std::sqrt(diagMean);
            std::cout << Form("  Mean value                : %.4f", meanVal) << std::endl;
            std::cout << Form("  Mean per-point sigma (diag of C) : %.4f (%.2f %% of mean)", meanSigma, meanSigma / meanVal * 100) << std::endl;

            for (int r = 0; r < R; ++r) {
                TVectorD d(N);
                for (int i = 0; i < N; ++i) d(i) = Y[r][i] - Ybar[i];
                TVectorD Cinv_d = d;
                chol.Solve(Cinv_d);
                double chi2 = 0; for (int i = 0; i < N; ++i) chi2 += d(i) * Cinv_d(i);
                double chi2_ndof = chi2 / N;
                std::cout << Form("  Repeat %d (%-20s) : chi2/N = %.2f  %s", r + 1, tags[r].Data(), chi2_ndof,
                        chi2_ndof <= 1.5 ? "-> consistent with the shared covariance structure"
                                          : "-> this repeat deviates beyond the ensemble's own covariance") << std::endl;
            }
        }
    }
    std::cout << "\n===================================================================\n" << std::endl;
    for (auto* f : files) f->Close();
}
