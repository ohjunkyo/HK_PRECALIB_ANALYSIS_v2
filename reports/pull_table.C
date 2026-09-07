#include <TFile.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <vector>
#include <cmath>
#include <cstdio>

void pull_table() {
    std::vector<TString> tags = {"20260815_0_45","20260816_0_45","20260816_100_145","20260817_0_45"};
    std::vector<TString> pmts = {"EM6400","EL5150"};
    const char* axes[2] = {"X","Y"};

    printf("%-10s %-8s %8s %8s\n", "PMT", "axis", "N", "RMS(pull)");
    for (auto& pmt : pmts) {
        for (int a = 0; a < 2; ++a) {
            TFile* fref = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tags[0].Data()));
            TGraphErrors* gref = (TGraphErrors*)fref->Get(Form("gr_%s_%s_MonNorm_CorrRawQE", pmt.Data(), axes[a]));
            std::vector<double> pulls;
            for (size_t t = 1; t < tags.size(); ++t) {
                TFile* ftgt = TFile::Open(Form("./Data/UNIFORMITY/Graphs_Uniformity_%s.root", tags[t].Data()));
                TGraphErrors* gtgt = (TGraphErrors*)ftgt->Get(Form("gr_%s_%s_MonNorm_CorrRawQE", pmt.Data(), axes[a]));
                if (!gref || !gtgt) { if(ftgt) ftgt->Close(); continue; }
                for (int i = 0; i < gref->GetN(); ++i) {
                    double aRef = gref->GetX()[i], vRef = gref->GetY()[i], eRef = gref->GetErrorY(i);
                    if (vRef <= 0) continue;
                    for (int j = 0; j < gtgt->GetN(); ++j) {
                        if (std::abs(gtgt->GetX()[j]-aRef) > 0.1) continue;
                        double vTgt = gtgt->GetY()[j], eTgt = gtgt->GetErrorY(j);
                        if (vTgt <= 0) break;
                        double den = std::sqrt(eRef*eRef+eTgt*eTgt);
                        if (den > 0) pulls.push_back((vTgt-vRef)/den);
                        break;
                    }
                }
                ftgt->Close();
            }
            fref->Close();
            int N = pulls.size();
            double mean=0; for(double v:pulls) mean+=v; mean/=N;
            double s2=0; for(double v:pulls) s2+=(v-mean)*(v-mean); double rms=std::sqrt(s2/N);
            double errRms = rms/std::sqrt(2.0*N);
            printf("%-10s %-8s %8d %6.2f+-%.2f\n", pmt.Data(), axes[a], N, rms, errRms);
        }
    }
}
