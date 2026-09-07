#include <TGraph.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

void acf_numbers(const char* csv) {
    std::ifstream in(csv);
    std::vector<double> ep, q0, q1, q2;
    std::string line;
    while (std::getline(in, line)) {
        std::stringstream ss(line); std::string f[5];
        for (int i=0;i<5&&std::getline(ss,f[i],',');++i){}
        try { ep.push_back(std::stod(f[1])); q0.push_back(std::stod(f[2]));
              q1.push_back(std::stod(f[3])); q2.push_back(std::stod(f[4])); } catch(...) {continue;}
    }
    int n = ep.size();
    std::vector<double> t(n);
    for (int i=0;i<n;++i) t[i] = (ep[i]-ep[0])/60.0;
    std::vector<double> dts;
    for (int i=1;i<n;++i) dts.push_back(t[i]-t[i-1]);
    std::sort(dts.begin(),dts.end());
    double dtMed = dts[dts.size()/2];
    std::vector<double> grid;
    for (double g=0; g<t.back(); g+=dtMed) grid.push_back(g);
    int N = grid.size();
    double bart = 1.96/std::sqrt((double)N);
    printf("N=%d grid, dt=%.2f min, total=%.0f min (%.1f h), Bartlett=+-%.3f\n",
           N, dtMed, t.back(), t.back()/60.0, bart);
    const char* nm[3] = {"Mon.","Rot#1","Rot#2"};
    std::vector<double>* qs[3] = {&q0,&q1,&q2};
    for (int c=0;c<3;++c) {
        TGraph g(n, t.data(), qs[c]->data());
        std::vector<double> yi(N);
        for (int i=0;i<N;++i) yi[i]=g.Eval(grid[i]);
        double mean=0; for(double v:yi) mean+=v; mean/=N;
        std::vector<double> ym(N); for(int i=0;i<N;++i) ym[i]=yi[i]-mean;
        double den=0; for(double v:ym) den+=v*v;
        printf("\n[%s] ACF (raw, mean-removed):\n", nm[c]);
        for (int k=1;k<=20&&k<N;++k) {
            double num=0; for(int i=0;i<N-k;++i) num+=ym[i]*ym[i+k];
            double r=num/den;
            printf("  lag=%6.1f min  ACF=%+.3f %s\n", k*dtMed, r,
                   std::abs(r)>bart ? "  <== outside 95% bound":"");
        }
    }
}
