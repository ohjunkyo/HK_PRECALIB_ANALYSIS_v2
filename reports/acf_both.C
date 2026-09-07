#include <TGraph.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

void acf_both(const char* csv, const char* label) {
    std::ifstream in(csv);
    std::vector<double> ep, q[3];
    std::string line;
    while (std::getline(in, line)) {
        std::stringstream ss(line); std::string f[5];
        for (int i=0;i<5&&std::getline(ss,f[i],',');++i){}
        try { ep.push_back(std::stod(f[1]));
              for(int c=0;c<3;++c) q[c].push_back(std::stod(f[2+c])); } catch(...) {continue;}
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
    int win = std::max(3,(int)std::round(180.0/dtMed));
    printf("\n########## %s ##########\n", label);
    printf("N=%d, dt=%.2f min, trend window=%d points (%.0f min), Bartlett=+-%.3f\n",
           N, dtMed, win, win*dtMed, bart);
    const char* nm[3] = {"Mon.","Rot#1","Rot#2"};
    for (int c=0;c<3;++c) {
        TGraph g(n, t.data(), q[c].data());
        std::vector<double> yi(N);
        for (int i=0;i<N;++i) yi[i]=g.Eval(grid[i]);
        double mean=0; for(double v:yi) mean+=v; mean/=N;
        std::vector<double> ym(N), tr(N), rs(N);
        for(int i=0;i<N;++i) ym[i]=yi[i]-mean;
        int half=win/2;
        for(int i=0;i<N;++i){int lo=std::max(0,i-half),hi=std::min(N,i+half+1);
            double s=0;for(int k=lo;k<hi;++k)s+=yi[k];tr[i]=s/(hi-lo);}
        for(int i=0;i<N;++i) rs[i]=yi[i]-tr[i];
        auto acf=[&](std::vector<double>&v,int k){double num=0,den=0;
            for(double x:v)den+=x*x; for(int i=0;i<N-k;++i)num+=v[i]*v[i+k];return num/den;};
        printf("\n[%s]  lag(min)   ACF_before   ACF_after\n", nm[c]);
        for (int k=1;k<=8;++k) {
            double a=acf(ym,k), b=acf(rs,k);
            printf("        %7.1f   %+7.3f%s  %+7.3f%s\n", k*dtMed,
                   a, std::abs(a)>bart?"*":" ", b, std::abs(b)>bart?"*":" ");
        }
    }
}
