// repro_common.h -- shared reproducibility statistics.
//
// Split out of Draw_Reproducibility.C so the pairwise macro and the
// multi-set macro (Draw_Reproducibility_Multi.C, first tag = reference)
// compute "systematic" the same way. Keeping two copies in sync by hand is
// exactly how the two would drift apart.
#ifndef REPRO_COMMON_H
#define REPRO_COMMON_H

#include <vector>
#include <cmath>
#include <algorithm>

struct MatchedPoint {
    double angle;
    double v1, e1, v2, e2;
    double ratio, ratioErr, pull;
};

struct ReproStats {
    int n = 0;
    double meanRatio = 0, obsSpread = 0, expSpread = 0, pullRMS = 0, sysComp = 0;
    double meanRelErr1 = 0, meanRelErr2 = 0;
    double sysUL = 0;        // 95% CL upper limit on the systematic
    bool   sysResolved = false;  // false -> quote sysUL as "< x", not sysComp
    bool   statOverestimated = false;  // obs scatter below stat expectation even at 95% UL
};

static ReproStats ComputeRepro(const std::vector<MatchedPoint>& pts) {
    ReproStats s;
    s.n = (int)pts.size();
    if (s.n < 2) return s;

    for (auto& p : pts) s.meanRatio += p.ratio;
    s.meanRatio /= s.n;

    double obs2 = 0, exp2 = 0, pull2 = 0;
    for (auto& p : pts) {
        obs2  += (p.ratio - s.meanRatio) * (p.ratio - s.meanRatio);
        exp2  += p.ratioErr * p.ratioErr;
        pull2 += p.pull * p.pull;
        s.meanRelErr1 += p.e1 / p.v1;
        s.meanRelErr2 += p.e2 / p.v2;
    }
    s.obsSpread = std::sqrt(obs2 / (s.n - 1));
    s.expSpread = std::sqrt(exp2 / s.n);
    s.pullRMS   = std::sqrt(pull2 / s.n);
    s.sysComp   = (s.obsSpread > s.expSpread)
                   ? std::sqrt(s.obsSpread * s.obsSpread - s.expSpread * s.expSpread) : 0.0;
    s.meanRelErr1 /= s.n; s.meanRelErr2 /= s.n;

    // "syst. = 0.00%" was never an honest number. It only means the observed
    // scatter came out at or below the statistical expectation, which with
    // N~23 points happens routinely by chance -- RMS(pull) itself carries an
    // uncertainty of ~1/sqrt(2N) (~15% here), so a true systematic of a
    // couple of percent is entirely compatible with measuring RMS(pull)<1.
    // Reporting 0 advertises a precision the data cannot support.
    //
    // Instead: propagate that uncertainty and quote a 95% CL upper limit.
    //   RMS(pull) = R +- R/sqrt(2N)
    //   R_UL      = R + 1.645*R/sqrt(2N)          (one-sided 95%)
    //   sys_UL    = sqrt(max(0, R_UL^2 - 1)) * expSpread
    // When the central value already exceeds statistics (R>1 by more than
    // its own error) the systematic IS resolved and sysComp is quoted as a
    // value; otherwise only the limit is meaningful.
    {
        double relErrR = 1.0 / std::sqrt(2.0 * s.n);
        double rUL = s.pullRMS * (1.0 + 1.645 * relErrR);
        s.sysUL = std::sqrt(std::max(0.0, rUL * rUL - 1.0)) * s.expSpread;
        // Resolved only if R is above 1 by >1.645 sigma (same one-sided 95%).
        s.sysResolved = (s.pullRMS - 1.0) > 1.645 * s.pullRMS * relErrR;

        // If even the UPPER bound on RMS(pull) sits below 1, the quoted
        // per-point statistical errors are larger than the scatter they are
        // supposed to describe -- i.e. they are overestimated, and a limit
        // derived from them ("syst < 0.00 %") would claim an impossibly tight
        // constraint. Fall back to the one bound that holds regardless of how
        // well the error bars are calibrated: a systematic can never exceed
        // the total observed scatter.
        s.statOverestimated = (rUL < 1.0);
        if (s.sysUL <= 0.0) s.sysUL = s.obsSpread;
    }
    return s;
}

#endif // REPRO_COMMON_H
