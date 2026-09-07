// angle_convert.h — single source of truth for the rotation-stage angle ->
// PMT incidence (Hamamatsu) angle conversion, shared by every analysis macro
// (prod_ntp_v7, read_ntp_v7, Draw_Contour_v3, Draw_Uniformity_Norm_v7, run_info).
//
// Geometry: cable pin A~H sits at a fixed angle around the PMT; the PMT's +X
// axis points along pin G, +Y along pin A. To scan the X (or Y) axis the stage
// rotates so that axis aligns with the fixed scan axis -> x_rot / y_rot below.
//
// Sign convention (R12860-22 bottom view): X+ = right, Y+ = down. With cable A
// mounted in the fixed hole the PMT frame == lab frame; measured zenith(+) probes
// the lab "left" side. Because the rotate stage only spans 0~135 deg, opposite
// cable pairs (A/E, B/F, C/G, D/H) share the same x_rot/y_rot (mod 180) yet face
// physically opposite cathode regions -- so the incidence-angle SIGN must come
// from the full cable azimuth, not from x_rot/y_rot alone. GetHamamatsuAngle below
// recovers that full azimuth (Delta) to assign the correct sign for all 8 cables.
#ifndef ANGLE_CONVERT_H
#define ANGLE_CONVERT_H

#include <TMath.h>
#include <cctype>
#include <cmath>

// Default cable direction per channel when RunInfo has no Direction branch
// (older data). CH0 = Monitor PMT (A); CH1/CH2 = test PMTs (B).
inline char DirForCh(int ch) { return (ch == 0) ? 'A' : 'B'; }

// Standard pin position (deg) around the PMT for each cable direction A~H.
inline int PosMapAngle(char dir) {
    switch (toupper(dir)) {
        case 'E': return 0;   case 'F': return 45;  case 'G': return 90;  case 'H': return 135;
        case 'A': return 180; case 'B': return 225; case 'C': return 270; case 'D': return 315;
    }
    return 180; // default = A
}

// Rotate-stage angles that align the PMT's X / Y axis with the fixed scan axis,
// derived from the cable direction. Verified: A -> 0/90, B -> 45/135.
inline void GetXYRotForDirection(char dir, int& x_rot, int& y_rot) {
    int pm = PosMapAngle(dir);
    int xm = (((pm - 90)  % 180) + 180) % 180 - 90;   // canonical (-90, 90]
    int ym = (((pm - 180) % 180) + 180) % 180 - 90;
    x_rot = (xm < 0) ? xm + 180 : xm;
    y_rot = (ym < 0) ? ym + 180 : ym;
}

// 2nd-order calibration polynomial: stage angle (KR) -> PMT incidence (Hamamatsu).
inline double ConvertKRtoHamamatsu(double kr) {
    return -0.0049 * TMath::Power(kr, 2) + 1.7515 * kr - 0.0402;
}

// Determine scan axis (X/Y) from cable direction + rotate angle, and return the
// signed Hamamatsu incidence angle. axisLabel is set to "X-axis"/"Y-axis"/"?-axis".
inline double GetHamamatsuAngle(char dir, double tilt_val, double rot_val, const char*& axisLabel) {
    int x_rot, y_rot;
    GetXYRotForDirection(dir, x_rot, y_rot);
    int rot = (int)rot_val;
    bool isX = (rot == x_rot);
    bool isY = (rot == y_rot);
    axisLabel = isX ? "X-axis" : (isY ? "Y-axis" : "?-axis");

    // Recover the full cable azimuth so opposite cables (A/E, B/F, C/G, D/H) get
    // distinct signs even though they collapse under x_rot/y_rot (mod 180).
    int delta = ((180 - PosMapAngle(dir)) % 360 + 360) % 360;
    // X scan: tilt(+) -> X-  (sign flip) for cables E,F,G,H; X+ for A,B,C,D.
    bool xflip = ((delta + x_rot) % 360 == 180);
    // Y scan: tilt(+) -> Y+  (sign flip) for cables C,D,E,F; Y- for A,B,G,H.
    bool yflip = ((delta + y_rot) % 360 == 270);

    // Anchor to cable A (validated): X-scan tilt(+) -> +X (+ham), Y-scan tilt(+) -> -Y.
    double sgn = (tilt_val < 0 ? -1.0 : 1.0);
    if (isX)      sgn *= xflip ? -1.0 :  1.0;
    else if (isY) sgn *= yflip ?  1.0 : -1.0;
    return sgn * ConvertKRtoHamamatsu(std::abs(tilt_val));
}

#endif // ANGLE_CONVERT_H
