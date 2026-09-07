#include <iostream>
#include <fstream>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TPolyLine3D.h>

// #define FLAG_OLD        // Old B-field code equivalent

// For HK PMT precalibration setup
#define HK_PRE_CALIB
#define N_PM_HK 6
#define N_COIL_HK 7

#define N_PM 11146
#define N_COIL_H 14
#define N_COIL_V 12
#define N_COIL_HH 6

// Output root file
char *outrootfile = nullptr;

// Length of coil segment for Biot-Savart calculation [mm]
Float_t l_seg = 100.;

// PMT position vectors
TVector3 v_pmt[N_PM];   // SK coordinate
TVector3 v_pmtC[N_PM];  // Coil coordinate (+X = south)

// For HK precalibration
TVector3 v_pmt_HK[N_PM_HK];

// B-field vectors
TVector3 vB[N_PM];
TVector3 vB_HK[N_PM_HK];

// Horizontal coil position parameters [mm]
// Float_t z_hor[N_COIL_H] = {-20700, -20700, -18200, -14200, -10200, -6200, -2200, 1800, 5800, 9800, 13800, 17800, 20700, 21200};
Float_t z_hor[N_COIL_H] = {21200, 20700, 17800, 13800, 9800, 5800, 1800, -2200, -6200, -10200, -14200, -18200, -20700, -20700};
Float_t r_hor[N_COIL_H] = {13000, 19000, 19000, 19000, 19000, 19000, 19000, 19000, 19000, 19000, 19000, 19000, 19000, 13000};
Float_t z_hh[N_COIL_HH] = {20500, 12400, 4100, -4100, -12600, -19600};

// Vertical coil position parameters in coil coordinate [mm]
Float_t x_ver_top[N_COIL_V] = {17639., 14420., 14900, 10700, 6400, 2130, -2130, -6400, -10700, -14900, -14382., -17639.};
Float_t x_ver_side[N_COIL_V] = {17639., 14420., 11193., 7850., 4735., 1483., -1717., -4927., -8142., -11371., -14382., -17639.};
Float_t z_ver_top[N_COIL_V] = {17800., 20700., 21000., 21000., 21000., 21000., 21000., 21000., 21000., 21000., 20700., 17800.};

// Rectangle coil position parameters for HK pre-calibration in coil coordinates[mm]
// coil 1~3 : x = constant 
// coil 4~5 : y = constant
// coil 6~7 : z = constant
// Need 4 line segments

// Original Geometry
Float_t x_HK_start[N_COIL_HK] = { -1700, 0, 1700, -1766, -1766, -1791, -1791};
Float_t y_HK_start[N_COIL_HK] = { -1598, -1598, -1598, 648, -1482, -1623, -1623};
Float_t z_HK_start[N_COIL_HK] = { 1020.5, 1020.5, 1020.5, 995.5, 995.5, 750, -750};

Float_t Coil_Length_x = 3532.;
Float_t Coil_Height_x = 2025.;
Float_t Coil_Length_y = 3532.;
Float_t Coil_Height_y = 1975.;
Float_t Coil_Length_z = 3582.;
Float_t Coil_Height_z = 3582.;

//Float_t x_HK_start[N_COIL_HK] = { -1350,  1050, -1800, -1800, -1800, -1800 }; // -1000, 1070 X-coil
//Float_t y_HK_start[N_COIL_HK] = {-1800, -1800,  1400,  -1400, -1800, -1800};
/////////// Move to 50 cm & -150 
//Float_t y_HK_start[N_COIL_HK] = { -1800, -1800, 500, -1500, -1800, -1800 };
//Float_t z_HK_start[N_COIL_HK] = { 1100,  1100,  1100,  1100,   550,  -550};
///////////// Increase the spacing  & Austrailia sys. height
//Float_t z_HK_start[N_COIL_HK] = { 1000, 1000, 1000, 1000, 750, -750 };

/*
   Float_t x_HK_top_end[N_COIL_HK] = { -900,  900,  1800,  1800, 1800, 1800};
   Float_t y_HK_top_end[N_COIL_HK] = { 1800, 1800,   900,  -900, -1800, -1800};
   Float_t z_HK_top_end[N_COIL_HK] = {};

   Float_t x_HK_side_start[N_COIL_HK] = { -900,   900, -1800, -1800, -1800, -1800};
   Float_t y_HK_side_start[N_COIL_HK] = {-1800,  -1800, 900, -900,   -1800, -1800};
   Float_t z_HK_side_start[N_COIL_HK] = { ,  ,};

   Float_t x_HK_side_end[N_COIL_HK] = { -900,  900,  1800,  1800, 1800, 1800};
   Float_t y_HK_side_end[N_COIL_HK] = { 1800, 1800,   900,  -900, -1800, -1800};
   Float_t z_HK_side_end[N_COIL_HK] = {};
   */


// Coil current [A]
Float_t curr_hor[N_COIL_H] = {
	31.35*4, 28.20*2*4, 31.35*4, 31.35*4, 31.35*4, 31.35*4, 
	31.35*4, 31.35*4,   31.35*4, 31.35*4, 31.35*4, 31.35*4, 28.20*2*4, 31.35*4
};

Float_t curr_ver[N_COIL_V] = {
	31.10*4, 31.10*4, 31.10*4, 31.10*4, 28.60*4, 28.60*4, 
	28.60*4, 28.60*4, 31.10*4, 31.10*4, 31.10*4, 31.10*4
};

Float_t curr_HK[N_COIL_HK] = {-40., -40., -40., -10., -10., 80., 80.};

// Conversion factor
Float_t am2mg=4*TMath::Pi()*1000;

// Radius of SK tank
#ifdef FLAG_OLD
Float_t r_tank = 19000.;
#else
Float_t r_tank = 19500.;
#endif

// Earth magnetic field (in coil coordinate)
/*
 * For SK,
 TVector3 b_earth_x(-303.7, 0.0, 0.0);
 TVector3 b_earth_z(0, 0.0, -351.8);
 TVector3 b_earth(-303.7, 0.0, -351.8);
 */ 
TVector3 b_earth(-288, 101, -386);

// Rotation angle to south
Float_t angle_south = -47.5*TMath::DegToRad();

// Magnetic field measure points [mm]
Float_t x_MF_measured[4] = {2000, 1000, 0.0, -1000};
Float_t y_MF_measured[4] = {1500, 500, -500, -1500}; 
//Float_t z_MF_measured[3] = {520, 900, 1300}; // 520-1000, 900-1000, 1300-1000
Float_t z_MF_measured[3] = {520-1000, 900-2000/2, 1300-2000/2}; // 520-1000, 900-1000, 1300-1000


Float_t earth_Bx[3][4][4] = {
	{{-300, -316, -322, -309}, {-305, -295, -301, -318}, {-252, -276, -282, -297}, {-204, -217, -238, -253}},
	{{-291, -315, -318, -309}, {-301, -299, -304, -316}, {-263, -280, -289, -300}, {-231, -245, -267, -286}},
	{{-290, -302, -302, -296}, {-283, -281, -291, -303}, {-260, -267, -280, -287}, {-238, -251, -277, -295}}};

Float_t earth_By[3][4][4] = {
	{{119, 116, 120, 144}, {166, 111, 105, 123}, {110, 116, 84, 127}, {93, 95, 81, 107}},
	{{131, 108, 101, 125}, {104, 106, 106, 114}, {105, 98, 84, 106}, {92, 87, 72, 90}},
	{{117, 110, 110, 131}, {111, 114, 107, 111}, {103, 105, 88, 110}, {99, 97, 73, 94}}};

Float_t earth_Bz[3][4][4] = {
	{{-351, -349, -347, -337}, {-391, -389, -386, -366}, {-462, -409, -399, -429}, {-409, -415, -406, -395}},
	{{-364, -356, -349, -338}, {-396, -390, -380, -365}, {-434, -409, -401, -399}, {-407, -409, -398, -382}},
	{{-354, -342, -335, -324}, {-375, -369, -359, -343}, {-401, -392, -384, -372}, {-388, -396, -386, -361}}};

// PolyLines for coil
TObjArray* poly_coil;

// Function to calculate B field at a specific grid point
void calc_segment_grid (TVector3 v_seg, TVector3 v_curr, Float_t d_seg, TVector3 point, TVector3& aB);
void calc_linear_grid (Float_t curr, TVector3 v_start, TVector3 v_end, TVector3 point, TVector3& aB);
void calc_coilHK_grid (TVector3 point, TVector3& aB);

// Function to calculate B field
void calc_segment (TVector3 v_seg, TVector3 v_curr, Float_t d_seg);
void calc_linear (Float_t curr, TVector3 v_start, TVector3 v_end);
void calc_radial (Float_t curr, TVector3 v_start, TVector3 v_end);

