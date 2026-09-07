#include <TError.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMarker.h>
#include <TPolyLine.h>
#include <TBox.h>
#include <calc_bfield.h>
#include <special_coil.h>

void calc_segment_grid(TVector3 v_seg, TVector3 v_curr, Float_t d_seg, TVector3 point, TVector3& aB) {
	// Position of grid point from coil segment
	TVector3 v_point_rel = point - v_seg;

	// Calculate vector product
	TVector3 v_prod = v_curr.Cross(v_point_rel);

	// Distance from coil segment to grid point
	Float_t r_this = v_point_rel.Mag();

	// Biot-Savart
	Float_t factor = 1 / (4*TMath::Pi()) / pow(r_this, 3) * am2mg * d_seg;
	TVector3 vB_seg = v_prod * factor;

	// Add to magnetic field vector for this grid point
	aB += vB_seg;
}

// Calculate B-field of all PMTs for linear coil part
void calc_linear_grid(Float_t curr, TVector3 v_start, TVector3 v_end, TVector3 point, TVector3& aB) {
	// Number of segments
	Float_t l_coil = (v_end - v_start).Mag();
	Int_t n_seg = (int)(l_coil / l_seg);
	Float_t d_seg = l_coil / n_seg;

	// Loop for segments
	for (Int_t iSeg=0; iSeg<n_seg; iSeg++) {

#ifdef FLAG_OLD
		Float_t factor = 1. / n_seg * (iSeg + 0.0);
#else
		Float_t factor = 1. / n_seg * (iSeg + 0.5);
#endif
		TVector3 v_seg = v_start + (v_end - v_start) * factor;

		// Current vector
		TVector3 v_curr = v_end - v_start;
		v_curr = v_curr.Unit() * curr;

		calc_segment_grid(v_seg, v_curr, d_seg, point, aB);
	}
}     

void calc_coilHK_grid(TVector3 point, TVector3& aB) { 
	for (Int_t iCoil = 0; iCoil < N_COIL_HK; iCoil++) { 
		if (iCoil <= 2) {                          
			TVector3 v_start_top (x_HK_start[iCoil], y_HK_start[iCoil], z_HK_start[iCoil]);
			TVector3 v_end_top   (x_HK_start[iCoil], y_HK_start[iCoil] + Coil_Length_x, z_HK_start[iCoil]);

			TVector3 v_start_side1 = v_end_top;
			TVector3 v_end_side1   = v_start_side1;

			v_end_side1[2] -= Coil_Height_x;               

			TVector3 v_start_bot = v_end_side1;
			TVector3 v_end_bot   = v_start_bot;

			v_end_bot[1] -= Coil_Length_x;

			TVector3 v_start_side2 = v_end_bot;
			TVector3 v_end_side2   = v_start_top;

			calc_linear_grid(curr_HK[iCoil], v_start_top, v_end_top, point, aB);
			calc_linear_grid(curr_HK[iCoil], v_start_side1, v_end_side1, point, aB);
			calc_linear_grid(curr_HK[iCoil], v_start_bot, v_end_bot, point, aB);
			calc_linear_grid(curr_HK[iCoil], v_start_side2, v_end_side2, point, aB);
		}
		else if (iCoil <= 4) {
			TVector3 v_start_top (x_HK_start[iCoil], y_HK_start[iCoil], z_HK_start[iCoil]);
			TVector3 v_end_top   (x_HK_start[iCoil] + Coil_Length_y, y_HK_start[iCoil], z_HK_start[iCoil]);

			TVector3 v_start_side1 = v_end_top;
			TVector3 v_end_side1   = v_start_side1;

			v_end_side1[2] -= Coil_Height_y;

			TVector3 v_start_bot = v_end_side1;
			TVector3 v_end_bot   = v_start_bot;

			v_end_bot[0] -= Coil_Length_y;

			TVector3 v_start_side2 = v_end_bot; 
			TVector3 v_end_side2   = v_start_top; 

			calc_linear_grid(curr_HK[iCoil], v_start_top, v_end_top, point, aB);
			calc_linear_grid(curr_HK[iCoil], v_start_side1, v_end_side1, point, aB); 
			calc_linear_grid(curr_HK[iCoil], v_start_bot, v_end_bot, point, aB);
			calc_linear_grid(curr_HK[iCoil], v_start_side2, v_end_side2, point, aB);
		}
		else { 
			TVector3 v_start1 (x_HK_start[iCoil], y_HK_start[iCoil], z_HK_start[iCoil]);
			TVector3 v_end1   (x_HK_start[iCoil] + Coil_Length_z, y_HK_start[iCoil], z_HK_start[iCoil]);

			TVector3 v_start2 = v_end1; 
			TVector3 v_end2   = v_start2;

			v_end2[1] += Coil_Height_z; 

			TVector3 v_start3 = v_end2;
			TVector3 v_end3   = v_start3;

			v_end3[0] -= Coil_Length_z;

			TVector3 v_start4 = v_end3;
			TVector3 v_end4   = v_start1;

			calc_linear_grid(curr_HK[iCoil], v_start1, v_end1, point, aB);
			calc_linear_grid(curr_HK[iCoil], v_start2, v_end2, point, aB);
			calc_linear_grid(curr_HK[iCoil], v_start3, v_end3, point, aB);
			calc_linear_grid(curr_HK[iCoil], v_start4, v_end4, point, aB);
		}                                                                                                                                                                                                                       
	}                                                                                                                                                                                                                                
}

TVector3 getAverageEarthField() {
	Float_t sum_Bx = 0.0, sum_By = 0.0, sum_Bz = 0.0;
	int count = 0;
	int max = 0;
	// Loop through all height levels (52cm, 90cm, 130cm) and all grid points
	for (int iz = 0; iz < 3; iz++) {
		for (int ix = 0; ix < 4; ix++) {
			for (int iy = 0; iy < 4; iy++) {
				// Accumulate the values of the magnetic field components
				sum_Bx += earth_Bx[iz][ix][iy];
				sum_By += earth_By[iz][ix][iy];
				sum_Bz += earth_Bz[iz][ix][iy];
				count++;  // Increment the count for averaging
			}
		}
	}

	// Compute the average for each component
	Float_t avg_Bx = sum_Bx / count;
	Float_t avg_By = sum_By / count;
	Float_t avg_Bz = sum_Bz / count;
	// Return the average values as a vector (Bx, By, Bz)
	return TVector3(avg_Bx, avg_By, avg_Bz);
}

void reset_magnetic_field(){
	for (Int_t iPM=0; iPM<N_PM; iPM++) {
		vB[iPM].SetXYZ(0,0,0);
	}
}

// Calculate B-field of all PMTs for one coil segment
void calc_segment (TVector3 v_seg, TVector3 v_curr, Float_t d_seg) {

	// Loop for PMTs
	for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {

		// Position of PMT from coil segment
		TVector3 v_pmt_rel = v_pmt_HK[iPM] - v_seg;

		// Calculate vector product
		TVector3 v_prod = v_curr.Cross(v_pmt_rel);

		// Distance from coil segment to PMT
		Float_t r_this = v_pmt_rel.Mag();

		// Biot-Savart
		Float_t factor = 1 / (4*TMath::Pi()) / pow(r_this, 3) * am2mg * d_seg;
		TVector3 vB_seg = v_prod * factor;

		// Add to B-field vector
		vB_HK[iPM] += vB_seg;

	}

}

// Calculate B-field of all PMTs for linear coil part
void calc_linear (Float_t curr, TVector3 v_start, TVector3 v_end) {

	// Number of segments
	Float_t l_coil = (v_end - v_start).Mag();
	Int_t n_seg = (int)(l_coil / l_seg);
	Float_t d_seg = l_coil / n_seg;

	// Initialize poly-line
	TPolyLine3D* pl = new TPolyLine3D(n_seg+1);
	pl->SetNextPoint(v_start.X(), v_start.Y(), v_start.Z());

	// Loop for segments
	for (Int_t iSeg=0; iSeg<n_seg; iSeg++) {

		// Position of this segment
#ifdef FLAG_OLD
		Float_t factor = 1. / n_seg * (iSeg + 0.0);
#else
		Float_t factor = 1. / n_seg * (iSeg + 0.5);
#endif      
		TVector3 v_seg = v_start + (v_end - v_start) * factor;

		// Current vector
		TVector3 v_curr = v_end - v_start;
		v_curr = v_curr.Unit() * curr;

		// Calculate B-field for this segment
		calc_segment(v_seg, v_curr, d_seg);

		// Add poly-line
		TVector3 v1 = v_seg + v_curr.Unit() * (d_seg / 2.);
		pl->SetNextPoint(v1.X(), v1.Y(), v1.Z());

	}

	// Add poly-line to TObjArray for poly-line 
	poly_coil->Add(pl);

}

// Calculate B-field of all PMTs for radial coil part
void calc_radial (Float_t curr, TVector3 v_start, TVector3 v_end) {

	// Radius of this coil
	TVector3 v_center(0, 0, v_start.Z());
	Float_t r_coil = (v_start - v_center).Mag();

	// Angle between start and end
	Float_t angle_this = TMath::Pi()*2;
	if (v_start != v_end) { // Part of vertical coils
		TVector3 v_start_rel = v_start - v_center;
		TVector3 v_end_rel   = v_end   - v_center;
		angle_this = v_start_rel.Angle(v_end_rel);
	}

	// Number of segments
	Float_t l_coil = r_coil * fabs(angle_this);
#ifdef FLAG_OLD   
	Int_t n_seg = (Int_t)(fabs(angle_this) * 180./ TMath::Pi());
#else
	Int_t n_seg = (int)(l_coil / l_seg);
#endif

	// One segment length
	Float_t d_seg = l_coil / n_seg;
#ifdef FLAG_OLD
	if (v_start != v_end) { // Part of vertical coils
		d_seg = 100.;
	}
#endif

	// Direction to go
	Float_t dir_sign = +1;
	if (v_start != v_end) {
		TVector3 v_start_rel = v_start - v_center;
		TVector3 v_end_rel   = v_end   - v_center;
		TVector3 v_prod = v_start_rel.Cross(v_end_rel);
		dir_sign = v_prod.Z() / fabs(v_prod.Z());
	}

	// Initialize poly-line
	TPolyLine3D* pl = new TPolyLine3D(n_seg+1);
	pl->SetNextPoint(v_start.X(), v_start.Y(), v_start.Z());

	// Loop for segments
	for (Int_t iSeg=0; iSeg<n_seg; iSeg++) {

		// Position of this segment
		TVector3 v_seg = v_start;
#ifdef FLAG_OLD      
		Float_t phi_this = angle_this / n_seg * (iSeg + 0.0);
#else
		Float_t phi_this = angle_this / n_seg * (iSeg + 0.5);
#endif      
		v_seg.RotateZ(dir_sign*phi_this);

		// Current vector
		TVector3 v_curr = v_seg - v_center;
		v_curr.RotateZ(dir_sign*TMath::Pi()/2);
		v_curr = v_curr.Unit() * curr;

		// Calculate B-field for this segment
		calc_segment(v_seg, v_curr, d_seg);

		// Add poly-line
		TVector3 v1 = v_seg + v_curr.Unit() * (d_seg / 2.);
		pl->SetNextPoint(v1.X(), v1.Y(), v1.Z());

	}

	// Add poly-line to TObjArray for poly-line 
	poly_coil->Add(pl);
}

// M.W Lee
// For HK pre-calibration simulation, 7 rectangle shaped-coils
void calc_coilHK()
{
	// Input : curr : current value
	// 		   v_start : coil position
	// 		   v_end : same

	std::cout <<"========= << COIL GEOMETRY >> ==========" << std::endl;
	for (Int_t iCoil = 0; iCoil < N_COIL_HK; iCoil++)
	{
		//std::cout << "HK PreCalibration Room coil " << iCoil << Form(" HK%02d", iCoil + 1) << std::endl;

		if (iCoil <= 2) // coil x position = constant --> Vertical
		{
			// Start/end position of each segment
			TVector3 v_start_top (x_HK_start[iCoil], y_HK_start[iCoil], z_HK_start[iCoil]);
			TVector3 v_end_top   (x_HK_start[iCoil], y_HK_start[iCoil] + Coil_Length_x, z_HK_start[iCoil]);
			std::cout<<'['<<N_COIL_HK<<']'<<std::endl;
			TVector3 v_start_side1 = v_end_top;
			TVector3 v_end_side1   = v_start_side1;

			v_end_side1[2] -= Coil_Height_x;

			TVector3 v_start_bot = v_end_side1;
			TVector3 v_end_bot   = v_start_bot;

			v_end_bot[1] -= Coil_Length_x;

			TVector3 v_start_side2 = v_end_bot;
			TVector3 v_end_side2   = v_start_top;

			double length_top = (v_end_top - v_start_top).Mag();
			double length_side1 = (v_end_side1 - v_start_side1).Mag();
			double length_bot = (v_end_bot - v_start_bot).Mag();
			double length_side2 = (v_end_side2 - v_start_side2).Mag();

			//			std::cout << Form("X-coil(%d) lengths: top=%.1f, height1=%.1f, bot=%.1f, height2=%.1f", iCoil, length_top, length_side1, length_bot, length_side2) << std::endl;
			std::cout << Form("X-coil(%d) top(%.0f, %.0f, %.0f)",
					iCoil, v_start_top.X(), v_start_top.Y(), v_start_top.Z()) << std::endl;

			// Calculate B-field for 4 segments
			calc_linear(curr_HK[iCoil], v_start_top, v_end_top); 
			calc_linear(curr_HK[iCoil], v_start_side1, v_end_side1); 
			calc_linear(curr_HK[iCoil], v_start_bot, v_end_bot); 
			calc_linear(curr_HK[iCoil], v_start_side2, v_end_side2); 
		}
		else if (iCoil <= 4) // coil y position = constant --> Vertical
		{
			TVector3 v_start_top (x_HK_start[iCoil], y_HK_start[iCoil], z_HK_start[iCoil]);
			TVector3 v_end_top   (x_HK_start[iCoil] + Coil_Length_y, y_HK_start[iCoil], z_HK_start[iCoil]);

			TVector3 v_start_side1 = v_end_top;
			TVector3 v_end_side1   = v_start_side1;

			v_end_side1[2] -= Coil_Height_y;

			TVector3 v_start_bot = v_end_side1;
			TVector3 v_end_bot   = v_start_bot;
			v_end_bot[0] -= Coil_Length_y;

			TVector3 v_start_side2 = v_end_bot;
			TVector3 v_end_side2   = v_start_top;

			double length_top = (v_end_top - v_start_top).Mag();
			double length_side1 = (v_end_side1 - v_start_side1).Mag();
			double length_bot = (v_end_bot - v_start_bot).Mag();
			double length_side2 = (v_end_side2 - v_start_side2).Mag();


			//			std::cout << Form("Y-coil(%d) lengths: top=%.1f, height1=%.1f, bot=%.1f, height2=%.1f", iCoil, length_top, length_side1, length_bot, length_side2) << std::endl;
			std::cout << Form("Y-coil(%d) top(%.0f, %.0f, %.0f)", 
					iCoil, v_start_top.X(), v_start_top.Y(), v_start_top.Z()) << std::endl;

			// Calculate B-field for 4 segments
			calc_linear(curr_HK[iCoil], v_start_top, v_end_top); 
			calc_linear(curr_HK[iCoil], v_start_side1, v_end_side1); 
			calc_linear(curr_HK[iCoil], v_start_bot, v_end_bot); 
			calc_linear(curr_HK[iCoil], v_start_side2, v_end_side2); 
		}
		else // coil z position = constant --> Horizontal
		{
			TVector3 v_start1 (x_HK_start[iCoil], y_HK_start[iCoil], z_HK_start[iCoil]);
			TVector3 v_end1   (x_HK_start[iCoil] + Coil_Length_z, y_HK_start[iCoil], z_HK_start[iCoil]);

			TVector3 v_start2 = v_end1;
			TVector3 v_end2   = v_start2;
			v_end2[1] += Coil_Height_z;

			TVector3 v_start3 = v_end2;
			TVector3 v_end3   = v_start3;
			v_end3[0] -= Coil_Length_z;

			TVector3 v_start4 = v_end3;
			TVector3 v_end4   = v_start1;

			double length1 = (v_end1 - v_start1).Mag();
			double length2 = (v_end2 - v_start2).Mag();
			double length3 = (v_end3 - v_start3).Mag();
			double length4 = (v_end4 - v_start4).Mag();

			//			std::cout << Form("Z-coil(%d) lengths: seg1=%.2f, seg2=%.2f, seg3=%.2f, seg4=%.2f", iCoil, length1, length2, length3, length4) << std::endl;
			std::cout << Form("Z-coil(%d): seg1(%.0f, %.0f, %.0f)",
					iCoil, v_start1.X(), v_start1.Y(), v_start1.Z()) << std::endl;

			// Calculate B-field for 4 segments
			calc_linear(curr_HK[iCoil], v_start1, v_end1); 
			calc_linear(curr_HK[iCoil], v_start2, v_end2); 
			calc_linear(curr_HK[iCoil], v_start3, v_end3); 
			calc_linear(curr_HK[iCoil], v_start4, v_end4); 
		}
	}// End of B-field calculation for each coils
}


// Main function
int main(int argc, char *argv[]) {

	// Read PMT position file
	std::ifstream ifs_pmt_pos("pmtpos.dat");
	Int_t cab_this; //std::cout<<"417!!!"<<std::endl;
	Float_t x_this, y_this, z_this;
	for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
		ifs_pmt_pos >> cab_this >> x_this >> y_this >> z_this;
		v_pmt_HK[cab_this-1].SetXYZ(x_this, y_this, z_this);
	}

	if (argc != 9)
	{
		std::cout << "No enough arguments(Need to 7 current values, and output file name)" << std::endl; 
		return -1;
	}

	// Get HK precalibration current
	for (int i = 1; i < argc-1; i++)
	{
		curr_HK[i-1] = std::stof(argv[i]);
	}

	outrootfile = argv[argc-1];

	// Initialize B-field vectors
	for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
		vB_HK[iPM].SetXYZ(0,0,0);
	}

	// Initialize TObjArray for poly-line 
	poly_coil = new TObjArray();

	// Draw poly-lines
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 1200, 1200);
	TCanvas* c2 = new TCanvas("th2d", "2D", 0, 0, 1200, 1200);
	TCanvas* c3 = new TCanvas("th3d", "3D", 0, 0, 1200, 1200);
	TCanvas* c4 = new TCanvas("th2d2", "2D", 0, 0, 1200, 1200);
	TCanvas* c5 = new TCanvas("th2d3", "2D", 0, 0, 1200, 1200);
	
	// Calculate B-field from HK precalibration coils
	calc_coilHK();

	c1->cd();
	// Draw HK precalibration coils
	Int_t nPoly = poly_coil->GetEntries();
	for (Int_t iPoly=0; iPoly<nPoly; iPoly++) {
		TPolyLine3D* poly_this = (TPolyLine3D*)poly_coil->At(iPoly);
		poly_this->SetLineWidth(2);
		if(iPoly<8)poly_this->SetLineColor(kRed);
		else if(iPoly<16)poly_this->SetLineColor(kBlue);
		else if(iPoly<24)poly_this->SetLineColor(kBlack);
		poly_this->Draw("same");
	}

	// Define tree to output results
	//TFile* ofile = new TFile("test.root", "RECREATE");
	TFile* ofile = new TFile(outrootfile, "RECREATE");
	if (!ofile || ofile->IsZombie()) {
		std::cerr << "Error: Could not create output file!" << std::endl;
		return -1;
	}

	// Variables for general data
	Float_t i1 = curr_HK[0];
	Float_t i2 = curr_HK[1];
	Float_t i3 = curr_HK[2];
	Float_t i4 = curr_HK[3];
	Float_t i5 = curr_HK[4];
	Float_t i6 = curr_HK[5];
	Float_t i7 = curr_HK[6];

	Int_t aPM;
	TVector3 aPos, aB;
	Float_t aBx, aBy, aBz, aBpara, aBperp, aBmag;

	//////////////////// PMT position //////////////////

	TTree* tree = new TTree("tree", "PMT Magnetic Field Data");

	tree->Branch("i1", &i1);
	tree->Branch("i2", &i2);
	tree->Branch("i3", &i3);
	tree->Branch("i4", &i4);
	tree->Branch("i5", &i5);
	tree->Branch("i6", &i6);
	tree->Branch("i7", &i7);
	tree->Branch("pm", &aPM);
	tree->Branch("pos", &aPos);
	tree->Branch("B", &aB);
	tree->Branch("Bx", &aBx);
	tree->Branch("By", &aBy);
	tree->Branch("Bz", &aBz);
	tree->Branch("Bmag", &aBmag);
	tree->Branch("Bperp", &aBperp);
	tree->Branch("Bpara", &aBpara);

	// Loop over PMTs
	TVector3 earthMF = getAverageEarthField();

	std::cout <<Form("Average EMF | Bx: %.2f, By: %.2f, Bz: %.2f, Bmag: %.2f"
			,earthMF.X(), earthMF.Y(), earthMF.Z(), earthMF.Mag()) <<std::endl;
	std::cout <<"========= << MAGNETIC FIELD CALCULATING >> ==========" << std::endl;
	std::cout <<"< PMT POSITION >" << std::endl;
	std::vector<float> aBmag_vec;

	for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {

		// PMT position
		aPos.SetXYZ(v_pmt_HK[iPM].X(), v_pmt_HK[iPM].Y(), v_pmt_HK[iPM].Z());

		// Calculated magnetic field for PMT
		aB.SetXYZ(vB_HK[iPM].X(), vB_HK[iPM].Y(), vB_HK[iPM].Z());
		// aB += b_earth_x;

		aB += earthMF;

		aBmag = aB.Mag();
		aBx = aB.X(); aBy = aB.Y(); aBz = aB.Z();

		aBmag_vec.push_back(aBmag);

		// Calculate perpendicular component
		TVector3 v_dir;
		v_dir = TVector3(0, 0, +1);
		aBperp = sqrt(aB.Mag2() - pow(aB * v_dir, 2)); 
		aBpara = aB*v_dir;

		// Fill tree
		aPM = iPM;

		std::cout <<Form("iPM: %d, Position (X): %.1f, (Y): %.1f, (Z): %.1f, Bx: %.2f, By: %.2f, Bz: %.2f, Bmag: %.2f"
				,iPM, v_pmt_HK[iPM].X()/10., v_pmt_HK[iPM].Y()/10., v_pmt_HK[iPM].Z()/10., aBx, aBy, aBz, aBmag) << std::endl;


		tree->Fill();


	} // End of loop over PMTs

	///////////////////// GRID ///////////////////
	std::cout <<"< GRID SYSTEM >" << std::endl;

	TTree* tree_Grid = new TTree("tree_Grid", "Grid Magnetic Field Data");

	Float_t gridaB = 0; 
	Float_t gridaBx = 0;
	Float_t gridaBy = 0;
	Float_t gridaBz = 0;
	Float_t gridaBmag = 0;

	tree_Grid->Branch("gridaB", &gridaB);
	tree_Grid->Branch("gridaBx", &gridaBx);
	tree_Grid->Branch("gridaBy", &gridaBy);
	tree_Grid->Branch("gridaBz", &gridaBz);
	tree_Grid->Branch("gridaBmag", &gridaBmag);

	Float_t x_min = -2100.;
	Float_t x_max = 2100.;
	Float_t y_min = -2100.;
	Float_t y_max = 2500.;
	Float_t z_min = -1300.;
	Float_t z_max = 1500.;

	Float_t dx = 100.0;
	Float_t dy = 100.0;
	Float_t dz = 100.0;

	Int_t gridX = (x_max - x_min)/dx;
	Int_t gridY = (y_max - y_min)/dy;
	Int_t gridZ = (z_max - z_min)/dz;

	Int_t total_size = gridX * gridY * gridZ;
	Int_t progress_counter = 0;

	gStyle->SetNumberContours(100);

	gStyle->SetOptStat(0);

	TH2F *H2F = new TH2F("H2F", "H2F", gridX, x_min, x_max, gridY, y_min, y_max); // Large area of PMT (-50cm)
	TH2F *H2F2 = new TH2F("H2F2", "H2F2", gridX, x_min, x_max, gridY, y_min, y_max); // PMT Photocathode (0m)
        TH2F *H2F3 = new TH2F("H2F3", "H2F3", gridX, x_min, x_max, gridY, y_min, y_max); // PMT Photocathode (50cm)
	TH3F *H3F = new TH3F("H3F", "H3F", gridX, x_min, x_max, gridY, y_min, y_max, gridZ, z_min, z_max); 

	H2F->SetMinimum(0);
        H2F->SetMaximum(100);
        H2F2->SetMinimum(0);
        H2F2->SetMaximum(100);
        H2F3->SetMinimum(0);
        H2F3->SetMaximum(100);

	Float_t target_Z = v_pmt_HK[0].Z(); // Z = 1000 mm
	std::cout << Form("Target Height: %.2f, 0 [mm]", target_Z) << std::endl;

	H2F->SetTitle(Form("Bmag at %.2fm = PMT Height",target_Z/1000));
	H3F->SetTitle("Coil Height & PMT Height");
	H2F2->SetTitle("Bmag at Center");
	H2F3->SetTitle(Form("Bmag at %.2fm = PMT Height",target_Z/1000));

	for (int ix = 0; ix < gridX; ix++){
		for (int iy = 0; iy < gridY; iy++){
			for (int iz = 0; iz < gridZ; iz++){

				Float_t varx = x_min + ix * dx;
				Float_t vary = y_min + iy * dy;
				Float_t varz = z_min + iz * dz;

				TVector3 gridaB_vec(0,0,0);
				TVector3 point(varx, vary, varz);

				calc_coilHK_grid(point, gridaB_vec);

				gridaB_vec += earthMF;

				gridaBx = gridaB_vec.X();
				gridaBy = gridaB_vec.Y();
				gridaBz = gridaB_vec.Z();

			        //gridaBx = 0;
                                //gridaBy = 0;
                                //gridaBz = 0;

				gridaBmag = gridaB_vec.Mag();  
				Float_t gridaB_fill = gridaBmag;
				if(gridaB_fill <= -100||gridaB_fill >= 100) gridaB_fill = -9999.;

				if (gridaB_fill) {
					if(varz >= -500 && varz <= -500){
						H2F->Fill(varx, vary, gridaB_fill);
					//	H3F->Fill(varx, vary, varz, gridaBmag);

					}
					if(varz >= 0 && varz <= 0){
						H2F2->Fill(varx, vary, gridaB_fill);
					}
					if(varz >= +500 && varz <= +500){
                                                H2F3->Fill(varx, vary, gridaB_fill);
                                        }
				}

				for (int iPM = 0; iPM < N_PM_HK; iPM++) {
					if (std::abs(varx - v_pmt_HK[iPM].X()) < dx / 2 &&
							std::abs(vary - v_pmt_HK[iPM].Y()) < dy / 2 &&
							std::abs(varz - v_pmt_HK[iPM].Z()) < dz / 2) {

						std::cout << Form("Grid Coordinates: (%.0f, %.0f, %.0f)\t >>> \t PMT %d: (%.0f, %.0f, %.0f) >>>  Bmag: %.2f ==> %.2f",
								varx, vary, varz, iPM + 1, 
								v_pmt_HK[iPM].X(), v_pmt_HK[iPM].Y(), v_pmt_HK[iPM].Z(), gridaBmag, aBmag_vec[iPM]) << std::endl;
					}
				}
				tree_Grid->Fill(); //std::cout<<"647!!!"<<std::endl;

			}
		}
	}

	c2->cd();
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	gStyle->SetPalette(55);
	gStyle->SetHistMinimumZero(kTRUE);
	H2F->Draw("ARR COLZ same");  // Large area of PMT (-19 cm)

	TMarker* mar[N_PM_HK];
	for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
		Float_t xPM = v_pmt_HK[iPM].X(); 
		Float_t yPM = v_pmt_HK[iPM].Y();
		Float_t zPM = v_pmt_HK[iPM].Z();
		mar[iPM] = new TMarker(xPM, yPM, kOpenCircle);
		mar[iPM]->SetMarkerSize(15);
		mar[iPM]->SetMarkerColor(kBlack);
		mar[iPM]->Draw("same");
	}

	for (int iCoil =0; iCoil < 5; iCoil++) {
		TPolyLine *PLine= new TPolyLine(2);
		PLine->SetLineColor(kBlack);
		PLine->SetLineWidth(3);

		if ( iCoil < 3 ) {
			PLine->SetPoint(0, x_HK_start[iCoil], y_HK_start[iCoil]);
			PLine->SetPoint(1, x_HK_start[iCoil], (y_HK_start[iCoil] + 3532));
		}
		else if ( iCoil < 5 ) {
			PLine->SetPoint(0, x_HK_start[iCoil], y_HK_start[iCoil]);
			PLine->SetPoint(1, (x_HK_start[iCoil] + 3532), y_HK_start[iCoil]);
		}
		PLine->Draw();
	}

	c3->cd(); // 3D area
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	H3F->SetMarkerStyle(20);
	H3F->SetFillStyle(3004);
	//H3F->Draw("BOX2");  
	H3F->Draw("LEGO");  

	c4->cd(); // PMT Photocathode height (0 m)
	gPad->SetRightMargin(0.15);
	gPad->SetLeftMargin(0.15);
	gStyle->SetPalette(55);
	H2F2->Draw("ARR COLZ same");  

	for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
		Float_t xPM = v_pmt_HK[iPM].X(); 
		Float_t yPM = v_pmt_HK[iPM].Y();
		Float_t zPM = v_pmt_HK[iPM].Z();
		mar[iPM] = new TMarker(xPM, yPM, kOpenCircle);
		mar[iPM]->SetMarkerSize(15);
		mar[iPM]->SetMarkerColor(kBlack);
		mar[iPM]->Draw("same");
	}

	for (int iCoil =0; iCoil < 5; iCoil++) {
		TPolyLine *PLine= new TPolyLine(2);
		PLine->SetLineColor(kBlack);
		PLine->SetLineWidth(3);

		if ( iCoil < 3 ) {
			PLine->SetPoint(0, x_HK_start[iCoil], y_HK_start[iCoil]);
			PLine->SetPoint(1, x_HK_start[iCoil], (y_HK_start[iCoil] + Coil_Length_x));
		}
		else if ( iCoil < 5 ) {
			PLine->SetPoint(0, x_HK_start[iCoil], y_HK_start[iCoil]);
			PLine->SetPoint(1, (x_HK_start[iCoil] + Coil_Length_y), y_HK_start[iCoil]);
		}
		PLine->Draw();
	}

        c5->cd(); // PMT Photocathode height (0 m)
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gStyle->SetPalette(55);
        H2F3->Draw("ARR COLZ same");

        for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
                Float_t xPM = v_pmt_HK[iPM].X();
                Float_t yPM = v_pmt_HK[iPM].Y();
                Float_t zPM = v_pmt_HK[iPM].Z();
                mar[iPM] = new TMarker(xPM, yPM, kOpenCircle);
                mar[iPM]->SetMarkerSize(15);
                mar[iPM]->SetMarkerColor(kBlack);
                mar[iPM]->Draw("same");
        }

        for (int iCoil =0; iCoil < 5; iCoil++) {
                TPolyLine *PLine= new TPolyLine(2);
                PLine->SetLineColor(kBlack);
                PLine->SetLineWidth(3);

                if ( iCoil < 3 ) {
                        PLine->SetPoint(0, x_HK_start[iCoil], y_HK_start[iCoil]);
                        PLine->SetPoint(1, x_HK_start[iCoil], (y_HK_start[iCoil] + Coil_Length_x));
                }
                else if ( iCoil < 5 ) {
                        PLine->SetPoint(0, x_HK_start[iCoil], y_HK_start[iCoil]);
                        PLine->SetPoint(1, (x_HK_start[iCoil] + Coil_Length_y), y_HK_start[iCoil]);
                }
                PLine->Draw();
        }

	ofile->cd();

	tree->Write();
	tree_Grid->Write();

	c1->Write();
	c2->Write(); 
	c3->Write();  
	c4->Write();  
	c5->Write();

	ofile->Close();

}
