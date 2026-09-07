#include <iostream>
#include <TCanvas.h>
#include <TH2F.h>
#include <TPolyLine.h>
#include <TLatex.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TMarker.h>
#include <TVector3.h>
#include <TTree.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include <TPaletteAxis.h>
#include <calc_bfield.h>
#define N_PM_HK 6

// Main function
int main(int argc, char *argv[]) {
	gStyle->SetPalette(1,0);
	gStyle->SetTextFont(132);

	gErrorIgnoreLevel = kFatal; // Suppressing Information message

	// Define arrays to store PMT positions and B-field vectors
	TVector3 v_pos[N_PM_HK];
	TVector3 v_dir[N_PM_HK];
	TVector3 v_mag[N_PM_HK];

	// Read B-field tree
	//
	TString fname;
	TString outname;
	if (argc == 3) {
		fname = argv[1];   // 입력받은 파일명
		outname = argv[2]; // 출력할 PDF 파일명
	} else {
		std::cout << "Usage: ./draw $B-field_root_file $output_draw_file(pdf)" << std::endl;
		return -1;
	}
	TFile* ifile = TFile::Open(fname);
	TTree* tree = (TTree*)ifile->Get("tree");
	Int_t aPM;
	TVector3* aPos = new TVector3();
	TVector3* aB = new TVector3();
	tree->SetBranchAddress("pm", &aPM);
	tree->SetBranchAddress("pos", &aPos);
	tree->SetBranchAddress("B", &aB);
	TCanvas* c2 = (TCanvas*)ifile->Get("c1");

	TTree *tree_Grid = (TTree*)ifile->Get("tree_Grid");
	TCanvas* grid2D= (TCanvas*)ifile->Get("th2d");
	TCanvas* grid2D2= (TCanvas*)ifile->Get("th2d2");
  //      TCanvas* grid2D3= (TCanvas*)ifile->Get("th2d3");
	TCanvas* grid3D= (TCanvas*)ifile->Get("th3d");
	Float_t gridaBmag = 0;
       	tree_Grid->SetBranchAddress("gridaBmag",&gridaBmag);

	// Magnetic south
	Float_t angle_south = -47.5*TMath::DegToRad();

	// Fill arrays
	for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
		tree->GetEntry(iPM);
		v_pos[aPM] = *aPos;
		v_mag[aPM] = *aB;
		v_dir[aPM] = TVector3(0, 0, +1);
	}

	// Prepare to draw
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 1600,1200);
	TH2F* hframe = new TH2F("hframe", "", 100, -180, 180, 100, -165, 195);
	hframe->SetStats(0);
	hframe->GetXaxis()->SetNdivisions(520);
	hframe->GetYaxis()->SetNdivisions(520);

	Int_t ncolors = gStyle->GetNumberOfColors();
	TH2F* h_dummy = new TH2F("h_dummy", "", 2, 0, 2, 1, 0, 1);
	//c1->SaveAs("test.pdf[");
	c1->SaveAs(outname+".pdf"+"[");
	Float_t B_draw[N_PM_HK];

	// Make histograms
	TH1F* h_perp = new TH1F("h_perp", "", 20, 0., 0.);
//	h_perp->SetCanExtend(TH1::kAllAxes);
//	h_perp->SetBins(100, 0, 10);

	// TBox objects creation
	TBox* outer = new TBox(-180, -165, 180, 195);
	TBox* aus_sys = new TBox(-149, -145.7, -3, 62.3);
	TBox* aus_mon = new TBox(-103, 67.3, -3, 167.3);                                                                     
	TBox* kor_sys1 = new TBox(3, -125.7, 103, -41.7);
	TBox* kor_sys2 = new TBox(3, -41.7, 103, 42.3); 
	TBox* kor_mon = new TBox(3, 67.3, 103, 167.3);

	aus_sys->SetFillStyle(0);
	aus_mon->SetFillStyle(0);
	kor_sys1->SetFillStyle(0);
	kor_sys2->SetFillStyle(0);
	kor_mon->SetFillStyle(0);

	TLatex *latex = new TLatex();
	latex->SetTextSize(0.05);
	latex->SetTextFont(132);

	double x_min_aus_sys = -149.;
	double x_max_aus_sys = -3.;
	double y_min_aus_sys = -145.7;
	double y_max_aus_sys = 62.3;

	double x_min_aus_mon = -103.;
	double x_max_aus_mon = -3.;
	double y_min_aus_mon = 67.3;
	double y_max_aus_mon = 167.3;

	double x_min_kor_sys1 = 3.;
	double x_max_kor_sys1 = 103.;
	double y_min_kor_sys1 = -125.7;
	double y_max_kor_sys1 = -41.7;

	double x_min_kor_sys2 = 3.;
	double x_max_kor_sys2 = 103.;
	double y_min_kor_sys2 = -41.7;
	double y_max_kor_sys2 = 42.3;

	double x_min_kor_mon = 3.;
	double x_max_kor_mon = 103.;
	double y_min_kor_mon = 67.3;
	double y_max_kor_mon = 167.3;

	// Make markers 
	TMarker* mar[N_PM_HK];
	for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
		Float_t xPM = v_pos[iPM].X()/10.; // [mm] -> [cm]
		Float_t yPM = v_pos[iPM].Y()/10.;
		Float_t zPM = v_pos[iPM].Z()/10.;

		mar[iPM] = new TMarker(xPM, yPM, kFullCircle);
		mar[iPM]->SetMarkerSize(15);

		double distance_x = 0.0;
		double distance_y = 0.0;

		if (iPM == 0) {
			// Distance from aus_sys for PMTs 0 and 1
			distance_x = fabs(xPM - x_min_aus_sys);  // 정면 거리
			distance_y = fabs(yPM - y_min_aus_sys);  // 측면 거리
		} else if (iPM == 1) {
			distance_x = fabs(xPM - x_max_aus_sys);  // 정면 거리
			distance_y = fabs(yPM - y_max_aus_sys);  // 측면 거리
		} else if (iPM == 2) {
			// Distance from aus_mon for PMT 2
			distance_x = fabs(xPM - x_max_aus_mon);
			distance_y = fabs(yPM - y_max_aus_mon);
		} else if (iPM == 3) {
			// Distance from kor_sys1 for PMT 3
			distance_x = fabs(xPM - x_max_kor_sys1);
			distance_y = fabs(yPM - y_max_kor_sys1);
		} else if (iPM == 4) {
			// Distance from kor_sys2 for PMT 4
			distance_x = fabs(xPM - x_max_kor_sys2);
			distance_y = fabs(yPM - y_max_kor_sys2);
		} else if (iPM == 5) {
			// Distance from kor_mon for PMT 5
			distance_x = fabs(xPM - x_max_kor_mon);
			distance_y = fabs(yPM - y_max_kor_mon);
		}

		// Print or store distances as needed
//		std::cout << "Distance for PMT " << iPM << ": Front = " << distance_x << " cm, Side = " << distance_y << " cm" << std::endl;

	}

	// Loop for draw patterns
	for (Int_t iDraw=0; iDraw<6; iDraw++) {

		// Find max B field
		Float_t B_max = -1e5, B_min = 1e5;
		for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
			if (iDraw ==0)     B_draw[iPM] = v_mag[iPM].X();
			else if (iDraw==1) B_draw[iPM] = v_mag[iPM].Y();
			else if (iDraw==2) B_draw[iPM] = v_mag[iPM].Z();
			else if (iDraw==3) B_draw[iPM] = v_mag[iPM].Mag();
			else if (iDraw==4) B_draw[iPM] = v_mag[iPM] * v_dir[iPM]; // Parralel to PMT direction
			else if (iDraw==5) B_draw[iPM] = sqrt(v_mag[iPM].Mag2() - pow(v_mag[iPM] * v_dir[iPM], 2)); // Perpendicular to PMT direction
			if (B_draw[iPM] > B_max) {
				B_max = B_draw[iPM];
			}
			if (B_draw[iPM] < B_min) {
				B_min = B_draw[iPM];
			}
		}

		// Make color palette
		h_dummy->SetBinContent(1, 1, B_min);
		h_dummy->SetBinContent(2, 1, B_max);
		h_dummy->Draw("COLZ");
		c1->Update();
		TPaletteAxis* palette = (TPaletteAxis*)h_dummy->GetListOfFunctions()->FindObject("palette");

		// Set title
		if (iDraw ==0)     hframe->SetTitle("X");
		else if (iDraw==1) hframe->SetTitle("Y");
		else if (iDraw==2) hframe->SetTitle("Z");
		else if (iDraw==3) hframe->SetTitle("Mag");
		else if (iDraw==4) hframe->SetTitle("Parallel");
		else if (iDraw==5) hframe->SetTitle("Perpendicular");

		// Draw frame
		hframe->Draw();
		palette->Draw("same");
		outer->Draw("L same");
		aus_sys->Draw("L same");
		aus_mon->Draw("L same");
		kor_sys1->Draw("L same");
		kor_sys2->Draw("L same");
		kor_mon->Draw("L same");

		latex->SetTextAlign(12);
		latex->DrawLatex(-100, -45, "Aus. System");
		latex->DrawLatex(-100, 80, "Aus. Mon.");
		latex->DrawLatex(110, -100, "Kor. Sys.1");
		latex->DrawLatex(110, 0, "Kor. Sys.2");
		latex->DrawLatex(110, 100, "Kor. Mon.");

		// Draw B-field of each PMT
		for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
			Int_t color_this = gStyle->GetColorPalette((int)((B_draw[iPM]-B_min)*(ncolors-1.)/(B_max-B_min)));
			mar[iPM]->SetMarkerColor(color_this);
			mar[iPM]->Draw("same");
		}

		for (int iCoil =0; iCoil < 5; iCoil++)
		{
			TPolyLine *PLine= new TPolyLine(2);

			if ( iCoil < 3 ){
				PLine->SetPoint(0, x_HK_start[iCoil]/10., y_HK_start[iCoil]/10.);
				PLine->SetPoint(1, x_HK_start[iCoil]/10., (y_HK_start[iCoil] + 3532)/10.);
				PLine->SetLineColor(kRed);
				PLine->SetLineWidth(2);
			}
			else if ( iCoil < 5 ){
				PLine->SetPoint(0, x_HK_start[iCoil]/10., y_HK_start[iCoil]/10.);
				PLine->SetPoint(1, (x_HK_start[iCoil] + 3532)/10., y_HK_start[iCoil]/10.);
				PLine->SetLineColor(kBlue);
				PLine->SetLineWidth(2);
			}
			PLine->Draw("same");
		}
		c1->Modified();
		// Fill histogram
		if (iDraw == 5) {
			for (Int_t iPM=0; iPM<N_PM_HK; iPM++) {
				h_perp->Fill(B_draw[iPM]);
			}
		}
		c1->Update();
		// Save canvas
		//c1->SaveAs("test.pdf");
		c1->SaveAs(outname+".pdf");

	}
	// Draw B-field histogram
	tree->Draw("B.Mag()");
	c1->SaveAs(outname+".pdf");
	h_perp->SetTitle(";B_{perp} [mG]");
	h_perp->Draw();
	c1->SaveAs(outname+".pdf");
	grid2D->SaveAs(outname+".pdf");
	grid2D2->SaveAs(outname+".pdf");
//        grid2D3->SaveAs(outname+".pdf");
	grid3D->SaveAs(outname+".pdf");
	tree_Grid->Draw("gridaBmag");
	c2->SaveAs(outname+".pdf");
	c1->SaveAs(outname + ".pdf" + "]");  // 출력 파일명으로 저장
}
