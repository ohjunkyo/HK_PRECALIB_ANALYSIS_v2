#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TColor.h"

void SetAxisStyle(TAxis* axis) {
	if (!axis) return;
	axis->SetTitleFont(132);   
	axis->SetLabelFont(132);   
	axis->SetTitleSize(0.04);
	axis->SetLabelSize(0.04);
	axis->CenterTitle(true);
//	axis->SetTitleOffset(1.0);
}

void SetHistStyle(TH1* hist) {
	if (!hist) return;

	hist->SetStats(0);
	gPad->SetLeftMargin(0.15);
	SetAxisStyle(hist->GetXaxis());
	SetAxisStyle(hist->GetYaxis());
}


void SetHist2DStyle(TH2* hist) {
	if (!hist) return;

	hist->SetStats(0);
	SetAxisStyle(hist->GetXaxis());
	SetAxisStyle(hist->GetYaxis());
	SetAxisStyle(hist->GetZaxis());
}
void SetGraphStyle(TGraph* graph) {
	if (!graph) return;

	graph->SetLineWidth(2);
	graph->SetLineColor(kBlack);
	graph->SetMarkerStyle(20);
	graph->SetMarkerSize(1.0);
	graph->SetMarkerColor(kBlack);

	SetAxisStyle(graph->GetXaxis());
	SetAxisStyle(graph->GetYaxis());
}

void SetErrGraphStyle(TGraphErrors* graph) {
	if (!graph) return;

	SetGraphStyle(graph);
	graph->SetFillColorAlpha(kGray, 0.5); 
	graph->SetFillStyle(1001);            
}
