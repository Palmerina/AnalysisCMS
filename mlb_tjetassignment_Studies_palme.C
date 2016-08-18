#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <fstream>
#include <iostream>




void mlb_tjetassignment_Studies_palme() {

	TFile *HistogramFile = TFile::Open("mlb_tjetassignmentHistos_palme.root");


	TH1D* h_mlb1_notb[2];

	TH1D* h_mlb2_notb[2];


	TH1D* h_mlb1true_notb[2];

	TH1D* h_mlb2true_notb[2];

	TString HistoName[2]={"_top", "_stop"};

	float mlb1_notb_Int[2];
	float mlb2_notb_Int[2];


	float mlb1true_notb_Int[2];
	float mlb2true_notb_Int[2];

	for (int i  = 0; i<2; i++) {


		h_mlb1_notb[i] = (TH1D*) HistogramFile->Get("h_mlb1_notb" + HistoName[i]);

		h_mlb2_notb[i] = (TH1D*) HistogramFile->Get("h_mlb2_notb" + HistoName[i]);


		h_mlb1true_notb[i] = (TH1D*) HistogramFile->Get("h_mlb1true_notb" + HistoName[i]);

		h_mlb2true_notb[i] = (TH1D*) HistogramFile->Get("h_mlb2true_notb" + HistoName[i]);

	}

	TCanvas *CC = new TCanvas("CC", "", 800, 800);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetGridx(); CC1->SetGridy(); 


	TString Option = "histosame";

	int HistoCol[2] = {4, 2};

	for (int dt=0; dt<2; dt++){

		mlb1_notb_Int[dt] = h_mlb1_notb[dt]->Integral(1, 3001);
		h_mlb1_notb[dt]->Scale(1./mlb1_notb_Int[dt]); // normalization of the histogram

		mlb2_notb_Int[dt] = h_mlb2_notb[dt]->Integral(1, 3001);
		h_mlb2_notb[dt]->Scale(1./mlb2_notb_Int[dt]); // normalization of the histogram


		h_mlb1_notb[dt] -> SetLineColor(4);
		h_mlb1_notb[dt] -> SetLineStyle(2);
		h_mlb1_notb[dt] -> SetLineWidth(2);

		h_mlb2_notb[dt] -> SetLineColor(2);
		h_mlb2_notb[dt] -> SetLineStyle(2);
		h_mlb2_notb[dt] -> SetLineWidth(2);



		mlb1true_notb_Int[dt] = h_mlb1true_notb[dt]->Integral(1, 3001);
		h_mlb1true_notb[dt]->Scale(1./mlb1true_notb_Int[dt]); // normalization of the histogram

		mlb2true_notb_Int[dt] = h_mlb2true_notb[dt]->Integral(1, 3001);
		h_mlb2true_notb[dt]->Scale(1./mlb2true_notb_Int[dt]); // normalization of the histogram

		h_mlb1true_notb[dt] -> SetLineColor(4);
		h_mlb1true_notb[dt] -> SetLineStyle(1);
		h_mlb1true_notb[dt] -> SetLineWidth(2);

		h_mlb2true_notb[dt] -> SetLineColor(2);
		h_mlb2true_notb[dt] -> SetLineStyle(1);
		h_mlb2true_notb[dt] -> SetLineWidth(2);


	}






	CC->cd(1); // se pone en el TPad 1 


	h_mlb1_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mlb1_notb[0]->GetXaxis()->SetTitle("mlb top");
	h_mlb1_notb[0]->DrawCopy("histo");

	h_mlb2_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mlb2_notb[0]->GetXaxis()->SetTitle("mlb top");
	h_mlb2_notb[0]->DrawCopy(Option);

	h_mlb1true_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mlb1true_notb[0]->GetXaxis()->SetTitle("mlbtrue top");
	h_mlb1true_notb[0]->DrawCopy(Option);

	h_mlb2true_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mlb2true_notb[0]->GetXaxis()->SetTitle("mlbtrue top");
	h_mlb2true_notb[0]->DrawCopy(Option);

	TLegend *leg1 = new TLegend(0.75,0.3,0.9,0.5);
	leg1->AddEntry(h_mlb1_notb[0],"mlb1_notb","l");
	leg1->AddEntry(h_mlb2_notb[0],"mlb2_notb","l");
	leg1->AddEntry(h_mlb1true_notb[0],"mlb1true_notb","l");
	leg1->AddEntry(h_mlb2true_notb[0],"mlb2true_notb","l");
	leg1->Draw();




	CC->Print("mlb_tjetassignment_palme.png");
}

