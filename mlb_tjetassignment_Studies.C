
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




void mlb_tjetassignment_Studies() {

	TFile *HistogramFile = TFile::Open("mlb_tjetassignmentHistos.root");


	TH1D* h_mlb1_notb[2];
	TH1D* h_mlb1_b[2];
	TH1D* h_mlb1_match[2];

	TH1D* h_mlb2_notb[2];
	TH1D* h_mlb2_b[2];
	TH1D* h_mlb2_match[2];

	TString HistoName[2]={"_top", "_stop"};



	for (int i  = 0; i<2; i++) {


		h_mlb1_b[i] = (TH1D*) HistogramFile->Get("h_mlb1_b" + HistoName[i]);
		h_mlb1_notb[i] = (TH1D*) HistogramFile->Get("h_mlb1_notb" + HistoName[i]);
		h_mlb1_match[i] = (TH1D*) HistogramFile->Get("h_mlb1_match" + HistoName[i]);

		h_mlb2_b[i] = (TH1D*) HistogramFile->Get("h_mlb2_b" + HistoName[i]);
		h_mlb2_notb[i] = (TH1D*) HistogramFile->Get("h_mlb2_notb" + HistoName[i]);
		h_mlb2_match[i] = (TH1D*) HistogramFile->Get("h_mlb2_match" + HistoName[i]);

	}

	TCanvas *CC = new TCanvas("CC", "", 1450, 800);
	CC->Divide(1, 2);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetLogy();CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2->SetLogy();CC2->SetGridx(); CC2->SetGridy(); 

	TString Option = "histo";

	int HistoCol[2] = {4, 2};



		h_mlb1_b[0] -> SetLineColor(HistoCol[0]);
		h_mlb1_b[0] -> SetLineStyle(1);
		h_mlb1_b[0] -> SetLineWidth(2);

		h_mlb2_b[0] -> SetLineColor(HistoCol[0]);
		h_mlb2_b[0] -> SetLineStyle(1);
		h_mlb2_b[0] -> SetLineWidth(2);

		h_mlb1_notb[0] -> SetLineColor(HistoCol[1]);
		h_mlb1_notb[0] -> SetLineStyle(1);
		h_mlb1_notb[0] -> SetLineWidth(2);
		
		h_mlb2_notb[0] -> SetLineColor(HistoCol[1]);
		h_mlb2_notb[0] -> SetLineStyle(1);
		h_mlb2_notb[0] -> SetLineWidth(2);
		
		h_mlb1_match[0] -> SetLineColor(1);
		h_mlb1_match[0] ->SetLineStyle(1);
		h_mlb1_match[0] ->SetLineWidth(2);

		h_mlb2_match[0] -> SetLineColor(1);
		h_mlb2_match[0] ->SetLineStyle(1);
		h_mlb2_match[0] ->SetLineWidth(2);


		h_mlb1_b[1] -> SetLineColor(HistoCol[0]);
		h_mlb1_b[1] -> SetLineStyle(1);
		h_mlb1_b[1] -> SetLineWidth(2);

		h_mlb2_b[1] -> SetLineColor(HistoCol[0]);
		h_mlb2_b[1] -> SetLineStyle(1);
		h_mlb2_b[1] -> SetLineWidth(2);

		h_mlb1_notb[1] -> SetLineColor(HistoCol[1]);
		h_mlb1_notb[1] -> SetLineStyle(1);
		h_mlb1_notb[1] -> SetLineWidth(2);
		
		h_mlb2_notb[1] -> SetLineColor(HistoCol[1]);
		h_mlb2_notb[1] -> SetLineStyle(1);
		h_mlb2_notb[1] -> SetLineWidth(2);
		
		h_mlb1_match[1] -> SetLineColor(1);
		h_mlb1_match[1] ->SetLineStyle(1);
		h_mlb1_match[1] ->SetLineWidth(2);

		h_mlb2_match[1] -> SetLineColor(1);
		h_mlb2_match[1] ->SetLineStyle(1);
		h_mlb2_match[1] ->SetLineWidth(2);




for (int dt = 0; dt<2; dt++) { // stop or top



	CC->cd(1); // se pone en el TPad 1 

	h_mlb1_b[0]->GetXaxis()->SetRange(1, 300);
	h_mlb1_b[0]->GetXaxis()->SetTitle("mlb top");
	h_mlb1_b[0]->DrawCopy(Option);

	h_mlb2_b[0]->GetXaxis()->SetRange(1, 300);
	h_mlb2_b[0]->GetXaxis()->SetTitle("mlb top");
	h_mlb2_b[0]->DrawCopy(Option);

	h_mlb1_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mlb1_notb[0]->GetXaxis()->SetTitle("mlb top");
	h_mlb1_notb[0]->DrawCopy(Option);

	h_mlb2_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mlb2_notb[0]->GetXaxis()->SetTitle("mlb top");
	h_mlb2_notb[0]->DrawCopy(Option);

	h_mlb1_match[0]->GetXaxis()->SetRange(1, 300);
	h_mlb1_match[0]->GetXaxis()->SetTitle("mlb top");
	h_mlb1_match[0]->DrawCopy(Option);

	h_mlb2_match[0]->GetXaxis()->SetRange(1, 300);
	h_mlb2_match[0]->GetXaxis()->SetTitle("mlb top");
	h_mlb2_match[0]->DrawCopy(Option);

	TLegend *leg1 = new TLegend(0.3,0.1,0.4,0.4);
	leg1->AddEntry(h_mlb1_b[0],"mlb1 tjet1assignment = 0","l");
	leg1->AddEntry(h_mlb1_notb[0],"mlb1 tjet1assignment = 1","l");
	leg1->AddEntry(h_mlb1_match[0],"mlb1 tjet1assignment = 2","l");
	leg1->AddEntry(h_mlb2_b[0],"mlb2 tjet2assignment = 0","l");
	leg1->AddEntry(h_mlb2_notb[0],"mlb2 tjet2assignment = 2","l");
	leg1->AddEntry(h_mlb2_match[0],"mlb2 tjet2assignment = 2","l");
	leg1->Draw();


	CC->cd(2); // se pone en el TPad 1 

	h_mlb1_b[1]->GetXaxis()->SetRange(1, 300);
	h_mlb1_b[1]->GetXaxis()->SetTitle("mlb stop");
	h_mlb1_b[1]->DrawCopy(Option);

	h_mlb2_b[1]->GetXaxis()->SetRange(1, 300);
	h_mlb2_b[1]->GetXaxis()->SetTitle("mlb stop");
	h_mlb2_b[1]->DrawCopy(Option);

	h_mlb1_notb[1]->GetXaxis()->SetRange(1, 300);
	h_mlb1_notb[1]->GetXaxis()->SetTitle("mlb stop");
	h_mlb1_notb[1]->DrawCopy(Option);

	h_mlb2_notb[1]->GetXaxis()->SetRange(1, 300);
	h_mlb2_notb[1]->GetXaxis()->SetTitle("mlb stop");
	h_mlb2_notb[1]->DrawCopy(Option);

	h_mlb1_match[1]->GetXaxis()->SetRange(1, 300);
	h_mlb1_match[1]->GetXaxis()->SetTitle("mlb stop");
	h_mlb1_match[1]->DrawCopy(Option);

	h_mlb2_match[1]->GetXaxis()->SetRange(1, 300);
	h_mlb2_match[1]->GetXaxis()->SetTitle("mlb stop");
	h_mlb2_match[1]->DrawCopy(Option);

	TLegend *leg2 = new TLegend(0.3,0.1,0.4,0.4);
	leg2->AddEntry(h_mlb1_b[1],"mlb1 tjet1assignment = 0","l");
	leg2->AddEntry(h_mlb1_notb[1],"mlb1 tjet1assignment = 1","l");
	leg2->AddEntry(h_mlb1_match[1],"mlb1 tjet1assignment = 2","l");
	leg2->AddEntry(h_mlb2_b[1],"mlb2 tjet2assignment = 0","l");
	leg2->AddEntry(h_mlb2_notb[1],"mlb2 tjet2assignment = 2","l");
	leg2->AddEntry(h_mlb2_match[1],"mlb2 tjet2assignment = 2","l");
	leg2->Draw();





	Option= "histosame";
}   


CC->Print("mlb_tjetassignment.png");
}

