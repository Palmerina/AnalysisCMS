#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
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




void mlb2D_tjetassignment_Studies() {

	TFile *HistogramFile = TFile::Open("mlb2D_tjetassignmentHistos.root");


	TH2D* h_mlb_notb[2];
	TH2D* h_mlb_b[2];
	TH2D* h_mlb_match[2];


	TH2D* h_mlbcomb_notb[2];
	TH2D* h_mlbcomb_b[2];
	TH2D* h_mlbcomb_match[2];

	TString HistoName[2]={"_top", "_stop"};




	h_mlb_b[0] = (TH2D*) HistogramFile->Get("h_mlb_b" + HistoName[0]);
	h_mlb_notb[0] = (TH2D*) HistogramFile->Get("h_mlb_notb" + HistoName[0]);
	h_mlb_match[0] = (TH2D*) HistogramFile->Get("h_mlb_match" + HistoName[0]);

	h_mlbcomb_b[0] = (TH2D*) HistogramFile->Get("h_mlbcomb_b" + HistoName[0]);
	h_mlbcomb_notb[0] = (TH2D*) HistogramFile->Get("h_mlbcomb_notb" + HistoName[0]);
	h_mlbcomb_match[0] = (TH2D*) HistogramFile->Get("h_mlbcomb_match" + HistoName[0]);


	TCanvas *CC = new TCanvas("CC", "", 1000, 800);
	CC->Divide(1, 3);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2->SetGridx(); CC2->SetGridy(); 
	TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mlblbtrue_cut
	CC3->SetGridx(); CC3->SetGridy(); 


	TCanvas *CComb = new TCanvas("CComb", "", 1000, 800);
	CComb->Divide(1, 3);
	TPad *CComb1 = (TPad*)CComb->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CComb1->SetGridx(); CComb1->SetGridy(); 
	TPad *CComb2 = (TPad*)CComb->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CComb2->SetGridx(); CComb2->SetGridy(); 
	TPad *CComb3 = (TPad*)CComb->GetPad(3); //Ahi va el h_mt2mlblbtrue_cut
	CComb3->SetGridx(); CComb3->SetGridy(); 





	CC->cd(1); // se pone en el TPad 1 

	h_mlb_match[0]->GetXaxis()->SetRange(1, 300);
	h_mlb_match[0]->GetXaxis()->SetTitle("mlb1_match");
	h_mlb_match[0]->GetYaxis()->SetRange(1, 300);
	h_mlb_match[0]->GetYaxis()->SetTitle("mlb2_match");
	h_mlb_match[0]->SetTitle("mlb_match 2D  top");
	h_mlb_match[0]->DrawCopy("box");




	CC->cd(2); // se pone en el TPad 1 

	h_mlb_b[0]->GetXaxis()->SetRange(1, 300);
	h_mlb_b[0]->GetXaxis()->SetTitle("mlb1_b");
	h_mlb_b[0]->GetYaxis()->SetRange(1, 300);
	h_mlb_b[0]->GetYaxis()->SetTitle("mlb2_b");
	h_mlb_b[0]->SetTitle("mlb_b 2D top");
	h_mlb_b[0]->DrawCopy("box");


	CC->cd(3); // se pone en el TPad 1 


	h_mlb_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mlb_notb[0]->GetXaxis()->SetTitle("mlb1_notb");
	h_mlb_notb[0]->GetYaxis()->SetRange(1, 300);
	h_mlb_notb[0]->GetYaxis()->SetTitle("mlb2_notb");
	h_mlb_notb[0]->SetTitle("mlb_notb top");
	h_mlb_notb[0]->DrawCopy("box");


	CComb->cd(1); // se pone en el TPad 1 

	h_mlbcomb_match[0]->GetXaxis()->SetRange(1, 300);
	h_mlbcomb_match[0]->GetXaxis()->SetTitle("mlb1comb_match");
	h_mlbcomb_match[0]->GetYaxis()->SetRange(1, 300);
	h_mlbcomb_match[0]->GetYaxis()->SetTitle("mlb2comb_match");
	h_mlbcomb_match[0]->SetTitle("mlbcomb_match top");
	h_mlbcomb_match[0]->DrawCopy("box");



	CComb->cd(2); // se pone en el TPad 1 

	h_mlbcomb_b[0]->GetXaxis()->SetRange(1, 300);
	h_mlbcomb_b[0]->GetXaxis()->SetTitle("mlb1comb_b");
	h_mlbcomb_b[0]->GetYaxis()->SetRange(1, 300);
	h_mlbcomb_b[0]->GetYaxis()->SetTitle("mlb1comb_b");
	h_mlbcomb_b[0]->SetTitle("mlbcomb_b top");
	h_mlbcomb_b[0]->DrawCopy("box");


	CComb->cd(3); // se pone en el TPad 1 

	h_mlbcomb_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mlbcomb_notb[0]->GetXaxis()->SetTitle("mlb1comb_notb");
	h_mlbcomb_notb[0]->GetYaxis()->SetRange(1, 300);
	h_mlbcomb_notb[0]->GetYaxis()->SetTitle("mlb1comb_notb");
	h_mlbcomb_notb[0]->SetTitle("mlbcomb_notb top");
	h_mlbcomb_notb[0]->DrawCopy("box");




	CC->Print("mlb2D_tjetassignment.png");
	CComb->Print("mlbcomb2D_tjetassignment.png");
}

