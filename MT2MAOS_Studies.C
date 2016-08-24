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




void MT2MAOS_Studies() {

	TFile *HistogramFile = TFile::Open("MT2MAOS_Histos.root");


	TH2D* h_mt2ll_mt2llMAOS;
	TH2D* h_mt2ll_mt2lblbMAOStrue;
	TH2D* h_mt2lblbtrue_mt2lblbMAOStrue;


	TH2D* h_mt2lblb_mt2lblbMAOS;
	TH2D* h_mt2ll_mt2lblbMAOS[2];
	TH2D* h_mt2ll_mt2lblbMAOS_cut[2];

	TString HistoName[2]={"_top", "_stop"};




	h_mt2ll_mt2lblbMAOStrue = (TH2D*) HistogramFile->Get("h_mt2ll_mt2lblbMAOStrue" + HistoName[0]);
	h_mt2ll_mt2llMAOS = (TH2D*) HistogramFile->Get("h_mt2ll_mt2llMAOS" + HistoName[0]);
	h_mt2lblbtrue_mt2lblbMAOStrue = (TH2D*) HistogramFile->Get("h_mt2lblbtrue_mt2lblbMAOStrue" + HistoName[0]);
	h_mt2lblb_mt2lblbMAOS = (TH2D*) HistogramFile->Get("h_mt2lblb_mt2lblbMAOS" + HistoName[0]);

	for (int dt = 0; dt < 2; dt++) {
		h_mt2ll_mt2lblbMAOS[dt] = (TH2D*) HistogramFile->Get("h_mt2ll_mt2lblbMAOS" + HistoName[dt]);
		h_mt2ll_mt2lblbMAOS_cut[dt] = (TH2D*) HistogramFile->Get("h_mt2ll_mt2lblbMAOS_cut" + HistoName[dt]);
	
	}


	TCanvas *CC = new TCanvas("CC", "", 1000, 800);
	CC->Divide(2, 2);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2->SetGridx(); CC2->SetGridy(); 
	TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mlblbtrue_cut
	CC3->SetGridx(); CC3->SetGridy(); 
	TPad *CC4 = (TPad*)CC->GetPad(4); //Ahi va el h_mt2mlblbtrue
	CC4->SetGridx(); CC4->SetGridy(); 


	TCanvas *CCstop = new TCanvas("CCstop", "", 1000, 800);
	CCstop->Divide(2, 2);
	TPad *CCstop1 = (TPad*)CCstop->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CCstop1->SetGridx(); CCstop1->SetGridy(); 
	TPad *CCstop2 = (TPad*)CCstop->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CCstop2->SetGridx(); CCstop2->SetGridy(); 
	TPad *CCstop3 = (TPad*)CCstop->GetPad(3); //Ahi va el h_mt2mlblbtrue_cut
	CCstop3->SetGridx(); CCstop3->SetGridy(); 
	TPad *CCstop4 = (TPad*)CCstop->GetPad(4); //Ahi va el h_mt2mlblbtrue
	CCstop4->SetGridx(); CCstop4->SetGridy(); 



	CC->cd(1); // se pone en el TPad 1 

	h_mt2lblbtrue_mt2lblbMAOStrue->GetXaxis()->SetRangeUser(1, 600);
	h_mt2lblbtrue_mt2lblbMAOStrue->GetXaxis()->SetTitle("mt2lblbtrue");
	h_mt2lblbtrue_mt2lblbMAOStrue->GetYaxis()->SetRangeUser(1, 600);
	h_mt2lblbtrue_mt2lblbMAOStrue->GetYaxis()->SetTitle("mt2lblbMAOStrue");
	h_mt2lblbtrue_mt2lblbMAOStrue->SetTitle("mt2lblbtrue_mt2lblbMAOStrue 2D  top");
	h_mt2lblbtrue_mt2lblbMAOStrue->DrawCopy("box");
	



	CC->cd(2); // se pone en el TPad 1 

	h_mt2ll_mt2lblbMAOStrue->GetXaxis()->SetRangeUser(1, 300);
	h_mt2ll_mt2lblbMAOStrue->GetXaxis()->SetTitle("mt2ll");
	h_mt2ll_mt2lblbMAOStrue->GetYaxis()->SetRangeUser(1, 600);
	h_mt2ll_mt2lblbMAOStrue->GetYaxis()->SetTitle("mt2lblbMAOStrue");
	h_mt2ll_mt2lblbMAOStrue->SetTitle("mt2ll_mt2lblbMAOStrue 2D top");
	h_mt2ll_mt2lblbMAOStrue->DrawCopy("box");


	CC->cd(3); // se pone en el TPad 1 


	h_mt2ll_mt2llMAOS->GetXaxis()->SetRangeUser(1, 300);
	h_mt2ll_mt2llMAOS->GetXaxis()->SetTitle("mt2ll");
	h_mt2ll_mt2llMAOS->GetYaxis()->SetRangeUser(1, 300);
	h_mt2ll_mt2llMAOS->GetYaxis()->SetTitle("mt2llMAOS");
	h_mt2ll_mt2llMAOS->SetTitle("mt2ll_mt2llMAOS top");
	h_mt2ll_mt2llMAOS->DrawCopy("box");


	CC->cd(4); // se pone en el TPad 1 

	h_mt2lblb_mt2lblbMAOS->GetXaxis()->SetRangeUser(1, 600);
	h_mt2lblb_mt2lblbMAOS->GetXaxis()->SetTitle("mt2lblb");
	h_mt2lblb_mt2lblbMAOS->GetYaxis()->SetRangeUser(1, 600);
	h_mt2lblb_mt2lblbMAOS->GetYaxis()->SetTitle("mt2lblbMAOS");
	h_mt2lblb_mt2lblbMAOS->SetTitle("mt2lblb_mt2lblbMAOS top");
	h_mt2lblb_mt2lblbMAOS->DrawCopy("box");

	for (int dt = 0; dt < 2; dt++) {

		CCstop->cd(dt+1); // se pone en el TPad 1 


		h_mt2ll_mt2lblbMAOS_cut[dt]->GetXaxis()->SetRangeUser(1, 300);
		h_mt2ll_mt2lblbMAOS_cut[dt]->GetXaxis()->SetTitle("mt2ll");
		h_mt2ll_mt2lblbMAOS_cut[dt]->GetYaxis()->SetRangeUser(1, 600);
		h_mt2ll_mt2lblbMAOS_cut[dt]->GetYaxis()->SetTitle("mt2lblbMAOS_cut");
		h_mt2ll_mt2lblbMAOS_cut[dt]->SetTitle("mt2ll_mt2lblbMAOS" + HistoName[dt] + " with mlb <= 160");
		h_mt2ll_mt2lblbMAOS_cut[dt]->DrawCopy("box");


		CCstop->cd(dt+3); // se pone en el TPad 1 

		h_mt2ll_mt2lblbMAOS[dt]->GetXaxis()->SetRangeUser(1, 300);
		h_mt2ll_mt2lblbMAOS[dt]->GetXaxis()->SetTitle("mt2ll");
		h_mt2ll_mt2lblbMAOS[dt]->GetYaxis()->SetRangeUser(1, 600);
		h_mt2ll_mt2lblbMAOS[dt]->GetYaxis()->SetTitle("mt2lblbMAOS");
		h_mt2ll_mt2lblbMAOS[dt]->SetTitle("mt2ll_mt2lblbMAOS" + HistoName[dt]);
		h_mt2ll_mt2lblbMAOS[dt]->DrawCopy("box");

	}


	CC->Print("MT2MAOS_Studies.png");
	CCstop->Print("MT2MAOS_Studies2.png");
}

