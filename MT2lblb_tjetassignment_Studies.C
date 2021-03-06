
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




void MT2lblb_tjetassignment_Studies() {

	TFile *HistogramFile = TFile::Open("MT2lblb_tjetassignment_Histos.root");


	TH1D* h_mt2lblb_notb[2];
	TH1D* h_mt2lblb_b[2];
	TH1D* h_mt2lblb_match[2];


	TH1D* h_mt2lblbcomb_notb[2];
	TH1D* h_mt2lblbcomb_b[2];
	TH1D* h_mt2lblbcomb_match[2];

	TString HistoName[2]={"_top", "_stop"};

	float mt2lblb_notb_Int[2];
	float mt2lblb_b_Int[2];
	float mt2lblb_match_Int[2];

	float mt2lblbcomb_notb_Int[2];
	float mt2lblbcomb_b_Int[2];
	float mt2lblbcomb_match_Int[2];

	for (int i  = 0; i<2; i++) {


		h_mt2lblb_b[i] = (TH1D*) HistogramFile->Get("h_mt2lblb_b" + HistoName[i]);
		h_mt2lblb_notb[i] = (TH1D*) HistogramFile->Get("h_mt2lblb_notb" + HistoName[i]);
		h_mt2lblb_match[i] = (TH1D*) HistogramFile->Get("h_mt2lblb_match" + HistoName[i]);


		h_mt2lblbcomb_b[i] = (TH1D*) HistogramFile->Get("h_mt2lblbcomb_b" + HistoName[i]);
		h_mt2lblbcomb_notb[i] = (TH1D*) HistogramFile->Get("h_mt2lblbcomb_notb" + HistoName[i]);
		h_mt2lblbcomb_match[i] = (TH1D*) HistogramFile->Get("h_mt2lblbcomb_match" + HistoName[i]);

	}

	TCanvas *CC = new TCanvas("CC", "", 1450, 800);
	CC->Divide(1, 2);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mt2lblblbtrue
	CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mt2lblblbtrue_cut
	CC2->SetGridx(); CC2->SetGridy(); 

	TCanvas *CComb = new TCanvas("CComb", "", 1450, 800);
	CComb->Divide(1, 2);
	TPad *CComb1 = (TPad*)CComb->GetPad(1); //Ahi va el h_mt2mt2lblblbtrue
	CComb1->SetGridx(); CComb1->SetGridy(); 
	TPad *CComb2 = (TPad*)CComb->GetPad(2); //Ahi va el h_mt2mt2lblblbtrue_cut
	CComb2->SetGridx(); CComb2->SetGridy(); 

	TString Option = "histosame";

	int HistoCol[2] = {4, 2};

	for (int dt=0; dt<2; dt++){

		mt2lblb_notb_Int[dt] = h_mt2lblb_notb[dt]->Integral(1, 3001);
		h_mt2lblb_notb[dt]->Scale(1./mt2lblb_notb_Int[dt]); // normalization of the histogram

		mt2lblb_b_Int[dt] = h_mt2lblb_b[dt]->Integral(1, 3001);
		h_mt2lblb_b[dt]->Scale(1./mt2lblb_b_Int[dt]); // normalization of the histogram

		mt2lblb_match_Int[dt] = h_mt2lblb_match[dt]->Integral(1, 3001);
		h_mt2lblb_match[dt]->Scale(1./mt2lblb_match_Int[dt]); // normalization of the histogram

		h_mt2lblb_b[dt] -> SetLineColor(2);
		h_mt2lblb_b[dt] -> SetLineStyle(1);
		h_mt2lblb_b[dt] -> SetLineWidth(2);

		h_mt2lblb_notb[dt] -> SetLineColor(4);
		h_mt2lblb_notb[dt] -> SetLineStyle(1);
		h_mt2lblb_notb[dt] -> SetLineWidth(2);

		h_mt2lblb_match[dt] -> SetLineColor(1);
		h_mt2lblb_match[dt] ->SetLineStyle(1);
		h_mt2lblb_match[dt] ->SetLineWidth(2);


		mt2lblbcomb_notb_Int[dt] = h_mt2lblbcomb_notb[dt]->Integral(1, 3001);
		h_mt2lblbcomb_notb[dt]->Scale(1./mt2lblbcomb_notb_Int[dt]); // normalization of the histogram

		mt2lblbcomb_b_Int[dt] = h_mt2lblbcomb_b[dt]->Integral(1, 3001);
		h_mt2lblbcomb_b[dt]->Scale(1./mt2lblbcomb_b_Int[dt]); // normalization of the histogram

		mt2lblbcomb_match_Int[dt] = h_mt2lblbcomb_match[dt]->Integral(1, 3001);
		h_mt2lblbcomb_match[dt]->Scale(1./mt2lblbcomb_match_Int[dt]); // normalization of the histogram

		h_mt2lblbcomb_b[dt] -> SetLineColor(2);
		h_mt2lblbcomb_b[dt] -> SetLineStyle(1);
		h_mt2lblbcomb_b[dt] -> SetLineWidth(2);

		h_mt2lblbcomb_notb[dt] -> SetLineColor(4);
		h_mt2lblbcomb_notb[dt] -> SetLineStyle(1);
		h_mt2lblbcomb_notb[dt] -> SetLineWidth(2);

		h_mt2lblbcomb_match[dt] -> SetLineColor(1);
		h_mt2lblbcomb_match[dt] ->SetLineStyle(1);
		h_mt2lblbcomb_match[dt] ->SetLineWidth(2);

	}






	CC->cd(1); // se pone en el TPad 1 

	h_mt2lblb_match[0]->GetXaxis()->SetRange(1, 300);
	h_mt2lblb_match[0]->GetXaxis()->SetTitle("mt2lblb top");
	h_mt2lblb_match[0]->SetTitle("mt2lblb top");
	h_mt2lblb_match[0]->DrawCopy("histo");

	h_mt2lblb_b[0]->GetXaxis()->SetRange(1, 300);
	h_mt2lblb_b[0]->GetXaxis()->SetTitle("mt2lblb top");
	h_mt2lblb_b[0]->DrawCopy(Option);


	h_mt2lblb_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mt2lblb_notb[0]->GetXaxis()->SetTitle("mt2lblb top");
	h_mt2lblb_notb[0]->DrawCopy(Option);

	TLegend *leg1 = new TLegend(0.7,0.3,0.9,0.5);
	leg1->AddEntry(h_mt2lblb_match[0],"mt2lblb_match (tjet1assignment = 2 and tjet2assignment = 2)","l");
	leg1->AddEntry(h_mt2lblb_b[0],"mt2lblb_b (tjet1assignment = 1 and tjet2assignment = 1)","l");
	leg1->AddEntry(h_mt2lblb_notb[0],"mt2lblb_notb (in any other case","l");
	leg1->Draw();


	CC->cd(2); // se pone en el TPad 1 

	h_mt2lblb_notb[1]->GetXaxis()->SetRange(1, 300);
	h_mt2lblb_notb[1]->GetXaxis()->SetTitle("mt2lblb stop");
	h_mt2lblb_notb[1]->SetTitle("mt2lblb stop");
	h_mt2lblb_notb[1]->DrawCopy("histo");

	h_mt2lblb_match[1]->GetXaxis()->SetRange(1, 300);
	h_mt2lblb_match[1]->GetXaxis()->SetTitle("mt2lblb stop");
	h_mt2lblb_match[1]->DrawCopy(Option);

	h_mt2lblb_b[1]->GetXaxis()->SetRange(1, 300);
	h_mt2lblb_b[1]->GetXaxis()->SetTitle("mt2lblb stop");
	h_mt2lblb_b[1]->DrawCopy(Option);

	TLegend *leg2 = new TLegend(0.7,0.3,0.9,0.5);
	leg2->AddEntry(h_mt2lblb_match[1],"mt2lblb_match (tjet1assignment = 2 and tjet2assignment = 2)","l");
	leg2->AddEntry(h_mt2lblb_b[1],"mt2lblb_b (tjet1assignment = 1 and tjet2assignment = 1)","l");
	leg2->AddEntry(h_mt2lblb_notb[1],"mt2lblb_notb (in any other case)","l");
	leg2->Draw();





	CComb->cd(1); // se pone en el TPad 1 


	h_mt2lblbcomb_b[0]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbcomb_b[0]->GetXaxis()->SetTitle("mt2lblbcomb top");
	h_mt2lblbcomb_b[0]->SetTitle("mt2lblbcomb top");
	h_mt2lblbcomb_b[0]->SetMaximum(0.022);
	h_mt2lblbcomb_b[0]->DrawCopy("histo");

	h_mt2lblbcomb_match[0]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbcomb_match[0]->GetXaxis()->SetTitle("mt2lblbcomb top");
	h_mt2lblbcomb_match[0]->DrawCopy(Option);


	h_mt2lblbcomb_notb[0]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbcomb_notb[0]->GetXaxis()->SetTitle("mt2lblbcomb top");
	h_mt2lblbcomb_notb[0]->DrawCopy(Option);

	TLegend *leg3 = new TLegend(0.7,0.3,0.9,0.5);
	leg3->AddEntry(h_mt2lblbcomb_match[0],"mt2lblb_match (tjet1assignment = 1 and tjet2assignment = 1)","l");
	leg3->AddEntry(h_mt2lblbcomb_b[0],"mt2lblb_b (tjet1assignment = 2 and tjet2assignment = 2)","l");
	leg3->AddEntry(h_mt2lblbcomb_notb[0],"mt2lblb_notb (in any other case","l");
	leg3->Draw();


	CComb->cd(2); // se pone en el TPad 1 

	h_mt2lblbcomb_b[1]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbcomb_b[1]->GetXaxis()->SetTitle("mt2lblbcomb stop");
	h_mt2lblbcomb_b[1]->SetTitle("mt2lblbcomb stop");
	h_mt2lblbcomb_b[1]->SetMaximum(0.022);
	h_mt2lblbcomb_b[1]->DrawCopy("histo");

	h_mt2lblbcomb_notb[1]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbcomb_notb[1]->GetXaxis()->SetTitle("mt2lblbcomb stop");
	h_mt2lblbcomb_notb[1]->DrawCopy(Option);

	h_mt2lblbcomb_match[1]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbcomb_match[1]->GetXaxis()->SetTitle("mt2lblbcomb stop");
	h_mt2lblbcomb_match[1]->DrawCopy(Option);

	TLegend *leg4 = new TLegend(0.7,0.3,0.9,0.5);
	leg4->AddEntry(h_mt2lblbcomb_match[1],"mt2lblb_match (tjet1assignment = 1 and tjet2assignment = 1)","l");
	leg4->AddEntry(h_mt2lblbcomb_b[1],"mt2lblb_b (tjet1assignment = 2 and tjet2assignment = 2)","l");
	leg4->AddEntry(h_mt2lblbcomb_notb[1],"mt2lblb_notb (in any other case","l");
	leg4->Draw();







	CC->Print("MT2lblb_tjetassignment.png");
	CComb->Print("MT2lblbcomb_tjetassignment.png");
}

