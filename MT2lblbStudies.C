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




void MT2lblbStudies() {

	TFile *HistogramFile = TFile::Open("MT2lblbHistos.root");


	TH1D* h_mt2lblbtrue[2];
	TH1D* h_mt2lblbtrue_cut[2];
	TH1D* h_mt2lblbtrueInteg[2];
	TH1D* h_mt2lblbtrueInteg_cut[2];
	TH1D* h_mt2lblbtrueSignif;
	TH1D* h_mt2lblbtrueSignif_cut;

	TH1D* h_mt2lblb[2];
	TH1D* h_mt2lblb_cut[2];
	TH1D* h_mt2lblbInteg[2];
	TH1D* h_mt2lblbInteg_cut[2];
	TH1D* h_mt2lblbSignif;
	TH1D* h_mt2lblbSignif_cut;

	TString HistoName[2]={"_top", "_stop"};


	float mt2lblbtrue_Int[2];
	float mt2lblbtrue_cut_Int[2];

	float mt2lblb_Int[2];
	float mt2lblb_cut_Int[2];

	for (int i  = 0; i<2; i++) {


		h_mt2lblbtrue[i] = (TH1D*) HistogramFile->Get("h_mt2lblbtrue" + HistoName[i]);
		h_mt2lblbtrue_cut[i] = (TH1D*) HistogramFile->Get("h_mt2lblbtrue_cut" + HistoName[i] + "_135");

		h_mt2lblbtrueInteg[i] = new TH1D("h_mt2lblbtrueInteg" + HistoName[i], "h_mt2lblbtrueInteg", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbtrueInteg_cut[i] = new TH1D("h_mt2lblbtrueInteg_cut" + HistoName[i] + "_135", "h_mt2lblbtrueInteg_cut", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbtrueSignif = new TH1D("h_mt2lblbtrueSignif" + HistoName[i], "h_mt2lblbtrueSignif", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbtrueSignif_cut = new TH1D("h_mt2lblbtrueSignif_cut" + HistoName[i] + "_135", "h_mt2lblbtrueSignif_cut", 3000, 0, 3000); //2 histos para top y stop


		h_mt2lblb[i] = (TH1D*) HistogramFile->Get("h_mt2lblb" + HistoName[i]);
		h_mt2lblb_cut[i] = (TH1D*) HistogramFile->Get("h_mt2lblb_cut" + HistoName[i] + "_135");

		h_mt2lblbInteg[i] = new TH1D("h_mt2lblbInteg" + HistoName[i], "h_mt2lblbInteg", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbInteg_cut[i] = new TH1D("h_mt2lblbInteg_cut" + HistoName[i] + "_135", "h_mt2lblbInteg_cut", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbSignif = new TH1D("h_mt2lblbSignif" + HistoName[i], "h_mt2lblbSignif", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbSignif_cut = new TH1D("h_mt2lblbSignif_cut" + HistoName[i] + "_135", "h_mt2lblbSignif_cut", 3000, 0, 3000); //2 histos para top y stop

	}

	TCanvas *CCtrue = new TCanvas("CCtrue", "", 1450, 800);
	CCtrue->Divide(1, 3);
	TPad *CC1true = (TPad*)CCtrue->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1true->SetLogy();CC1true->SetGridx(); CC1true->SetGridy(); 
	TPad *CC2true = (TPad*)CCtrue->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2true->SetLogy();CC2true->SetGridx(); CC2true->SetGridy(); 
	TPad *CC3true = (TPad*)CCtrue->GetPad(3); //Ahi va el h_mt2mlblbtrue_cut
	CC3true->SetGridx(); CC3true->SetGridy(); 


	TCanvas *CC = new TCanvas("CC", "", 1450, 800);
	CC->Divide(1, 3);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetLogy();CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2->SetLogy();CC2->SetGridx(); CC2->SetGridy(); 
	TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mlblbtrue_cut
	CC3->SetGridx(); CC3->SetGridy(); 

	TString Option = "histo";

	int HistoCol[2] = {4, 2};



	for (int dt  = 0; dt<2; dt++) {


		mt2lblbtrue_Int[dt] = h_mt2lblbtrue[dt]->Integral(1, 3001);
		h_mt2lblbtrue[dt]->Scale(1./mt2lblbtrue_Int[dt]); // normalization of the histogram

		mt2lblbtrue_cut_Int[dt] = h_mt2lblbtrue_cut[dt]->Integral(1, 3001);
		h_mt2lblbtrue_cut[dt]->Scale(1./mt2lblbtrue_Int[dt]); // normalization of the histogram


		h_mt2lblbtrue[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbtrue_cut[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbtrueInteg[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbtrueInteg_cut[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbtrueSignif -> SetLineColor(1);
		h_mt2lblbtrueSignif_cut -> SetLineColor(1);

		h_mt2lblbtrue[dt] -> SetLineStyle(1);
		h_mt2lblbtrue_cut[dt] ->SetLineStyle(2);
		h_mt2lblbtrueInteg[dt] -> SetLineStyle(1);
		h_mt2lblbtrueInteg_cut[dt] ->SetLineStyle(2);
		h_mt2lblbtrueSignif -> SetLineStyle(1);
		h_mt2lblbtrueSignif_cut -> SetLineStyle(2);

		h_mt2lblbtrue[dt] -> SetLineWidth(2);
		h_mt2lblbtrue_cut[dt] ->SetLineWidth(2);
		h_mt2lblbtrueInteg[dt] -> SetLineWidth(2);
		h_mt2lblbtrueInteg_cut[dt] ->SetLineWidth(2);
		h_mt2lblbtrueSignif -> SetLineWidth(2);
		h_mt2lblbtrueSignif_cut -> SetLineWidth(2);


		mt2lblb_Int[dt] = h_mt2lblb[dt]->Integral(1, 3001);
		h_mt2lblb[dt]->Scale(1./mt2lblb_Int[dt]); // normalization of the histogram

		mt2lblb_cut_Int[dt] = h_mt2lblb_cut[dt]->Integral(1, 3001);
		h_mt2lblb_cut[dt]->Scale(1./mt2lblb_Int[dt]); // normalization of the histogram


		h_mt2lblb[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblb_cut[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbInteg[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbInteg_cut[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbSignif -> SetLineColor(1);
		h_mt2lblbSignif_cut -> SetLineColor(1);

		h_mt2lblb[dt] -> SetLineStyle(1);
		h_mt2lblb_cut[dt] ->SetLineStyle(2);
		h_mt2lblbInteg[dt] -> SetLineStyle(1);
		h_mt2lblbInteg_cut[dt] ->SetLineStyle(2);
		h_mt2lblbSignif -> SetLineStyle(1);
		h_mt2lblbSignif_cut -> SetLineStyle(2);

		h_mt2lblb[dt] -> SetLineWidth(2);
		h_mt2lblb_cut[dt] ->SetLineWidth(2);
		h_mt2lblbInteg[dt] -> SetLineWidth(2);
		h_mt2lblbInteg_cut[dt] ->SetLineWidth(2);
		h_mt2lblbSignif -> SetLineWidth(2);
		h_mt2lblbSignif_cut -> SetLineWidth(2);

	}      

	for (int hf = 0; hf<2; hf++) { // stop or top


		int nBinsX = h_mt2lblbtrue[hf]->GetNbinsX();

		for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

			float recursive_integral_true = h_mt2lblbtrue[hf]->Integral(ib, 3001);
			h_mt2lblbtrueInteg[hf]->SetBinContent(ib, recursive_integral_true); // asigna el valor ThisBinContent al bin ib

			float recursive_integral = h_mt2lblb[hf]->Integral(ib, 3001);
			h_mt2lblbInteg[hf]->SetBinContent(ib, recursive_integral); // asigna el valor ThisBinContent al bin ib


			float recursive_integral_true_cut = h_mt2lblbtrue_cut[hf]->Integral(ib, 3001);
			h_mt2lblbtrueInteg_cut[hf]->SetBinContent(ib, recursive_integral_true_cut/(mt2lblbtrue_cut_Int[hf]/mt2lblbtrue_Int[hf])); 		

			float recursive_integral_cut = h_mt2lblb_cut[hf]->Integral(ib, 3001);
			h_mt2lblbInteg_cut[hf]->SetBinContent(ib, recursive_integral_cut/(mt2lblb_cut_Int[hf]/mt2lblb_Int[hf])); 		

		}

	}



int nBinsX = h_mt2lblbtrue[0]->GetNbinsX();

for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

	float	topBackground_true = mt2lblbtrue_Int[0] * h_mt2lblbtrueInteg[0]->GetBinContent(ib); 
	float	stopEvents_true = mt2lblbtrue_Int[1] * h_mt2lblbtrueInteg[1]->GetBinContent(ib);
	if (topBackground_true + stopEvents_true <= 0.) continue;
	float	significance_true = stopEvents_true/(std::sqrt(stopEvents_true+topBackground_true));
	h_mt2lblbtrueSignif->SetBinContent(ib,significance_true); 

	float	topBackground = mt2lblb_Int[0] * h_mt2lblbInteg[0]->GetBinContent(ib); 
	float	stopEvents = mt2lblb_Int[1] * h_mt2lblbInteg[1]->GetBinContent(ib);
	if (topBackground + stopEvents <= 0.) continue;
	float	significance = stopEvents/(std::sqrt(stopEvents+topBackground));
	h_mt2lblbSignif->SetBinContent(ib,significance); 


	float	topBackground_true_cut = mt2lblbtrue_Int[0] * h_mt2lblbtrueInteg_cut[0]->GetBinContent(ib); 
	float	stopEvents_true_cut =  mt2lblbtrue_Int[1] * h_mt2lblbtrueInteg_cut[1]->GetBinContent(ib);
	if (topBackground_true_cut + stopEvents_true_cut <= 0.) continue;
	float	significance_true_cut = stopEvents_true_cut/(std::sqrt(stopEvents_true_cut+topBackground_true_cut));
	h_mt2lblbtrueSignif_cut->SetBinContent(ib,significance_true_cut); 

	float	topBackground_cut = mt2lblb_Int[0] * h_mt2lblbInteg_cut[0]->GetBinContent(ib); 
	float	stopEvents_cut =  mt2lblb_Int[1] * h_mt2lblbInteg_cut[1]->GetBinContent(ib);
	if (topBackground_cut + stopEvents_cut <= 0.) continue;
	float	significance_cut = stopEvents_cut/(std::sqrt(stopEvents_cut+topBackground_cut));
	h_mt2lblbSignif_cut->SetBinContent(ib,significance_cut); 

}




for (int dt = 0; dt<2; dt++) { // stop or top

	CCtrue->cd(1); // se pone en el TPad 1 
	h_mt2lblbtrue[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbtrue[dt]->GetXaxis()->SetTitle("MT2lblb");
	h_mt2lblbtrue_cut[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbtrue[dt]->DrawCopy(Option);
	h_mt2lblbtrue_cut[dt]->DrawCopy("histosame");

	TLegend *leg1 = new TLegend(0.3,0.1,0.4,0.4);
	leg1->AddEntry(h_mt2lblbtrue[0],"top","l");
	leg1->AddEntry(h_mt2lblbtrue[1],"stop","l");
	leg1->AddEntry(h_mt2lblbtrue_cut[0],"top with cut","l");
	leg1->AddEntry(h_mt2lblbtrue_cut[1],"stop with cut","l");
	leg1->Draw();


	CCtrue->cd(2); // se pone en el TPad 1 
	h_mt2lblbtrueInteg[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbtrueInteg[dt]->GetXaxis()->SetTitle("MT2lblb integral");
	h_mt2lblbtrueInteg_cut[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbtrueInteg[dt]->DrawCopy(Option);
	h_mt2lblbtrueInteg_cut[dt]->DrawCopy("histosame");



	CCtrue->cd(3); // se pone en el TPad 1 
	h_mt2lblbtrueSignif->GetXaxis()->SetRange(1, 300);
	h_mt2lblbtrueSignif_cut->GetXaxis()->SetRange(1, 300);
	h_mt2lblbtrueSignif->GetXaxis()->SetTitle("Significance");
	h_mt2lblbtrueSignif_cut->DrawCopy(Option);
	h_mt2lblbtrueSignif->DrawCopy("histosame");


	CC->cd(1); // se pone en el TPad 1 
	h_mt2lblb[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2lblb[dt]->GetXaxis()->SetTitle("MT2lblb");
	h_mt2lblb_cut[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2lblb[dt]->DrawCopy(Option);
	h_mt2lblb_cut[dt]->DrawCopy("histosame");

	TLegend *leg2 = new TLegend(0.3,0.1,0.4,0.4);
	leg2->AddEntry(h_mt2lblb[0],"top","l");
	leg2->AddEntry(h_mt2lblb[1],"stop","l");
	leg2->AddEntry(h_mt2lblb_cut[0],"top with cut","l");
	leg2->AddEntry(h_mt2lblb_cut[1],"stop with cut","l");
	leg2->Draw();


	CC->cd(2); // se pone en el TPad 1 
	h_mt2lblbInteg[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbInteg[dt]->GetXaxis()->SetTitle("MT2lblb integral");
	h_mt2lblbInteg_cut[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2lblbInteg[dt]->DrawCopy(Option);
	h_mt2lblbInteg_cut[dt]->DrawCopy("histosame");



	CC->cd(3); // se pone en el TPad 1 
	h_mt2lblbSignif->GetXaxis()->SetRange(1, 300);
	h_mt2lblbSignif_cut->GetXaxis()->SetRange(1, 300);
	h_mt2lblbSignif->GetXaxis()->SetTitle("Significance");
	h_mt2lblbSignif_cut->DrawCopy(Option);
	h_mt2lblbSignif->DrawCopy("histosame");

	Option= "histosame";
}   


CCtrue->Print("MT2lblbtrueStudies_cut135.png");
CC->Print("MT2lblbStudies_cut135.png");
}

