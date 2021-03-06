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




void MT2bbStudies() {

	TFile *HistogramFile = TFile::Open("MT2bbHistos.root");


	TH1D* h_mt2bbtrue[2];
	TH1D* h_mt2bbtrue_cut[2];
	TH1D* h_mt2bbtrueInteg[2];
	TH1D* h_mt2bbtrueInteg_cut[2];
	TH1D* h_mt2bbtrueSignif;
	TH1D* h_mt2bbtrueSignif_cut;

	TH1D* h_mt2bb[2];
	TH1D* h_mt2bb_cut[2];
	TH1D* h_mt2bbInteg[2];
	TH1D* h_mt2bbInteg_cut[2];
	TH1D* h_mt2bbSignif;
	TH1D* h_mt2bbSignif_cut;

	TString HistoName[2]={"_top", "_stop"};


	float mt2bbtrue_Int[2];
	float mt2bbtrue_cut_Int[2];

	float mt2bb_Int[2];
	float mt2bb_cut_Int[2];

	for (int i  = 0; i<2; i++) {


		h_mt2bbtrue[i] = (TH1D*) HistogramFile->Get("h_mt2bbtrue" + HistoName[i]);
		h_mt2bbtrue_cut[i] = (TH1D*) HistogramFile->Get("h_mt2bbtrue_cut" + HistoName[i]);

		h_mt2bbtrueInteg[i] = new TH1D("h_mt2bbtrueInteg" + HistoName[i], "h_mt2bbtrueInteg", 3000, 0, 3000); //2 histos para top y stop
		h_mt2bbtrueInteg_cut[i] = new TH1D("h_mt2bbtrueInteg_cut" + HistoName[i], "h_mt2bbtrueInteg_cut", 3000, 0, 3000); //2 histos para top y stop
		h_mt2bbtrueSignif = new TH1D("h_mt2bbtrueSignif" + HistoName[i], "h_mt2bbtrueSignif", 3000, 0, 3000); //2 histos para top y stop
		h_mt2bbtrueSignif_cut = new TH1D("h_mt2bbtrueSignif_cut" + HistoName[i], "h_mt2bbtrueSignif_cut", 3000, 0, 3000); //2 histos para top y stop


		h_mt2bb[i] = (TH1D*) HistogramFile->Get("h_mt2bb" + HistoName[i]);
		h_mt2bb_cut[i] = (TH1D*) HistogramFile->Get("h_mt2bb_cut" + HistoName[i]);

		h_mt2bbInteg[i] = new TH1D("h_mt2bbInteg" + HistoName[i], "h_mt2bbInteg", 3000, 0, 3000); //2 histos para top y stop
		h_mt2bbInteg_cut[i] = new TH1D("h_mt2bbInteg_cut" + HistoName[i], "h_mt2bbInteg_cut", 3000, 0, 3000); //2 histos para top y stop
		h_mt2bbSignif = new TH1D("h_mt2bbSignif" + HistoName[i], "h_mt2bbSignif", 3000, 0, 3000); //2 histos para top y stop
		h_mt2bbSignif_cut = new TH1D("h_mt2bbSignif_cut" + HistoName[i], "h_mt2bbSignif_cut", 3000, 0, 3000); //2 histos para top y stop

	}

	TCanvas *CCtrue = new TCanvas("CCtrue", "", 1450, 800);
	CCtrue->Divide(1, 3);
	TPad *CC1true = (TPad*)CCtrue->GetPad(1); //Ahi va el h_mt2mbbtrue
	CC1true->SetLogy();CC1true->SetGridx(); CC1true->SetGridy(); 
	TPad *CC2true = (TPad*)CCtrue->GetPad(2); //Ahi va el h_mt2mbbtrue_cut
	CC2true->SetLogy();CC2true->SetGridx(); CC2true->SetGridy(); 
	TPad *CC3true = (TPad*)CCtrue->GetPad(3); //Ahi va el h_mt2mbbtrue_cut
	CC3true->SetGridx(); CC3true->SetGridy(); 


	TCanvas *CC = new TCanvas("CC", "", 1450, 800);
	CC->Divide(1, 3);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mbbtrue
	CC1->SetLogy();CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mbbtrue_cut
	CC2->SetLogy();CC2->SetGridx(); CC2->SetGridy(); 
	TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mbbtrue_cut
	CC3->SetGridx(); CC3->SetGridy(); 

	TString Option = "histo";

	int HistoCol[2] = {4, 2};



	for (int dt  = 0; dt<2; dt++) {


		mt2bbtrue_Int[dt] = h_mt2bbtrue[dt]->Integral(1, 3001);
		h_mt2bbtrue[dt]->Scale(1./mt2bbtrue_Int[dt]); // normalization of the histogram

		mt2bbtrue_cut_Int[dt] = h_mt2bbtrue_cut[dt]->Integral(1, 3001);
		h_mt2bbtrue_cut[dt]->Scale(1./mt2bbtrue_Int[dt]); // normalization of the histogram


		h_mt2bbtrue[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbtrue_cut[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbtrueInteg[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbtrueInteg_cut[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbtrueSignif -> SetLineColor(1);
		h_mt2bbtrueSignif_cut -> SetLineColor(1);

		h_mt2bbtrue[dt] -> SetLineStyle(1);
		h_mt2bbtrue_cut[dt] ->SetLineStyle(2);
		h_mt2bbtrueInteg[dt] -> SetLineStyle(1);
		h_mt2bbtrueInteg_cut[dt] ->SetLineStyle(2);
		h_mt2bbtrueSignif -> SetLineStyle(1);
		h_mt2bbtrueSignif_cut -> SetLineStyle(2);

		h_mt2bbtrue[dt] -> SetLineWidth(2);
		h_mt2bbtrue_cut[dt] ->SetLineWidth(2);
		h_mt2bbtrueInteg[dt] -> SetLineWidth(2);
		h_mt2bbtrueInteg_cut[dt] ->SetLineWidth(2);
		h_mt2bbtrueSignif -> SetLineWidth(2);
		h_mt2bbtrueSignif_cut -> SetLineWidth(2);


		mt2bb_Int[dt] = h_mt2bb[dt]->Integral(1, 3001);
		h_mt2bb[dt]->Scale(1./mt2bb_Int[dt]); // normalization of the histogram

		mt2bb_cut_Int[dt] = h_mt2bb_cut[dt]->Integral(1, 3001);
		h_mt2bb_cut[dt]->Scale(1./mt2bb_Int[dt]); // normalization of the histogram


		h_mt2bb[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bb_cut[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbInteg[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbInteg_cut[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbSignif -> SetLineColor(1);
		h_mt2bbSignif_cut -> SetLineColor(1);

		h_mt2bb[dt] -> SetLineStyle(1);
		h_mt2bb_cut[dt] ->SetLineStyle(2);
		h_mt2bbInteg[dt] -> SetLineStyle(1);
		h_mt2bbInteg_cut[dt] ->SetLineStyle(2);
		h_mt2bbSignif -> SetLineStyle(1);
		h_mt2bbSignif_cut -> SetLineStyle(2);

		h_mt2bb[dt] -> SetLineWidth(2);
		h_mt2bb_cut[dt] ->SetLineWidth(2);
		h_mt2bbInteg[dt] -> SetLineWidth(2);
		h_mt2bbInteg_cut[dt] ->SetLineWidth(2);
		h_mt2bbSignif -> SetLineWidth(2);
		h_mt2bbSignif_cut -> SetLineWidth(2);

	}      

	for (int hf = 0; hf<2; hf++) { // stop or top


		int nBinsX = h_mt2bbtrue[hf]->GetNbinsX();

		for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

			float recursive_integral_true = h_mt2bbtrue[hf]->Integral(ib, 3001);
			h_mt2bbtrueInteg[hf]->SetBinContent(ib, recursive_integral_true); // asigna el valor ThisBinContent al bin ib

			float recursive_integral = h_mt2bb[hf]->Integral(ib, 3001);
			h_mt2bbInteg[hf]->SetBinContent(ib, recursive_integral); // asigna el valor ThisBinContent al bin ib


			float recursive_integral_true_cut = h_mt2bbtrue_cut[hf]->Integral(ib, 3001);
			h_mt2bbtrueInteg_cut[hf]->SetBinContent(ib, recursive_integral_true_cut/(mt2bbtrue_cut_Int[hf]/mt2bbtrue_Int[hf])); 		

			float recursive_integral_cut = h_mt2bb_cut[hf]->Integral(ib, 3001);
			h_mt2bbInteg_cut[hf]->SetBinContent(ib, recursive_integral_cut/(mt2bb_cut_Int[hf]/mt2bb_Int[hf])); 		

		}

	}



int nBinsX = h_mt2bbtrue[0]->GetNbinsX();

for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

	float	topBackground_true = mt2bbtrue_Int[0] * h_mt2bbtrueInteg[0]->GetBinContent(ib); 
	float	stopEvents_true = mt2bbtrue_Int[1] * h_mt2bbtrueInteg[1]->GetBinContent(ib);
	if (topBackground_true + stopEvents_true <= 0.) continue;
	float	significance_true = stopEvents_true/(std::sqrt(stopEvents_true+topBackground_true));
	h_mt2bbtrueSignif->SetBinContent(ib,significance_true); 

	float	topBackground = mt2bb_Int[0] * h_mt2bbInteg[0]->GetBinContent(ib); 
	float	stopEvents = mt2bb_Int[1] * h_mt2bbInteg[1]->GetBinContent(ib);
	if (topBackground + stopEvents <= 0.) continue;
	float	significance = stopEvents/(std::sqrt(stopEvents+topBackground));
	h_mt2bbSignif->SetBinContent(ib,significance); 


	float	topBackground_true_cut = mt2bbtrue_Int[0] * h_mt2bbtrueInteg_cut[0]->GetBinContent(ib); 
	float	stopEvents_true_cut =  mt2bbtrue_Int[1] * h_mt2bbtrueInteg_cut[1]->GetBinContent(ib);
	if (topBackground_true_cut + stopEvents_true_cut <= 0.) continue;
	float	significance_true_cut = stopEvents_true_cut/(std::sqrt(stopEvents_true_cut+topBackground_true_cut));
	h_mt2bbtrueSignif_cut->SetBinContent(ib,significance_true_cut); 

	float	topBackground_cut = mt2bb_Int[0] * h_mt2bbInteg_cut[0]->GetBinContent(ib); 
	float	stopEvents_cut =  mt2bb_Int[1] * h_mt2bbInteg_cut[1]->GetBinContent(ib);
	if (topBackground_cut + stopEvents_cut <= 0.) continue;
	float	significance_cut = stopEvents_cut/(std::sqrt(stopEvents_cut+topBackground_cut));
	h_mt2bbSignif_cut->SetBinContent(ib,significance_cut); 

}




for (int dt = 0; dt<2; dt++) { // stop or top

	CCtrue->cd(1); // se pone en el TPad 1 
	h_mt2bbtrue[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bbtrue[dt]->GetXaxis()->SetTitle("MT2bb");
	h_mt2bbtrue_cut[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bbtrue[dt]->DrawCopy(Option);
	h_mt2bbtrue_cut[dt]->DrawCopy("histosame");

	TLegend *leg1 = new TLegend(0.3,0.1,0.4,0.4);
	leg1->AddEntry(h_mt2bbtrue[0],"top","l");
	leg1->AddEntry(h_mt2bbtrue[1],"stop","l");
	leg1->AddEntry(h_mt2bbtrue_cut[0],"top with cut","l");
	leg1->AddEntry(h_mt2bbtrue_cut[1],"stop with cut","l");
	leg1->Draw();


	CCtrue->cd(2); // se pone en el TPad 1 
	h_mt2bbtrueInteg[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bbtrueInteg[dt]->GetXaxis()->SetTitle("MT2bb integral");
	h_mt2bbtrueInteg_cut[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bbtrueInteg[dt]->DrawCopy(Option);
	h_mt2bbtrueInteg_cut[dt]->DrawCopy("histosame");



	CCtrue->cd(3); // se pone en el TPad 1 
	h_mt2bbtrueSignif->GetXaxis()->SetRange(1, 300);
	h_mt2bbtrueSignif_cut->GetXaxis()->SetRange(1, 300);
	h_mt2bbtrueSignif->GetXaxis()->SetTitle("Significance");
	h_mt2bbtrueSignif_cut->DrawCopy(Option);
	h_mt2bbtrueSignif->DrawCopy("histosame");


	CC->cd(1); // se pone en el TPad 1 
	h_mt2bb[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bb[dt]->GetXaxis()->SetTitle("MT2bb");
	h_mt2bb_cut[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bb[dt]->DrawCopy(Option);
	h_mt2bb_cut[dt]->DrawCopy("histosame");

	TLegend *leg2 = new TLegend(0.3,0.1,0.4,0.4);
	leg2->AddEntry(h_mt2bb[0],"top","l");
	leg2->AddEntry(h_mt2bb[1],"stop","l");
	leg2->AddEntry(h_mt2bb_cut[0],"top with cut","l");
	leg2->AddEntry(h_mt2bb_cut[1],"stop with cut","l");
	leg2->Draw();


	CC->cd(2); // se pone en el TPad 1 
	h_mt2bbInteg[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bbInteg[dt]->GetXaxis()->SetTitle("MT2bb integral");
	h_mt2bbInteg_cut[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bbInteg[dt]->DrawCopy(Option);
	h_mt2bbInteg_cut[dt]->DrawCopy("histosame");



	CC->cd(3); // se pone en el TPad 1 
	h_mt2bbSignif->GetXaxis()->SetRange(1, 300);
	h_mt2bbSignif_cut->GetXaxis()->SetRange(1, 300);
	h_mt2bbSignif->GetXaxis()->SetTitle("Significance");
	h_mt2bbSignif_cut->DrawCopy(Option);
	h_mt2bbSignif->DrawCopy("histosame");

	Option= "histosame";
}   


CCtrue->Print("MT2bbtrueStudies_cut135.png");
CC->Print("MT2bbStudies_cut135.png");
}

