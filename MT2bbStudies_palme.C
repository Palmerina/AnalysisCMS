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




void MT2bbStudies_palme() {

	TFile *HistogramFile = TFile::Open("MT2bbHistos_palme.root");

	TH1D* h_mt2bb[2];
	TH1D* h_mt2bbInteg[2];
	TH1D* h_mt2bbSignif;

	TString HistoName[2]={"_top", "_stop"};


	float mt2bb_Int[2];

	for (int i  = 0; i<2; i++) {


		h_mt2bb[i] = (TH1D*) HistogramFile->Get("h_mt2bb" + HistoName[i]);

		h_mt2bbInteg[i] = new TH1D("h_mt2bbInteg" + HistoName[i], "h_mt2bbInteg", 3000, 0, 3000); //2 histos para top y stop
		h_mt2bbSignif = new TH1D("h_mt2bbSignif" + HistoName[i], "h_mt2bbSignif", 3000, 0, 3000); //2 histos para top y stop

	}


	TCanvas *CC = new TCanvas("CC", "", 600, 800);
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

		mt2bb_Int[dt] = h_mt2bb[dt]->Integral(1, 3001);
		h_mt2bb[dt]->Scale(1./mt2bb_Int[dt]); // normalization of the histogram



		h_mt2bb[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbInteg[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2bbSignif -> SetLineColor(1);

		h_mt2bb[dt] -> SetLineStyle(1);
		h_mt2bbInteg[dt] -> SetLineStyle(1);
		h_mt2bbSignif -> SetLineStyle(1);

		h_mt2bb[dt] -> SetLineWidth(2);
		h_mt2bbInteg[dt] -> SetLineWidth(2);
		h_mt2bbSignif -> SetLineWidth(2);

	}      

	for (int hf = 0; hf<2; hf++) { // stop or top


		int nBinsX = h_mt2bb[hf]->GetNbinsX();

		for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

			float recursive_integral = h_mt2bb[hf]->Integral(ib, 3001);
			h_mt2bbInteg[hf]->SetBinContent(ib, recursive_integral); // asigna el valor ThisBinContent al bin ib

		}

	}



int nBinsX = h_mt2bb[0]->GetNbinsX();

for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

	float	topBackground = mt2bb_Int[0] * h_mt2bbInteg[0]->GetBinContent(ib); 
	float	stopEvents = mt2bb_Int[1] * h_mt2bbInteg[1]->GetBinContent(ib);
	if (topBackground + stopEvents <= 0.) continue;
	float	significance = stopEvents/(std::sqrt(stopEvents+topBackground));
	h_mt2bbSignif->SetBinContent(ib,significance); 


}




for (int dt = 0; dt<2; dt++) { // stop or top

	CC->cd(1); // se pone en el TPad 1 
	h_mt2bb[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bb[dt]->GetXaxis()->SetTitle("MT2bb");
	h_mt2bb[dt]->DrawCopy(Option);

	TLegend *leg2 = new TLegend(0.3,0.1,0.4,0.4);
	leg2->AddEntry(h_mt2bb[0],"top","l");
	leg2->AddEntry(h_mt2bb[1],"stop","l");
	leg2->Draw();


	CC->cd(2); // se pone en el TPad 1 
	h_mt2bbInteg[dt]->GetXaxis()->SetRange(1, 300);
	h_mt2bbInteg[dt]->GetXaxis()->SetTitle("MT2bb integral");
	h_mt2bbInteg[dt]->DrawCopy(Option);



	CC->cd(3); // se pone en el TPad 1 
	h_mt2bbSignif->GetXaxis()->SetRange(1, 300);
	h_mt2bbSignif->GetXaxis()->SetTitle("Significance");
	h_mt2bbSignif->DrawCopy(Option);

	Option= "histosame";
}   


CC->Print("MT2bbStudies_palme.png");
}
