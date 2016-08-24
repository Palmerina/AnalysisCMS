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




void MT2lblbMAOS_Studies() {

	TFile *HistogramFile = TFile::Open("MT2lblbMAOS_Histos.root");


	TH1D* h_mt2lblb[2];
	TH1D* h_mt2lblbInteg[2];
	TH1D* h_mt2lblbSignif;

	TH1D* h_mt2lblbMAOS[2];
	TH1D* h_mt2lblbMAOSInteg[2];
	TH1D* h_mt2lblbMAOSSignif;

	TString HistoName[2]={"_top", "_stop"};


	float mt2lblb_Int[2];

	float mt2lblbMAOS_Int[2];

	for (int i  = 0; i<2; i++) {


		h_mt2lblb[i] = (TH1D*) HistogramFile->Get("h_mt2lblb" + HistoName[i]);
		h_mt2lblbMAOS[i] = (TH1D*) HistogramFile->Get("h_mt2lblbMAOS" + HistoName[i]);

		h_mt2lblbInteg[i] = new TH1D("h_mt2lblbInteg" + HistoName[i], "h_mt2lblbInteg", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbMAOSInteg[i] = new TH1D("h_mt2lblbMAOSInteg" + HistoName[i], "h_mt2lblbMAOSInteg", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbSignif = new TH1D("h_mt2lblbSignif" + HistoName[i], "h_mt2lblbSignif", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbMAOSSignif = new TH1D("h_mt2lblbMAOSSignif" + HistoName[i], "h_mt2lblbMAOSSignif", 3000, 0, 3000); //2 histos para top y stop


	}

	TCanvas *CC = new TCanvas("CC", "", 1200, 850);
	CC->Divide(2, 3);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblb
	CC1->SetLogy();CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblb_cut
	CC2->SetLogy();CC2->SetGridx(); CC2->SetGridy(); 
	TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mlblb_cut
	CC3->SetLogy(); CC3->SetGridx(); CC3->SetGridy(); 
	TPad *CC4 = (TPad*)CC->GetPad(4); //Ahi va el h_mt2mlblb
	CC4->SetLogy();CC4->SetGridx(); CC4->SetGridy(); 
	TPad *CC5 = (TPad*)CC->GetPad(5); //Ahi va el h_mt2mlblb
	CC5->SetGridx(); CC5->SetGridy(); 
	TPad *CC6 = (TPad*)CC->GetPad(6); //Ahi va el h_mt2mlblb
	CC6->SetGridx(); CC6->SetGridy(); 


	TString Option = "histo";

	int HistoCol[2] = {4, 2};



	for (int dt  = 0; dt<2; dt++) {


		mt2lblb_Int[dt] = h_mt2lblb[dt]->Integral(1, 3001);
		h_mt2lblb[dt]->Scale(1./mt2lblb_Int[dt]); // normalization of the histogram

		mt2lblbMAOS_Int[dt] = h_mt2lblbMAOS[dt]->Integral(1, 3001);
		h_mt2lblbMAOS[dt]->Scale(1./mt2lblb_Int[dt]); // normalization of the histogram


		h_mt2lblb[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbMAOS[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbInteg[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbMAOSInteg[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbSignif -> SetLineColor(1);
		h_mt2lblbMAOSSignif -> SetLineColor(1);

		h_mt2lblb[dt] -> SetLineStyle(1);
		h_mt2lblbMAOS[dt] ->SetLineStyle(1);
		h_mt2lblbInteg[dt] -> SetLineStyle(1);
		h_mt2lblbMAOSInteg[dt] ->SetLineStyle(1);
		h_mt2lblbSignif -> SetLineStyle(1);
		h_mt2lblbMAOSSignif -> SetLineStyle(1);

		h_mt2lblb[dt] -> SetLineWidth(2);
		h_mt2lblbMAOS[dt] ->SetLineWidth(2);
		h_mt2lblbInteg[dt] -> SetLineWidth(2);
		h_mt2lblbMAOSInteg[dt] ->SetLineWidth(2);
		h_mt2lblbSignif -> SetLineWidth(2);
		h_mt2lblbMAOSSignif -> SetLineWidth(2);


	}      

	for (int hf = 0; hf<2; hf++) { // stop or top

		int nBinsX = h_mt2lblbMAOS[hf]->GetNbinsX();

		for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

			float recursive_integral = h_mt2lblb[hf]->Integral(ib, 3001);
			h_mt2lblbInteg[hf]->SetBinContent(ib, recursive_integral); // asigna el valor ThisBinContent al bin ib

			float recursive_integral_MAOS = h_mt2lblbMAOS[hf]->Integral(ib, 3001);
			h_mt2lblbMAOSInteg[hf]->SetBinContent(ib, recursive_integral_MAOS); // asigna el valor ThisBinContent al bin ib

		}

	}



	int nBinsX = h_mt2lblbMAOS[0]->GetNbinsX();

	for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

		float	topBackground = mt2lblb_Int[0] * h_mt2lblbInteg[0]->GetBinContent(ib); 
		float	stopEvents = mt2lblb_Int[1] * h_mt2lblbInteg[1]->GetBinContent(ib);
		if (topBackground + stopEvents <= 0.) continue;
		float	significance = stopEvents/(std::sqrt(stopEvents+topBackground));
		h_mt2lblbSignif->SetBinContent(ib,significance); 

		float	topBackground_MAOS = mt2lblbMAOS_Int[0] * h_mt2lblbMAOSInteg[0]->GetBinContent(ib); 
		float	stopEvents_MAOS = mt2lblbMAOS_Int[1] * h_mt2lblbMAOSInteg[1]->GetBinContent(ib);
		if (topBackground_MAOS + stopEvents_MAOS <= 0.) continue;
		float	significance_MAOS = stopEvents_MAOS/(std::sqrt(stopEvents_MAOS+topBackground_MAOS));
		h_mt2lblbMAOSSignif->SetBinContent(ib,significance_MAOS); 

	}




	for (int dt = 0; dt<2; dt++) { // stop or top

		CC->cd(1); // se pone en el TPad 1 
		h_mt2lblb[dt]->GetXaxis()->SetRange(1, 300);
		h_mt2lblb[dt]->GetXaxis()->SetTitle("mt2lblb");
		h_mt2lblb[dt]->DrawCopy(Option);

		TLegend *leg1 = new TLegend(0.3,0.1,0.4,0.3);
		leg1->AddEntry(h_mt2lblb[0],"top","l");
		leg1->AddEntry(h_mt2lblb[1],"stop","l");
		leg1->Draw();

		CC->cd(2); // se pone en el TPad 1 
		h_mt2lblbMAOS[dt]->GetXaxis()->SetRange(1, 300);
		h_mt2lblbMAOS[dt]->GetXaxis()->SetTitle("mt2lblbMAOS");
		h_mt2lblbMAOS[dt]->DrawCopy(Option);

		TLegend *leg2 = new TLegend(0.3,0.1,0.4,0.3);
		leg2->AddEntry(h_mt2lblbMAOS[0],"top","l");
		leg2->AddEntry(h_mt2lblbMAOS[1],"stop","l");
		leg2->Draw();


		CC->cd(3); // se pone en el TPad 1 
		h_mt2lblbInteg[dt]->GetXaxis()->SetRange(1, 300);
		h_mt2lblbInteg[dt]->GetXaxis()->SetTitle("mt2lblb integral");
		h_mt2lblbInteg[dt]->DrawCopy(Option);

		TLegend *leg3 = new TLegend(0.3,0.1,0.4,0.3);
		leg3->AddEntry(h_mt2lblbInteg[0],"top","l");
		leg3->AddEntry(h_mt2lblbInteg[1],"stop","l");
		leg3->Draw();


		CC->cd(4); // se pone en el TPad 1 
		h_mt2lblbMAOSInteg[dt]->GetXaxis()->SetRange(1, 300);
		h_mt2lblbMAOSInteg[dt]->GetXaxis()->SetTitle("mt2lblbMAOS integral");
		h_mt2lblbMAOSInteg[dt]->DrawCopy(Option);

		TLegend *leg4 = new TLegend(0.3,0.1,0.4,0.3);
		leg4->AddEntry(h_mt2lblbMAOSInteg[0],"top","l");
		leg4->AddEntry(h_mt2lblbMAOSInteg[1],"stop","l");
		leg4->Draw();


		CC->cd(5); // se pone en el TPad 1 
		h_mt2lblbSignif->GetXaxis()->SetRange(1, 400);
		h_mt2lblbSignif->GetXaxis()->SetTitle("Significance");
		h_mt2lblbSignif->DrawCopy(Option);


		CC->cd(6); // se pone en el TPad 1 
		h_mt2lblbMAOSSignif->GetXaxis()->SetRange(1, 400);
		h_mt2lblbMAOSSignif->GetXaxis()->SetTitle("mt2lblbMAOS Significance");
		h_mt2lblbMAOSSignif->DrawCopy(Option);




		Option= "histosame";
	}   


	CC->Print("MT2lblbMAOS.png");
}
