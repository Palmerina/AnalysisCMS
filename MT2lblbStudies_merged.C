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




void MT2lblbStudies_merged() {

	TFile *HistogramFile = TFile::Open("MT2lblbHistos.root");


	TH1D* h_mt2lblbtrue[2];
	TH1D* h_mt2lblbtrue_cut[2][50];
	TH1D* h_mt2lblbInteg[2];
	TH1D* h_mt2lblbInteg_cut[2][50];
	TH1D* h_mt2lblbSignif;
	TH1D* h_mt2lblbSignif_cut[50];

	TString HistoName[2]={"_top", "_stop"};
	int HistoCol[2] = {4, 2};
	int HistoColSignif[50];
	TString Legends[50];
	int LineStyle[50];

	for (int i=0; i<50; i++) {
		LineStyle[i] = i+2;
		HistoColSignif[i] = i+2;
	}

	for (float cut = 150.0; cut <= 170.0; cut += 5.0) {

		int icut = (cut-150)/5;
		char scut [50];
		sprintf(scut, "_%.f", cut);
		Legends[icut] = "mlbtrue <=  ";
		Legends[icut].Append(scut);

	}

	TCanvas *CC = new TCanvas("CC", "", 1200, 830);
	CC->Divide(1, 3);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetLogy();CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2->SetLogy();CC2->SetGridx(); CC2->SetGridy(); 
	TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mlblbtrue_cut
	CC3->SetGridx(); CC3->SetGridy(); 

	TString Option = "histo";


	float mt2lblbtrue_Int[2];
	float mt2lblbtrue_cut_Int[2][50];

	for (int i  = 0; i<2; i++) {


		h_mt2lblbtrue[i] = (TH1D*) HistogramFile->Get("h_mt2lblbtrue" + HistoName[i]);
		h_mt2lblbInteg[i] = new TH1D("h_mt2lblbInteg" + HistoName[i], "h_mt2lblbInteg", 3000, 0, 3000); //2 histos para top y stop
		h_mt2lblbSignif= new TH1D("h_mt2lblbSignif", "h_mt2lblbSignif", 3000, 0, 3000); //2 histos para top y stop
		for (float cut = 150.0; cut <= 170.0; cut += 5.0) {

			int icut = (cut-150)/5;
			char scut [50];
			sprintf(scut, "_%.f", cut);
			h_mt2lblbtrue_cut[i][icut] = (TH1D*) HistogramFile -> Get("h_mt2lblbtrue_cut" + HistoName[i] + scut); //2 histos para top y stop
			h_mt2lblbInteg_cut[i][icut] = new TH1D("h_mt2lblbInteg_cut" + HistoName[i] + scut, "h_mt2lblbInteg_cut" + HistoName[i] + scut, 3000, 0, 3000); //2 histos para top y stop
			h_mt2lblbSignif_cut[icut] = new TH1D("h_mt2lblbSignif_cut" + HistoName[i] + scut, "h_mt2lblbSignif_cut" + HistoName[i] + scut, 3000, 0, 3000); //2 histos para top y stop
		}
	}

	for (int dt  = 0; dt<2; dt++) {


		mt2lblbtrue_Int[dt] = h_mt2lblbtrue[dt]->Integral();
		h_mt2lblbtrue[dt]->Scale(1./mt2lblbtrue_Int[dt]); // normalization of the histogram


		h_mt2lblbtrue[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbInteg[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblbSignif -> SetLineColor(1);

		h_mt2lblbtrue[dt] -> SetLineStyle(1);
		h_mt2lblbInteg[dt] -> SetLineStyle(1);
		h_mt2lblbSignif -> SetLineStyle(1);

		h_mt2lblbtrue[dt] -> SetLineWidth(2);
		h_mt2lblbInteg[dt] -> SetLineWidth(2);
		h_mt2lblbSignif -> SetLineWidth(2);


		for (float cut = 150.0; cut <= 170.0; cut += 5.0) {

			int icut = (cut-150)/5;

			mt2lblbtrue_cut_Int[dt][icut] = h_mt2lblbtrue_cut[dt][icut]->Integral();
			h_mt2lblbtrue_cut[dt][icut]->Scale(1./mt2lblbtrue_Int[dt]); // normalization of the histogram


			h_mt2lblbtrue_cut[dt][icut] -> SetLineColor(HistoCol[dt]);
			h_mt2lblbInteg_cut[dt][icut] -> SetLineColor(HistoCol[dt]);
			h_mt2lblbSignif_cut[icut] -> SetLineColor(HistoColSignif[icut]);

			h_mt2lblbtrue_cut[dt][icut] ->SetLineStyle(LineStyle[icut]);
			h_mt2lblbInteg_cut[dt][icut] ->SetLineStyle(LineStyle[icut]);
			h_mt2lblbSignif_cut[icut] -> SetLineStyle(1);

			h_mt2lblbtrue_cut[dt][icut] ->SetLineWidth(2);
			h_mt2lblbInteg_cut[dt][icut] ->SetLineWidth(2);
			h_mt2lblbSignif_cut[icut] -> SetLineWidth(2);

		}
	}     


	for (int hf = 0; hf<2; hf++) { // stop or top


		int nBinsX = h_mt2lblbtrue[hf]->GetNbinsX();

		for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

			float recursive_integral = h_mt2lblbtrue[hf]->Integral(ib, 3001);
			h_mt2lblbInteg[hf]->SetBinContent(ib, recursive_integral); // asigna el valor ThisBinContent al bin ib


			for (float cut = 150.0; cut <= 170.0; cut += 5.0) {

				int icut = (cut-150)/5;
				float recursive_integral_cut = h_mt2lblbtrue_cut[hf][icut]->Integral(ib, 3001);
				h_mt2lblbInteg_cut[hf][icut]->SetBinContent(ib, recursive_integral_cut/(mt2lblbtrue_cut_Int[hf][icut]/mt2lblbtrue_Int[hf])); 			}

		}
	}





	int nBinsX = h_mt2lblbtrue[0]->GetNbinsX();

	for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

		float	topBackground = mt2lblbtrue_Int[0] * h_mt2lblbInteg[0]->GetBinContent(ib); 
		float	stopEvents = mt2lblbtrue_Int[1] * h_mt2lblbInteg[1]->GetBinContent(ib);
		if (topBackground + stopEvents <= 0.) continue;
		float	significance = stopEvents/(std::sqrt(stopEvents+topBackground));
		h_mt2lblbSignif->SetBinContent(ib,significance); 


		for (float cut = 150.0; cut <= 170.0; cut += 5.0) {

			int icut = (cut-150)/5;

			float	topBackground_cut = mt2lblbtrue_cut_Int[0][icut] * h_mt2lblbInteg_cut[0][icut]->GetBinContent(ib); 
			float	stopEvents_cut =  mt2lblbtrue_cut_Int[1][icut] * h_mt2lblbInteg_cut[1][icut]->GetBinContent(ib);
			if (topBackground_cut + stopEvents_cut <= 0.) continue;
			float	significance_cut = stopEvents_cut/(std::sqrt(stopEvents_cut+topBackground_cut));
			h_mt2lblbSignif_cut[icut]->SetBinContent(ib,significance_cut); 

		}

	}




	for (int dt = 0; dt<2; dt++) { // stop or top

		CC->cd(1); // se pone en el TPad 1 
		h_mt2lblbtrue[dt]->GetXaxis()->SetRange(0, 300);
		h_mt2lblbtrue[dt]->GetXaxis()->SetTitle("MT2lblb");
		h_mt2lblbtrue[dt]->SetMaximum(1);
		h_mt2lblbtrue[dt]->DrawCopy(Option);

		TLegend *leg1 = new TLegend(0.3,0.1,0.6,0.5);
		leg1->AddEntry(h_mt2lblbtrue[0],"top","l");
		leg1->AddEntry(h_mt2lblbtrue[1],"stop","l");

		for (float cut = 150.0; cut <= 170.0; cut += 5.0) {

			int icut = (cut-150)/5;
			char scut [50];
			sprintf(scut, "_%.f", cut);

			h_mt2lblbtrue_cut[dt][icut]->GetXaxis()->SetRange(0, 300);
			h_mt2lblbtrue_cut[dt][icut]->DrawCopy("histosame");

			leg1->AddEntry(h_mt2lblbtrue_cut[0][icut],Legends[icut], "l");

		}

		leg1->Draw();


		CC->cd(2); // se pone en el TPad 1 
		h_mt2lblbInteg[dt]->GetXaxis()->SetRange(1, 300);
		h_mt2lblbInteg[dt]->GetXaxis()->SetTitle("MT2lblb integral");
		h_mt2lblbInteg[dt]->DrawCopy(Option);

		for (float cut = 150.0; cut <= 170.0; cut += 5.0) {

			int icut = (cut-150)/5;
			h_mt2lblbInteg_cut[dt][icut]->GetXaxis()->SetRange(1, 300);
			h_mt2lblbInteg_cut[dt][icut]->DrawCopy("histosame");

		}

		CC->cd(3); // se pone en el TPad 1 
		h_mt2lblbSignif->GetXaxis()->SetRange(1, 300);
		h_mt2lblbSignif->GetXaxis()->SetTitle("Significance");
		h_mt2lblbSignif->DrawCopy("histosame");

		TLegend *leg2 = new TLegend(0.1,0.5,0.4,0.9);
		leg2->AddEntry(h_mt2lblbSignif,"no cut","l");

		for (float cut = 150.0; cut <= 170.0; cut += 5.0) {

			int icut = (cut-150)/5;
			h_mt2lblbSignif_cut[icut]->GetXaxis()->SetRange(1, 300);
			h_mt2lblbSignif_cut[icut]->DrawCopy(Option);


			leg2->AddEntry(h_mt2lblbSignif_cut[icut],Legends[icut], "l");
		}

		leg2->Draw();

		Option= "histosame";
	}   


	CC->Print("MT2lblbStudies_merged.png");
}
