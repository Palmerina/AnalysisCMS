
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




void neutrinoMET_Studies_palme() {

	TFile *HistogramFile = TFile::Open("neutrinoMET_Histos_palme.root");


	TH2D* h_neutrinoMET[2];

	TString HistoName[2]={"_top", "_stop"};

	for (int dt = 0; dt < 2; dt++) {

		h_neutrinoMET[dt] = (TH2D*) HistogramFile->Get("h_neutrinoMET" + HistoName[dt]);

	}


	TCanvas *CC = new TCanvas("CC", "", 1000, 800);
	CC->Divide(1, 2);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2->SetGridx(); CC2->SetGridy(); 



	for (int dt = 0; dt < 2; dt++) {

		CC->cd(dt+1); // se pone en el TPad 1 

		h_neutrinoMET[dt]->GetXaxis()->SetRangeUser(1, 300);
		h_neutrinoMET[dt]->GetXaxis()->SetTitle("neutrino MET");
		h_neutrinoMET[dt]->GetYaxis()->SetRangeUser(1, 300);
		h_neutrinoMET[dt]->GetYaxis()->SetTitle("MET");
		h_neutrinoMET[dt]->SetTitle("MET" + HistoName[dt]);
		h_neutrinoMET[dt]->DrawCopy("box");

	}






	CC->Print("neutrinoMET.png");
}

