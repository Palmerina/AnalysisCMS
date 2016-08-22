
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
	TH2D* h_neutrinoPhi[2];

	TString HistoName[2]={"_top", "_stop"};

	for (int dt = 0; dt < 2; dt++) {

		h_neutrinoMET[dt] = (TH2D*) HistogramFile->Get("h_neutrinoMET" + HistoName[dt]);
		h_neutrinoPhi[dt] = (TH2D*) HistogramFile->Get("h_neutrinoPhi" + HistoName[dt]);
		cout << "bucle for" << endl;
	}


	TCanvas *CC = new TCanvas("CC", "", 1300, 800);
	CC->Divide(2, 2);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2->SetGridx(); CC2->SetGridy(); 
	TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mlblbtrue
	CC3->SetGridx(); CC3->SetGridy(); 
	TPad *CC4 = (TPad*)CC->GetPad(4); //Ahi va el h_mt2mlblbtrue_cut
	CC4->SetGridx(); CC4->SetGridy(); 
	cout << "TCanvas"<< endl;


	for (int dt = 0; dt < 2; dt++) {

		CC->cd(dt+1); // se pone en el TPad 1 

		h_neutrinoMET[dt]->GetXaxis()->SetRangeUser(1, 300);
		h_neutrinoMET[dt]->GetXaxis()->SetTitle("neutrino MET");
		h_neutrinoMET[dt]->GetYaxis()->SetRangeUser(1, 400);
		h_neutrinoMET[dt]->GetYaxis()->SetTitle("MET");
		h_neutrinoMET[dt]->SetTitle("MET" + HistoName[dt]);
		h_neutrinoMET[dt]->DrawCopy("box");
	
		cout << "cd(dt+1)"<<endl;
	}	
	for (int dt = 0; dt < 2; dt++) {
		CC->cd(dt+3); // se pone en el TPad 1 

		h_neutrinoPhi[dt]->GetXaxis()->SetRangeUser(1, 300);
		h_neutrinoPhi[dt]->GetXaxis()->SetTitle("neutrino Phi");
		h_neutrinoPhi[dt]->GetYaxis()->SetRangeUser(1, 400);
		h_neutrinoPhi[dt]->GetYaxis()->SetTitle("Phi");
		h_neutrinoPhi[dt]->SetTitle("Phi" + HistoName[dt]);
		h_neutrinoPhi[dt]->DrawCopy("box");
		
		cout << "cd(dt+3)"<<endl;
	}	





	CC->Print("neutrinoMET.png");
}

