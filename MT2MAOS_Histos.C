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

float eventW, metPfType1;
float njet, channel, nbjet30csvv2l, nbjet30csvv2m, nbjet30csvv2t;
float mt2ll, mt2llMAOS, mt2lltrue, mt2lblb, mt2lblbMAOS, mt2lblbcomb, mt2lblbMAOScomb, mt2lblbtrue, mt2lblbMAOStrue;
float mlb1, mlb1true, mlb1comb, mlb1truecomb;
float mlb2, mlb2true, mlb2comb, mlb2truecomb;
float tjet1pt, tjet1phi, tjet1eta, tjet1mass, tjet1csvv2ivf, tjet1assignment;
float tjet2pt, tjet2phi, tjet2eta, tjet2mass, tjet2csvv2ivf, tjet2assignment;
float bjet1pt, bjet1phi, bjet1eta, bjet1mass, bjet1csvv2ivf;
float bjet2pt, bjet2phi, bjet2eta, bjet2mass, bjet2csvv2ivf;

//top
TH2D *h_mt2ll_mt2llMAOS;
TH2D *h_mt2ll_mt2lblbMAOStrue; 
TH2D *h_mt2lblbtrue_mt2lblbMAOStrue; 
TH2D *h_mt2lblb_mt2lblbMAOS; 
//top y stop
TH2D *h_mt2ll_mt2lblbMAOS[2]; 
TH2D *h_mt2ll_mt2lblbMAOS_cut[2]; 


TTree *GetMiniTree(TFile *MiniTreeFile) {

	TTree *MiniTree = (TTree*) MiniTreeFile->Get("latino"); 

	MiniTree->SetBranchAddress("eventW",          &eventW);
	MiniTree->SetBranchAddress("metPfType1",      &metPfType1);
	MiniTree->SetBranchAddress("njet",            &njet);
	MiniTree->SetBranchAddress("channel",         &channel);

	MiniTree->SetBranchAddress("nbjet30csvv2l",   &nbjet30csvv2l);
	MiniTree->SetBranchAddress("nbjet30csvv2m",   &nbjet30csvv2m);
	MiniTree->SetBranchAddress("nbjet30csvv2t",   &nbjet30csvv2t);

	MiniTree->SetBranchAddress("mt2ll",           &mt2ll);
	MiniTree->SetBranchAddress("mt2llMAOS",       &mt2llMAOS);
	MiniTree->SetBranchAddress("mt2lltrue",       &mt2lltrue);
	MiniTree->SetBranchAddress("mt2lblb",         &mt2lblb);
	MiniTree->SetBranchAddress("mt2lblbMAOS",     &mt2lblbMAOS);
	MiniTree->SetBranchAddress("mt2lblbtrue",     &mt2lblbtrue);
	MiniTree->SetBranchAddress("mt2lblbMAOStrue", &mt2lblbMAOStrue);
	MiniTree->SetBranchAddress("mt2lblbcomb",     &mt2lblbcomb);
	MiniTree->SetBranchAddress("mt2lblbMAOScomb", &mt2lblbMAOScomb);

	MiniTree->SetBranchAddress("mlb1",            &mlb1);
	MiniTree->SetBranchAddress("mlb2",            &mlb2);
	MiniTree->SetBranchAddress("mlb1comb",        &mlb1comb);
	MiniTree->SetBranchAddress("mlb2comb",        &mlb2comb);
	MiniTree->SetBranchAddress("mlb1true",        &mlb1true);
	MiniTree->SetBranchAddress("mlb2true",        &mlb2true);
	MiniTree->SetBranchAddress("mlb1truecomb",    &mlb1truecomb);
	MiniTree->SetBranchAddress("mlb2truecomb",    &mlb2truecomb);

	MiniTree->SetBranchAddress("bjet1pt",         &bjet1pt);        
	MiniTree->SetBranchAddress("bjet1eta",        &bjet1eta);       
	MiniTree->SetBranchAddress("bjet1phi",        &bjet1phi);       
	MiniTree->SetBranchAddress("bjet1mass",       &bjet1mass);      
	MiniTree->SetBranchAddress("bjet1csvv2ivf",   &bjet1csvv2ivf);  
	MiniTree->SetBranchAddress("bjet2pt",         &bjet2pt);        
	MiniTree->SetBranchAddress("bjet2eta",        &bjet2eta);       
	MiniTree->SetBranchAddress("bjet2phi",        &bjet2phi);       
	MiniTree->SetBranchAddress("bjet2mass",       &bjet2mass);      
	MiniTree->SetBranchAddress("bjet2csvv2ivf",   &bjet2csvv2ivf);  
	MiniTree->SetBranchAddress("tjet1pt",         &tjet1pt);        
	MiniTree->SetBranchAddress("tjet1eta",        &tjet1eta);       
	MiniTree->SetBranchAddress("tjet1phi",        &tjet1phi);       
	MiniTree->SetBranchAddress("tjet1mass",       &tjet1mass);      
	MiniTree->SetBranchAddress("tjet1csvv2ivf",   &tjet1csvv2ivf);  
	MiniTree->SetBranchAddress("tjet1assignment", &tjet1assignment);
	MiniTree->SetBranchAddress("tjet2pt",         &tjet2pt);        
	MiniTree->SetBranchAddress("tjet2eta",        &tjet2eta);       
	MiniTree->SetBranchAddress("tjet2phi",        &tjet2phi);
	MiniTree->SetBranchAddress("tjet2mass",       &tjet2mass);      
	MiniTree->SetBranchAddress("tjet2csvv2ivf",   &tjet2csvv2ivf);  
	MiniTree->SetBranchAddress("tjet2assignment", &tjet2assignment);

	return MiniTree;

}

void MT2MAOS_Histos() {

	TString FileName[2] = {"./minitrees/minitreesLuca/TTTo2L2Nu.root",
		"./minitrees/minitreesLuca/T2tt_mStop500-525-550_mLSP1to425-325to450-1to475.root"};

	TString Option = "histo";
	TString Histoname[2] = {"_top", "_stop"};

	int HistoCol[2] = {4, 2};


	for (int dt = 0; dt<2; dt++) {

		TFile *MiniTreeFile = TFile::Open(FileName[dt]);

		TTree *MiniTree = GetMiniTree(MiniTreeFile);

		Int_t nentries = (Int_t) MiniTree->GetEntries();

		if (dt == 0){

			h_mt2ll_mt2llMAOS = new TH2D("h_mt2ll_mt2llMAOS" + Histoname[dt],"h_mt2ll_mt2llMAOS", 300, 0, 3000, 300, 0, 3000); //2 histos para top y stop
			h_mt2ll_mt2lblbMAOStrue = new TH2D("h_mt2ll_mt2lblbMAOStrue" + Histoname[dt], "h_mt2ll_mt2llMAOS", 300, 0, 3000, 300, 0, 3000); //2 histos para top y stop
			h_mt2lblbtrue_mt2lblbMAOStrue = new TH2D("h_mt2lblbtrue_mt2lblbMAOStrue" + Histoname[dt],"h_mt2lblbtrue_mt2lblbMAOStrue",   300, 0, 3000, 300, 0, 3000); //2 histos para top y stop
			h_mt2lblb_mt2lblbMAOS = new TH2D("h_mt2lblb_mt2lblbMAOS" + Histoname[dt],"h_mt2lblb_mt2lblbMAOS",   300, 0, 3000, 300, 0, 3000); //2 histos para top y stop

		}

		h_mt2ll_mt2lblbMAOS[dt] = new TH2D("h_mt2ll_mt2lblbMAOS" + Histoname[dt], "h_mt2lblb_mt2lblbMAOS", 300, 0, 3000, 300, 0, 3000); //2 histos para top y stop
		h_mt2ll_mt2lblbMAOS_cut[dt] = new TH2D("h_mt2ll_mt2lblbMAOS_cut" + Histoname[dt],"h_mt2ll_mt2lblbMAOS_cut",   300, 0, 3000, 300, 0, 3000); //2 histos para top y stop




		for (Int_t i = 0; i<nentries; i++) {

			MiniTree->GetEntry(i);
			if (njet<2) continue;
			if (nbjet30csvv2m < 1) continue;

			if (dt == 0){

				h_mt2ll_mt2llMAOS -> Fill(mt2ll, mt2llMAOS, eventW);
				h_mt2ll_mt2lblbMAOStrue -> Fill(mt2ll, mt2lblbMAOStrue, eventW);
				h_mt2lblbtrue_mt2lblbMAOStrue -> Fill(mt2lblbtrue, mt2lblbMAOStrue, eventW);
				h_mt2lblb_mt2lblbMAOS -> Fill(mt2lblb, mt2lblbMAOS, eventW);

			}

			h_mt2ll_mt2lblbMAOS[dt] -> Fill(mt2ll, mt2lblbMAOS, eventW);

			if (mlb1 <= 160 && mlb2 <= 160) {

				h_mt2ll_mt2lblbMAOS_cut[dt] -> Fill(mt2ll, mt2lblbMAOS, eventW);

			}
		}
	} 
     
	TFile *OutFile = new TFile("MT2MAOS_Histos.root", "recreate");

	h_mt2ll_mt2llMAOS->Write();
	h_mt2ll_mt2lblbMAOStrue->Write();
	h_mt2lblbtrue_mt2lblbMAOStrue->Write();
	h_mt2lblb_mt2lblbMAOS->Write();

	for (int hf = 0; hf<2; hf++) { // stop or top

		h_mt2ll_mt2lblbMAOS[hf]->Write();
		h_mt2ll_mt2lblbMAOS_cut[hf]->Write();

	}

	OutFile->Close(); 
}
