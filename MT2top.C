#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
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
#include "TLorentzVector.h"

#include <fstream>
#include <iostream>

#include "Mt2/Basic_Mt2_332_Calculator.h"
#include "Mt2/Basic_Mt2_Top_Calculator.h"

float eventW, metPfType1, pfType1Met, pfType1Metphi; 
float njet, channel, nbjet30csvv2l, nbjet30csvv2m, nbjet30csvv2t;
float mt2ll, mt2bb, mt2bbtrue, mt2lblb, mt2lblbcomb, mt2lblbtrue;
float mlb1, mlb1true, mlb1comb, mlb1truecomb;
float mlb2, mlb2true, mlb2comb, mlb2truecomb;
float tjet1pt, tjet1phi, tjet1eta, tjet1mass, tjet1csvv2ivf, tjet1assignment;
float tjet2pt, tjet2phi, tjet2eta, tjet2mass, tjet2csvv2ivf, tjet2assignment;
float bjet1pt, bjet1phi, bjet1eta, bjet1mass, bjet1csvv2ivf;
float bjet2pt, bjet2phi, bjet2eta, bjet2mass, bjet2csvv2ivf;
float lep1pt, lep1phi, lep1eta, lep2pt, lep2phi, lep2eta;
float neutrino1px, neutrino1py, neutrino1pz, neutrino2px, neutrino2py, neutrino2pz;

TTree *GetMiniTree(TFile *MiniTreeFile) {

	TTree *MiniTree = (TTree*) MiniTreeFile->Get("latino"); 

	MiniTree->SetBranchAddress("eventW",          &eventW);
	MiniTree->SetBranchAddress("metPfType1",      &metPfType1);
	MiniTree->SetBranchAddress("pfType1Met",      &pfType1Met);
	MiniTree->SetBranchAddress("pfType1Metphi",   &pfType1Metphi);
	MiniTree->SetBranchAddress("njet",            &njet);
	MiniTree->SetBranchAddress("channel",         &channel);

	MiniTree->SetBranchAddress("nbjet30csvv2l",   &nbjet30csvv2l);
	MiniTree->SetBranchAddress("nbjet30csvv2m",   &nbjet30csvv2m);
	MiniTree->SetBranchAddress("nbjet30csvv2t",   &nbjet30csvv2t);

	MiniTree->SetBranchAddress("mt2ll",           &mt2ll);
	MiniTree->SetBranchAddress("mt2bb",           &mt2bb);
	MiniTree->SetBranchAddress("mt2lblb",         &mt2lblb);
	MiniTree->SetBranchAddress("mt2bbtrue",       &mt2bbtrue);
	MiniTree->SetBranchAddress("mt2lblbcomb",     &mt2lblbcomb);
	MiniTree->SetBranchAddress("mt2lblbtrue",     &mt2lblbtrue);

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

	MiniTree->SetBranchAddress("lep1pt",          &lep1pt);
	MiniTree->SetBranchAddress("lep1phi",         &lep1phi);
	MiniTree->SetBranchAddress("lep1eta",         &lep1eta);
	MiniTree->SetBranchAddress("lep2pt",          &lep2pt);
	MiniTree->SetBranchAddress("lep2phi",         &lep2phi);
	MiniTree->SetBranchAddress("lep2eta",         &lep2eta);

	MiniTree->SetBranchAddress("neutrino1px",     &neutrino1px);
	MiniTree->SetBranchAddress("neutrino1py",     &neutrino1py);
	MiniTree->SetBranchAddress("neutrino1pz",     &neutrino1pz);
	MiniTree->SetBranchAddress("neutrino2px",     &neutrino2px);
	MiniTree->SetBranchAddress("neutrino2py",     &neutrino2py);
	MiniTree->SetBranchAddress("neutrino2pz",     &neutrino2pz);

	return MiniTree;

}

Mt2::Basic_Mt2_332_Calculator Mt2Calcolator;
Mt2::Basic_Mt2_Top_Calculator Mt2TopCalcolator;

struct Mt2Result{float Mt2ll; float Mt2lblb; float Mlb1; float Mlb2; float Neu1Px; float Neu1Py; float Neu2Px; float Neu2Py;};

TLorentzVector Lepton[2], Bottom[2], LepB[2][2];


TH2D *h_neutrino1_pxy[2];
TH2D *h_neutrino2_pxy[2];

TH1D *h_mt2ll_minitrees[2];
TH1D *h_mt2ll_minitrees_Integ[2];
TH1D *h_mt2ll_minitrees_Signif;
TH1D *h_mt2ll_S0_TW05[2];
TH1D *h_mt2ll_S0_TW05_Integ[2];
TH1D *h_mt2ll_S0_TW05_Signif;
TH1D *h_mt2lblb_minitrees[2];
TH1D *h_mt2lblb_minitrees_Integ[2];
TH1D *h_mt2lblb_minitrees_Signif;
TH1D *h_mt2lblb_minitrees_cut[2][50]; // cut en mt2ll
TH1D *h_mt2lblb_minitrees_Integ_cut[2][50];
TH1D *h_mt2lblb_minitrees_Signif_cut[50];
TH1D *h_mt2lblb_S0_TW05[2];
TH1D *h_mt2lblb_S0_TW05_Integ[2];
TH1D *h_mt2lblb_S0_TW05_Signif;
TH1D *h_mt2lblb_S0_TW05_cut[2][50];
TH1D *h_mt2lblb_S0_TW05_Integ_cut[2][50];
TH1D *h_mt2lblb_S0_TW05_Signif_cut[50];

float mt2ll_minitrees_Int[2];
float mt2lblb_minitrees_Int[2];
float mt2lblb_minitrees_Int_cut[2][50];
float mt2ll_S0_TW05_Int[2];
float mt2lblb_S0_TW05_Int[2];
float mt2lblb_S0_TW05_Int_cut[2][50];



int low_icut = 200;
int high_icut = 300;
int istep = 30;

int HistoColSignif[50];
TString Legends_mini[2][50];
TString Legends_S0_TW05[2][50];
int LineStyle[2]={1,2};


Mt2Result ComputeMt2Top(Mt2::LorentzVector& Lepton1,  Mt2::LorentzVector& Lepton2,
		Mt2::LorentzVector& Bottom1,  Mt2::LorentzVector& Bottom2,
		Mt2::TwoVector& ptmiss, double mEachInvisible, 
		int Mt2TopStrategy, float Mt2TopWeight) {

	cout << "ComputeMt2Top"<< endl;

	double Mt2W[2], Mt2Top[2], NeuAPx[2], NeuAPy[2], NeuBPx[2], NeuBPy[2];

	const double mT2Total1 = Mt2TopCalcolator.mt2_Top(Lepton1, Lepton2, Bottom1, Bottom2, ptmiss, mEachInvisible, Mt2TopStrategy, Mt2TopWeight);

	Mt2W[0] = Mt2TopCalcolator.getMt2W_atMt2Solution();
	Mt2Top[0] = Mt2TopCalcolator.getMt2Top_atMt2Solution();

	NeuAPx[0] = Mt2TopCalcolator.getPXInvisA_atMt2Solution();
	NeuAPy[0] = Mt2TopCalcolator.getPYInvisA_atMt2Solution();
	NeuBPx[0] = Mt2TopCalcolator.getPXInvisB_atMt2Solution();
	NeuBPy[0] = Mt2TopCalcolator.getPYInvisB_atMt2Solution();

	const double mT2Total2 = Mt2TopCalcolator.mt2_Top(Lepton1, Lepton2, Bottom2, Bottom1, ptmiss, mEachInvisible, Mt2TopStrategy, Mt2TopWeight);

	Mt2W[1] = Mt2TopCalcolator.getMt2W_atMt2Solution();
	Mt2Top[1] = Mt2TopCalcolator.getMt2Top_atMt2Solution();

	NeuAPx[1] = Mt2TopCalcolator.getPXInvisA_atMt2Solution();
	NeuAPy[1] = Mt2TopCalcolator.getPYInvisA_atMt2Solution();
	NeuBPx[1] = Mt2TopCalcolator.getPXInvisB_atMt2Solution();
	NeuBPy[1] = Mt2TopCalcolator.getPYInvisB_atMt2Solution();

	int BestPair = (mT2Total1<=mT2Total2) ? 0 : 1;

	Mt2Result ThisMt2Result;

	ThisMt2Result.Mt2ll = Mt2W[BestPair];
	ThisMt2Result.Mt2lblb = Mt2Top[BestPair];

	ThisMt2Result.Mlb1 = (Lepton[0] + Bottom[BestPair]).M();
	ThisMt2Result.Mlb2 = (Lepton[1] + Bottom[1-BestPair]).M();

	ThisMt2Result.Neu1Px = NeuAPx[BestPair];
	ThisMt2Result.Neu1Py = NeuAPy[BestPair];
	ThisMt2Result.Neu2Px = NeuBPx[BestPair];
	ThisMt2Result.Neu2Py = NeuBPy[BestPair];

	return ThisMt2Result;

}

void MT2top(bool TestStandardMt2 = false) {

	cout << "MT2top" << endl;

	//TString FileName[2] = {"/afs/cern.ch/user/p/palmerin/public/TTTo2L2Nu.root",
	//"/afs/cern.ch/user/p/palmerin/public/T2tt_mStop500-525-550_mLSP1to425-325to450-1to475.root"};

	TString FileName[2] = {"minitrees/minitreesLuca/TTTo2L2Nu.root",
		"minitrees/minitreesLuca/T2tt_mStop600-950_mLSP1to450.root"};

	cout << "FileName" << endl;

	TString Histoname[2]={"_top", "_stop"};

	for (int i=0; i<50; i++) {
		HistoColSignif[i] = i+2;
		
	}

	for (int dt = 0; dt<2; dt++) {
		for (int cut = low_icut; cut <= high_icut; cut += istep) {

			int icut = (cut-low_icut)/istep;
			char scut [50];
			sprintf(scut, "%.d", cut);
			Legends_mini[dt][icut] = "Mt2lblb (minitrees) with mt2ll<=";
			Legends_mini[dt][icut].Append(scut);
			Legends_mini[dt][icut].Append(Histoname[dt]);

		}
	}
	
	for (int dt = 0; dt<2; dt++) {
		for (int cut = low_icut; cut <= high_icut; cut += istep) {

			int icut = (cut-low_icut)/istep;
			char scut [50];
			sprintf(scut, "%.d", cut);
			Legends_S0_TW05[dt][icut] = "Mt2lblb (Strategy 0, TopWeight 0.5) with mt2ll<=";
			Legends_S0_TW05[dt][icut].Append(scut);
			Legends_S0_TW05[dt][icut].Append(Histoname[dt]);

		}
	}

	for (int dt = 0; dt<2; dt++) {

		h_neutrino1_pxy[dt] = new TH2D("h_neutrino1_pxy" + Histoname[dt],"h_neutrino1_pxy" + Histoname[dt], 300, 0, 3000, 300, 0, 3000); //2 histos para top y stop
		h_neutrino2_pxy[dt] = new TH2D("h_neutrino2_pxy" + Histoname[dt],"h_neutrino2_pxy" + Histoname[dt], 300, 0, 3000, 300, 0, 3000); //2 histos para top y stop

		h_mt2ll_minitrees[dt] = new TH1D("h_mt2ll_minitrees" + Histoname[dt],"h_mt2ll_minitrees" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2ll_S0_TW05[dt] = new TH1D("h_mt2ll_S0_TW05" + Histoname[dt],"h_mt2ll_S0_TW05" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2lblb_minitrees[dt] = new TH1D("h_mt2lblb_minitrees" + Histoname[dt],"h_mt2lblb_minitrees" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2lblb_S0_TW05[dt] = new TH1D("h_mt2lblb_S0_TW05" + Histoname[dt],"h_mt2lblb_S0_TW05" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop

		h_mt2ll_minitrees_Integ[dt] = new TH1D("h_mt2ll_minitrees_Integ" + Histoname[dt],"h_mt2ll_minitrees_Integ" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2ll_S0_TW05_Integ[dt] = new TH1D("h_mt2ll_S0_TW05_Integ" + Histoname[dt],"h_mt2ll_S0_TW05_Integ" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2lblb_minitrees_Integ[dt] = new TH1D("h_mt2lblb_minitrees_Integ" + Histoname[dt],"h_mt2lblb_minitrees_Integ" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2lblb_S0_TW05_Integ[dt] = new TH1D("h_mt2lblb_S0_TW05_Integ" + Histoname[dt],"h_mt2lblb_S0_TW05_Integ" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop

		h_mt2ll_minitrees_Signif = new TH1D("h_mt2ll_minitrees_Signif" + Histoname[dt],"h_mt2ll_minitrees_Signif" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2ll_S0_TW05_Signif = new TH1D("h_mt2ll_S0_TW05_Signif" + Histoname[dt],"h_mt2ll_S0_TW05_Signif" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2lblb_minitrees_Signif = new TH1D("h_mt2lblb_minitrees_Signif" + Histoname[dt],"h_mt2lblb_minitrees_Signif" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop
		h_mt2lblb_S0_TW05_Signif = new TH1D("h_mt2lblb_S0_TW05_Signif" + Histoname[dt],"h_mt2lblb_S0_TW05_Signif" + Histoname[dt], 300, 0, 3000); //2 histos para top y stop

		for (int cut = low_icut; cut <= high_icut; cut += istep) {

			int icut = (cut-low_icut)/istep;
			char scut [50];
			sprintf(scut, "_%.f", cut);

			h_mt2lblb_minitrees_cut[dt][icut] = new TH1D("h_mt2lblb_minitrees_cut" + Histoname[dt] + scut,"h_mt2lblb_minitrees_cut" + Histoname[dt] + scut, 300, 0, 3000); //2 histos para top y stop

			h_mt2lblb_S0_TW05_cut[dt][icut] = new TH1D("h_mt2lblb_S0_TW05_cut" + Histoname[dt] + scut,"h_mt2lblb_S0_TW05_cut" + Histoname[dt] + scut, 300, 0, 3000); //2 histos para top y stop

			h_mt2lblb_minitrees_Integ_cut[dt][icut] = new TH1D("h_mt2lblb_minitrees_Integ_cut" + Histoname[dt] + scut,"h_mt2lblb_minitrees_Integ_cut" + Histoname[dt] + scut, 300, 0, 3000); //2 histos para top y stop
			h_mt2lblb_S0_TW05_Integ_cut[dt][icut] = new TH1D("h_mt2lblb_S0_TW05_Integ_cut" + Histoname[dt] + scut,"h_mt2lblb_S0_TW05_Integ_cut" + Histoname[dt] + scut, 300, 0, 3000); //2 histos para top y stop
			h_mt2lblb_minitrees_Signif_cut[icut] = new TH1D("h_mt2lblb_minitrees_Signif_cut" + Histoname[dt] + scut,"h_mt2lblb_minitrees_Signif_cut" + Histoname[dt] + scut, 300, 0, 3000); //2 histos para top y stop
			h_mt2lblb_S0_TW05_Signif_cut[icut] = new TH1D("h_mt2lblb_S0_TW05_Signif_cut" + Histoname[dt] + scut,"h_mt2lblb_S0_TW05_Signif_cut" + Histoname[dt] + scut, 300, 0, 3000); //2 histos para top y stop
			cout << "new TH1D" << endl;

		}



		TFile *MiniTreeFile = TFile::Open(FileName[dt]);

		TTree *MiniTree = GetMiniTree(MiniTreeFile);

		Int_t nentries = (Int_t) MiniTree->GetEntries();


		for (Int_t i = 0; i<nentries; i++) {

			cout << "eventW" <<Histoname[dt] << eventW << endl;
			MiniTree->GetEntry(i);

			// Apply ttbar selection
			if (njet<2) continue;
			if (nbjet30csvv2m<1) continue;

			// Fill info for mt2
			Lepton[0].SetPtEtaPhiM(lep1pt, lep1eta, lep1phi,  0.1);
			Lepton[1].SetPtEtaPhiM(lep2pt, lep2eta, lep2phi,  0.1);

			Bottom[0].SetPtEtaPhiM(tjet1pt, tjet1eta, tjet1phi,  5.);
			Bottom[1].SetPtEtaPhiM(tjet2pt, tjet2eta, tjet2phi,  5.);

			LepB[0][0] = Lepton[0] + Bottom[0]; LepB[1][0] = Lepton[1] + Bottom[1];
			LepB[0][1] = Lepton[0] + Bottom[1]; LepB[1][1] = Lepton[1] + Bottom[0];

			double m_invis_mass=0.;
			Mt2::TwoVector pT_Miss(pfType1Met*cos(pfType1Metphi), pfType1Met*sin(pfType1Metphi));

			Mt2::LorentzVector Lepton1; Lepton1.setVectM(Lepton[0].Px(), Lepton[0].Py(), Lepton[0].Pz(), Lepton[0].M());
			Mt2::LorentzVector Lepton2; Lepton2.setVectM(Lepton[1].Px(), Lepton[1].Py(), Lepton[1].Pz(), Lepton[1].M());

			Mt2::LorentzVector Bottom1; Bottom1.setVectM(Bottom[0].Px(), Bottom[0].Py(), Bottom[0].Pz(), Bottom[0].M());
			Mt2::LorentzVector Bottom2; Bottom2.setVectM(Bottom[1].Px(), Bottom[1].Py(), Bottom[1].Pz(), Bottom[1].M());

			//if (TestStandardMt2) {

			//Mt2::LorentzTransverseVector visA = Lepton1.getLorentzTransverseVector();
			//Mt2::LorentzTransverseVector visB = Lepton2.getLorentzTransverseVector();
			Mt2::LorentzTransverseVector visA(Mt2::TwoVector(lep1pt*cos(lep1phi), lep1pt*sin(lep1phi)), 0.1);
			Mt2::LorentzTransverseVector visB(Mt2::TwoVector(lep2pt*cos(lep2phi), lep2pt*sin(lep2phi)), 0.1);

			// Let´s cross check mt2ll first  
			const double mT2ll = Mt2Calcolator.mt2_332(visA, visB, pT_Miss, m_invis_mass);	
			Mt2Result Mt2Strategy0p0 = ComputeMt2Top(Lepton1, Lepton2, Bottom1, Bottom2, pT_Miss, m_invis_mass, 0, 0.);


			//std::cout << "Test mt2ll " << mt2ll << " " << mT2ll << " " << Mt2Strategy0p0.Mt2ll << std::endl;

			// Then cross check mt2lblb
			Mt2::LorentzTransverseVector visUA(Mt2::TwoVector(LepB[0][0].Px(), LepB[0][0].Py()), LepB[0][0].M());
			Mt2::LorentzTransverseVector visUB(Mt2::TwoVector(LepB[1][0].Px(), LepB[1][0].Py()), LepB[1][0].M());

			const double mT2lblb = Mt2Calcolator.mt2_332(visUA, visUB, pT_Miss, m_invis_mass);

			Mt2Result Mt2Strategy0p999 = ComputeMt2Top(Lepton1, Lepton2, Bottom1, Bottom2, pT_Miss, m_invis_mass, 0, 999.);
			Mt2Result Mt2Strategy0p05 = ComputeMt2Top(Lepton1, Lepton2, Bottom1, Bottom2, pT_Miss, m_invis_mass, 0, 500.0);

			//cout << "peso 999: " << Mt2Strategy0p999.Mt2lblb << endl;
			//std::cout << "Test mt2lblb " << mt2lblb << " " << mt2lblbcomb << " " << mT2lblb << " " << Mt2Strategy0p999.Mt2lblb << std::endl;
			//std::cout << "Test mlb     " << mlb1 << " " << mlb2 << " " << Mt2Strategy0p999.Mlb1 << " " << Mt2Strategy0p999.Mlb2 << std::endl;

			h_neutrino1_pxy[dt] -> Fill(neutrino1px-Mt2Strategy0p999.Neu1Px, neutrino1py-Mt2Strategy0p999.Neu1Py, eventW);
			h_neutrino2_pxy[dt] -> Fill(fabs(neutrino2px-Mt2Strategy0p999.Neu2Px), fabs(neutrino2py-Mt2Strategy0p999.Neu2Py), eventW);

			if (mlb1 <= 160.0 && mlb2 <= 160.0) {

				h_mt2ll_minitrees[dt] -> Fill(mt2ll, eventW);
				h_mt2ll_S0_TW05[dt] -> Fill(Mt2Strategy0p05.Mt2ll, eventW);
				h_mt2lblb_minitrees[dt] -> Fill(mt2lblb, eventW);
				h_mt2lblb_S0_TW05[dt] -> Fill(Mt2Strategy0p05.Mt2lblb, eventW);
				//cout << "mt2ll minitrees: " << mt2ll << endl;
				//cout << "mt2lblb minitrees: " << mt2lblb << endl;
				//cout << "Mt2ll S0_TW05: " <<Mt2Strategy0p05.Mt2ll << endl;
				//cout << "Mt2lblb S0_TW05: " << Mt2Strategy0p05.Mt2lblb << endl;
			}

			for (int cut = low_icut; cut <= high_icut; cut += istep) {

				int icut = (cut-low_icut)/istep;
				if (mt2ll <= float(cut)) {

					h_mt2lblb_minitrees_cut[dt][icut] -> Fill(mt2lblb, eventW);
					h_mt2lblb_S0_TW05_cut[dt][icut] -> Fill(Mt2Strategy0p05.Mt2lblb, eventW);

				}


			}
		}


	}

	for (int dt  = 0; dt<2; dt++) {


		mt2ll_minitrees_Int[dt] = h_mt2ll_minitrees[dt]->Integral(1, 3001);
		h_mt2ll_minitrees[dt]->Scale(1./mt2ll_minitrees_Int[dt]); // normalization of the histogram

		mt2ll_S0_TW05_Int[dt] = h_mt2ll_S0_TW05[dt]->Integral(1, 3001);
		h_mt2ll_S0_TW05[dt]->Scale(1./mt2ll_S0_TW05_Int[dt]); // normalization of the histogram

		mt2lblb_minitrees_Int[dt] = h_mt2lblb_minitrees[dt]->Integral(1, 3001);
		h_mt2lblb_minitrees[dt]->Scale(1./mt2lblb_minitrees_Int[dt]); // normalization of the histogram

		mt2lblb_S0_TW05_Int[dt] = h_mt2lblb_S0_TW05[dt]->Integral(1, 3001);
		h_mt2lblb_S0_TW05[dt]->Scale(1./mt2lblb_S0_TW05_Int[dt]); // normalization of the histogram

		for (int cut = low_icut; cut <= high_icut; cut += istep) {

			int icut = (cut-low_icut)/istep;

			mt2lblb_minitrees_Int_cut[dt][icut] = h_mt2lblb_minitrees_cut[dt][icut]->Integral(1, 3001);
			h_mt2lblb_minitrees_cut[dt][icut]->Scale(1./mt2lblb_minitrees_Int_cut[dt][icut]); // normalization of the histogram

			mt2lblb_S0_TW05_Int_cut[dt][icut] = h_mt2lblb_S0_TW05_cut[dt][icut]->Integral(1, 3001);
			h_mt2lblb_S0_TW05_cut[dt][icut]->Scale(1./mt2lblb_S0_TW05_Int_cut[dt][icut]); // normalization of the histogram


		}


	}


	int nBinsX = h_mt2ll_minitrees[0]->GetNbinsX();

	for (int hf = 0; hf<2; hf++) { // stop or top


		for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

			float recursive_integral_mt2ll_minitrees = h_mt2ll_minitrees[hf]->Integral(ib, 3001);
			h_mt2ll_minitrees_Integ[hf]->SetBinContent(ib, recursive_integral_mt2ll_minitrees); // asigna el valor ThisBinContent al bin ib


			float recursive_integral_mt2ll_S0_TW05 = h_mt2ll_S0_TW05[hf]->Integral(ib, 3001);
			h_mt2ll_S0_TW05_Integ[hf]->SetBinContent(ib, recursive_integral_mt2ll_S0_TW05); // asigna el valor ThisBinContent al bin ib



			float recursive_integral_mt2lblb_minitrees = h_mt2lblb_minitrees[hf]->Integral(ib, 3001);
			h_mt2lblb_minitrees_Integ[hf]->SetBinContent(ib, recursive_integral_mt2lblb_minitrees); // asigna el valor ThisBinContent al bin ib


			float recursive_integral_mt2lblb_S0_TW05 = h_mt2lblb_S0_TW05[hf]->Integral(ib, 3001);
			h_mt2lblb_S0_TW05_Integ[hf]->SetBinContent(ib, recursive_integral_mt2lblb_S0_TW05); // asigna el valor ThisBinContent al bin ib

			for (int cut = low_icut; cut <= high_icut; cut += istep) {

				int icut = (cut-low_icut)/istep;

				float recursive_integral_mt2lblb_minitrees_cut = h_mt2lblb_minitrees_cut[hf][icut]->Integral(ib, 3001);
				h_mt2lblb_minitrees_Integ_cut[hf][icut]->SetBinContent(ib, recursive_integral_mt2lblb_minitrees_cut); // asigna el valor ThisBinContent al bin ib

				float recursive_integral_mt2lblb_S0_TW05_cut = h_mt2lblb_S0_TW05_cut[hf][icut]->Integral(ib, 3001);
				h_mt2lblb_S0_TW05_Integ_cut[hf][icut]->SetBinContent(ib, recursive_integral_mt2lblb_S0_TW05_cut); // asigna el valor ThisBinContent al bin ib


			}
		}
	}



	for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins

		float	topBackground_mt2ll_minitrees = mt2ll_minitrees_Int[0] * h_mt2ll_minitrees_Integ[0]->GetBinContent(ib); 
		float	stopEvents_mt2ll_minitrees = mt2ll_minitrees_Int[1] * h_mt2ll_minitrees_Integ[1]->GetBinContent(ib);
		if (topBackground_mt2ll_minitrees + stopEvents_mt2ll_minitrees <= 0.) continue;
		float	significance_mt2ll_minitrees = stopEvents_mt2ll_minitrees/(std::sqrt(stopEvents_mt2ll_minitrees+topBackground_mt2ll_minitrees));
		h_mt2ll_minitrees_Signif->SetBinContent(ib,significance_mt2ll_minitrees); 

		float	topBackground_mt2lblb_minitrees = mt2lblb_minitrees_Int[0] * h_mt2lblb_minitrees_Integ[0]->GetBinContent(ib); 
		float	stopEvents_mt2lblb_minitrees = mt2lblb_minitrees_Int[1] * h_mt2lblb_minitrees_Integ[1]->GetBinContent(ib);
		if (topBackground_mt2lblb_minitrees + stopEvents_mt2lblb_minitrees <= 0.) continue;
		float	significance_mt2lblb_minitrees = stopEvents_mt2lblb_minitrees/(std::sqrt(stopEvents_mt2lblb_minitrees+topBackground_mt2lblb_minitrees));
		h_mt2lblb_minitrees_Signif->SetBinContent(ib,significance_mt2lblb_minitrees); 

		float	topBackground_mt2ll_S0_TW05 = mt2ll_S0_TW05_Int[0] * h_mt2ll_S0_TW05_Integ[0]->GetBinContent(ib); 
		float	stopEvents_mt2ll_S0_TW05 = mt2ll_S0_TW05_Int[1] * h_mt2ll_S0_TW05_Integ[1]->GetBinContent(ib);
		if (topBackground_mt2ll_S0_TW05 + stopEvents_mt2ll_S0_TW05 <= 0.) continue;
		float	significance_mt2ll_S0_TW05 = stopEvents_mt2ll_S0_TW05/(std::sqrt(stopEvents_mt2ll_S0_TW05+topBackground_mt2ll_S0_TW05));
		h_mt2ll_S0_TW05_Signif->SetBinContent(ib,significance_mt2ll_S0_TW05); 

		float	topBackground_mt2lblb_S0_TW05 = mt2lblb_S0_TW05_Int[0] * h_mt2lblb_S0_TW05_Integ[0]->GetBinContent(ib); 
		float	stopEvents_mt2lblb_S0_TW05 = mt2lblb_S0_TW05_Int[1] * h_mt2lblb_S0_TW05_Integ[1]->GetBinContent(ib);
		if (topBackground_mt2lblb_S0_TW05 + stopEvents_mt2lblb_S0_TW05 <= 0.) continue;
		float	significance_mt2lblb_S0_TW05 = stopEvents_mt2lblb_S0_TW05/(std::sqrt(stopEvents_mt2lblb_S0_TW05+topBackground_mt2lblb_S0_TW05));
		h_mt2lblb_S0_TW05_Signif->SetBinContent(ib,significance_mt2lblb_S0_TW05); 

		for (int cut = low_icut; cut <= high_icut; cut += istep) {

			int icut = (cut-low_icut)/istep;

			float	topBackground_mt2lblb_minitrees_cut = mt2lblb_minitrees_Int_cut[0][icut] * h_mt2lblb_minitrees_Integ_cut[0][icut]->GetBinContent(ib); 
			float	stopEvents_mt2lblb_minitrees_cut = mt2lblb_minitrees_Int_cut[1][icut] * h_mt2lblb_minitrees_Integ_cut[1][icut]->GetBinContent(ib);
			if (topBackground_mt2lblb_minitrees_cut + stopEvents_mt2lblb_minitrees_cut <= 0.) continue;
			float	significance_mt2lblb_minitrees_cut = stopEvents_mt2lblb_minitrees_cut/(std::sqrt(stopEvents_mt2lblb_minitrees_cut+topBackground_mt2lblb_minitrees_cut));
			h_mt2lblb_minitrees_Signif_cut[icut]->SetBinContent(ib,significance_mt2lblb_minitrees_cut); 


			float	topBackground_mt2lblb_S0_TW05_cut = mt2lblb_S0_TW05_Int_cut[0][icut] * h_mt2lblb_S0_TW05_Integ_cut[0][icut]->GetBinContent(ib); 
			float	stopEvents_mt2lblb_S0_TW05_cut = mt2lblb_S0_TW05_Int_cut[1][icut] * h_mt2lblb_S0_TW05_Integ_cut[1][icut]->GetBinContent(ib);
			if (topBackground_mt2lblb_S0_TW05_cut + stopEvents_mt2lblb_S0_TW05_cut <= 0.) continue;
			float	significance_mt2lblb_S0_TW05_cut = stopEvents_mt2lblb_S0_TW05_cut/(std::sqrt(stopEvents_mt2lblb_S0_TW05_cut+topBackground_mt2lblb_S0_TW05_cut));
			h_mt2lblb_S0_TW05_Signif_cut[icut]->SetBinContent(ib,significance_mt2lblb_S0_TW05_cut); 

		}

	}




	// Here we can study different mt2ll-mt2lblb minimisations

	TCanvas *CC = new TCanvas("CC", "", 1000, 800);
	CC->Divide(2, 2);
	TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1->SetGridx(); CC1->SetGridy(); 
	TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2->SetGridx(); CC2->SetGridy(); 
	TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mlblbtrue
	CC3->SetGridx(); CC3->SetGridy(); 
	TPad *CC4 = (TPad*)CC->GetPad(4); //Ahi va el h_mt2mlblbtrue
	CC4->SetGridx(); CC4->SetGridy(); 

	TCanvas *CCmt2ll = new TCanvas("CCmt2ll", "", 1200, 800);
	CCmt2ll->Divide(2, 3);
	TPad *CC1mt2ll = (TPad*)CCmt2ll->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CC1mt2ll->SetLogy(); CC1mt2ll->SetGridx(); CC1mt2ll->SetGridy(); 
	TPad *CC2mt2ll = (TPad*)CCmt2ll->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CC2mt2ll->SetLogy(); CC2mt2ll->SetGridx(); CC2mt2ll->SetGridy(); 
	TPad *CC3mt2ll = (TPad*)CCmt2ll->GetPad(3); //Ahi va el h_mt2mlblbtrue
	CC3mt2ll->SetLogy(); CC3mt2ll->SetGridx(); CC3mt2ll->SetGridy(); 
	TPad *CC4mt2ll = (TPad*)CCmt2ll->GetPad(4); //Ahi va el h_mt2mlblbtrue
	CC4mt2ll->SetLogy(); CC4mt2ll->SetGridx(); CC4mt2ll->SetGridy(); 
	TPad *CC5mt2ll = (TPad*)CCmt2ll->GetPad(5); //Ahi va el h_mt2mlblbtrue
	CC5mt2ll->SetGridx(); CC5mt2ll->SetGridy(); 
	TPad *CC6mt2ll = (TPad*)CCmt2ll->GetPad(6); //Ahi va el h_mt2mlblbtrue
	CC6mt2ll->SetGridx(); CC6mt2ll->SetGridy(); 


	TCanvas *CCcut = new TCanvas("CCcut", "", 1500, 800);
	CCcut->Divide(2, 2);
	TPad *CCcut1 = (TPad*)CCcut->GetPad(1); //Ahi va el h_mt2mlblbtrue
	CCcut1->SetLogy(); CCcut1->SetGridx(); CCcut1->SetGridy(); 
	TPad *CCcut2 = (TPad*)CCcut->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
	CCcut2->SetLogy(); CCcut2->SetGridx(); CCcut2->SetGridy(); 
	TPad *CCcut3 = (TPad*)CCcut->GetPad(3); //Ahi va el h_mt2mlblbtrue
	CCcut3->SetGridx(); CCcut3->SetGridy(); 

	int HistoCol[2] = {4, 2};

	TString Option = "histo";

	TLegend *leg1 = new TLegend(0.4,0.1,0.7,0.4);
	leg1->AddEntry(h_mt2ll_minitrees[0],"minitree top","l");
	leg1->AddEntry(h_mt2ll_minitrees[1],"minitree stop","l");
	leg1->AddEntry(h_mt2ll_S0_TW05[0],"Strategy=0, TopWeight=0.5 top","l");
	leg1->AddEntry(h_mt2ll_S0_TW05[1],"Strategy=0, TopWeight=0.5 stop","l");

	TLegend *leg2 = new TLegend(0.6,0.1,0.9,0.4);
	leg2->AddEntry(h_mt2lblb_minitrees[0],"minitree top","l");
	leg2->AddEntry(h_mt2lblb_minitrees[1],"minitree stop","l");
	leg2->AddEntry(h_mt2lblb_S0_TW05[0],"Strategy=0, TopWeight=0.5 top","l");
	leg2->AddEntry(h_mt2lblb_S0_TW05[1],"Strategy=0, TopWeight=0.5 stop","l");

	TLegend *leg3 = new TLegend(0.6,0.1,0.9,0.4);
	leg3->AddEntry(h_mt2ll_minitrees_Integ[0],"minitree top","l");
	leg3->AddEntry(h_mt2ll_minitrees_Integ[1],"minitree stop","l");
	leg3->AddEntry(h_mt2ll_S0_TW05_Integ[0],"Strategy=0, TopWeight=0.5 top","l");
	leg3->AddEntry(h_mt2ll_S0_TW05_Integ[1],"Strategy=0, TopWeight=0.5 stop","l");

	TLegend *leg4 = new TLegend(0.6,0.1,0.9,0.4);
	leg4->AddEntry(h_mt2lblb_minitrees_Integ[0],"minitree top","l");
	leg4->AddEntry(h_mt2lblb_minitrees_Integ[1],"minitree stop","l");
	leg4->AddEntry(h_mt2lblb_S0_TW05_Integ[0],"Strategy=0, TopWeight=0.5 top","l");
	leg4->AddEntry(h_mt2lblb_S0_TW05_Integ[1],"Strategy=0, TopWeight=0.5 stop","l");

	TLegend *leg5 = new TLegend(0.6,0.3,0.9,0.6);
	leg5->AddEntry(h_mt2ll_minitrees_Signif,"minitree","l");
	leg5->AddEntry(h_mt2ll_S0_TW05_Signif,"Strategy=0, TopWeight=0.5","l");

	TLegend *leg6 = new TLegend(0.4,0.1,0.7,0.4);
	leg6->AddEntry(h_mt2lblb_minitrees_Signif,"minitree","l");
	leg6->AddEntry(h_mt2lblb_S0_TW05_Signif,"Strategy=0, TopWeight=0.5","l");


	TLegend *legcut1 = new TLegend(0.1,0.1,0.5,0.5);
	TLegend *legcut2 = new TLegend(0.1,0.1,0.5,0.5);
	TLegend *legcut3 = new TLegend(0.4,0.5,0.9,0.9);

	legcut1->AddEntry(h_mt2lblb_minitrees_cut[1][0], Legends_mini[1][0],"l");
	legcut1->AddEntry(h_mt2lblb_S0_TW05_cut[1][0], Legends_S0_TW05[1][0],"l");
	legcut2->AddEntry(h_mt2lblb_minitrees_Integ_cut[1][0], Legends_mini[1][0],"l");
	legcut2->AddEntry(h_mt2lblb_S0_TW05_Integ_cut[1][0], Legends_S0_TW05[1][0],"l");

	for (int cut = low_icut; cut <= high_icut; cut += istep) {

		int icut = (cut-low_icut)/istep;

		legcut1->AddEntry(h_mt2lblb_minitrees_cut[0][icut], Legends_mini[0][icut],"l");
		legcut1->AddEntry(h_mt2lblb_S0_TW05_cut[0][icut], Legends_S0_TW05[0][icut],"l");


		legcut2->AddEntry(h_mt2lblb_minitrees_Integ_cut[0][icut], Legends_mini[0][icut],"l");
		legcut2->AddEntry(h_mt2lblb_S0_TW05_Integ_cut[0][icut], Legends_S0_TW05[0][icut],"l");


		legcut3->AddEntry(h_mt2lblb_minitrees_Signif_cut[icut], Legends_mini[0][icut],"l");
		legcut3->AddEntry(h_mt2lblb_S0_TW05_Signif_cut[icut], Legends_S0_TW05[0][icut],"l");

	}


	for (int dt = 0; dt < 2; dt++) {

		h_mt2ll_minitrees[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2ll_minitrees[dt] -> SetLineStyle(1);
		h_mt2ll_minitrees[dt] -> SetLineWidth(2);
		h_mt2ll_S0_TW05[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2ll_S0_TW05[dt] -> SetLineStyle(2);
		h_mt2ll_S0_TW05[dt] -> SetLineWidth(2);

		h_mt2lblb_minitrees[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblb_minitrees[dt] -> SetLineStyle(1);
		h_mt2lblb_minitrees[dt] -> SetLineWidth(2);
		h_mt2lblb_S0_TW05[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblb_S0_TW05[dt] -> SetLineStyle(2);
		h_mt2lblb_S0_TW05[dt] -> SetLineWidth(2);

		h_mt2ll_minitrees_Integ[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2ll_minitrees_Integ[dt] -> SetLineStyle(1);
		h_mt2ll_minitrees_Integ[dt] -> SetLineWidth(2);
		h_mt2ll_S0_TW05_Integ[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2ll_S0_TW05_Integ[dt] -> SetLineStyle(2);
		h_mt2ll_S0_TW05_Integ[dt] -> SetLineWidth(2);

		h_mt2lblb_minitrees_Integ[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblb_minitrees_Integ[dt] -> SetLineStyle(1);
		h_mt2lblb_minitrees_Integ[dt] -> SetLineWidth(2);
		h_mt2lblb_S0_TW05_Integ[dt] -> SetLineColor(HistoCol[dt]);
		h_mt2lblb_S0_TW05_Integ[dt] -> SetLineStyle(2);
		h_mt2lblb_S0_TW05_Integ[dt] -> SetLineWidth(2);

		h_mt2ll_minitrees_Signif -> SetLineColor(1);
		h_mt2ll_minitrees_Signif -> SetLineStyle(1);
		h_mt2ll_minitrees_Signif -> SetLineWidth(2);
		h_mt2ll_S0_TW05_Signif -> SetLineColor(1);
		h_mt2ll_S0_TW05_Signif -> SetLineStyle(2);
		h_mt2ll_S0_TW05_Signif -> SetLineWidth(2);

		h_mt2lblb_minitrees_Signif -> SetLineColor(1);
		h_mt2lblb_minitrees_Signif -> SetLineStyle(1);
		h_mt2lblb_minitrees_Signif -> SetLineWidth(2);
		h_mt2lblb_S0_TW05_Signif -> SetLineColor(1);
		h_mt2lblb_S0_TW05_Signif -> SetLineStyle(2);
		h_mt2lblb_S0_TW05_Signif -> SetLineWidth(2);


		for (int cut = low_icut; cut <= high_icut; cut += istep) {

			int icut = (cut-low_icut)/istep;


			h_mt2lblb_minitrees_cut[dt][icut] -> SetLineColor(HistoColSignif[icut]);
			h_mt2lblb_minitrees_Integ_cut[dt][icut] -> SetLineColor(HistoColSignif[icut]);
			h_mt2lblb_minitrees_Signif_cut[icut] -> SetLineColor(HistoColSignif[icut]);

			h_mt2lblb_S0_TW05_cut[dt][icut] ->SetLineStyle(LineStyle[dt]);
			h_mt2lblb_S0_TW05_Integ_cut[dt][icut] ->SetLineStyle(LineStyle[dt]);
			h_mt2lblb_S0_TW05_Signif_cut[icut] -> SetLineStyle(1);


		}


		CC->cd(dt+1); // se pone en el TPad 1 

		h_neutrino1_pxy[dt]->GetXaxis()->SetRangeUser(1, 300);
		h_neutrino1_pxy[dt]->GetXaxis()->SetTitle("neutrino1 px_minitree - px_calc");
		h_neutrino1_pxy[dt]->GetYaxis()->SetRangeUser(1, 300);
		h_neutrino1_pxy[dt]->GetYaxis()->SetTitle("neutrino1 py_minitree - py_calc");
		h_neutrino1_pxy[dt]->SetTitle("neutrino1" + Histoname[dt]);
		h_neutrino1_pxy[dt]->DrawCopy("box");

		CC->cd(dt+3); // se pone en el TPad 1 

		h_neutrino2_pxy[dt]->GetXaxis()->SetRangeUser(1, 300);
		h_neutrino2_pxy[dt]->GetXaxis()->SetTitle("neutrino2 px_minitree - px_calc");
		h_neutrino2_pxy[dt]->GetYaxis()->SetRangeUser(1, 300);
		h_neutrino2_pxy[dt]->GetYaxis()->SetTitle("neutrino2 py_minitree - py_calc");
		h_neutrino2_pxy[dt]->SetTitle("neutrino2" + Histoname[dt]);
		h_neutrino2_pxy[dt]->DrawCopy("box");






		CCmt2ll->cd(1); // se pone en el TPad 1 

		h_mt2ll_minitrees[dt]->GetXaxis()->SetRangeUser(1, 600);
		h_mt2ll_minitrees[dt]->GetXaxis()->SetTitle("Mt2ll");
		h_mt2ll_minitrees[dt]->SetTitle("Mt2ll");
		h_mt2ll_minitrees[dt]->DrawCopy(Option);
		h_mt2ll_S0_TW05[dt]->GetXaxis()->SetRangeUser(1, 600);
		h_mt2ll_S0_TW05[dt]->DrawCopy("histosame");

		leg1->Draw();

		CCmt2ll->cd(2); // se pone en el TPad 1 

		h_mt2lblb_minitrees[dt]->GetXaxis()->SetRangeUser(1, 600);
		h_mt2lblb_minitrees[dt]->GetXaxis()->SetTitle("Mt2lblb");
		h_mt2lblb_minitrees[dt]->SetTitle("Mt2lblb");
		h_mt2lblb_minitrees[dt]->DrawCopy(Option);
		h_mt2lblb_S0_TW05[dt]->GetXaxis()->SetRangeUser(1, 600);
		h_mt2lblb_S0_TW05[dt]->DrawCopy("histosame");

		leg2->Draw();

		CCmt2ll->cd(3); // se pone en el TPad 1 

		h_mt2ll_minitrees_Integ[dt]->GetXaxis()->SetRangeUser(1, 600);
		h_mt2ll_minitrees_Integ[dt]->GetXaxis()->SetTitle("Mt2ll Integ");
		h_mt2ll_minitrees_Integ[dt]->SetTitle("Mt2ll Integ ");
		h_mt2ll_minitrees_Integ[dt]->DrawCopy(Option);
		h_mt2ll_S0_TW05_Integ[dt]->GetXaxis()->SetRangeUser(1, 600);
		h_mt2ll_S0_TW05_Integ[dt]->DrawCopy("histosame");

		leg3->Draw();

		CCmt2ll->cd(4); // se pone en el TPad 1 

		h_mt2lblb_minitrees_Integ[dt]->GetXaxis()->SetRangeUser(1, 600);
		h_mt2lblb_minitrees_Integ[dt]->GetXaxis()->SetTitle("Mt2lblb Integ");
		h_mt2lblb_minitrees_Integ[dt]->SetTitle("Mt2lblb Integ");
		h_mt2lblb_minitrees_Integ[dt]->DrawCopy(Option);
		h_mt2lblb_S0_TW05_Integ[dt]->GetXaxis()->SetRangeUser(1, 600);
		h_mt2lblb_S0_TW05_Integ[dt]->DrawCopy("histosame");

		leg4->Draw();

		CCmt2ll->cd(5); // se pone en el TPad 1 

		h_mt2ll_minitrees_Signif->GetXaxis()->SetRangeUser(1, 600);
		h_mt2ll_minitrees_Signif->GetXaxis()->SetTitle("Mt2ll Significance");
		h_mt2ll_minitrees_Signif->SetTitle("Mt2ll Significance");
		h_mt2ll_minitrees_Signif->DrawCopy(Option);
		h_mt2ll_S0_TW05_Signif->GetXaxis()->SetRangeUser(1, 600);
		h_mt2ll_S0_TW05_Signif->DrawCopy("histosame");

		leg5->Draw();

		CCmt2ll->cd(6); // se pone en el TPad 1 

		h_mt2lblb_minitrees_Signif->GetXaxis()->SetRangeUser(1, 600);
		h_mt2lblb_minitrees_Signif->GetXaxis()->SetTitle("Mt2lblb Significance");
		h_mt2lblb_minitrees_Signif->SetTitle("Mt2lblb Significance");
		h_mt2lblb_minitrees_Signif->DrawCopy(Option);
		h_mt2lblb_S0_TW05_Signif->GetXaxis()->SetRangeUser(1, 600);
		h_mt2lblb_S0_TW05_Signif->DrawCopy("histosame");

		leg6->Draw();


		for (int cut = low_icut; cut <= high_icut; cut += istep) {

			int icut = (cut-low_icut)/istep;

			CCcut->cd(1); // se pone en el TPad 1 

			h_mt2lblb_minitrees_cut[dt][icut]->GetXaxis()->SetRangeUser(1, 900);
			h_mt2lblb_minitrees_cut[dt][icut]->GetXaxis()->SetTitle("Mt2lblb");
			h_mt2lblb_minitrees_cut[dt][icut]->SetTitle("Mt2lblb with cut in Mt2ll");
			h_mt2lblb_minitrees_cut[dt][icut]->DrawCopy(Option);
			h_mt2lblb_S0_TW05_cut[dt][icut]->GetXaxis()->SetRangeUser(1, 900);
			h_mt2lblb_S0_TW05_cut[dt][icut]->DrawCopy("histosame");

			legcut1->Draw();


			CCcut->cd(2); // se pone en el TPad 1 

			h_mt2lblb_minitrees_Integ_cut[dt][icut]->GetXaxis()->SetRangeUser(1, 900);
			h_mt2lblb_minitrees_Integ_cut[dt][icut]->GetXaxis()->SetTitle("Mt2lblb Integ");
			h_mt2lblb_minitrees_Integ_cut[dt][icut]->SetTitle("Mt2lblb Integ with cut in Mt2ll");
			h_mt2lblb_minitrees_Integ_cut[dt][icut]->DrawCopy(Option);
			h_mt2lblb_S0_TW05_Integ_cut[dt][icut]->GetXaxis()->SetRangeUser(1, 900);
			h_mt2lblb_S0_TW05_Integ_cut[dt][icut]->DrawCopy("histosame");

			legcut2->Draw();


			CCcut->cd(3); // se pone en el TPad 1 

			h_mt2lblb_minitrees_Signif_cut[icut]->GetXaxis()->SetRangeUser(1, 900);
			h_mt2lblb_minitrees_Signif_cut[icut]->GetYaxis()->SetRangeUser(0, 0.02);
			h_mt2lblb_minitrees_Signif_cut[icut]->GetXaxis()->SetTitle("Mt2lblb Significance");
			h_mt2lblb_minitrees_Signif_cut[icut]->SetTitle("Mt2lblb Significance with cut in Mt2ll");
			h_mt2lblb_minitrees_Signif_cut[icut]->DrawCopy(Option);
			h_mt2lblb_S0_TW05_Signif_cut[icut]->GetXaxis()->SetRangeUser(1, 900);
			h_mt2lblb_S0_TW05_Signif_cut[icut]->DrawCopy("histosame");

			legcut3->Draw();

			Option = "histosame";

		}


	}	

	CC->Print("MT2top_neutrino_pxy.png");
	CCmt2ll->Print("MT2top_S0_TW05.png");
	CCcut->Print("MT2top_mt2llcut.png");

}
