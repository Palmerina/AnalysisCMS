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

double ComputeMt() {


}

Mt2Result ComputeMt2Top(Mt2::LorentzVector& Lepton1,  Mt2::LorentzVector& Lepton2,
			Mt2::LorentzVector& Bottom1,  Mt2::LorentzVector& Bottom2,
			Mt2::TwoVector& ptmiss, double mEachInvisible, 
			int Mt2TopStrategy, float Mt2TopWeight) {

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

  //TString FileName[2] = {"/afs/cern.ch/user/p/palmerin/public/TTTo2L2Nu.root",
  //"/afs/cern.ch/user/p/palmerin/public/T2tt_mStop500-525-550_mLSP1to425-325to450-1to475.root"};

  TString FileName[2] = {"minitrees/minitreesLuca/TTTo2L2Nu_test.root",
			 "minitrees/minitreesLuca/T2tt_mStop500-525-550_mLSP1to425-325to450-1to475.root"};

  for (int dt = 0; dt<1; dt++) {

    TFile *MiniTreeFile = TFile::Open(FileName[dt]);
    
    TTree *MiniTree = GetMiniTree(MiniTreeFile);
    
    Int_t nentries = (Int_t) MiniTree->GetEntries();
    
    for (Int_t i = 0; i<nentries; i++) {
      
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

      if (TestStandardMt2) {

	//Mt2::LorentzTransverseVector visA = Lepton1.getLorentzTransverseVector();
	//Mt2::LorentzTransverseVector visB = Lepton2.getLorentzTransverseVector();
	Mt2::LorentzTransverseVector visA(Mt2::TwoVector(lep1pt*cos(lep1phi), lep1pt*sin(lep1phi)), 0.1);
	Mt2::LorentzTransverseVector visB(Mt2::TwoVector(lep2pt*cos(lep2phi), lep2pt*sin(lep2phi)), 0.1);
	
	// Let´s cross check mt2ll first  
	const double mT2ll = Mt2Calcolator.mt2_332(visA, visB, pT_Miss, m_invis_mass);	
	Mt2Result Mt2Strategy0p0 = ComputeMt2Top(Lepton1, Lepton2, Bottom1, Bottom2, pT_Miss, m_invis_mass, 0, 0.);


	std::cout << "Test mt2ll " << mt2ll << " " << mT2ll << " " << Mt2Strategy0p0.Mt2ll << std::endl;
	
	// Then cross check mt2lblb
	Mt2::LorentzTransverseVector visUA(Mt2::TwoVector(LepB[0][0].Px(), LepB[0][0].Py()), LepB[0][0].M());
	Mt2::LorentzTransverseVector visUB(Mt2::TwoVector(LepB[1][0].Px(), LepB[1][0].Py()), LepB[1][0].M());
	
	const double mT2lblb = Mt2Calcolator.mt2_332(visUA, visUB, pT_Miss, m_invis_mass);
	
	Mt2Result Mt2Strategy0p999 = ComputeMt2Top(Lepton1, Lepton2, Bottom1, Bottom2, pT_Miss, m_invis_mass, 0, 999.);

	std::cout << "Test mt2lblb " << mt2lblb << " " << mt2lblbcomb << " " << mT2lblb << " " << Mt2Strategy0p999.Mt2lblb << std::endl;
	std::cout << "Test mlb     " << mlb1 << " " << mlb2 << " " << Mt2Strategy0p999.Mlb1 << " " << Mt2Strategy0p999.Mlb2 << std::endl;
     
      }

      // Here we can study different mt2ll-mt2lblb minimisations
      
    }

  }

}
