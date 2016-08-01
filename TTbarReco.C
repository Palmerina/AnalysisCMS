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

float eventW, metPfType1;
float njet, channel, nbjet30csvv2l, nbjet30csvv2m, nbjet30csvv2t;
float mt2ll, mt2bb, mt2bbtrue, mt2lblb, mt2lblbcomb, mt2lblbtrue;
float mlb1, mlb1true, mlb1comb, mlb1truecomb;
float mlb2, mlb2true, mlb2comb, mlb2truecomb;
float tjet1pt, tjet1phi, tjet1eta, tjet1mass, tjet1csvv2ivf, tjet1assignment;
float tjet2pt, tjet2phi, tjet2eta, tjet2mass, tjet2csvv2ivf, tjet2assignment;
float bjet1pt, bjet1phi, bjet1eta, bjet1mass, bjet1csvv2ivf;
float bjet2pt, bjet2phi, bjet2eta, bjet2mass, bjet2csvv2ivf;


TH1D *h_mt2lblbtrue_top[2]; 
TH1D *h_mt2lblbtrue_stop[2]; 


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

  return MiniTree;

}

void TTbarReco() {

  TString FileName[2] = {"./minitrees/nominal/Stop/TTTo2L2Nu.root",
			 "./minitrees/nominal/Stop/T2tt_mStop500-525-550_mLSP1to425-325to450-1to475.root"};


    TString Option = "histo";
    TString Histoname[2] = {"h_mt2lblb_top", "h_mt2lblb_stop"};

    int HistoCol[2] = {4, 2};
    float mt2lblb_top_Int[2];
    float mt2lblb_stop_Int[2];


  for (int dt = 0; dt<2; dt++) {

    TFile *MiniTreeFile = TFile::Open(FileName[dt]);
    
    TTree *MiniTree = GetMiniTree(MiniTreeFile);
    
    Int_t nentries = (Int_t) MiniTree->GetEntries();
    

    //Histogramas mt2lblb con corte y sin corte
    h_mt2lblbtrue_top[dt] = new TH1D(Histoname[0],"h_mt2lblbtrue_top",   3000, 0, 3000); //2 histos para top y stop
    h_mt2lblbtrue_stop[dt] = new TH1D(Histoname[1], "h_mt2lblbtrue_stop", 3000, 0, 3000); //2 histos para top y stop



    for (Int_t i = 0; i<nentries; i++) {
      
      MiniTree->GetEntry(i);
      
      // Apply ttbar selection
      if (njet<2) continue;
      if (nbjet30csvv2m<1) continue;


      //Histogramas mt2lblbtrue para el top y el stop
       if (dt == 0) {
		h_mt2lblbtrue_top[0] -> Fill(mt2lblbtrue, eventW);	
      
      		if (mlb1true <= 160 && mlb2true <= 160) {
			h_mt2lblbtrue_top[1] -> Fill(mt2lblbtrue, eventW);
      		}		
       	
      		else 
			continue;
       }
       else {
		h_mt2lblbtrue_stop[0] -> Fill(mt2lblbtrue, eventW);	
      
      		if (mlb1true <= 160 && mlb2true <= 160) {
			h_mt2lblbtrue_stop[1] -> Fill(mt2lblbtrue, eventW);
      		}	
      		else 
			continue;
       }   
     }
   }

    TCanvas *CC = new TCanvas("CC", "", 1200, 700);
    CC->Divide(1, 2);
    TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
    CC1->SetGridx(); CC1->SetGridy(); 
    TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut

   for (int dt = 0; dt<2; dt++) {

    mt2lblb_top_Int[dt] = h_mt2lblbtrue_top[dt]->Integral();
    h_mt2lblbtrue_top[dt]->Scale(1./mt2lblb_top_Int[dt]); // normalization of the histogram

    mt2lblb_stop_Int[dt] = h_mt2lblbtrue_stop[dt]->Integral();
    h_mt2lblbtrue_stop[dt]->Scale(1./mt2lblb_stop_Int[dt]); // normalization of the histogram

    h_mt2lblbtrue_top[dt] -> SetLineColor(HistoCol[dt]);
    h_mt2lblbtrue_stop[dt] -> SetLineColor(HistoCol[dt]);

    h_mt2lblbtrue_top[dt] -> SetLineStyle(1);
    h_mt2lblbtrue_stop[dt] -> SetLineStyle(1);

    h_mt2lblbtrue_top[dt] -> SetLineWidth(2);
    h_mt2lblbtrue_stop[dt] -> SetLineWidth(2);
    CC->cd(1); // se pone en el TPad 1 
    h_mt2lblbtrue_top[dt]->GetXaxis()->SetRange(1, 300);
    h_mt2lblbtrue_top[dt]->GetXaxis()->SetTitle("MT2lblb top");
    h_mt2lblbtrue_stop[dt]->DrawCopy(Option);
    

    CC->cd(2); // se pone en el TPad 2 
    h_mt2lblbtrue_stop[dt]->GetXaxis()->SetRange(1, 300);
    h_mt2lblbtrue_stop[dt]->GetXaxis()->SetTitle("MT2lblb stop");
    h_mt2lblbtrue_stop[dt]->DrawCopy(Option);
    

    Option= "histosame";
     
 }
 
    TLegend *leg1 = new TLegend(0.5,0.75,0.7,0.9);
    leg1->AddEntry(h_mt2lblbtrue_top[0],"top","l");
    leg1->AddEntry(h_mt2lblbtrue_top[1],"top (mlbtrue < 160)","l");
    leg1->Draw();

    TLegend *leg2 = new TLegend(0.5,0.75,0.7,0.9);
    leg1->AddEntry(h_mt2lblbtrue_stop[0],"stop","l");
    leg1->AddEntry(h_mt2lblbtrue_stop[1],"stop (mlbtrue < 160)","l");
    leg1->Draw();

  CC->Print("mt2lblb.png");
}

