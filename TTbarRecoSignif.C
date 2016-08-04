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


TH1D *h_mt2lblbtrue[2]; 
TH1D *h_mt2lblbtrue_cut[2]; 
TH1D* h_mt2lblbInteg[2];
TH1D* h_mt2lblbInteg_cut[2];
TH1D* h_mt2lblbSignif;
TH1D* h_mt2lblbSignif_cut;


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

void TTbarRecoSignif() {

  TString FileName[2] = {"./minitrees/nominal/Stop/TTTo2L2Nu.root",
			 "./minitrees/nominal/Stop/T2tt_mStop500-525-550_mLSP1to425-325to450-1to475.root"};

    TCanvas *CC = new TCanvas("CC", "", 1500, 750);
    CC->Divide(2, 3);
    TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el h_mt2mlblbtrue
    CC1->SetLogy();CC1->SetGridx(); CC1->SetGridy(); 
    TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el h_mt2mlblbtrue_cut
    CC2->SetLogy();CC2->SetGridx(); CC2->SetGridy(); 
    TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el h_mt2mlblbtrue_cut
    CC3->SetLogy();CC3->SetGridx(); CC3->SetGridy(); 
    TPad *CC4 = (TPad*)CC->GetPad(4);//Ahi va el h_mt2mlblbtrue_cut
    CC4->SetLogy();CC4->SetGridx(); CC4->SetGridy(); 
    TPad *CC5 = (TPad*)CC->GetPad(5); //Ahi va el h_mt2mlblbtrue_cut
    CC5->SetGridx(); CC5->SetGridy(); 
    TPad *CC6 = (TPad*)CC->GetPad(6); //Ahi va el h_mt2mlblbtrue_cut
    CC6->SetGridx(); CC6->SetGridy(); 

    TString Option = "histo";
    TString Histoname[2] = {"h_mt2lblb", "h_mt2lblb_cut"};

    int HistoCol[2] = {4, 2};
    float mt2lblbtrue_Int[2];
    float mt2lblbtrue_cut_Int[2];


  for (int dt = 0; dt<2; dt++) {

    TFile *MiniTreeFile = TFile::Open(FileName[dt]);
    
    TTree *MiniTree = GetMiniTree(MiniTreeFile);
    
    Int_t nentries = (Int_t) MiniTree->GetEntries();
    
    h_mt2lblbtrue[dt] = new TH1D(Histoname[dt],"h_mt2lblbtrue",   3000, 0, 3000); //2 histos para top y stop
    h_mt2lblbtrue_cut[dt] = new TH1D(Histoname[dt], "h_mt2lblbtrue_cut", 3000, 0, 3000); //2 histos para top y stop
    h_mt2lblbInteg[dt] = new TH1D("h_mt2lblbInteg", "h_mt2lblbInteg", 3000, 0, 3000); //2 histos para top y stop
    h_mt2lblbInteg_cut[dt] = new TH1D("h_mt2lblbInteg_cut", "h_mt2lblbInteg_cut", 3000, 0, 3000); //2 histos para top y stop
    h_mt2lblbSignif = new TH1D("h_mt2lblbSignif", "h_mt2lblbSignif", 3000, 0, 3000); //2 histos para top y stop
    h_mt2lblbSignif_cut = new TH1D("h_mt2lblbSignif_cut", "h_mt2lblbSignif_cut", 3000, 0, 3000); //2 histos para top y stop



    for (Int_t i = 0; i<nentries; i++) {
      
      MiniTree->GetEntry(i);
             h_mt2lblbtrue[dt] -> Fill(mt2lblbtrue, eventW);	
      
      if (mlb1true <= 160 && mlb2true <= 160) {
	h_mt2lblbtrue_cut[dt] -> Fill(mt2lblbtrue, eventW);
      }
      else 
	continue;
    }

    mt2lblbtrue_Int[dt] = h_mt2lblbtrue[dt]->Integral();
    h_mt2lblbtrue[dt]->Scale(1./mt2lblbtrue_Int[dt]); // normalization of the histogram

    mt2lblbtrue_cut_Int[dt] = h_mt2lblbtrue_cut[dt]->Integral();
    h_mt2lblbtrue_cut[dt]->Scale(1./mt2lblbtrue_cut_Int[dt]); // normalization of the histogram

    h_mt2lblbtrue[dt] -> SetLineColor(HistoCol[dt]);
    h_mt2lblbtrue_cut[dt] -> SetLineColor(HistoCol[dt]);
    h_mt2lblbInteg[dt] -> SetLineColor(HistoCol[dt]);
    h_mt2lblbInteg_cut[dt] -> SetLineColor(HistoCol[dt]);
    h_mt2lblbSignif -> SetLineColor(1);
    h_mt2lblbSignif_cut -> SetLineColor(1);

    h_mt2lblbtrue[dt] -> SetLineStyle(1);
    h_mt2lblbtrue_cut[dt] -> SetLineStyle(1);
    h_mt2lblbInteg[dt] -> SetLineStyle(1);
    h_mt2lblbInteg_cut[dt] ->SetLineStyle(1);
    h_mt2lblbSignif -> SetLineStyle(1);
    h_mt2lblbSignif_cut -> SetLineStyle(1);

    h_mt2lblbtrue[dt] -> SetLineWidth(2);
    h_mt2lblbtrue_cut[dt] -> SetLineWidth(2);
    h_mt2lblbInteg[dt] -> SetLineWidth(2);
    h_mt2lblbInteg_cut[dt] ->SetLineWidth(2);
    h_mt2lblbSignif -> SetLineWidth(2);
    h_mt2lblbSignif_cut -> SetLineWidth(2);

}      

    for (int hf = 0; hf<2; hf++) { // stop or top

  
 	int nBinsX = h_mt2lblbtrue[hf]->GetNbinsX();
 	int nBinsX_cut = h_mt2lblbtrue_cut[hf]->GetNbinsX();
    	float ThisBinContent = 1.;
    	float ThisBinContent_cut = 1.;

	for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins
    	// Asigna el valor 1 al primer bin y va restando el contenido del bin anterior a los siguientes
      		h_mt2lblbInteg[hf]->SetBinContent(ib, ThisBinContent); // asigna el valor ThisBinContent al bin ib
      		ThisBinContent -= h_mt2lblbtrue[hf]->GetBinContent(ib); // le resta el valor del bin ib del histograma MT2Histo a ThisBinContent. 
    
  
   	

      		h_mt2lblbInteg_cut[hf]->SetBinContent(ib, ThisBinContent_cut); // asigna el valor ThisBinContent al bin ib
      		ThisBinContent_cut -= h_mt2lblbtrue_cut[hf]->GetBinContent(ib); // le resta el valor del bin ib del histograma MT2Histo a ThisBinContent. 
	}


    }


 	int nBinsX = h_mt2lblbtrue[0]->GetNbinsX();

	for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins
  
     		float	topBackground = mt2lblbtrue_Int[0] * h_mt2lblbInteg[0]->GetBinContent(ib); 
     		float	stopEvents = mt2lblbtrue_Int[1] * h_mt2lblbInteg[1]->GetBinContent(ib);
		if (topBackground + stopEvents <= 0.) continue;
		float	significance = stopEvents/(std::sqrt(stopEvents+topBackground));
     		h_mt2lblbSignif->SetBinContent(ib,significance); 
  
   	
     		float	topBackground_cut = mt2lblbtrue_cut_Int[0] * h_mt2lblbInteg_cut[0]->GetBinContent(ib); 
     		float	stopEvents_cut =  mt2lblbtrue_cut_Int[1] * h_mt2lblbInteg_cut[1]->GetBinContent(ib);
		if (topBackground_cut + stopEvents_cut <= 0.) continue;
		float	significance_cut = stopEvents_cut/(std::sqrt(stopEvents_cut+topBackground_cut));
     		h_mt2lblbSignif_cut->SetBinContent(ib,significance_cut); 
      
   	}

  

    
    for (int dt = 0; dt<2; dt++) { // stop or top

    CC->cd(1); // se pone en el TPad 1 
    h_mt2lblbtrue[dt]->GetXaxis()->SetRange(1, 300);
    h_mt2lblbtrue[dt]->GetXaxis()->SetTitle("MT2lblb");
    h_mt2lblbtrue[dt]->DrawCopy(Option);
    
    TLegend *leg1 = new TLegend(0.5,0.75,0.7,0.9);
    leg1->AddEntry(h_mt2lblbtrue[0],"top","l");
    leg1->AddEntry(h_mt2lblbtrue[1],"stop","l");
    leg1->Draw();

    CC->cd(2); // se pone en el TPad 2 
    h_mt2lblbtrue_cut[dt]->GetXaxis()->SetRange(1, 300);
    h_mt2lblbtrue_cut[dt]->GetXaxis()->SetTitle("MT2lblb with mlb<160");
    h_mt2lblbtrue_cut[dt]->DrawCopy(Option);
    
    CC->cd(3); // se pone en el TPad 1 
    h_mt2lblbInteg[dt]->GetXaxis()->SetRange(1, 300);
    h_mt2lblbInteg[dt]->GetXaxis()->SetTitle("MT2lblb integral");
    h_mt2lblbInteg[dt]->DrawCopy(Option);
    

    CC->cd(4); // se pone en el TPad 1 
    h_mt2lblbInteg_cut[dt]->GetXaxis()->SetRange(1, 300);
    h_mt2lblbInteg_cut[dt]->GetXaxis()->SetTitle("MT2lblb integral with cut");
    h_mt2lblbInteg_cut[dt]->DrawCopy(Option);
    
    CC->cd(5); // se pone en el TPad 1 
    h_mt2lblbSignif->GetXaxis()->SetRange(1, 300);
    h_mt2lblbSignif->GetXaxis()->SetTitle("Significance");
    h_mt2lblbSignif->DrawCopy(Option);
    
    CC->cd(6); // se pone en el TPad 1 
    h_mt2lblbSignif_cut->GetXaxis()->SetRange(1, 300);
    h_mt2lblbSignif_cut->GetXaxis()->SetTitle("Significance with cut");
    h_mt2lblbSignif_cut->DrawCopy(Option);
    

    Option= "histosame";
  }   

 

  CC->Print("TTbarRecoSignif.png");
}
