#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <TLegend.h>
#include <iostream>


void mlbStudies() {
  
  TFile *HistoFile;
  HistoFile = TFile::Open("rootfiles/nominal/Stop/TTTo2L2Nu.root"); 

  TH1D* mlbHisto[2];
  TH1D* mlbcombHisto[2];
  TH1D* mlbtrueHisto[2];
  TH1D* mlbtruecombHisto[2];

  int HistoCol[2] = {4, 2};

  TString BTag = "02_Has2BJet"; // "02_Has2BJet","02_Has1BJet"
  TString mlb[2]= {"mlb1","mlb2"}; 
  TString mlbcomb[2]= {"mlb1comb","mlb2comb"}; 
  TString mlbtrue[2]= {"mlb1true","mlb2true"}; 
  TString mlbtruecomb[2]= {"mlb1truecomb","mlb2truecomb"}; 
  TString Reco = ""; //      "",  "true"

  float HistoInt[2];

  for (int hf=0; hf<2; hf++) { 

    TString HistoName_mlb = "Stop/" + BTag + "/h_" + mlb[hf]  + Reco + "_ll";
    TString HistoName_mlbcomb = "Stop/" + BTag + "/h_" + mlbcomb[hf]  + Reco + "_ll";
    TString HistoName_mlbtrue = "Stop/" + BTag + "/h_" + mlbtrue[hf]  + Reco + "_ll";
    TString HistoName_mlbtruecomb = "Stop/" + BTag + "/h_" + mlbtruecomb[hf]  + Reco + "_ll";

    mlbHisto[hf] = (TH1D*) HistoFile->Get(HistoName_mlb);
    mlbcombHisto[hf] = (TH1D*) HistoFile->Get(HistoName_mlbcomb);
    mlbtrueHisto[hf] = (TH1D*) HistoFile->Get(HistoName_mlbtrue);
    mlbtruecombHisto[hf] = (TH1D*) HistoFile->Get(HistoName_mlbtruecomb);
    
    //MT2Histo[hf]->Rebin(10);
    //HistoInt[hf] = MT2Histo[hf]->Integral();
    //MT2Histo[hf]->Scale(1./HistoInt[hf]); // normalization of the histogram
    mlbHisto[hf]->SetLineColor(HistoCol[hf]);
    mlbcombHisto[hf]->SetLineColor(HistoCol[hf]);
    mlbtrueHisto[hf]->SetLineColor(HistoCol[hf]);
    mlbtruecombHisto[hf]->SetLineColor(HistoCol[hf]);

    mlbHisto[hf]->SetLineStyle(1);
    mlbcombHisto[hf]->SetLineStyle(1);
    mlbtrueHisto[hf]->SetLineStyle(1);
    mlbtruecombHisto[hf]->SetLineStyle(1);

    mlbHisto[hf]->SetLineWidth(2);
    mlbcombHisto[hf]->SetLineWidth(2);
    mlbtrueHisto[hf]->SetLineWidth(2);
    mlbtruecombHisto[hf]->SetLineWidth(2);
  
      
    }

  
  TCanvas *CC = new TCanvas("CC", "", 1200, 400);
  CC->Divide(2, 2);
  TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el mlbHisto
  CC1->SetGridx(); CC1->SetGridy(); 
  TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el mlbcombHisto
  CC2->SetGridx(); CC2->SetGridy(); 
  TPad *CC3 = (TPad*)CC->GetPad(3); //Ahi va el mlbtrueHisto
  CC3->SetGridx(); CC3->SetGridy(); 
  TPad *CC4 = (TPad*)CC->GetPad(4); //Ahi va el mlbtruecombHisto
  CC4->SetGridx(); CC4->SetGridy(); 

    
  TString Option = "histo";

  for (int hf = 0; hf<2; hf++) {

    CC->cd(1); // se pone en el TPad 1 (mlbHisto)
    mlbHisto[hf]->GetXaxis()->SetRange(1, 300);
    mlbHisto[hf]->SetMaximum(30);

    mlbHisto[hf]->GetXaxis()->SetTitle("mlb");
    mlbHisto[hf]->DrawCopy(Option);
    
    TLegend *leg1 = new TLegend(0.5,0.75,0.7,0.9);
    leg1->AddEntry(mlbHisto[0],"mlb1","l");
    leg1->AddEntry(mlbHisto[1],"mlb2","l");
    leg1->Draw();


    CC->cd(2); // se pone en el TPad 2 (mlbcombHisto)
    mlbcombHisto[hf]->GetXaxis()->SetRange(1, 300);
    mlbcombHisto[hf]->SetMaximum(30);
    mlbcombHisto[hf]->GetXaxis()->SetTitle("mlbcomb");
    mlbcombHisto[hf]->DrawCopy(Option);
    
    TLegend *leg2 = new TLegend(0.5,0.75,0.7,0.9);
    leg2->AddEntry(mlbHisto[0],"mlb1comb","l");
    leg2->AddEntry(mlbHisto[1],"mlb2comb","l");
    leg2->Draw();
 

    CC->cd(3); // se pone en el TPad 1 (mlbtrueHisto)
    mlbtrueHisto[hf]->GetXaxis()->SetRange(1, 300);
    mlbtrueHisto[hf]->SetMaximum(30);
    mlbtrueHisto[hf]->GetXaxis()->SetTitle("mlbtrue");
    mlbtrueHisto[hf]->DrawCopy(Option);
    
    TLegend *leg3 = new TLegend(0.5,0.75,0.7,0.9);
    leg3->AddEntry(mlbHisto[0],"mlb1true","l");
    leg3->AddEntry(mlbHisto[1],"mlb2true","l");
    leg3->Draw();


    CC->cd(4); // se pone en el TPad 4 (mlbtruecombHisto)
    mlbtruecombHisto[hf]->GetXaxis()->SetRange(1, 300);
    mlbtruecombHisto[hf]->SetMaximum(30);
    mlbtruecombHisto[hf]->GetXaxis()->SetTitle("mlbtruecomb");
    mlbtruecombHisto[hf]->DrawCopy(Option);
    
    TLegend *leg4 = new TLegend(0.5,0.75,0.7,0.9);
    leg4->AddEntry(mlbHisto[0],"mlb1truecomb","l");
    leg4->AddEntry(mlbHisto[1],"mlb2truecomb","l");
    leg4->Draw();


    Option = "histosame";


  }   

  //CC->Print(BTag + "_" + mlb + "_" + Reco + ".png");
  
}
