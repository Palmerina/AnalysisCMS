#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <TLegend.h>
#include <iostream>


void Mt2RecoStudies() {
  
  TFile *HistoFile[2];
  HistoFile[1] = TFile::Open("rootfiles/nominal/Stop/T2tt_mStop500-525-550_mLSP1to425-325to450-1to475.root"); 
  HistoFile[0] = TFile::Open("rootfiles/nominal/Stop/TTTo2L2Nu.root"); 

  TH1D* MT2Histo[2];
  TH1D* MT2Integ[2];
  TH1D* MT2Signif;

  int HistoCol[2] = {4, 2};

  TString BTag = "02_Has1BJet"; // "02_Has2BJet"
  TString Mt2T = "mt2lblb"; //    "mt2bb" "mt2lblb", "mt2ll"
  TString Reco = "true"; //      "",  "true"

  float HistoInt[2];

  for (int hf = 0; hf<2; hf++) { // stop or top

    TString HistoName = "Stop/" + BTag + "/h_" + Mt2T + Reco + "_ll";
    MT2Histo[hf] = (TH1D*) HistoFile[hf]->Get(HistoName);
    
    MT2Histo[hf]->Rebin(10);
    HistoInt[hf] = MT2Histo[hf]->Integral();
    MT2Histo[hf]->Scale(1./HistoInt[hf]); // normalization of the histogram
    MT2Histo[hf]->SetLineColor(HistoCol[hf]);
    MT2Histo[hf]->SetLineStyle(1);
    MT2Histo[hf]->SetLineWidth(2);
  
    int nBinsX = MT2Histo[hf]->GetNbinsX();
    HistoName = "h_" + Mt2T + Reco + "_ll_integ";
    MT2Integ[hf] = new TH1D(HistoName, "", nBinsX, 0, 3000);
    float ThisBinContent = 1.;
    for (int ib = 1; ib<=nBinsX; ib++) { // loop through the bins
    // Asigna el valor 1 al primer bin y va restando el contenido del bin anterior a los siguientes
      MT2Integ[hf]->SetBinContent(ib, ThisBinContent); // asigna el valor ThisBinContent al bin ib
      ThisBinContent -= MT2Histo[hf]->GetBinContent(ib); // le resta el valor del bin ib del histograma MT2Histo a ThisBinContent. 
      
    }
    MT2Integ[hf]->SetLineColor(HistoCol[hf]);
    MT2Integ[hf]->SetLineStyle(1);
    MT2Integ[hf]->SetLineWidth(2);



    
  }
  
  int nBinsStop = MT2Integ[1]->GetNbinsX();
  TString HistoName = "h_" + Mt2T + Reco + "_ll_signif";
  MT2Signif = new TH1D(HistoName, "", nBinsStop, 0, 3000);

  for (int ib = 1; ib<=nBinsStop; ib++) {
     	float	topBackground = HistoInt[0]*MT2Integ[0]->GetBinContent(ib); 
     	float	stopEvents = HistoInt[1]*MT2Integ[1]->GetBinContent(ib);
	if (topBackground + stopEvents <= 0.) continue;
	//cout << topBackground + stopEvents << " " << HistoInt[0] << " " << HistoInt[1] << " " << MT2Integ[0]->GetBinContent(ib) << " " << MT2Integ[1]->GetBinContent(ib) << endl;
	float	significance = stopEvents/(std::sqrt(stopEvents+topBackground));
     		MT2Signif->SetBinContent(ib,significance); 
  }	 
  
  MT2Signif->SetLineColor(1);
  MT2Signif->SetLineStyle(1);
  MT2Signif->SetLineWidth(2);

  TCanvas *CC = new TCanvas("CC", "", 1200, 400);
  CC->Divide(2, 2);
  TPad *CC1 = (TPad*)CC->GetPad(1); //Ahi va el MT2Histo
  CC1->SetLogy(); CC1->SetGridx(); CC1->SetGridy(); 
  TPad *CC2 = (TPad*)CC->GetPad(2); //Ahi va el MT2Integ
  CC2->SetLogy(); CC2->SetGridx(); CC2->SetGridy(); 

    
  TString Option = "histo";

  for (int hf = 0; hf<2; hf++) {

    CC->cd(1); // se pone en el TPad 1 (MT2Histo)
    MT2Histo[hf]->GetXaxis()->SetRange(1, 80);
    MT2Histo[hf]->GetXaxis()->SetTitle("M_T2");
    MT2Histo[hf]->DrawCopy(Option);
  
    CC->cd(2); // se pone en el TPad 2 (MT2InteMT2Integ)
    MT2Integ[hf]->GetXaxis()->SetRange(1, 80);
    MT2Integ[hf]->GetXaxis()->SetTitle("Integrated M_T2");
    MT2Integ[hf]->DrawCopy(Option);
   
 
    Option = "histosame";

  }   
 
    TLegend *leg = new TLegend(0.5,0.75,0.7,0.9);
    leg->AddEntry(MT2Histo[0],"top","l");
    leg->AddEntry(MT2Histo[1],"stop","l");

    leg->Draw();

    //TLegend *leg2 = new TLegend(0.5,0.75,3.5,0.9);
    //leg2->AddEntry(MT2Integ[0],"top","l");
    //leg2->AddEntry(MT2Integ[1],"stop","l");
    //leg2->Draw();

    TCanvas *c= new TCanvas("c","",1200,400);
    MT2Signif->GetXaxis()->SetTitle("Significance");
    MT2Signif->Draw();

  CC->Print(BTag + "_" + Mt2T + "_" + Reco + ".png");
  
}
