#include <iostream>

#include "TFile.h"
#include "TMVA/Config.h"
#include "TMVA/Factory.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"


// Constants
//------------------------------------------------------------------------------
const TString inputdir       = "../minitrees/nominal/TTDM/";
const TString trainingdir    = "output/training/";
const TString weightsdir     = "output/weights/";
const TString applicationdir = "output/application/";


// Functions
//------------------------------------------------------------------------------
void MVATrain  (TString signal);

void MVARead   (TString signal,
		TString filename);

void AddProcess(TString kind,
		TString filename);


// Data members
//------------------------------------------------------------------------------
TTree*              _signaltree;
std::vector<TTree*> _mctree;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// MVA
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void MVA(TString signal     = "ttDM0001scalar00010",
	 bool    doMVATrain = true,
	 bool    doMVARead  = true)
{
  if (!doMVATrain && !doMVARead) return;

  gInterpreter->ExecuteMacro("../test/PaperStyle.C");

  gSystem->mkdir(trainingdir,    kTRUE);
  gSystem->mkdir(applicationdir, kTRUE);

  (TMVA::gConfig().GetIONames()).fWeightFileDir = weightsdir;


  // Training
  //----------------------------------------------------------------------------
  if (doMVATrain) MVATrain(signal);


  // Reading
  //----------------------------------------------------------------------------
  if (doMVARead)
    {
      MVARead(signal, signal);

      MVARead(signal, "01_Data");
      MVARead(signal, "14_HZ");
      MVARead(signal, "10_HWW");
      MVARead(signal, "06_WW");
      MVARead(signal, "02_WZTo3LNu");
      MVARead(signal, "03_ZZ");
      MVARead(signal, "11_Wg");
      MVARead(signal, "07_ZJets");
      MVARead(signal, "09_TTV");
      MVARead(signal, "04_TTTo2L2Nu");
      MVARead(signal, "05_ST");
      MVARead(signal, "00_Fakes");
    }
}


//------------------------------------------------------------------------------
// MVATrain
//------------------------------------------------------------------------------
void MVATrain(TString signal)
{
  TFile* outputfile = TFile::Open(trainingdir + signal + ".root", "recreate");


  // Factory
  //----------------------------------------------------------------------------
  TMVA::Factory* factory = new TMVA::Factory(signal, outputfile,    
					     "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");


  // Get the trees
  //----------------------------------------------------------------------------
  _mctree.clear();

  AddProcess("signal", signal);
  AddProcess("background", "TTTo2L2Nu");

  //  AddProcess("background", "14_HZ");
  //  AddProcess("background", "10_HWW");
  //  AddProcess("background", "06_WW");
  //  AddProcess("background", "02_WZTo3LNu");
  //  AddProcess("background", "03_ZZ");
  //  AddProcess("background", "11_Wg");
  //  AddProcess("background", "07_ZJets");
  //  AddProcess("background", "09_TTV");
  //  AddProcess("background", "05_ST");
  //  AddProcess("background", "00_Fakes");

  Double_t weight = 1.0;

  factory->AddSignalTree(_signaltree, weight);

  for (UInt_t i=0; i<_mctree.size(); i++) factory->AddBackgroundTree(_mctree[i], weight);

  factory->SetWeightExpression("eventW");


  // Add variables
  //----------------------------------------------------------------------------
  // Be careful with the order: it must be respected at the reading step
  // factory->AddVariable("<var1>+<var2>", "pretty title", "unit", 'F');

  factory->AddVariable("channel",        "", "", 'F');
  factory->AddVariable("metPfType1",     "", "", 'F');
  factory->AddVariable("m2l",            "", "", 'F');
  //factory->AddVariable("njet",           "", "", 'F');
  //factory->AddVariable("nbjet30csvv2m",  "", "", 'F');
  factory->AddVariable("lep1pt",         "", "", 'F');
  factory->AddVariable("lep2pt",         "", "", 'F');
  factory->AddVariable("jet1pt",         "", "", 'F');
  factory->AddVariable("jet2pt",         "", "", 'F');
  factory->AddVariable("ht",             "", "", 'F');
  factory->AddVariable("mtw1",           "", "", 'F');
  factory->AddVariable("mtw2",           "", "", 'F');
  factory->AddVariable("dphill",         "", "", 'F');
  factory->AddVariable("dphilep1jet1",   "", "", 'F');
  factory->AddVariable("dphilep1jet2",   "", "", 'F');
  factory->AddVariable("dphilmet1",      "", "", 'F');
  factory->AddVariable("dphilep2jet1",   "", "", 'F');
  factory->AddVariable("dphilep2jet2",   "", "", 'F');
  factory->AddVariable("dphilmet2",      "", "", 'F');
  factory->AddVariable("dphijj",         "", "", 'F');
  factory->AddVariable("dphijet1met",    "", "", 'F');
  factory->AddVariable("dphijet2met",    "", "", 'F');
  factory->AddVariable("dphillmet",      "", "", 'F');


  // Preselection cuts and preparation
  //----------------------------------------------------------------------------
  factory->PrepareTrainingAndTestTree("", "NormMode=EqualNumEvents:nTrain_Signal=1000:nTest_Signal=1000:nTrain_Background=1000:nTest_Background=1000:SplitMode=Alternate:!V");


  // Book MVA
  //----------------------------------------------------------------------------
  factory->BookMethod(TMVA::Types::kMLP, "MLP01",
		      "H:!V:NeuronType=sigmoid:NCycles=200:VarTransform=Norm:HiddenLayers=25,10:TestRate=2");
  //factory->BookMethod(TMVA::Types::kMLP, "MLP02",
  //		      "H:!V:NeuronType=sigmoid:NCycles=200:VarTransform=Norm:HiddenLayers=15,10:TestRate=2");
  //factory->BookMethod(TMVA::Types::kMLP, "MLP03",
  // 		      "H:!V:NeuronType=sigmoid:NCycles=200:VarTransform=Norm:HiddenLayers=10,9:TestRate=2");
  //factory->BookMethod(TMVA::Types::kMLP, "MLP04",
  //		      "H:!V:NeuronType=sigmoid:NCycles=200:VarTransform=Norm:HiddenLayers=10,5:TestRate=2");
  factory->BookMethod(TMVA::Types::kMLP, "MLP05",
		      "H:!V:NeuronType=sigmoid:NCycles=400:VarTransform=Norm:HiddenLayers=6,5:TestRate=2");
  factory->BookMethod(TMVA::Types::kMLP, "MLP06",
  		      "H:!V:NeuronType=sigmoid:NCycles=200:VarTransform=Norm:HiddenLayers=20,10:TestRate=2");
  factory->BookMethod(TMVA::Types::kMLP, "MLP07",
  		      "H:!V:NeuronType=sigmoid:NCycles=200:VarTransform=Norm:HiddenLayers=20,20:TestRate=2");
  //factory->BookMethod(TMVA::Types::kMLP, "MLP08",
  //		      "H:!V:NeuronType=sigmoid:NCycles=200:VarTransform=Norm:HiddenLayers=10,10:TestRate=2");


  // Train, test and evaluate MVA
  //----------------------------------------------------------------------------
  factory->TrainAllMethods();     // Train using the set of training events
  factory->TestAllMethods();      // Evaluate using the set of test events
  factory->EvaluateAllMethods();  // Evaluate and compare performance


  // Save the output
  //----------------------------------------------------------------------------
  outputfile->Close();

  delete factory;
}


//------------------------------------------------------------------------------
// MVARead
//------------------------------------------------------------------------------
void MVARead(TString signal, TString filename)
{

  TMVA::Reader* reader   = new TMVA::Reader("!Color:!Silent");    

  float channel;
  float metPfType1;
  float m2l;
  float njet;
  float nbjet20cmvav2l;
  float lep1pt;
  float lep2pt;
  float jet1pt;
  float jet2pt;
  float ht;
  float mtw1;
  float mtw2;
  float dphill;
  float dphilep1jet1;
  float dphilep1jet2;
  float dphilmet1;
  float dphilep2jet1;
  float dphilep2jet2;
  float dphilmet2;
  float dphijj;
  float dphijet1met;
  float dphijet2met;
  float dphillmet;
  float eventW;
  float mva01; float mva05; float mva06; float mva07; 


  reader->AddVariable("channel",        &channel);
  reader->AddVariable("metPfType1",     &metPfType1);
  reader->AddVariable("m2l",            &m2l);
  //reader->AddVariable("njet",           &njet);
  //reader->AddVariable("nbjet20cmvav2l", &nbjet20cmvav2l);
  reader->AddVariable("lep1pt",         &lep1pt);
  reader->AddVariable("lep2pt",         &lep2pt);
  reader->AddVariable("jet1pt",         &jet1pt);
  reader->AddVariable("jet2pt",         &jet2pt);
  reader->AddVariable("ht",             &ht);
  reader->AddVariable("mtw1",           &mtw1);
  reader->AddVariable("mtw2",           &mtw2);
  reader->AddVariable("dphill",         &dphill);
  reader->AddVariable("dphilep1jet1",   &dphilep1jet1);
  reader->AddVariable("dphilep1jet2",   &dphilep1jet2);
  reader->AddVariable("dphilmet1",      &dphilmet1);
  reader->AddVariable("dphilep2jet1",   &dphilep2jet1);
  reader->AddVariable("dphilep2jet2",   &dphilep2jet2);
  reader->AddVariable("dphilmet2",      &dphilmet2);
  reader->AddVariable("dphijj",         &dphijj);
  reader->AddVariable("dphijet1met",    &dphijet1met);
  reader->AddVariable("dphijet2met",    &dphijet2met);
  reader->AddVariable("dphillmet",      &dphillmet);

 
  // Book MVA methods
  //----------------------------------------------------------------------------
  reader->BookMVA("01", weightsdir + signal + "_MLP01.weights.xml");
  reader->BookMVA("05", weightsdir + signal + "_MLP05.weights.xml");
  reader->BookMVA("06", weightsdir + signal + "_MLP06.weights.xml");
  reader->BookMVA("07", weightsdir + signal + "_MLP07.weights.xml");

  // Get MVA response
  //----------------------------------------------------------------------------
  TH1D* h_mva = new TH1D("h_mva_" + signal, "", 100, -0.05, 1.05);

  TFile* input = TFile::Open(inputdir + filename + ".root", "update");

  TTree* theTree = (TTree*)input->Get("latino");

  TBranch* b_mva01 = theTree->Branch("mva01_" + signal, &mva01, "mva/F" );
  TBranch* b_mva05 = theTree->Branch("mva05_" + signal, &mva05, "mva/F" );
  TBranch* b_mva06 = theTree->Branch("mva06_" + signal, &mva06, "mva/F" );
  TBranch* b_mva07 = theTree->Branch("mva07_" + signal, &mva07, "mva/F" );

  theTree->SetBranchAddress("channel",        &channel);
  theTree->SetBranchAddress("metPfType1",     &metPfType1);
  theTree->SetBranchAddress("m2l",            &m2l);
  //theTree->SetBranchAddress("njet",           &njet);
  //theTree->SetBranchAddress("nbjet20cmvav2l", &nbjet20cmvav2l);
  theTree->SetBranchAddress("lep1pt",         &lep1pt);
  theTree->SetBranchAddress("lep2pt",         &lep2pt);
  theTree->SetBranchAddress("jet1pt",         &jet1pt);
  theTree->SetBranchAddress("jet2pt",         &jet2pt);
  theTree->SetBranchAddress("ht",             &ht);
  theTree->SetBranchAddress("mtw1",           &mtw1);
  theTree->SetBranchAddress("mtw2",           &mtw2);
  theTree->SetBranchAddress("dphill",         &dphill);
  theTree->SetBranchAddress("dphilep1jet1",   &dphilep1jet1);
  theTree->SetBranchAddress("dphilep1jet2",   &dphilep1jet2);
  theTree->SetBranchAddress("dphilmet1",      &dphilmet1);
  theTree->SetBranchAddress("dphilep2jet1",   &dphilep2jet1);
  theTree->SetBranchAddress("dphilep2jet2",   &dphilep2jet2);
  theTree->SetBranchAddress("dphilmet2",      &dphilmet2);
  theTree->SetBranchAddress("dphijj",         &dphijj);
  theTree->SetBranchAddress("dphijet1met",    &dphijet1met);
  theTree->SetBranchAddress("dphijet2met",    &dphijet2met);
  theTree->SetBranchAddress("dphillmet",      &dphillmet);
  theTree->SetBranchAddress("eventW",         &eventW);

  Long64_t nentries = theTree->GetEntries();

  for (Long64_t ievt=0; ievt<nentries; ievt++) {

    theTree->GetEntry(ievt);

    mva01 = reader->EvaluateMVA("01");
    mva05 = reader->EvaluateMVA("05");  
    mva06 = reader->EvaluateMVA("06");
    mva07 = reader->EvaluateMVA("07");

    b_mva01->Fill();
    b_mva05->Fill();
    b_mva06->Fill();
    b_mva07->Fill();

    //if (njet > 1) h_mva->Fill(mva, eventW);
  }


  // Save
  //----------------------------------------------------------------------------
  theTree->Write("", TObject::kOverwrite);

  input->Close();

  TFile* target = TFile::Open(applicationdir + signal + "__" + filename + ".root", "recreate");

  h_mva->Write();
  
  target->Close();

  h_mva->Delete();

  delete reader;
}


//------------------------------------------------------------------------------
// AddProcess
//------------------------------------------------------------------------------
void AddProcess(TString kind, TString filename)
{
  TString fullname = inputdir + filename + ".root";

  if (gSystem->AccessPathName(fullname))
    {
      printf(" [MVA::AddProcess] Cannot access %s\n", fullname.Data());
      return;
    }

  TFile* file = new TFile(fullname, "read");

  TTree* tree = (TTree*)file->Get("latino");

  if (kind.EqualTo("signal"))
    _signaltree = tree;
  else
    _mctree.push_back(tree);
}
