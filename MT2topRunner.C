void MT2topRunner() {
  
  // Get Mt2 library, and its pre-requisite minuit2 library:
  gSystem->Load("libMinuit2.so");
  gSystem->Load("/afs/cern.ch/user/p/palmerin/test/CMSSW_7_6_3/lib/liboxbridgekinetics-1.0.so");
  
  // tell the interpreter where the include files are for the MT2 library:
  gInterpreter->AddIncludePath("/afs/cern.ch/user/p/palmerin/test/CMSSW_7_6_3/include/oxbridgekinetics-1.0");
  
  gROOT->ProcessLine(".L MT2top.C++");

  gROOT->ProcessLine(".x MT2topExec.C");

}
