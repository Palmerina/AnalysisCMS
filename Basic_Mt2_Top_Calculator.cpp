// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "Mt2/Basic_Mt2_Top_Calculator.h"
#include "Mt2/Mt2Units.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/FunctionMinimum.h"
//#include <fstream>

using ROOT::Minuit2::MnUserParameters;
using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnSimplex;
using ROOT::Minuit2::MnScan;
using ROOT::Minuit2::FunctionMinimum;

namespace Mt2 {

  /*
    
  stealing parts of files in
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/?cvsroot=atlas
  
  like
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/SusyUserExampleHistTool.cxx?rev=1.4;content-type=text%2Fplain;cvsroot=atlas
  
  and
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/~checkout~/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/mT2Fcn.cxx?rev=1.1;content-type=text%2Fplain;cvsroot=atlas
  
  BUT CRUCIALLY:

  then modifying them to remove the assumption about massless visible particles!
  */



  double Basic_Mt2_Top_Calculator::mt2_Top(const LorentzVector& visA, 
					   const LorentzVector& visB,
					   const LorentzVector& visb1, 
					   const LorentzVector& visb2,
					   const TwoVector& ptmiss, 
					   const double mEachInvisible,
					   const int Mt2TopStrategy, const float Mt2TopWeight) {
    const double mSq = mt2_Top_Sq(visA, visB, visb1, visb2, ptmiss, mEachInvisible, Mt2TopStrategy, Mt2TopWeight);

    if (mSq>=0) {
      return sqrt(mSq);
    } else {
      // something went wrong -- maybe just rounding errors;
      return sqrt(fabs(mSq));
    }
  }


  double Basic_Mt2_Top_Calculator::mt2_Top_Sq(const LorentzVector& visA,
					      const LorentzVector& visB,
					      const LorentzVector& visb1, 
					      const LorentzVector& visb2,
					      const TwoVector& ptmiss,
					      const double mEachInvisible,
					      const int Mt2TopStrategy, const float Mt2TopWeight){
    
    LorentzVector LB1 = visA + visb1, LB2 = visB + visb2;

    mT2Fcn theFCN(ptmiss.px(),
		  ptmiss.py(),
		  mEachInvisible,
		  visA.getLorentzTransverseVector().px(),
		  visA.getLorentzTransverseVector().py(),
		  visA.m(),
		  visB.getLorentzTransverseVector().px(),
		  visB.getLorentzTransverseVector().py(),
		  visB.m(),
		  visb1.getLorentzTransverseVector().px(),
		  visb1.getLorentzTransverseVector().py(),
		  visb1.m(),
		  visb2.getLorentzTransverseVector().px(),
		  visb2.getLorentzTransverseVector().py(),
		  visb2.m(),
		  LB1.m(),
		  LB2.m(),
		  Mt2TopStrategy, 
		  Mt2TopWeight);

    const double massScale = (
			      ptmiss.pt() +
			      mEachInvisible +
			      visA.getLorentzTransverseVector().pt()  +
			      visA.m() +
			      visB.getLorentzTransverseVector().pt() +
			      visB.m() +
			      visb1.getLorentzTransverseVector().pt() +
			      visb1.m() +
			      visb2.getLorentzTransverseVector().pt() +
			      visb2.m() 
			      )/10.0; //6.0; ???
    // DANG! Try to get rid of Minuit output:
    //std::ofstream    DANG_log("/dev/null");
    //std::streambuf * DANG_save = std::cerr.rdbuf();
    //if (DANG_log.is_open()) {
    //  std::cerr.rdbuf(DANG_log.rdbuf());
    // }

    double guessx = 0.5*(ptmiss.px());
    double guessy = 0.5*(ptmiss.py());
    
    MnUserParameters upar;
    upar.Add("etx", guessx, 0.02*massScale); 
    upar.Add("ety", guessy, 0.02*massScale);
    const int highQuality=2;    

    // Usually migrad produces the best minumum.
    // But when the minimum is in a fold, migrad can fail badly.
    // On the fold, simplex does well.  We therefore do both separately
    // and record the answer of the one that did best.
    
    // Further to the above notes, it now seems that by choosing the massScale sensibly, and by making the "tolerance" (the second argument to the call that extracts the FunctionMinimum below) much smaller, the simplex algorithm seems to work so well that we don't need the migrad algorithm.  This is good news, as the migrad algorithm produces lots of error output that we can't get rid of.

    MnSimplex simplex(theFCN, upar, highQuality);
    FunctionMinimum minS = simplex(0,massScale*0.000001);

    etxAt=minS.UserState().Value("etx");
    etyAt=minS.UserState().Value("ety");
    etxBt=ptmiss.px()-etxAt;
    etyBt=ptmiss.py()-etyAt;
    //MnMigrad migrad(theFCN, upar, highQuality);
    //FunctionMinimum minM = migrad(0,massScale*0.000001);
    //const double best = fmin(minS.Fval(), minM.Fval());
    const double best = minS.Fval();

    // DANG! Undoing our attempt to get rid of Minuit output:
    //if (DANG_log.is_open()) {
    //  std::cerr.rdbuf(DANG_save);
    //}

    // Set output for Mt2
    double ETchiA = sqrt(etxAt*etxAt + etyAt*etyAt + mEachInvisible*mEachInvisible);
    double AuxMt2_A = visA.m()*visA.m() + mEachInvisible*mEachInvisible + 
      2.0*(visA.ET()*ETchiA - (visA.getLorentzTransverseVector().px()*etxAt + visA.getLorentzTransverseVector().py()*etyAt));
    double ETchiB = sqrt(etxBt*etxBt + etyBt*etyBt + mEachInvisible*mEachInvisible);
    double AuxMt2_B = visB.m()*visB.m() + mEachInvisible*mEachInvisible + 
      2.0*(visB.ET()*ETchiB - (visB.getLorentzTransverseVector().px()*etxBt + visB.getLorentzTransverseVector().py()*etyBt));
    Mt2W = sqrt(fmax(AuxMt2_A, AuxMt2_B));

    if (Mt2TopWeight==999.) Mt2Top = sqrt(best);
    //else Mt2Top = sqrt(sqrt(fabs(best*best - Mt2W*Mt2W)))/Mt2TopWeight;
    else Mt2Top = sqrt(fabs(best - Mt2W*Mt2W))/Mt2TopWeight;

    std::cout << "Mt2Top: " << Mt2Top << std::endl;

    return best;    

  }
  
  double  Basic_Mt2_Top_Calculator::mT2Fcn::operator()(const std::vector<double>& par) const {
    
    double qT1x = par[0];
    double qT1y = par[1];
    
    double qT2x = theExmiss - qT1x;
    double qT2y = theEymiss - qT1y;
    
    double ETj1 = sqrt(theMass1*theMass1 +  thePT1x*thePT1x + thePT1y*thePT1y);    // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
    double ETj2 = sqrt(theMass2*theMass2 +  thePT2x*thePT2x + thePT2y*thePT2y);     // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless! 
    
    double ETchi1 = sqrt(qT1x*qT1x + qT1y*qT1y + theMchi*theMchi);
    double ETchi2 = sqrt(qT2x*qT2x + qT2y*qT2y + theMchi*theMchi);
    
    double mTsq1 
      = theMass1*theMass1 // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
      + theMchi*theMchi + 2.0*(ETj1*ETchi1 - (thePT1x*qT1x + thePT1y*qT1y));
    double mTsq2 
      = theMass2*theMass2  // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
      + theMchi*theMchi + 2.0*(ETj2*ETchi2 - (thePT2x*qT2x + thePT2y*qT2y));
    
    double mTtsq1 = -1., mTtsq2 = -1.;

    if (theTopStrategy==0) { // "classic" mt2lblb with Mlb (as in AnalysisCMS)       

      float thelbPT1x = thebPT1x + thePT1x;
      float thelbPT1y = thebPT1y + thePT1y;
      float thelbPT2x = thebPT2x + thePT2x;
      float thelbPT2y = thebPT2y + thePT2y;

      double ETlb1 = sqrt(theMlb1*theMlb1 +  thelbPT1x*thelbPT1x + thelbPT1y*thelbPT1y);
      double ETlb2 = sqrt(theMlb2*theMlb2 +  thelbPT2x*thelbPT2x + thelbPT2y*thelbPT2y);

      mTtsq1 = theMlb1*theMlb1 + theMchi*theMchi + 2.0*(ETlb1*ETchi1 - (thelbPT1x*qT1x + thelbPT1y*qT1y) );
      mTtsq2 = theMlb2*theMlb2 + theMchi*theMchi + 2.0*(ETlb2*ETchi2 - (thelbPT2x*qT2x + thelbPT2y*qT2y) );

    } else if (theTopStrategy==1) { // mt2lblb with no Mlb (someone does like this...)
    
      double ETb1 = sqrt(thebMass1*thebMass1 +  thebPT1x*thebPT1x + thebPT1y*thebPT1y);
      double ETb2 = sqrt(thebMass2*thebMass2 +  thebPT2x*thebPT2x + thebPT2y*thebPT2y);

      mTtsq1 = theMass1*theMass1 + thebMass1*thebMass1 + theMchi*theMchi + 
	2.0*(ETj1*ETchi1 + ETb1*ETchi1 + ETb1*ETj1 - 
	     (thePT1x*qT1x + thePT1y*qT1y + thebPT1x*qT1x + thebPT1y*qT1y + thebPT1x*thePT1x + thebPT1y*thePT1y)
	     );
    
      mTtsq2 = theMass2*theMass2 + thebMass2*thebMass2 + theMchi*theMchi + 
	2.0*(ETj2*ETchi2 + ETb2*ETchi2 + ETb2*ETj2 - 
	     (thePT2x*qT2x + thePT2y*qT2y + thebPT2x*qT2x + thebPT2y*qT2y + thebPT2x*thePT2x + thebPT2y*thePT2y)
	     );

    } else if (theTopStrategy==2) { // mt2bb (the leptons are added to the MET; it's physically sounded, but lot of information loss...)

      float theWMass = 80.385;

      double ETb1 = sqrt(thebMass1*thebMass1 +  thebPT1x*thebPT1x + thebPT1y*thebPT1y);
      double ETb2 = sqrt(thebMass2*thebMass2 +  thebPT2x*thebPT2x + thebPT2y*thebPT2y);

      // This is just wrong, this strategy is incompatible with a simultaneous minimisation with mt2ll ...
      double invPT1x = qT1x;
      double invPT1y = qT1y;
      double ETinv1 = sqrt(theWMass*theWMass + invPT1x*invPT1x + invPT1y*invPT1y);

      double invPT2x = qT2x;
      double invPT2y = qT2y;
      double ETinv2 = sqrt(theWMass*theWMass + invPT2x*invPT2x + invPT2y*invPT2y);

      // ... so will use mt2ll MAOS values instead
      mTtsq1 = thebMass1*thebMass1 + theWMass*theWMass + 2.0*(ETb1*ETinv1 - (thebPT1x*invPT1x + thebPT1y*invPT1y) );
      mTtsq2 = thebMass2*thebMass2 + theWMass*theWMass + 2.0*(ETb2*ETinv2 - (thebPT2x*invPT2x + thebPT2y*invPT2y) );
    
    } else if (theTopStrategy==3) { // Same as above, but without adding the leptons to the missing ET (???)

      float theWMass = 80.385;

      double ETb1 = sqrt(thebMass1*thebMass1 +  thebPT1x*thebPT1x + thebPT1y*thebPT1y);
      double ETb2 = sqrt(thebMass2*thebMass2 +  thebPT2x*thebPT2x + thebPT2y*thebPT2y);

      double invPT1x = qT1x;
      double invPT1y = qT1y;
      double ETinv1 = sqrt(theWMass*theWMass + invPT1x*invPT1x + invPT1y*invPT1y);

      double invPT2x = qT2x;
      double invPT2y = qT2y;
      double ETinv2 = sqrt(theWMass*theWMass + invPT2x*invPT2x + invPT2y*invPT2y);

      mTtsq1 = thebMass1*thebMass1 + theWMass*theWMass + 2.0*(ETb1*ETinv1 - (thebPT1x*invPT1x + thebPT1y*invPT1y) );
      mTtsq2 = thebMass2*thebMass2 + theWMass*theWMass + 2.0*(ETb2*ETinv2 - (thebPT2x*invPT2x + thebPT2y*invPT2y) );

    } else if (theTopStrategy==4) { // This should be more correct

      float theWMass = 80.385;

      double ETb1 = sqrt(thebMass1*thebMass1 +  thebPT1x*thebPT1x + thebPT1y*thebPT1y);
      double ETb2 = sqrt(thebMass2*thebMass2 +  thebPT2x*thebPT2x + thebPT2y*thebPT2y);

      double invPT1x = thePT1x + qT1x;
      double invPT1y = thePT1y + qT1y;
      double ETinv1 = sqrt(theWMass*theWMass + invPT1x*invPT1x + invPT1y*invPT1y);

      double invPT2x = thePT2x + qT2x;
      double invPT2y = thePT2y + qT2y;
      double ETinv2 = sqrt(theWMass*theWMass + invPT2x*invPT2x + invPT2y*invPT2y);

      mTtsq1 = thebMass1*thebMass1 + theWMass*theWMass + 2.0*(ETb1*ETinv1 - (thebPT1x*invPT1x + thebPT1y*invPT1y) );
      mTtsq2 = thebMass2*thebMass2 + theWMass*theWMass + 2.0*(ETb2*ETinv2 - (thebPT2x*invPT2x + thebPT2y*invPT2y) );

    }

    // std::cout << "MOO " << sqrt(fmax(mTsq1,mTsq2)) << "\t" << qT1x << "\t" << qT1y << "\t" << 0.5*(qT1x-qT2x) << "\t" << 0.5*(qT1y-qT2y ) << std::endl;
 
    double WMt2 = fmax(mTsq1, mTsq2);
    double TopMt2 = fmax(mTtsq1, mTtsq2);
     
    if (theTopWeight==999.) WMt2 = 0.;
    else TopMt2 *= theTopWeight*theTopWeight;

    // See comment above
    if (theTopStrategy==2) TopMt2 = 0.;

    //return sqrt(WMt2*WMt2 + TopMt2*TopMt2);
    return WMt2 + TopMt2;

  }

} // end of Mt2 Namespace
