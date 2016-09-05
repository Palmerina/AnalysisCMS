// Header file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#ifndef Basic_Mt2_Top_Calculator_H
#define Basic_Mt2_Top_Calculator_H

#include "Mt2/Mt2_Top_Calculator.h"
#include "Minuit2/FCNBase.h"
#include <vector>

namespace Mt2 {

  class Basic_Mt2_Top_Calculator : public Mt2_Top_Calculator {
  public:
    double mt2_Top(const LorentzVector& visibleA,  // 3 d.o.f. 
		   const LorentzVector& visibleB,  // 3 d.o.f.
		   const LorentzVector& visibleb1, // 3 d.o.f. 
		   const LorentzVector& visibleb2, // 3 d.o.f.
		   const TwoVector& ptmiss,        // 2 d.o.f.
		   const double mInvisible,
		   const int Mt2TopStrategy, const float Mt2TopWeight);
    double mt2_Top_Sq(const LorentzVector& visibleA,  // 3 d.o.f. 
		      const LorentzVector& visibleB,  // 3 d.o.f.
		      const LorentzVector& visibleb1, // 3 d.o.f. 
		      const LorentzVector& visibleb2, // 3 d.o.f.
		      const TwoVector& ptmiss,        // 2 d.o.f.
		      const double mInvisible,
		      const int Mt2TopStrategy, const float Mt2TopWeight);
    Basic_Mt2_Top_Calculator() : Mt2_Top_Calculator("Basic_Mt2_Top") {};

      double getPXInvisA_atMt2Solution() const {return etxAt;}
      double getPYInvisA_atMt2Solution() const {return etyAt;}
      double getPXInvisB_atMt2Solution() const {return etxBt;}
      double getPYInvisB_atMt2Solution() const {return etyBt;}
      double getMt2W_atMt2Solution() const {return Mt2W;}
      double getMt2Top_atMt2Solution() const {return Mt2Top;}

  private:
      double etxAt;
      double etyAt;
      double etxBt;
      double etyBt;
      double Mt2W;
      double Mt2Top;

    class mT2Fcn : public ROOT::Minuit2::FCNBase {
    public:
      mT2Fcn(const double exmiss, 
	     const double eymiss, 
	     const double mchi,
	     const double pT1x,
	     const double pT1y,
	     const double mass1,
	     const double pT2x,
	     const double pT2y,
	     const double mass2,
	     const double bpT1x,
	     const double bpT1y,
	     const double bmass1,
	     const double bpT2x,
	     const double bpT2y,
	     const double bmass2,
	     const double Mlb1,
	     const double Mlb2,
	     const int Mt2TopStrategy, const float Mt2TopWeight) : 
	theExmiss(exmiss),
	theEymiss(eymiss),
        theMchi(mchi),
	thePT1x(pT1x),
	thePT1y(pT1y),
	theMass1(mass1),
	thePT2x(pT2x),
	thePT2y(pT2y),
	theMass2(mass2),
	thebPT1x(bpT1x),
	thebPT1y(bpT1y),
	thebMass1(bmass1),
	thebPT2x(bpT2x),
	thebPT2y(bpT2y),
        thebMass2(bmass2),
	theMlb1(Mlb1), 
        theMlb2(Mlb2),
	theTopStrategy(Mt2TopStrategy),
	theTopWeight(Mt2TopWeight),
	theErrorDef(1.) {}
	
      ~mT2Fcn() {}
      
      virtual double Up() const {return theErrorDef;}
      virtual double operator()(const std::vector<double>&) const;
      
      const double exMiss() const {return theExmiss;}
      const double eyMiss() const {return theEymiss;}
      const double mChi() const {return theMchi;}
      const double pT1x() const {return thePT1x;}
      const double pT1y() const {return thePT1y;}
      const double pT2x() const {return thePT2x;}
      const double pT2y() const {return thePT2y;}
      const double bpT1x() const {return thebPT1x;}
      const double bpT1y() const {return thebPT1y;}
      const double bpT2x() const {return thebPT2x;}
      const double bpT2y() const {return thebPT2y;}
      
      void setErrorDef(double def) {theErrorDef = def;}
      
    private:

      double theExmiss;
      double theEymiss;
      double theMchi;
      
      double thePT1x;
      double thePT1y;
      double theMass1;
      
      double thePT2x;
      double thePT2y;
      double theMass2;
      
      double thebPT1x;
      double thebPT1y;
      double thebMass1;
      
      double thebPT2x;
      double thebPT2y;
      double thebMass2;
      
      double theMlb1;
      double theMlb2;

      int    theTopStrategy;
      double    theTopWeight;

      double theErrorDef;


    };
  };

} //end namespace Mt2

#endif
