#ifndef KPIRESONANCEFLATTE_H 
#define KPIRESONANCEFLATTE_H

#include "KPiResonanceBase.h" 

#include <complex>
#include <string>

//Derived class to implement the Flatte resonance lineshape
//This is useful for channels where there are more than one set of
//contributions, e.g. pi-pi, pi-K etc.
//https://cds.cern.ch/record/2719235/files/scoap3-fulltext.pdf

class KPiResonanceFlatte : public KPiResonanceBase { 
public: 
  
  KPiResonanceFlatte( const double mean,
		      const double gkk,
		      const double gpipi,
		      const double alpha,
		      const double mB,
		      const double mpsi,
		      const double mkaon,
		      const double mpion,
		      const unsigned int LB,
		      const double radius );
  
  KPiResonanceFlatte(const KPiResonanceFlatte& other);
  
  virtual ~KPiResonanceFlatte() {};
  
  std::complex<double> lineshape( const double m, 
				  const double pB,
				  const double pR) const;
  
  KPiResonanceBase* clone() const{
    return new KPiResonanceFlatte(*this);
  }
  
  unsigned int spin() const;
  
  std::string getType() const { return "Flatte"; }
  
 private:
  double mean_;
  double gkk_;
  double gpipi_;
  double alpha_;
  
  unsigned int LB_;

  double radius_;

};

inline unsigned int KPiResonanceFlatte::spin() const {
  return 0;
}


#endif
