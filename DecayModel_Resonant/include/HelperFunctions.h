#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

// ROOT
#include "TLorentzVector.h"

// ROOT Math
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

// std
#include <complex> 

/*
 * Useful helper functions to perform useful conversions etc
 */

namespace Helper {
    
    typedef ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<double> > LorentzVector;
    
    void BelleAngles( const int BID,
                     const double mkpi,
                     const double coskpi,
                     const double cospsi,
                     const double phi,
                     double& mz,
                     double& cosz,
                     double& cospsiz,
                     double& phiz,
                     const double mB ,
                     const double mpsi ,
                     const double mkaon,
                     const double mpion,
                     const double mmuon  );
    
    void BelleAngles( const int BID,
                     const TLorentzVector& muplus,
                     const TLorentzVector& muminus,
                     const TLorentzVector& kaon,
                     const TLorentzVector& pion,
                     double& mz,
                     double& cosz,
                     double& cospsiz,
                     double& phiz  );
    
    void HelicityAngles( const int BID,
                        const TLorentzVector& muplus,
                        const TLorentzVector& muminus,
                        const TLorentzVector& kaon,
                        const TLorentzVector& pion,
                        double& mkpi,
                        double& coskpi,
                        double& cospsi,
                        double& phi  );
    
    double momentum( const double m, const double m1, const double m2 ) ;
    
    double cosangle( const LorentzVector& vecA,
                     const LorentzVector& vecB );
    
    double cosangle( const LorentzVector& vecAB,
                     const LorentzVector& vecA,
                     const LorentzVector& vecC );
 
    double BellePlaneAngle( const LorentzVector& vecA,
                            const LorentzVector& vecB,
			    const LorentzVector& vecC );

    void convertToHelicityBasis( const std::complex<double>& Azero, 
				 const std::complex<double>& Apara,
				 const std::complex<double>& Aperp, 
				 std::complex<double>& H0, 
				 std::complex<double>& Hp, 
				 std::complex<double>& Hm );

    void convertToTransversityBasis( const std::complex<double>& H0, 
				     const std::complex<double>& Hp, 
				     const std::complex<double>& Hm, 
				     std::complex<double>& Azero, 
				     std::complex<double>& Apara, 
				     std::complex<double>& Aperp );

    void cpBasisSwap( const std::complex<double>& H0,
                      const std::complex<double>& Hp,
                      const std::complex<double>& Hm,
		       std::complex<double> &H0CP,
		       std::complex<double> &HpCP,
		       std::complex<double> &HmCP, 
		      const int spin );

    void transFracScale( std::complex<double> &Azero,
		       std::complex<double> &Apara,
		       std::complex<double> &Aperp, 
		       const double fpara,
		       const double fperp );

    void fitFracScale(double &scalefitfrac,
			  std::complex<double> &Azero,
			  std::complex<double> &Apara,
			  std::complex<double> &Aperp,
			  const double &fitfrac);

}
#endif 
