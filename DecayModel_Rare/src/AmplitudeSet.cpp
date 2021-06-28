// Local
#include "AmplitudeSet.h"

// std
#include <iostream>


//Class to store the amplitude polarisations
AmplitudeSet::AmplitudeSet( const AmplitudeSet& other ){ 
  AZEROL = other.AZEROL; 
  AZEROR = other.AZEROR; 
  APARAL = other.APARAL;
  APARAR = other.APARAR;
  APERPL = other.APERPL;
  APERPR = other.APERPR;
  ATIME  = other.ATIME ;
}

AmplitudeSet::AmplitudeSet( const std::complex<double> aZEROL,
		            const std::complex<double> aZEROR,
			    const std::complex<double> aPARAL,
			    const std::complex<double> aPARAR,
			    const std::complex<double> aPERPL,
			    const std::complex<double> aPERPR,
			    const std::complex<double> aTIME  ){
  
  AZEROL = aZEROL;
  AZEROR = aZEROR;
  APARAL = aPARAL;
  APARAR = aPARAR;
  APERPL = aPERPL;
  APERPR = aPERPR;
  ATIME  = aTIME ;
}

AmplitudeSet::AmplitudeSet( ){}  

void AmplitudeSet::print() const { 
  std::cout 
    << "  A_{zero}^{L} = " << AZEROL 
    << ", A_{zero}^{R} = " << AZEROR 
    << ", A_{para}^{L} = " << APARAL 
    << ", A_{para}^{R} = " << APARAR 
    << ", A_{perp}^{L} = " << APERPR
    << ", A_{perp}^{R} = " << APERPL 
    << ", A_{time}     = " << ATIME 
    << std::endl;    
}

void AmplitudeSet::flipCP(const int spin){

	if(spin % 2 == 0){

		AZEROR = -AZEROR;
		AZEROL = -AZEROL;
		APARAR = -APARAR;
		APARAL = -APARAL;
	}
	else{
		APERPL = -APERPL;
		APERPR = -APERPR;
	}
}
