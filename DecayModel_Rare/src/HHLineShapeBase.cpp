#include "HHLineShapeBase.h"
#include "Parameters.h"
#include "LineShapes.h"

#include "TRandom3.h"

#include <iostream>

//Base class to describe the particle invariant mass lineshape
double HHLineShapeBase::normalise(const double mother, const double daug1, const double daug2) const {

	double integral = 0.0;

	const unsigned int numiter = 10000000;

	double mkpimax = mother - 2*Parameters::mMuon;
	//Currently hard coded for kpi
	double mkpimin = daug1 + daug2; 
	double volume = mkpimax - mkpimin;
	double mkpi{0};

	for( unsigned int i=0; i < numiter; i++ ){
		mkpi = gRandom->Uniform(mkpimin,mkpimax);
		integral += std::norm( this->lineshape(mkpi) );
	}

	integral *= (volume/double(numiter));

	double norm = std::sqrt( 1.0 / integral );

	return norm;

}

double HHLineShapeBase::shape( const double mkpi ) const{
	
	return std::norm( this->lineshape(mkpi) );
}
