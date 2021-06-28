#include "Event.h"
#include "PsiKPiAmplitudeModel.h"
#include "Variables.h" 
#include "OptionsParser.h" 


#include <iostream>
#include <string> 

#include "TGraph.h"
#include "TFile.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

//Determine the fit fraction via Monte Carlo Integration
double integral(double &sum, int num, double &volume);

int main( int argc, char* argv[] ) {

	OptionsParser parser( argc, argv );
	std::cout << parser << std::endl;

	PsiKPiAmplitudeModel model;
	model.initialise( parser.modelname );

	double mkpi;
	double cospsi{0};
	double coskpi{0};
	double phi{0};
	double allresresult{0};
	std::vector<std::string> resonancevector = {"K*(892)","K0(1400)","K0(800)","K2(1430)","K1(1410)","K1(1680)","K0(1950)","K2(1980)","K3(1780)","K4(2045)","Z(4430)","Z(4200)"};
	std::vector<double> resintvector(12,0.0);

	double mkpi_max = Variables::mB - Variables::mPsi;
	double mkpi_min = Variables::mKaon + Variables::mPion;
	double cos_max = 1.0, cos_min = -1.0, phi_max = TMath::Pi(), phi_min = -TMath::Pi();
	double volume = (mkpi_max-mkpi_min)*(cos_max-cos_min)*(cos_max-cos_min)*(phi_max-phi_min);

	for ( unsigned int i = 0; i < parser.nsample ; ++i ){

		mkpi   = gRandom->Uniform( mkpi_min, mkpi_max);
		coskpi = gRandom->Uniform( cos_min, cos_max );
		cospsi = gRandom->Uniform( cos_min, cos_max );
		phi    = gRandom->Uniform( phi_min, phi_max );

		Event evt(511, mkpi,coskpi,cospsi,phi);

		allresresult += model.evaluate(evt);

		for(unsigned int j=0; j<resonancevector.size(); j++){
			resintvector[j] += (model.evaluate(resonancevector[j],evt));
		} 
	}

	double denom = integral(allresresult,parser.nsample,volume);

	for(unsigned int i=0; i<resonancevector.size();i++){
		double x = (integral(resintvector[i],parser.nsample,volume) / denom);
		std::cout << "Fit fraction for resonance: " << resonancevector[i] << " is: " << x << std::endl;
	}

	return 0;
}

double integral(double &sum, int num, double &volume){

	return (volume/num)*sum;

}
