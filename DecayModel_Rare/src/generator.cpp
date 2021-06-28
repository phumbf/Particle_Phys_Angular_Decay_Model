#include "random"
#include "generator.h"

#include "TVectorD.h"
#include "TRandom3.h"

#include <iostream>

generator::generator(){}

generator::generator(const UInt_t seed){
	m_rand = TRandom(seed);
}

generator::~generator(){}

void generator::setSeed(const UInt_t seed){

	m_rand.SetSeed(seed);
}

double generator::fluctuate(const double mean, const double error) {

	return m_rand.Gaus(mean,error);
}

void generator::fluctuate(const std::vector<double> &mean,
	       const std::vector<double> &error,
	       std::vector<double> &vec) {

	for(unsigned int i=0; i<mean.size(); i++){
	
		vec.push_back(m_rand.Gaus(mean[i],error[i]));
	}
}

void generator::fluctuate(const TVectorD &mean,
	       const TVectorD &error,
	       std::vector<double> &vec) {

	for(Int_t i=0; i<mean.GetNoElements(); i++){
		
		vec.push_back(m_rand.Gaus(mean[i],error[i]));
	}
}

TVectorD generator::fluctuate(const TVectorD &mean,
		          const TVectorD &error){

	TVectorD vec(mean.GetNoElements());

	for(Int_t i=0; i<mean.GetNoElements(); i++){
		
		vec[i] = m_rand.Gaus(mean[i],error[i]);
	}

	return vec;
}

