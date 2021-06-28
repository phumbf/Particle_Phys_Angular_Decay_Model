#ifndef FLUCTUATION_H  
#define FLUCTUATION_H 1

#include "TVectorD.h"
#include "TRandom3.h"

//Generator class for random number generation
class generator{ 

	public: 
		TRandom m_rand;

		generator(const UInt_t seed);
		generator();
		~generator();

		void setSeed(const UInt_t seed);

		double fluctuate(const double mean, const double error) ;

		void fluctuate(const std::vector<double> &mean,
				const std::vector<double> &error,
				std::vector<double> &vec );

		void fluctuate(const TVectorD &mean,
				const TVectorD &error,
				std:: vector<double> &vec);

		TVectorD fluctuate(const TVectorD &mean,
				const TVectorD &error);
}; 

#endif
