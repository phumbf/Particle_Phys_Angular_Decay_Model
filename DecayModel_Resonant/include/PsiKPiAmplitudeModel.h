#ifndef PSIKPIAMPLITUDEMODEL_H
#define PSIKPIAMPLITUDEMODEL_H

#include "Event.h" 
#include "KPiResonanceBase.h"
#include "PsiPiResonance.h" 

#include <vector>
#include <string>

// ROOT 
#include "TTree.h" 
#include "TRandom3.h"
#include "Minuit2/MnUserParameters.h" 


/*
 * Amplitude model for B --> psi K pi, including
 * sum of conventional (K pi) resonances and exotic
 * (psi pi) states. 
 *
 * Model is configured using an external json file in 
 * PsiKPIAmplitudeModel::initialise 
 * 
 */


class PsiKPiAmplitudeModel {
public:
    PsiKPiAmplitudeModel() ;

    PsiKPiAmplitudeModel( const std::string& json ) ; 
    
    ~PsiKPiAmplitudeModel() ;
    
    PsiKPiAmplitudeModel( const PsiKPiAmplitudeModel& other );
    
    double operator()( const double* x ) const ;
    
    double evaluate( const Event& event )  const ;
    
    double evaluateCP( const Event& event ) const ;

    double evaluateCP( const std::string& name, const Event& event ) const ;

    double evaluate( const unsigned int j, const Event& event ) const ;
    
    double evaluate( const std::string& name, const Event& event ) const ;

    void initialise( const std::string& json );

    unsigned getParameters( ROOT::Minuit2::MnUserParameters& par ) ; 
    
    void setParameters( const std::vector< double >& par ) ;
    
    const KPiResonanceBase* getResonance( const std::string& name ) const ;
    
 public:
    std::complex<double>  getAmplitude( const int lambda, const Event& event, const bool isCP) const;
     
    std::complex<double>  getAmplitude( const std::string& name, const int lambda, const Event& event, const bool isCP) const;

    std::complex<double>  getAmplitude( const unsigned int j, const int lambda, const Event& event, const bool isCP ) const;

    std::complex<double>  getAmplitude( const int lambda, const Event& event) const;
     
    std::complex<double>  getAmplitude( const std::string& name, const int lambda, const Event& event) const;

    std::complex<double>  getAmplitude( const unsigned int j, const int lambda, const Event& event) const;

    std::map<std::string,double> fitfractions(const bool isCP);

    void rescaleres(const std::string& refres, std::map<std::string,double>& fitfracs, std::map<std::string,double>& targetfracs);

    double integrate(std::string name, const unsigned int numevts, const bool isCP);

    std::vector<std::string> resonancelist();

    void print();

private:
    
    void clearResonanceList() ; 
    
    // Summed conventional and exotic resonances
    std::vector< KPiResonanceBase* > kpi_resonances_ ;
    std::vector< PsiPiResonance*   > psipi_resonances_; 
    
};

#endif
