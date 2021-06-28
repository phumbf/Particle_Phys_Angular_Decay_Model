#ifndef KPIRESONANCEBASE_H 
#define KPIRESONANCEBASE_H

#include <complex>
#include <string>

#include "Event.h"

// ROOT 
#include "Minuit2/MnUserParameters.h" 


//Base Class to store particle resonance information such as masses, spins etc. 
class KPiResonanceBase {
    
public:
    
    KPiResonanceBase( const double mB,
                      const double mpsi,
                      const double mkaon,
                      const double mpion );
    
    KPiResonanceBase( const KPiResonanceBase& other ) ;
    
    virtual ~KPiResonanceBase() { } ;
    
    virtual std::complex<double> lineshape( const double m, const double pB, const double pR ) const = 0 ;

    virtual std::string getType() const = 0;
    
    virtual unsigned int spin() const = 0;
    
    virtual KPiResonanceBase* clone() const = 0 ;

    virtual void setParameters( const std::vector< double >& par ) ;

 public:
    void setAmplitudes( const std::complex<double>& A0, const std::complex<double>& Ap, const std::complex<double>& Am );

    void setAmplitudesPolar( const double magA0, const double phaseA0 ) ;
    
    void setAmplitudesPolar( const double magA0  , const double magAp, const double magAm,
			     const double phaseA0, const double phaseAp, const double phaseAm ) ;
    
    void setAmplitudesCartesian( const double reA0, const double imA0 );

    void setAmplitudesCartesian( const double reA0, const double imA0, const double reAp, const double imAp, const double reAm, const double imAm  );
    
    std::complex<double> evaluate( const Event& event, const int lambdapsi, const int lambda, const bool isCP ) const ;
    
    std::complex<double> evaluate( const Event& event, const int lambda, const bool isCP ) const ;
    
    std::complex<double> evaluate( const Event& event, const int lambdapsi, const int lambda) const ;
    
    std::complex<double> evaluate( const Event& event, const int lambda) const ;

    bool forbidden( const double m ) const ;
    
 public:
    
    std::string name() const ;
    
    void setName( const std::string& name ) ;
    
    bool isNamed() const ;
    
    bool isNamed( const std::string& name ) const ;

    bool isNameIn( const std::string& name ) const ;

    void usePolar();  
    
    void useCartesian();
    
    double project( const double m ) const ;
    
    unsigned getParameters( ROOT::Minuit2::MnUserParameters& par ) ; 

    void Rescale(const double scalefactor);

public:
    std::complex<double> getA0() const;
    
    std::complex<double> getAp() const;
    
    std::complex<double> getAm() const;

    double getA0Fraction() const ; 
    
    double getApFraction() const ;

    double getAmFraction() const ;
    
    double getAparaFraction() const ;
    
    double getAperpFraction() const ;

 private: 

    double getNorm() const ;
    
protected:
    
    double mB_;
    double mpsi_;
    double mkaon_;
    double mpion_;
    double min_;
    double max_;
    
 private:
    
    bool polar_;
    
    unsigned int offset_; 
    
    std::complex<double> A0_;
    std::complex<double> Ap_;
    std::complex<double> Am_;

    std::complex<double> A0CP_;
    std::complex<double> ApCP_;
    std::complex<double> AmCP_;
    
    std::string name_;
    
    friend std::ostream& operator<<( std::ostream& stream, const KPiResonanceBase& resonance ) ;
};


#endif
