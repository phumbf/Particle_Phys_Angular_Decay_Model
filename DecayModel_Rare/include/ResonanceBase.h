#ifndef RESONANCEBASE_H
#define RESONANCEBASE_H 

#include <string>
#include <vector>
#include <iostream>
#include "AmplitudeSet.h"
#include "HHLineShapeBase.h"
#include "FormFactorBase.h"
#include "Wilson.h"
#include "Parameters.h"
#include "TH1D.h"

//The base class for a resonant object
class ResonanceBase {

 public:
  
  virtual ~ResonanceBase() {} ;
  
  void print() const;
  
  void setMotherMass( const double mothermass );
  
  void setName( const std::string name );

  void setBFfac( const double BFfac );
  
  void setisCP( const bool isCP );

  bool isCP();
  
  void setisEvtGenWC( const bool isEvtGenWC );

  bool isEvtGenWC();

  void setisbtod( const bool isbtod );

  bool isbtod();

  std::string getName() const;
  
  bool hasName( const std::string name ) const ;

  void setMass( const double mass);
  
  double getMotherMass() const;
  
  double getMass() const;
  
  int getL() const;
  
  void setL( const int l );

  void setDaughterMasses( const double mone, const double mtwo );

  void setShape( std::unique_ptr<HHLineShapeBase>& shape );
  
  virtual std::complex<double> APERP( const double s,
				      const double m,
				      const double LR ) const = 0;
  
  virtual std::complex<double> APARA( const double s,
				      const double m,
				      const double LR ) const = 0;
  
  
  virtual std::complex<double> AZERO( const double s,
				      const double m,
				      const double LR ) const = 0;

  virtual std::complex<double> ATIME( const double s,
				      const double m ) const = 0;

  
  void getAmplitudes( const double s,
		      const double mhh, 
		      AmplitudeSet& amp) const;

  void getAmplitudesCPswap( const double s,
                            const double mhh,
                            AmplitudeSet& amp) const;
  
  double normalisation( const double s, const double m ) const ;

  void setCorrections(  const std::complex<double>& cZERO, 
		        const std::complex<double>& cPARA,
		        const std::complex<double>& cPERP ) ;

  void setCorrections( const std::complex<double>& cZERO );

  void setWilsons( Wilson &wc,
		   Wilson &wp);

  double prefactor( const int mproj,
		    const double mass,
		    const double lambdaM,
		    const double lambdaH) const;
 
  void setPolarResFactor( const double alpha,
                          const double phi);

  void setResFactor( const double alpha,
                     const double phi);

  Wilson m_WC;
  Wilson m_WP;

  private: 
 
  double m_mass{0};
  
  int m_L{0};
  
  std::string m_name{""};

 protected:
  
  double m_mothermass{Parameters::mB};
  double m_mone{Parameters::mPion};
  double m_mtwo{Parameters::mPion};

  double m_branching{1.0};
  bool m_isCP{false};
  bool m_isEvtGenWC{false};
  bool m_isbtod{false};

  std::complex<double> m_resfactor{1.0,0.0};
  
  std::unique_ptr<HHLineShapeBase> m_lineshape;
  // non-factorisable corrections
  std::complex<double> m_cPARA = std::complex<double>(0,0); 
  std::complex<double> m_cPERP = std::complex<double>(0,0); 
  std::complex<double> m_cZERO = std::complex<double>(0,0); 
};


#endif
