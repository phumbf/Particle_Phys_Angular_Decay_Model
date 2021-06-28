#ifndef SCALARRESONANCE_H
#define SCALARRESONANCE_H

#include "ResonanceBase.h"
#include "FormFactorBase.h"

#include <memory>

//Class for a scalar (i.e. spin = 0) resonance 
class ScalarResonanceBase : virtual public ResonanceBase {

 public:
  
  ScalarResonanceBase();
  
  ~ScalarResonanceBase();
  
  void setFFs( std::unique_ptr<ScalarFormFactorBase>& ff );
  
  ScalarFormFactorBase* getFFs() const;
  
  std::complex<double> APERP( const double s,
			      const double m, 
			      const double LR ) const override;
  
  std::complex<double> APARA( const double s,
			      const double m, 
			      const double LR ) const override;
  
  
  std::complex<double> AZERO( const double s,
			      const double m,
			      const double LR ) const override;
  
  std::complex<double> ATIME( const double s,
			      const double m ) const override;


  std::unique_ptr<ScalarFormFactorBase>  m_ff ;
};

#endif
