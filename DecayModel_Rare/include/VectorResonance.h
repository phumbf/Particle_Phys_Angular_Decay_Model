#ifndef VECTORRESONANCE_H
#define VECTORRESONANCE_H 

#include "ResonanceBase.h"
#include "FormFactorBase.h"
#include <memory>
#include "TH1D.h"

//Class for a vector (i.e. spin > 0) resonance 
class VectorResonanceBase : virtual public ResonanceBase {
  
 public:
  VectorResonanceBase();
  
  ~VectorResonanceBase();
  
  void setFFs( std::unique_ptr<VectorFormFactorBase>& ff );

  VectorFormFactorBase* getFFs() const ;
  
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

 private:
  
  std::unique_ptr<VectorFormFactorBase> m_ff ;
};

#endif
