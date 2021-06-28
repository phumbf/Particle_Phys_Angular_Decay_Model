#ifndef RESONANCEPARSER_H
#define RESONANCEPARSER_H

#include "ResonanceBase.h"
#include "FormFactorBase.h"
#include "LineShapes.h"

#include <vector>
#include <memory>
#include <string>

//Class to allow the resonance parser to be read in 
class Resparser{

 public:
  Resparser();
  
  ~Resparser();
  
  std::vector<std::unique_ptr<ResonanceBase>> getResonances( const std::string filename ) const;
  
  std::unique_ptr<VectorFormFactorBase> getVectorFFs( const std::string filename ) const;
  
  std::unique_ptr<ScalarFormFactorBase> getScalarFFs(const std::string filename ) const;
 
  void setDebug( const bool flag ){ m_debug = flag; }

 private:
  bool m_debug{false};
 
};

#endif

