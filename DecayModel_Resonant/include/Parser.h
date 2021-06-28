#ifndef PARSER_H 
#define PARSER_H 

#include "KPiResonanceBase.h"
#include "PsiPiResonance.h" 

#include <string>
#include <map> 

#include <boost/property_tree/json_parser.hpp>

//Class to Parse in resonances from a json list
extern void parse_resonance_list( const std::string json,
                                  std::vector<KPiResonanceBase*>& kpi_resonances,
				  std::vector<PsiPiResonance*  >& psipi_resonances,
				  const bool debug = false ) ;


extern void update_fraction( const std::string json, 
			     std::map< std::string, double >& values );

extern std::map<std::string,double> get_target_fractions( const std::string json); 

#endif
