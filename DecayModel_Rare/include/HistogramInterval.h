#include "TH1.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <utility>

extern std::pair<double,double> centralInterval( std::vector<double>& values, const double CL = 0.68 );

extern std::pair<double,double> shortestInterval( TH1* hist, const double CL = 0.68 ) ;

extern std::pair<double,double> centralInterval( TH1* hist, const double CL = 0.68 ) ;

extern void printHistogramProperties( TH1* hist );
