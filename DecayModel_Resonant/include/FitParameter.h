#ifndef FITPARAMETER_H 
#define FITPARAMETER_H

#include <limits>   
#include <algorithm>
#include <string>

#include "Minuit2/MnUserParameters.h"

/*
 * A class to hold parameters for likelihood fits
 * */

class FitParameter {
public:
    
    
    
    FitParameter() :
        name_     ( "" ),
        floating_ ( false ),
        min_      ( std::numeric_limits<double>::min() ),
        max_      ( std::numeric_limits<double>::max() ),
        err_      ( std::numeric_limits<double>::max() ),
        val_      ( 0 ) {} ;
    
    FitParameter( const std::string& name, const double val ) :
        name_     ( name ),
        floating_ ( false ),
        min_      ( std::numeric_limits<double>::min() ),
        max_      ( std::numeric_limits<double>::max() ),
        err_      ( std::numeric_limits<double>::max() ),
        val_      ( val ) {} ;
  
    FitParameter( const std::string& name, const double val, const double min, const double max, const double err ) :
        name_     ( name ),
        floating_ ( true ),
        min_      ( min  ),
        max_      ( max  ),
        err_      ( err  ),
        val_      ( val ){} ;
    
    FitParameter( const FitParameter& other ) :
        name_     ( other.name_ ),
        floating_ ( other.floating_ ),
        min_      ( other.min_ ),
        max_      ( other.max_ ),
        err_      ( other.err_ ),
        val_      ( other.val_ ){} ;
    
    ~FitParameter() {} ;
    
    bool isFloating() const {
        return floating_;
    }
    
    void floatInRange( const double min, const double max, const double err ) {
        min_      = min;
        max_      = max;
        err_      = err;
        floating_ = true;
    }
    
    double min() const {
        return min_;
    }
    
    double max() const {
        return max_;
    }
    
    double err() const {
        return err_;
    }
    
    bool hasChanged( const double val ) const {
        return ( fabs( val_  - val ) < std::numeric_limits<double>::epsilon()  ) ;
    }
    
    void set( double val ) {
        val_ = val;
    }
    
    int parameter( const std::string& parent, ROOT::Minuit2::MnUserParameters& par ) const ;
    
    
private:
    std::string name_ ;
    
    bool floating_;
    
    double min_;
    double max_;
    double err_;
    double val_;
};

#endif
