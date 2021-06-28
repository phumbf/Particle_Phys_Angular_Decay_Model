#include "FitParameter.h"

/*
 * A class to hold parameters for likelihood fits
 * */

int FitParameter::parameter( const std::string& parent, ROOT::Minuit2::MnUserParameters& par ) const {
    if ( floating_ ){
        par.Add( parent + "_" + name_ , val_ , err_ , min_ , max_ ) ;
        return 1;
    }
    return 0;
}
