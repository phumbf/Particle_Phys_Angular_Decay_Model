
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>

#include "TH1.h"

//#include "TGraph.h"
//#include "TGraphAsymmErrors.h"

std::pair<double,double> centralInterval( TH1* hist, const double CL = 0.68 ){
    
    int low = 0;
    int upp = hist->GetNbinsX() + 1;
    
    const double limit = 0.5*(1.0-CL)*hist->Integral();
    
    double sumlow  = 0;
    double sumupp  = 0;
    
    for ( int i = 1; i <= hist->GetNbinsX(); ++i ){
        sumlow += hist->GetBinContent( i );
        sumupp += hist->GetBinContent( hist->GetNbinsX() + 1 - i );
        
        if ( (low == 0) && (sumlow >limit) ) low = i;
        if ( (upp == hist->GetNbinsX() + 1) && (sumupp > limit) ) upp = hist->GetNbinsX() + 1 - i;
    }
    
    return std::make_pair( hist->GetXaxis()->GetBinLowEdge(low),
                          hist->GetXaxis()->GetBinUpEdge(upp) );
}

std::pair<double,double> shortestInterval( TH1* hist, const double CL = 0.68 ){
    
    const double limit = CL*hist->Integral();
    
    const int mpv   = hist->GetMaximumBin();
    const int nbins = hist->GetNbinsX() ;
    
    int low = mpv;
    int upp = mpv;
    
    double integral = hist->GetBinContent( mpv );
    double lowval   = 0;
    double uppval   = 0;
    
    while ( integral < limit ){
        lowval = ( low - 1 >= 1 ) ? hist->GetBinContent( low - 1 ) : 0 ;
        uppval = ( upp + 1 <= nbins ) ? hist->GetBinContent( upp + 1 ) : 0 ;
        
        if ( lowval > uppval ){
            integral += lowval;
            low--;
        }
        else {
            integral += uppval;
            upp++;
        }
    }
    
    
    return std::make_pair( hist->GetXaxis()->GetBinLowEdge(low),
                          hist->GetXaxis()->GetBinUpEdge(upp) );
}


std::pair<double,double> centralInterval( std::vector<double>& values, const double CL = 0.68 ){
    
    std::sort( values.begin(), values.end() );
    
    unsigned int lowval = (0.5*(1.0-CL))*values.size();
    unsigned int uppval = (1.0 - 0.5*(1.0-CL))*values.size();
    
    return std::make_pair( values[ lowval ], values[ uppval ] ) ;
}



void printHistogramProperties( TH1* hist ) {
    
    const double mpv  = hist->GetXaxis()->GetBinCenter( hist->GetMaximumBin() );
    const double mean = hist->GetMean() ;
    const double rms  = hist->GetRMS() ;
    
    std::cout << " Mean = " << mean << "\n"
              << " MPV  = " << mpv  << "\n"
              << " RMS  = " << rms  << std::endl;
    
    std::pair<double,double> limits  = centralInterval( hist, 0.68 );
    
    std::cout << " Limits (68%) = ( "
              << limits.first  << " , "
              << limits.second << " )"
              << std::endl;
    
   limits = shortestInterval( hist, 0.68 );
    
    std::cout << " Shorest interval (68%) = ( "
              << limits.first  << " , "
              << limits.second << " )"
              << std::endl;
   
    
    return;
}

