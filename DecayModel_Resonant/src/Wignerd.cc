#include "Wignerd.h"
#include "TMath.h"

#include <iostream>
#include <cassert>

#include "TGraph.h"

#include "TLatex.h"
#include "TCanvas.h"
#include "TAxis.h"

//Class to calculate Wigner D function contributions

double Wignerd::djmpm( const double costheta, const unsigned int j, int mp, int m ){
   
    
    if ( (unsigned)abs( mp ) > j || (unsigned)abs( m ) > j ) {
        return 0;
    }
    
    if ( abs(m) > abs(mp) ){
        return  std::pow( -1, m-mp )*Wignerd::djmpm( costheta, j, m, mp ) ;
    }
    
    
    if ( mp <= 0 ) {
        if ( m < 0 ){
            return Wignerd::djmpm( costheta, j, -m, -mp ) ;
        }
        else if ( m > 0 ){
            return std::pow( -1, m-mp )*Wignerd::djmpm( costheta, j, m, mp ) ;
        }
        else if ( m == 0 && mp < 0 ){
            return Wignerd::djmpm( costheta, j, -m, -mp ) ;
        }
    }
    
    if      ( j == 0 ){
        if  ( mp == 0 && m == 0 ) return 1.0;
    }
    else if ( j == 1 ){
        if      ( mp == 0 && m ==  0 ) return d100 ( costheta );
        else if ( mp == 1 && m == -1 ) return d11m1( costheta );
        else if ( mp == 1 && m ==  0 ) return d110 ( costheta );
        else if ( mp == 1 && m ==  1 ) return d11p1( costheta );
    }
    else if ( j == 2 ){
        if      ( mp == 0 && m ==  0 ) return d200 ( costheta );
        else if ( mp == 1 && m == -1 ) return d21m1( costheta );
        else if ( mp == 1 && m ==  0 ) return d210 ( costheta );
        else if ( mp == 1 && m ==  1 ) return d21p1( costheta );
        else if ( mp == 2 && m == -2 ) return d22m2( costheta );
        else if ( mp == 2 && m == -1 ) return d22m1( costheta );
        else if ( mp == 2 && m ==  0 ) return d220 ( costheta );
        else if ( mp == 2 && m ==  1 ) return d22p1( costheta );
        else if ( mp == 2 && m ==  2 ) return d22p2( costheta );
    }
    else if ( j == 3 ){
        if      ( mp == 0 && m ==  0 ) return d300 ( costheta );
        else if ( mp == 1 && m == -1 ) return d31m1( costheta );
        else if ( mp == 1 && m ==  0 ) return d310 ( costheta );
        else if ( mp == 1 && m ==  1 ) return d31p1( costheta );
        else if ( mp == 2 && m == -2 ) return d32m2( costheta );
        else if ( mp == 2 && m == -1 ) return d32m1( costheta );
        else if ( mp == 2 && m ==  0 ) return d320 ( costheta );
        else if ( mp == 2 && m ==  1 ) return d32p1( costheta );
        else if ( mp == 2 && m ==  2 ) return d32p2( costheta );
        else if ( mp == 3 && m == -3 ) return d33m3( costheta );
        else if ( mp == 3 && m == -2 ) return d33m2( costheta );
        else if ( mp == 3 && m == -1 ) return d33m1( costheta );
        else if ( mp == 3 && m ==  0 ) return d330 ( costheta );
        else if ( mp == 3 && m ==  1 ) return d33p1( costheta );
        else if ( mp == 3 && m ==  2 ) return d33p2( costheta );
        else if ( mp == 3 && m ==  3 ) return d33p3( costheta );
    }
    else if ( j == 4 ){
        if      ( mp == 0 && m ==  0 ) return d400 ( costheta );
        else if ( mp == 1 && m == -1 ) return d41m1( costheta );
        else if ( mp == 1 && m ==  0 ) return d410 ( costheta );
        else if ( mp == 1 && m ==  1 ) return d41p1( costheta );
        else if ( mp == 2 && m == -2 ) return d42m2( costheta );
        else if ( mp == 2 && m == -1 ) return d42m1( costheta );
        else if ( mp == 2 && m ==  0 ) return d420 ( costheta );
        else if ( mp == 2 && m ==  1 ) return d42p1( costheta );
        else if ( mp == 2 && m ==  2 ) return d42p2( costheta );
        else if ( mp == 3 && m == -3 ) return d43m3( costheta );
        else if ( mp == 3 && m == -2 ) return d43m2( costheta );
        else if ( mp == 3 && m == -1 ) return d43m1( costheta );
        else if ( mp == 3 && m ==  0 ) return d430 ( costheta );
        else if ( mp == 3 && m ==  1 ) return d43p1( costheta );
        else if ( mp == 3 && m ==  2 ) return d43p2( costheta );
        else if ( mp == 3 && m ==  3 ) return d43p3( costheta );
        else if ( mp == 4 && m == -4 ) return d44m3( costheta );
        else if ( mp == 4 && m == -3 ) return d44m3( costheta );
        else if ( mp == 4 && m == -2 ) return d44m2( costheta );
        else if ( mp == 4 && m == -1 ) return d44m1( costheta );
        else if ( mp == 4 && m ==  0 ) return d440 ( costheta );
        else if ( mp == 4 && m ==  1 ) return d44p1( costheta );
        else if ( mp == 4 && m ==  2 ) return d44p2( costheta );
        else if ( mp == 4 && m ==  3 ) return d44p3( costheta );
        else if ( mp == 4 && m ==  4 ) return d44p4( costheta );
    }
    else {
        return Wignerd::generic( acos(costheta), j, mp, m );
    }
    return 0;
}




// J = 0

double Wignerd::d000( const double /*costheta*/ ) {
    return 1.0;
}

// J = 1

double Wignerd::d100( const double costheta ){
    return costheta;
}

double Wignerd::d11p1( const double costheta ){
    return ( 1. + costheta )/2.;
}

double Wignerd::d110( const double costheta ){
    return ( -std::sqrt( (1. - costheta*costheta )/2. ) );
}

double Wignerd::d11m1( const double costheta ){
    return ( 1. - costheta )/2.;
}


// J = 2

double Wignerd::d22p2( const double costheta ){
    return std::pow( (1. + costheta)/2. , 2 );
}

double Wignerd::d22p1( const double costheta ){
    return -(1. + costheta)*std::sqrt( 1. - costheta*costheta )/2. ;
}

double Wignerd::d220( const double costheta ){
    return std::sqrt(3./8.)*( 1. - costheta*costheta );
    
}

double Wignerd::d22m1( const double costheta ){
    return -(1. - costheta)*std::sqrt( 1. - costheta*costheta )/2. ;
}

double Wignerd::d22m2( const double costheta ){
    return std::pow( (1. - costheta)/2. , 2 );
}

double Wignerd::d21p1( const double costheta ){
    return ( 1. + costheta )*(2.*costheta - 1. )/2.;
}

double Wignerd::d21m1( const double costheta ){
    return ( 1. - costheta )*(2.*costheta + 1. )/2.;
}

double Wignerd::d210( const double costheta ){
    return -std::sqrt( 1. - costheta*costheta )*costheta*std::sqrt(3./2.);
}

double Wignerd::d200( const double costheta ){
    return (3./2)*costheta*costheta - (1./2.);
}


// J = 3

double Wignerd::d33m3( const double costheta ){
    return -std::pow(-1. + costheta, 3)/8. ;
}

double Wignerd::d33m2( const double costheta ){
    return -(std::pow(-1 + costheta,2)*std::sqrt(6 - 6* std::pow(costheta,2)))/8.;
}

double Wignerd::d33m1( const double costheta ){
    return (std::sqrt(15.)*std::pow(-1. + costheta,2)*(1. + costheta))/8.;
}

double Wignerd::d330( const double costheta ){
    return -(std::sqrt(5.)*std::pow(1. - std::pow(costheta,2),1.5))/4.;
}

double Wignerd::d33p1( const double costheta ){
    return -(std::sqrt(15.)*(-1. + costheta)*std::pow(1. + costheta,2))/8.;
}

double Wignerd::d33p2( const double costheta ){
    return -(std::pow(1. + costheta,2)*std::sqrt(6. - 6.*std::pow(costheta,2)))/8.;
}

double Wignerd::d33p3( const double costheta ){
    return std::pow(1. + costheta,3)/8.;
}


double Wignerd::d32m2( const double costheta ){
    return (std::pow(-1 + costheta,2)*(2 + 3*costheta))/4.;
}

double Wignerd::d32m1( const double costheta ){
    return -(std::sqrt(2.5)*(1. + (2. - 3.*costheta)*costheta)*std::sqrt(1 - std::pow(costheta,2)))/4.;
}

double Wignerd::d320( const double costheta ){
    return -(std::sqrt(7.5)*costheta*(-1 + std::pow(costheta,2)))/2.;
}

double Wignerd::d32p1( const double costheta ){
    return -(std::sqrt(2.5)*(1 + costheta)*(-1 + 3*costheta)*std::sqrt(1 - std::pow(costheta,2)))/4.;
}

double Wignerd::d32p2( const double costheta ){
    return (std::pow(1 + costheta,2)*(-2 + 3*costheta))/4.;
}

double Wignerd::d31m1( const double costheta ){
    return (-1 + costheta*(11 + 5*(1 - 3*costheta)*costheta))/8.;
}

double Wignerd::d310( const double costheta ){
    return -(std::sqrt(3 - 3*std::pow(costheta,2))*(-1 + 5*std::pow(costheta,2)))/4.;
}

double Wignerd::d31p1( const double costheta ){
    return ((1 + costheta)*(-1 + 5*costheta*(-2 + 3*costheta)))/8.;
}

double Wignerd::d300( const double costheta ){
    return (costheta*(-3 + 5*std::pow(costheta,2)))/2.;
}
// J = 4

double Wignerd::d44m4( const double costheta ){
    return std::pow(-1 + costheta,4)/16.;
}

double Wignerd::d44m3( const double costheta ){
    return (std::pow(-1 + costheta,3)*std::sqrt(2 - 2*std::pow(costheta,2)))/8.;
}

double Wignerd::d44m2( const double costheta ){
    return -(std::sqrt(7)*std::pow(-1 + costheta,3)*(1 + costheta))/8.;
}

double Wignerd::d44m1( const double costheta ){
    return -(std::sqrt(3.5)*std::pow(-1 + costheta,2)*(1 + costheta)*std::sqrt(1 - std::pow(costheta,2)))/4.;
}

double Wignerd::d440( const double costheta ){
    return (std::sqrt(17.5)*std::pow(-1 + std::pow(costheta,2),2))/8.;
}

double Wignerd::d44p1( const double costheta ){
    return (std::sqrt(3.5)*(-1 + costheta)*std::pow(1 + costheta,2)*std::sqrt(1 - std::pow(costheta,2)))/4.;
}

double Wignerd::d44p2( const double costheta ){
    return -(std::sqrt(7)*(-1 + costheta)*std::pow(1 + costheta,3))/8.;
}

double Wignerd::d44p3( const double costheta ){
    return -(std::pow(1 + costheta,3)*std::sqrt(2 - 2*std::pow(costheta,2)))/8.;
}

double Wignerd::d44p4( const double costheta ){
    return std::pow(1 + costheta,4)/16.;
}

// 3

double Wignerd::d43m3( const double costheta ){
    return -(std::pow(-1 + costheta,3)*(3 + 4*costheta))/8.;
}

double Wignerd::d43m2( const double costheta ){
    return -(std::sqrt(3.5)*std::pow(-1 + costheta,2)*(1 + 2*costheta)*std::sqrt(1 - std::pow(costheta,2)))/4.;
}

double Wignerd::d43m1( const double costheta ){
    return (std::sqrt(7)*std::pow(-1 + costheta,2)*(1 + costheta)*(1 + 4*costheta))/8.;
}

double Wignerd::d430( const double costheta ){
    return -(std::sqrt(35)*costheta*std::pow(1 - std::pow(costheta,2),1.5))/4.;
}

double Wignerd::d43p1( const double costheta ){
    return -(std::sqrt(7)*(-1 + costheta)*std::pow(1 + costheta,2)*(-1 + 4*costheta))/8.;
}

double Wignerd::d43p2( const double costheta ){
    return -(std::sqrt(3.5)*std::pow(1 + costheta,2)*(-1 + 2*costheta)*std::sqrt(1 - std::pow(costheta,2)))/4.;
}

double Wignerd::d43p3( const double costheta ){
    return (std::pow(1 + costheta,3)*(-3 + 4*costheta))/8.;
}

// 2

double Wignerd::d42m2( const double costheta ){
    return (std::pow(-1 + costheta,2)*(1 + 7*costheta*(1 + costheta)))/4.;
}

double Wignerd::d42m1( const double costheta ){
    return -(std::sqrt(2 - 2*std::pow(costheta,2))*(-1 + costheta*(8 + 7*(1 - 2*costheta)*costheta)))/8.;
}

double Wignerd::d420( const double costheta ){
    return -(std::sqrt(2.5)*(1 - 8*std::pow(costheta,2) + 7*std::pow(costheta,4)))/4.;
}

double Wignerd::d42p1( const double costheta ){
    return -((1 + costheta)*std::sqrt(2 - 2*std::pow(costheta,2))*(-1 + 7*costheta*(-1 + 2*costheta)))/8.;
}

double Wignerd::d42p2( const double costheta ){
    return (std::pow(1 + costheta,2)*(1 + 7*(-1 + costheta)*costheta))/4.;
}

// 1

double Wignerd::d41m1( const double costheta ){
    return (-3 + costheta*(-3 + costheta*(27 + 7*(1 - 4*costheta)*costheta)))/8.;
}

double Wignerd::d410( const double costheta ){
    return -(costheta*std::sqrt(5 - 5*std::pow(costheta,2))*(-3 + 7*std::pow(costheta,2)))/4.;
}

double Wignerd::d41p1( const double costheta ){
    return ((1 + costheta)*(3 + costheta*(-6 + 7*costheta*(-3 + 4*costheta))))/8.;
}

double Wignerd::d400( const double costheta ){
    return (3 - 30*std::pow(costheta,2) + 35*std::pow(costheta,4))/8.;
}





double Wignerd::generic( const double beta, const int j, int mp, int m ){

    if ( j - abs(mp) < 0 ||
         j - abs(m)  < 0 ) {
        return 0;
    }
    
    int smin = std::max( 0, m - mp );
    int smax = std::min( j + m , j - mp );
    
    double func = 0;
    
    for ( int s = smin; s <= smax; ++s ){
        double numer = std::pow( -1., mp - m + s );
        
        double denom = ( TMath::Factorial( j + m - s )*
                         TMath::Factorial( s )*
                         TMath::Factorial( mp - m + s )*
                         TMath::Factorial( j - mp - s ) );
    
        
        func +=  ( (numer/denom)*
                    std::pow( cos( 0.5*beta ), 2*j + m - mp - 2*s )*
                    std::pow( sin( 0.5*beta ), mp - m + 2*s ) );
    }
    
    
    return std::sqrt( TMath::Factorial( j + mp )*
                      TMath::Factorial( j - mp )*
                      TMath::Factorial( j + m )*
                      TMath::Factorial( j - m ) )*func;
}

TGraph* plot( int j,  int mp,  int m ){
    
    const double step = 0.01;
    
    TGraph* gr_fixedfn = new TGraph();
    TGraph* gr_generic = new TGraph();
    
    double costheta = -1.0 + step;
    
    while ( costheta < 1.0 ){
        gr_fixedfn->SetPoint( gr_fixedfn->GetN(), costheta, Wignerd::djmpm(costheta, j, mp, m ) );
        gr_generic->SetPoint( gr_generic->GetN(), costheta, Wignerd::generic( acos(costheta), j, mp, m ) );
        costheta += step;
    }
    
    gr_fixedfn->SetMaximum(+1.1);
    gr_fixedfn->SetMinimum(-1.1);
    
    gr_fixedfn->Draw("apl");
    gr_generic->Draw("l+");
    gr_generic->SetLineColor( kRed );
    
    gr_fixedfn->GetXaxis()->SetLimits(-1,1);
    gr_fixedfn->GetXaxis()->SetTitle("cos  #theta");
    
    char label[20];
    sprintf(label,"d^{%i}_{%i,%i}",j,mp,m);
    
    TLatex* text = new TLatex(0.25,0.15,label);
    text->SetNDC();
    text->Draw();
    
    return gr_fixedfn;
}
