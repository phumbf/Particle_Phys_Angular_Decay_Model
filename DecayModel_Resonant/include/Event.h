#ifndef RESONANCEEVENT_H
#define RESONANCEEVENT_H

#include "TLorentzVector.h"

#include <iostream>

/*
 * Simple class that knows about the mass and angular 
 * properties of a given event.
 * 
 */

class Event {
    
public:
    
    Event( const double mkpi,
           const double mz,
           const double coskpi,
           const double cospsi,
           const double phi,
           const double cosz,
           const double cospsiz,
           const double phiz );
    
    Event( const int BID,
           const double mkpi,
           const double coskpi,
           const double cospsi,
           const double phi );

    Event( const int BID,
           const TLorentzVector& muplus,
           const TLorentzVector& muminus,
           const TLorentzVector& kaon,
           const TLorentzVector& pion ) ;
    
    Event( const Event& other );
    
    ~Event() {};
    
    friend std::ostream& operator<<( std::ostream& stream, const Event& evt );
    
    // standard variables
    double mkpi_;
    double mz_;
    
    double coskpi_;
    double cospsi_;
    double phi_;
    
    double cosz_;
    double cospsiz_;
    double phiz_;
    double deltaphi_;
    
    // particle momentum for lineshape
    double pB_ ;
    double pR_ ;
    double pZ_ ;
    double pBz_;
    
    // efficiency (currently unused)
    double eff_ ;
    
    // double angular( const unsigned int j , const int lambdapsi, const int lambda ) const ;
    
    unsigned int ID_;

    static double mB;
    static double mBs;
    static double mPsi;
    static double mKaon;
    static double mPion;
    static double D;

protected:
    static unsigned int sID;
    
    
private:
    
};


#endif
