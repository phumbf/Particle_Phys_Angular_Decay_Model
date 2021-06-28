#include "Event.h"
#include <cmath>

#include "HelperFunctions.h"
#include "Variables.h"

/*
 * Simple class that knows about the mass and angular 
 * properties of a given event.
 * 
 */

// Static ID that is incremented every time a new event is created
unsigned int Event::sID = 0;

double Event::mB    = 5.279;

double Event::mBs    = 5.367;

double Event::mPsi  = 3.097;

double Event::mKaon = 0.494;

double Event::mPion = 0.140;

double Event::D = 0.01;


Event::Event( const Event& other ) :
  mkpi_   ( other.mkpi_ ),
  mz_     ( other.mz_ ),
  coskpi_ ( other.coskpi_ ),
  cospsi_ ( other.cospsi_ ),
  phi_    ( other.phi_ ),
  cosz_   ( other.cosz_ ),
  cospsiz_( other.cospsiz_ ),
  pB_     ( other.pB_ ),
  pR_     ( other.pR_ ),
  pZ_     ( other.pZ_ ),
  ID_     ( other.ID_ ) {};


Event::Event( const double mkpi,
              const double mz,
              const double coskpi,
              const double cospsi,
              const double phi,
              const double cosz,
              const double cospsiz,
              const double phiz ) :
mkpi_   ( mkpi ),
mz_     ( mz ),
coskpi_ ( coskpi ),
cospsi_ ( cospsi ),
phi_    ( phi ),
cosz_   ( cosz ),
cospsiz_( cospsiz ),
phiz_   ( phiz )
{
  ID_ = sID++;
    
    pB_   = Helper::momentum( mB, mkpi_, mPsi );
    pR_   = Helper::momentum( mkpi_, mKaon, mPion );
    pZ_   = Helper::momentum( mz_  , mPsi, mPion );
    pBz_  = Helper::momentum( mB, mz_  , mPion );
} ;

Event::Event( const int BID,
              const double mkpi,
              const double coskpi,
              const double cospsi,
              const double phi ) :
  mkpi_   ( mkpi ),
  coskpi_ ( coskpi ),
  cospsi_ ( cospsi ),
  phi_    ( phi )
{
  /*
   * Constructor that computes the Z properties 
   * 
   */
  
  ID_ = sID++;
  
  // compute z angles
  Helper::BelleAngles( BID, mkpi, coskpi, cospsi, phi,
		       mz_, cosz_, cospsiz_, phiz_ ,
		       mB, mPsi, mKaon, mPion, Variables::mMuon );
  
  
  pB_   = Helper::momentum( mB, mkpi_, mPsi );
  pR_   = Helper::momentum( mkpi_, mKaon, mPion );
  pZ_   = Helper::momentum( mz_  , mPsi, mPion );
  pBz_  = Helper::momentum( mB, mz_  , mPion );
};

Event::Event( const int BID,
              const TLorentzVector& muplus,
              const TLorentzVector& muminus,
              const TLorentzVector& kaon,
              const TLorentzVector& pion )
{
  /*
   * Constructor from particle Lorentz vectors
   * 
   */
  
  ID_ = sID++;
  
  // compute all angles from 4-vectors
  Helper::HelicityAngles( BID, muplus, muminus, kaon, pion,
			  mkpi_, coskpi_, cospsi_, phi_ );
  
  Helper::BelleAngles( BID, muplus, muminus, kaon, pion,
		       mz_, cosz_, cospsiz_, phi_ );
  
  
  pB_   = Helper::momentum( mB, mkpi_, mPsi );
  pR_   = Helper::momentum( mkpi_, mKaon, mPion );
  pZ_   = Helper::momentum( mz_  , mPsi, mPion );
  pBz_  = Helper::momentum( mB, mz_  , mPion );
} ;

std::ostream& operator<<( std::ostream& stream, const Event& evt ){
  // streamer operator
  stream << " ==> INFO (Event): Event \n " ;
  stream << "     m(K#pi)                = " << evt.mkpi_    << "\n" ;
  stream << "     cos(#theta_{#psi})     = " << evt.cospsi_  << "\n" ;
  stream << "     cos(#theta_{K#pi})     = " << evt.coskpi_  << "\n" ;
  stream << "     #phi                   = " << evt.phi_     << "\n" ;
  stream << "     m(#psi#pi)             = " << evt.mz_      << "\n" ;
  stream << "     cos(#theta_{Z})        = " << evt.cosz_    << "\n" ;
  stream << "     cos(#theta_{#psi}^{z}) = " << evt.cospsiz_ << "\n" ;
  stream << "     #phi_{z}               = " << evt.phiz_    << "\n" ;
  return stream;
}
