#include "HelperFunctions.h"

#include "Math/Vector3D.h"
#include "Math/Boost.h"
#include<math.h>


/*
 * Helper functions that are used to derive the angles 
 * from 4-vectors or the Z K angles from the conventional
 * K pi psi ones. 
 */ 

double Helper::BellePlaneAngle( const Helper::LorentzVector& vecA,
                                const Helper::LorentzVector& vecB,
                                const Helper::LorentzVector& vecC ) {
    
    const ROOT::Math::XYZVector dirA = ROOT::Math::XYZVector( vecA ).unit();
    const ROOT::Math::XYZVector momB = ROOT::Math::XYZVector( vecB );
    const ROOT::Math::XYZVector momC = ROOT::Math::XYZVector( vecC );
    
    
    const ROOT::Math::XYZVector aB = (momB - (momB.Dot( dirA ))*dirA).unit();
    const ROOT::Math::XYZVector aC = (momC - (momC.Dot( dirA ))*dirA).unit();
    
    double cosphi = aB.Dot( aC ) ;
    double sinphi = (dirA.Cross( aB )).Dot( aC );
    
    return  ( sinphi > 0.0 ? acos( cosphi ) : -acos( cosphi ) );
}

void Helper::BelleAngles( const int BID ,
                         const double mkpi,
                         const double coskpi,
                         const double cospsi,
                         const double phi,
                         double& mz,
                         double& cosz,
                         double& cospsiz,
                         double& phiz,
                         const double mB ,
                         const double mpsi ,
                         const double mkaon,
                         const double mpion,
                         const double mmuon ){
    
    // create 4 vectors from angles
    Helper::LorentzVector mother(0,0,0,mB);
 
    const double mommuon = Helper::momentum( mpsi, mmuon, mmuon );
    const double momkaon = Helper::momentum( mkpi, mkaon, mpion );
    const double mompsi  = Helper::momentum( mB  , mpsi , mkpi  );
    
    const double emuon = sqrt( mmuon*mmuon + mommuon*mommuon );
    const double ekaon = sqrt( mkaon*mkaon + momkaon*momkaon );
    const double epion = sqrt( mpion*mpion + momkaon*momkaon );
    
    Helper::LorentzVector psi( 0, 0, -mompsi, sqrt( mompsi*mompsi + mpsi*mpsi ) );
    Helper::LorentzVector kpi( 0, 0,  mompsi, sqrt( mompsi*mompsi + mkpi*mkpi ) );
    
    const double sinpsi = sqrt( 1. - cospsi*cospsi );
    
    Helper::LorentzVector mupinpsi(  mommuon*sinpsi,0,  mommuon*cospsi, emuon );
    Helper::LorentzVector muninpsi( -mommuon*sinpsi,0, -mommuon*cospsi, emuon );
    
    const double sinkpi = sqrt( 1. - coskpi*coskpi );
    const double kaonpx = momkaon*sinkpi*cos(phi);
    const double kaonpy = momkaon*sinkpi*sin(phi);
    const double kaonpz = momkaon*coskpi;
    Helper::LorentzVector kaoninkpi( kaonpx, kaonpy, kaonpz , ekaon );
    Helper::LorentzVector pioninkpi(-kaonpx,-kaonpy,-kaonpz , epion );
    
    ROOT::Math::Boost boostToPsi( ROOT::Math::BoostZ( psi.Pz()/psi.E() ) );
    ROOT::Math::Boost boostToKpi( ROOT::Math::BoostZ( kpi.Pz()/kpi.E() ) );
    Helper::LorentzVector kaon = boostToKpi( kaoninkpi );
    Helper::LorentzVector pion = boostToKpi( pioninkpi );
    Helper::LorentzVector mup  = boostToPsi( mupinpsi );
    Helper::LorentzVector mun  = boostToPsi( muninpsi );

    TLorentzVector vkaon( kaon.Px(), kaon.Py(), kaon.Pz(), kaon.E() );
    TLorentzVector vpion( pion.Px(), pion.Py(), pion.Pz(), pion.E() );
    TLorentzVector vmup ( mup.Px(), mup.Py(), mup.Pz(), mup.E() );
    TLorentzVector vmun ( mun.Px(), mun.Py(), mun.Pz(), mun.E() );
    
    BelleAngles( BID, vmup, vmun, vkaon, vpion, mz, cosz, cospsiz, phiz );
    
    return ;
}

double Helper::cosangle( const Helper::LorentzVector& vecA,
                        const Helper::LorentzVector& vecB ){
    
    // mother --> (AB) C in mother rest frame
    const Helper::LorentzVector vecAB = vecA + vecB;
    
    ROOT::Math::Boost boostToRest( vecAB.BoostToCM() );
    
    const Helper::LorentzVector vecAinAB = boostToRest( vecA );
    const ROOT::Math::XYZVector dirAB = ROOT::Math::XYZVector( vecAB ).unit();
    const ROOT::Math::XYZVector dirA  = ROOT::Math::XYZVector( vecAinAB ).unit();
    
    return dirA.Dot( dirAB );
}

double Helper::cosangle( const Helper::LorentzVector& vecAB,
                        const Helper::LorentzVector& vecA,
                        const Helper::LorentzVector& vecC ){
    
    // mother --> (AB) C in mother rest frame
    ROOT::Math::Boost boostToRest( vecAB.BoostToCM() );
    
    const Helper::LorentzVector vecAinAB = boostToRest( vecA );
    const Helper::LorentzVector vecCinAB = boostToRest( vecC );
    
    const ROOT::Math::XYZVector dirA = ROOT::Math::XYZVector( vecAinAB ).unit();
    const ROOT::Math::XYZVector dirC = ROOT::Math::XYZVector( vecCinAB ).unit();
    
    return dirA.Dot( dirC );
}

void Helper::BelleAngles( const int BID,
                         const TLorentzVector& mup,
                         const TLorentzVector& mun,
                         const TLorentzVector& kaon,
                         const TLorentzVector& pion,
                         double& mz,
                         double& cosz,
                         double& cospsiz,
                         double& phiz ) {
    
  Helper::LorentzVector pmup ( mup.Px() , mup.Py() , mup.Pz() , mup.E() );
  Helper::LorentzVector pmun ( mun.Px() , mun.Py() , mun.Pz() , mun.E() );
  Helper::LorentzVector pkaon( kaon.Px(), kaon.Py(), kaon.Pz(), kaon.E() );
  Helper::LorentzVector ppion( pion.Px(), pion.Py(), pion.Pz(), pion.E() );
  
  const Helper::LorentzVector mother = ( pmup + pmun + pkaon + ppion );

  
  // Boost to B frame
  ROOT::Math::Boost boostToMother( mother.BoostToCM() );


  pmup  = boostToMother( pmup );
  pmun  = boostToMother( pmun );
  pkaon = boostToMother( pkaon );
  ppion = boostToMother( ppion );

  const Helper::LorentzVector psi = pmup  + pmun ;
  const Helper::LorentzVector pz  = psi   + ppion ;
  const Helper::LorentzVector kpi = pkaon + ppion ;
  
  mz = pz.M();
  
  ROOT::Math::Boost boostToZ( pz.BoostToCM() );
  
  // Boost to Z frame and calculate angle
  const Helper::LorentzVector pioninz = boostToZ( ppion );
  const Helper::LorentzVector kaoninz = boostToZ( pkaon );
  const Helper::LorentzVector mupinz  = boostToZ( pmup );
  const Helper::LorentzVector muninz  = boostToZ( pmun );
  const Helper::LorentzVector psiinz  = boostToZ( psi ) ;
  const Helper::LorentzVector kpiinz  = boostToZ( kpi );
  
  const ROOT::Math::XYZVector dirpioninz = ROOT::Math::XYZVector( pioninz ).unit();
  const ROOT::Math::XYZVector dirkaoninz = ROOT::Math::XYZVector( kaoninz ).unit();
  const ROOT::Math::XYZVector dirmupinz  = ROOT::Math::XYZVector( mupinz ).unit();
  const ROOT::Math::XYZVector dirmuninz  = ROOT::Math::XYZVector( muninz ).unit();
  
  cosz = -dirkaoninz.Dot( dirpioninz );
  
  // Boost to Z frame, then psi frame to compute angles
  ROOT::Math::Boost boostToPsi( psiinz.BoostToCM() );
  
  const Helper::LorentzVector pioninpsi = boostToPsi( pioninz );
  const Helper::LorentzVector mupinpsi  = boostToPsi( mupinz );
  
  const ROOT::Math::XYZVector dirpioninpsi = ROOT::Math::XYZVector( pioninpsi ).unit();
  const ROOT::Math::XYZVector dirmupinpsi  = ROOT::Math::XYZVector( mupinpsi ).unit();
  
  cospsiz = -dirpioninpsi.Dot( dirmupinpsi );
  
  
  // double check this definition
  const ROOT::Math::XYZVector ek = ( dirpioninz.Cross( dirkaoninz ) ).unit();
  const ROOT::Math::XYZVector el = ( dirmupinz.Cross( dirmuninz ) ).unit();
  
  const ROOT::Math::XYZVector ez = ROOT::Math::XYZVector(kpiinz).unit();
  
  double cosphi = ( ek.Dot( el ) );
  double sinphi = ( el.Cross( ek ) ).Dot( ez );
  
  phiz    = ( sinphi > 0.0 ? acos( cosphi ) : -acos( cosphi ) );
  
  if ( BID < 0 ){
    cospsiz = -cospsiz;
    
    if ( phiz > 0 ) {
      phiz =  TMath::Pi() - phiz;
    } else {
      phiz = -TMath::Pi() - phiz;
    }
  }
  return ;
}

void Helper::HelicityAngles( const int BID,
                            const TLorentzVector& mup,
                            const TLorentzVector& mun,
                            const TLorentzVector& kaon,
                            const TLorentzVector& pion,
                            double& mkpi,
                            double& coskpi,
                            double& cospsi,
                            double& phi ) {
    
  Helper::LorentzVector pmup ( mup.Px() , mup.Py() , mup.Pz() , mup.E() );
  Helper::LorentzVector pmun ( mun.Px() , mun.Py() , mun.Pz() , mun.E() );
  Helper::LorentzVector pkaon( kaon.Px(), kaon.Py(), kaon.Pz(), kaon.E() );
  Helper::LorentzVector ppion( pion.Px(), pion.Py(), pion.Pz(), pion.E() );
  
  const Helper::LorentzVector mother = ( pmup + pmun + pkaon + ppion );
  
  ROOT::Math::Boost boostToMother( mother.BoostToCM() );
  
  pmup  = boostToMother( pmup );
  pmun  = boostToMother( pmun );
  pkaon = boostToMother( pkaon );
  ppion = boostToMother( ppion );
  
  // Compute K pi mass
  const Helper::LorentzVector kpi = pkaon + ppion ;
  mkpi = kpi.M();
  
  // Compute coskpi
  coskpi = Helper::cosangle( pkaon, ppion );
  
  // Compute cospsi
  cospsi = Helper::cosangle( pmup, pmun );
  
  if ( BID < 0 ){
    cospsi = -cospsi;
  }
  
  const ROOT::Math::XYZVector dirmup  = ROOT::Math::XYZVector( pmup ).unit();
  const ROOT::Math::XYZVector dirmun  = ROOT::Math::XYZVector( pmun ).unit();
  const ROOT::Math::XYZVector dirkaon = ROOT::Math::XYZVector( pkaon ).unit();
  const ROOT::Math::XYZVector dirpion = ROOT::Math::XYZVector( pkaon ).unit();
  
  const ROOT::Math::XYZVector ek = ( dirkaon.Cross( dirpion ) ).unit();
  const ROOT::Math::XYZVector ez = ROOT::Math::XYZVector(kpi).unit();
  
      std::cout << "check" << std::endl;

  
  if ( BID > 0 ){
    // Compute phi
    const ROOT::Math::XYZVector el = ( dirmup.Cross( dirmun ) ).unit();
    
    double cosphi = ( ek.Dot( el ) );
    double sinphi = ( el.Cross( ek ) ).Dot( ez );
    
    phi    = ( sinphi > 0.0 ? acos( cosphi ) : -acos( cosphi ) );
    double testphi{0};
    testphi = TMath::ATan2(sinphi,cosphi);
    std::cout << "check2" << std::endl;
    if(testphi != phi){
	std::cout << "ATan2 method is: " << testphi << " Normal method is: " << phi << std::endl;
    }

  }
  else {
    
    // Compute phi
    const ROOT::Math::XYZVector el = ( dirmun.Cross( dirmup ) ).unit();
    
    double cosphi = ( ek.Dot( el ) );
    double sinphi = ( el.Cross( ek ) ).Dot( ez );
    
    phi    = ( sinphi > 0.0 ? acos( cosphi ) : -acos( cosphi ) );
    phi    = -phi;
  }
  
  return ;
}


double Helper::momentum( const double m, const double m1, const double m2 )  {
  double msq = m*m;
  
  double psq = 0.25*(msq - (m1 + m2)*(m1 + m2))*(msq - (m1 - m2)*(m1 - m2))/msq;

 if(psq < 0 ){
	std::cout << "ERROR: psq variable is negative in Helper::momentum. Probably using the wrong kinematical variables in Event. Please check." << std::endl;
  	std::cout << "DEBUG HELPER: msq, m1, m2 " << msq << " " << m1 << " " << m2 << std::endl;
	return 0;
 }

 else{
  return std::sqrt( psq );
 }
}

void Helper::convertToTransversityBasis( const std::complex<double>& H0, 
					 const std::complex<double>& Hp, 
					 const std::complex<double>& Hm, 
					 std::complex<double>& Azero, 
					 std::complex<double>& Apara, 
					 std::complex<double>& Aperp ) {

  Azero = H0;
  Apara = ( Hp + Hm )/std::sqrt(2.);
  Aperp = ( Hp - Hm )/std::sqrt(2.);

  return; 
}

void Helper::convertToHelicityBasis(  const std::complex<double>& Azero, 
				      const std::complex<double>& Apara,
				      const std::complex<double>& Aperp, 
				      std::complex<double>& H0, 
				      std::complex<double>& Hp, 
				      std::complex<double>& Hm ){ 
  
  H0 = Azero;
  Hp = ( Apara + Aperp )/std::sqrt(2.);
  Hm = ( Apara - Aperp )/std::sqrt(2.);
  
  return;
}

void Helper::cpBasisSwap( const std::complex<double>& H0,
			  const std::complex<double>& Hp,
			  const std::complex<double>& Hm,
			  std::complex<double>& H0CP, 
			  std::complex<double>& HpCP, 
			  std::complex<double>& HmCP, 
			  const int spin ){

  std::complex<double> Azero, Apara, Aperp;
  
  convertToTransversityBasis( H0 ,Hp, Hm, Azero, Apara, Aperp );

  if(spin % 2 == 0 ){
    Azero = -Azero;
    Apara = -Apara;
  }
  else{
    Aperp = -Aperp;
  }
  
  convertToHelicityBasis( Azero, Apara, Aperp, H0CP, HpCP, HmCP );

  return;
}

void Helper::transFracScale( std::complex<double> &Azero,
			   std::complex<double> &Apara,
			   std::complex<double> &Aperp, 
			   const double fpara,
			   const double fperp ){
  
  const double fzero = 1. - fpara - fperp; 

  Azero = std::polar(sqrt(fzero),std::arg(Azero));
  Apara = std::polar(sqrt(fpara),std::arg(Apara));
  Aperp = std::polar(sqrt(fperp),std::arg(Aperp));

  return;
}

void Helper::fitFracScale(double &scalefitfrac,
			  std::complex<double> &Azero,
			  std::complex<double> &Apara,
			  std::complex<double> &Aperp,
			  const double &fitfrac){

	Azero = Azero*sqrt(fitfrac / scalefitfrac);
	Apara = Apara*sqrt(fitfrac / scalefitfrac);
	Aperp = Aperp*sqrt(fitfrac / scalefitfrac);

	return;
}
