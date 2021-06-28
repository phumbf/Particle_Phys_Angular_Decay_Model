
#include "FFcalculation.h"
#include "Parameters.h"
#include "Kinematic.h" 
#include "Parameters.h" 

#include <iostream>
#include <cmath>
#include <cassert>

#include "TVectorD.h"

//Classes to determine the different form factors as a function of kinematic terms. The basis are taken from the 
//different papers seen in the form_factors folder
//

//=============================================================================
// BSZ15 Constructors/setting
//=============================================================================


BSZ15FFCalculation::BSZ15FFCalculation( ) { } 

BSZ15FFCalculation::BSZ15FFCalculation( const BSZ15FFCalculation& other ) : 
  params_( other.params_ ) {} 

BSZ15FFCalculation::BSZ15FFCalculation( TVectorD& params ) {
  setParams( params );
  print();
}

BSZ15FFCalculation::BSZ15FFCalculation( const std::vector<double>& params ) { 
  setParams( params );
  print();
}


BSZ15FFCalculation::~BSZ15FFCalculation() {}


void BSZ15FFCalculation::setParams( const TVectorD& params ) {
    
  assert( params.GetNrows() == 21 );
  
  for ( int i = 0; i < 21; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void BSZ15FFCalculation::setParams( const std::vector< double >& params ) {
    
  assert( params.size() == 21 );
  
  for ( int i = 0; i < 21; ++i ){
    params_[i] = params[i];
  }
  
  return;
}


void BSZ15FFCalculation::print() const {

  std::cout << "Printing BSZ15 parameters:" << std::endl;
  
  for(int i=0; i < 21; i++){
    std::cout << params_[i] << std::endl;
  }
  
  return;
}

TVectorD BSZ15FFCalculation::getParams() const{
  
  TVectorD vec(21);
  
  for(int i=0; i < 21; i++){
    vec[i] = params_[i];
  }
  
  return vec;
}

//=============================================================================
// BSZ15 Form-factors
//=============================================================================

double BSZ15FFCalculation::zseries( const double qsq,
				    const double a0,
				    const double a1,
				    const double a2 ) const {
  
  const double zdiff = ( Kinematic::z( qsq, m_mass, m_motherMass ) - 
			 Kinematic::z( 0  , m_mass, m_motherMass ) );
  
  double result = a0;
  result += a1*( zdiff );
  result += a2*std::pow( zdiff, 2 ); 
  
  return result;
}




double BSZ15FFCalculation::pole( const double qsq, const double mass ) const {
  return 1./(1. - qsq/(mass*mass));
}



double BSZ15FFCalculation::A0( const double qsq ) const {
  return zseries( qsq, params_[0], params_[1], params_[2] )*pole( qsq, m_mB0m );
}

double BSZ15FFCalculation::A1( const double qsq ) const {
  return zseries( qsq, params_[3], params_[4], params_[5] )*pole( qsq, m_mB1p );
}

// A12?
double BSZ15FFCalculation::A12( const double qsq ) const {
  return zseries( qsq, params_[6], params_[7], params_[8] )*pole( qsq, m_mB1p );
}

double BSZ15FFCalculation::V( const double qsq )  const {
  return zseries( qsq, params_[9],params_[10],params_[11] )*pole( qsq, m_mB1m );
}

double BSZ15FFCalculation::T1( const double qsq ) const {
  return zseries( qsq, params_[12],params_[13],params_[14] )*pole( qsq, m_mB1m );
}

double BSZ15FFCalculation::T2( const double qsq ) const {
  return zseries( qsq, params_[15],params_[16],params_[17] )*pole( qsq, m_mB1p );
}

// T23?
double BSZ15FFCalculation::T23( const double qsq ) const {
  return zseries( qsq, params_[18],params_[19],params_[20] )*pole( qsq, m_mB1p );
}

double BSZ15FFCalculation::A2( const double qsq ) const {
  double result = 0;

  double mBsq = std::pow( m_motherMass, 2 );
  double mVsq = std::pow( m_mass, 2 );
  
  double mSum = m_motherMass + m_mass ;
  
  result -= 16.*(m_motherMass)*mVsq*( mSum )*A12( qsq );
  result += std::pow( mSum, 2 )*( mBsq - mVsq - qsq )*A1( qsq );
  
  result /= Kinematic::lambda( qsq, m_mass, m_motherMass );
  
  return result;
}

double BSZ15FFCalculation::T3( const double qsq ) const {
  double result = 0;
  
  double mBsq = std::pow( m_motherMass, 2 );
  double mVsq = std::pow( m_mass, 2 );
  
  result -= 8.*m_motherMass*mVsq*( m_motherMass - m_mass )*T23( qsq );
  result += ( mBsq - mVsq )*( mBsq + 3.*mVsq - qsq )*T2( qsq );
  
  result /= Kinematic::lambda( qsq, m_mass, m_motherMass ); 
  
  return result;
}

double BSZ15FFCalculation::T3tilde( const double qsq) const{
  
  double mBsq = std::pow( m_motherMass, 2 );
  double mVsq = std::pow( m_mass, 2 );

  return T2(qsq) + ( T3(qsq) * ( qsq / ( mBsq - mVsq ) ) );

}


//=============================================================================
// BZ04 Constructors/setting 
//=============================================================================



BZ04FFCalculation::BZ04FFCalculation( ){ } 

BZ04FFCalculation::BZ04FFCalculation( const BZ04FFCalculation& other ) : 
  params_( other.params_ ) {} 

BZ04FFCalculation::BZ04FFCalculation( TVectorD& params ) {
  setParams( params );
  print();
}

BZ04FFCalculation::BZ04FFCalculation( const std::vector<double>& params ) { 


  setParams( params );
  print();
}


BZ04FFCalculation::~BZ04FFCalculation() {}


void BZ04FFCalculation::setParams( const TVectorD& params ) {
    
  assert( params.GetNrows() == 19 );
  
  for ( int i = 0; i < 19; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void BZ04FFCalculation::setParams( const std::vector< double >& params ) {
    
  assert( params.size() == 19 );
  
  for ( int i = 0; i < 19; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void BZ04FFCalculation::print() const {
  
  std::cout << "Printing BZ04 parameters" << std::endl;
  
  for ( int i = 0; i < 19; ++i ){
    std::cout << params_[i] << std::endl;
  }

  return;
}

TVectorD BZ04FFCalculation::getParams() const{
  
  TVectorD vec(19);
  
  for(int i=0; i < 19; i++){
    vec[i] = params_[i];
  }
  
  return vec;
  
}

//=============================================================================
// BZ04 Form-factors
//=============================================================================

double BZ04FFCalculation::V( const double qsq ) const { 
  return ( pole( qsq, params_[0], std::pow( m_mB1m, 2 ) ) + 
	   pole( qsq, params_[1], params_[2] ) );
}

double BZ04FFCalculation::A0( const double qsq ) const { 
  return  ( pole( qsq, params_[3], std::pow( m_mB0m, 2 ) ) + 
	    pole( qsq, params_[4], params_[5] ) );
}

double BZ04FFCalculation::A1( const double qsq ) const { 
  return pole( qsq, params_[6], params_[7] );
}

double BZ04FFCalculation::A2( const double qsq ) const { 
  return ( pole( qsq  , params_[8], params_[10] ) +
	   sqpole( qsq, params_[9], params_[10] ) );
}

double BZ04FFCalculation::T1( const double qsq ) const { 
  return ( pole( qsq, params_[11], std::pow( m_mB1m, 2 ) ) + 
	   pole( qsq, params_[12], params_[13] ) );
}

double BZ04FFCalculation::T2( const double qsq ) const { 
  return pole( qsq, params_[14], params_[15] );
}

double BZ04FFCalculation::T3tilde( const double qsq ) const { 
  return ( pole( qsq  , params_[16], params_[18] ) + 
	   sqpole( qsq, params_[17], params_[18] ) );
}

double BZ04FFCalculation::T3( const double qsq ) const { 
  
  const double t3tilde = T3tilde( qsq ); 
  const double t2      = T2( qsq );
  
  double t3 = 0;
  
  if ( std::abs( qsq ) > 1e-10 ){
    t3 = ( std::pow(m_motherMass,2) - 
	   std::pow(m_mass,2) )*(t3tilde - t2)/qsq;
  }
  
  return t3; 
}

double BZ04FFCalculation::A12( const double qsq ) const{

  double mbmk = m_motherMass + m_mass;
  
  double mkmkqsq = m_motherMass*m_motherMass -
  		 m_mass*m_mass - qsq;
  
  double mbmkmk = 16.0 * m_motherMass * m_mass*m_mass;
  
  
  double numer = mbmk*mbmk*mkmkqsq*A1(qsq) -
  	       Kinematic::lambda(qsq, m_mass, m_motherMass)*A2(qsq);
  
  double denom = mbmkmk * mbmk;

	return numer / denom;
}

double BZ04FFCalculation::T23( const double qsq ) const{

  double mbmk = m_motherMass*m_motherMass - m_mass*m_mass;
  
  double mbmkqsq = m_motherMass*m_motherMass +
  	         3.0*m_mass*m_mass - qsq;
  
  double denom = 8.0*m_motherMass*m_mass*m_mass*
                 	   ( m_motherMass - m_mass );
  
  double numer = mbmk * mbmkqsq * T2(qsq) -
  	       Kinematic::lambda(qsq, m_mass, m_motherMass) * T3(qsq);
  
  return numer / denom ;
}

double BZ04FFCalculation::pole( const double qsq, const double r, const double msq ) const { 
  return r/( 1. - ( qsq/msq ) );
}


double BZ04FFCalculation::sqpole( const double qsq, const double r, const double msq ) const { 
  return r/std::pow( 1. - ( qsq/msq ), 2 );
}

//Scalar form factors add
AAS07FFCalculation::AAS07FFCalculation( ) { } 

AAS07FFCalculation::AAS07FFCalculation( const AAS07FFCalculation& other ):
  params_( other.params_ ) {} 

AAS07FFCalculation::AAS07FFCalculation( TVectorD& params ) {
  setParams( params );
  print();
}

AAS07FFCalculation::AAS07FFCalculation( const std::vector<double>& params ) { 
  setParams( params );
  print();
}

AAS07FFCalculation::~AAS07FFCalculation() {}

void AAS07FFCalculation::setParams( const TVectorD& params ) {
    
  assert( params.GetNrows() == 9 );
  
  for ( int i = 0; i < 9; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void AAS07FFCalculation::setParams( const std::vector< double >& params ) {
    
  assert( params.size() == 9 );
  
  for ( int i = 0; i < 9; ++i ){
    params_[i] = params[i];
  }
  
  return;
}


void AAS07FFCalculation::print() const {

  std::cout << "Printing Scalar parameters:" << std::endl;
  
  for(int i=0; i < 9; i++){
    std::cout << params_[i] << std::endl;
  }
  
  return;
}

TVectorD AAS07FFCalculation::getParams() const{
  
  TVectorD vec(9);
  
  for(int i=0; i < 9; i++){
    vec[i] = params_[i];
  }
  
  return vec;
}

double calcFabFF( const double qsq,
		  const double mB , 
		  const double f0,
		  const double a,
		  const double b ) {

  const double shat  = qsq/( mB*mB );

  return  f0 / ( 1.0 - a*shat + b*shat*shat  ); 
}

double calcFabFFPV( const double qsq,
		    const double m,
		    const double f0,
		    const double a,
		    const double b ){
  
	const double shat  = qsq/( m*m );

	return f0 / ( (1 - shat )*(1 - a*shat + b*shat*shat) );

}

double calcCCCFF( const double qsq,
		  const double mB,
		  const double f0,
		  const double c1,
		  const double c2,
		  const double c3){

	const double shat = qsq/ ( mB*mB );

	return f0*std::exp(c1*shat + c2*shat*shat + c3*shat*shat*shat);
}

//--------------------------------------------
//ABHH99
//--------------------------------------------
ABHH99FFCalculation::ABHH99FFCalculation( ){ } 

ABHH99FFCalculation::ABHH99FFCalculation( const ABHH99FFCalculation& other ) : 
  params_( other.params_ ) {} 

ABHH99FFCalculation::ABHH99FFCalculation( TVectorD& params ) {
  setParams( params );
  print();
}

ABHH99FFCalculation::ABHH99FFCalculation( const std::vector<double>& params ) { 
  setParams( params );
  print();
}

ABHH99FFCalculation::~ABHH99FFCalculation() {}

void ABHH99FFCalculation::setParams( const TVectorD& params ) {
    
  assert( params.GetNrows() == 28 );
  
  for ( int i = 0; i < 28; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void ABHH99FFCalculation::setParams( const std::vector< double >& params ) {
    
  assert( params.size() == 28 );
  
  for ( int i = 0; i < 28; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void ABHH99FFCalculation::print() const {

  std::cout << "Printing ABHH99 parameters:" << std::endl;
  
  for(int i=0; i < 28; i++){
    std::cout << params_[i] << std::endl;
  }
  
  return;
}

TVectorD ABHH99FFCalculation::getParams() const{
  
  TVectorD vec(28);
  
  for(int i=0; i < 28; i++){
    vec[i] = params_[i];
  }
  
  return vec;
}

double ABHH99FFCalculation::A1( const double qsq ) const { 
  return calcCCCFF(qsq, m_motherMass, params_[0], params_[1], params_[2], params_[3]);
}

double ABHH99FFCalculation::A2( const double qsq ) const { 
  return calcCCCFF(qsq, m_motherMass, params_[4], params_[5], params_[6], params_[7]);
}

double ABHH99FFCalculation::A0( const double qsq ) const { 
  return calcCCCFF(qsq, m_motherMass, params_[8], params_[9], params_[10], params_[11]);
}

double ABHH99FFCalculation::V( const double qsq ) const { 
  return calcCCCFF(qsq, m_motherMass, params_[12], params_[13], params_[14], params_[15]);
}

double ABHH99FFCalculation::T1( const double qsq ) const { 
  return calcCCCFF(qsq, m_motherMass, params_[16], params_[17], params_[18], params_[19]);
}

double ABHH99FFCalculation::T2( const double qsq ) const { 
  return calcCCCFF(qsq, m_motherMass, params_[20], params_[21], params_[22], params_[23]);
}

double ABHH99FFCalculation::T3( const double qsq ) const { 
  return calcCCCFF(qsq, m_motherMass, params_[24], params_[25], params_[26], params_[27]);
}

double ABHH99FFCalculation::T3tilde( const double qsq) const{

  double mBsq = std::pow( m_motherMass, 2 );
  double mVsq = std::pow( m_mass, 2 );

  return T2(qsq) + ( T3(qsq) * ( qsq / ( mBsq - mVsq ) ) );

}

double ABHH99FFCalculation::A12( const double qsq ) const{

  double mbmk = m_motherMass + m_mass;
 
  double mkmkqsq = m_motherMass*m_motherMass -
               	  m_mass*m_mass - qsq;

  double mbmkmk = 16.0 * m_motherMass * m_mass*m_mass;


  double numer = mbmk*mbmk*mkmkqsq*A1(qsq) -
	  Kinematic::lambda(qsq, m_motherMass, m_mass)*A2(qsq);

  double denom = mbmkmk * mbmk;

	return numer / denom;
}

double ABHH99FFCalculation::T23( const double qsq ) const{

  double mbmk = m_motherMass*m_motherMass - m_mass*m_mass;
  
  double mbmkqsq = m_motherMass*m_motherMass +
        	3.0*m_mass*m_mass - qsq;
  
  double denom = 8.0*m_motherMass*m_mass*m_mass*
  	( m_motherMass - m_mass );
  
  double numer = mbmk * mbmkqsq * T2(qsq) -
		Kinematic::lambda(qsq, m_motherMass, m_mass) * T3(qsq);

  return numer / denom ;
}
//--------------------------------------------
//MS00
//--------------------------------------------
MS00FFCalculation::MS00FFCalculation( ){ } 

MS00FFCalculation::MS00FFCalculation( const MS00FFCalculation& other ) : 
  params_( other.params_ ) {} 

MS00FFCalculation::MS00FFCalculation( TVectorD& params ) {
  setParams( params );
  print();
}

MS00FFCalculation::MS00FFCalculation( const std::vector<double>& params ) { 
  setParams( params );
  print();
}

MS00FFCalculation::~MS00FFCalculation() {}

void MS00FFCalculation::setParams( const TVectorD& params ) {
    
  assert( params.GetNrows() == 21 );
  
  for ( int i = 0; i < 21; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void MS00FFCalculation::setParams( const std::vector< double >& params ) {
    
  assert( params.size() == 21 );
  
  for ( int i = 0; i < 21; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void MS00FFCalculation::print() const {

  std::cout << "Printing MS00 parameters:" << std::endl;
  
  for(int i=0; i < 21; i++){
    std::cout << params_[i] << std::endl;
  }
  
  return;
}

TVectorD MS00FFCalculation::getParams() const{
  
  TVectorD vec(21);
  
  for(int i=0; i < 21; i++){
    vec[i] = params_[i];
  }
  
  return vec;
}

double MS00FFCalculation::A1( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass,  params_[0], params_[1], params_[2]);
}

double MS00FFCalculation::A2( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[3],params_[4],params_[5]);
}

double MS00FFCalculation::A0( const double qsq ) const { 
  return calcFabFFPV(qsq, m_motherMass, params_[6], params_[7], params_[8]);
}

double MS00FFCalculation::V( const double qsq ) const { 
  return calcFabFFPV(qsq, m_mBst, params_[9],params_[10], params_[11]);
}

double MS00FFCalculation::T1( const double qsq ) const { 
  return calcFabFF(qsq, m_mBst, params_[12], params_[13], params_[14]);
}

double MS00FFCalculation::T2( const double qsq ) const { 
  return calcFabFF(qsq,m_motherMass, params_[15], params_[16], params_[17]);
}

double MS00FFCalculation::T3( const double qsq ) const { 
  return calcFabFF(qsq,m_motherMass, params_[18], params_[19], params_[20]);
}

double MS00FFCalculation::T3tilde( const double qsq) const{

  double mBsq = std::pow( m_motherMass, 2 );
  double mVsq = std::pow( m_mass, 2 );

  return T2(qsq) + ( T3(qsq) * ( qsq / ( mBsq - mVsq ) ) );

}

double MS00FFCalculation::A12( const double qsq ) const{

  double mbmk = m_motherMass + m_mass;
 
  double mkmkqsq = m_motherMass*m_motherMass -
               	  m_mass*m_mass - qsq;

  double mbmkmk = 16.0 * m_motherMass * m_mass*m_mass;


  double numer = mbmk*mbmk*mkmkqsq*A1(qsq) -
	  Kinematic::lambda(qsq, m_motherMass, m_mass)*A2(qsq);

  double denom = mbmkmk * mbmk;

	return numer / denom;
}

double MS00FFCalculation::T23( const double qsq ) const{

  double mbmk = m_motherMass*m_motherMass - m_mass*m_mass;
  
  double mbmkqsq = m_motherMass*m_motherMass +
        	3.0*m_mass*m_mass - qsq;
  
  double denom = 8.0*m_motherMass*m_mass*m_mass*
  	( m_motherMass - m_mass );
  
  double numer = mbmk * mbmkqsq * T2(qsq) -
		Kinematic::lambda(qsq, m_motherMass, m_mass) * T3(qsq);

  return numer / denom ;
}
//---------------------------------------
//ESS18
//---------------------------------------
ESS18FFCalculation::ESS18FFCalculation( ){ } 

ESS18FFCalculation::ESS18FFCalculation( const ESS18FFCalculation& other ) : 
  params_( other.params_ ) {} 

ESS18FFCalculation::ESS18FFCalculation( TVectorD& params ) {
  setParams( params );
  print();
}

ESS18FFCalculation::ESS18FFCalculation( const std::vector<double>& params ) { 
  setParams( params );
  print();
}

ESS18FFCalculation::~ESS18FFCalculation() {}

void ESS18FFCalculation::setParams( const TVectorD& params ) {
    
  assert( params.GetNrows() == 21 );
  
  for ( int i = 0; i < 21; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void ESS18FFCalculation::setParams( const std::vector< double >& params ) {
    
  assert( params.size() == 21 );
  
  for ( int i = 0; i < 21; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void ESS18FFCalculation::print() const {

  std::cout << "Printing Scalar parameters:" << std::endl;
  
  for(int i=0; i < 21; i++){
    std::cout << params_[i] << std::endl;
  }
  
  return;
}

TVectorD ESS18FFCalculation::getParams() const{
  
  TVectorD vec(21);
  
  for(int i=0; i < 21; i++){
    vec[i] = params_[i];
  }
  
  return vec;
}

double ESS18FFCalculation::A1( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass,  params_[0], params_[1], params_[2]);
}

double ESS18FFCalculation::A2( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[3],params_[4],params_[5]);
}

double ESS18FFCalculation::A0( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[6], params_[7], params_[8]);
}

double ESS18FFCalculation::V( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[9],params_[10], params_[11]);
}

double ESS18FFCalculation::T1( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[12], params_[13], params_[14]);
}

double ESS18FFCalculation::T2( const double qsq ) const { 
  return calcFabFF(qsq,m_motherMass, params_[15], params_[16], params_[17]);
}

double ESS18FFCalculation::T3( const double qsq ) const { 
  return calcFabFF(qsq,m_motherMass, params_[18], params_[19], params_[20]);
}

double ESS18FFCalculation::T3tilde( const double qsq) const{

  double mBsq = std::pow( m_motherMass, 2 );
  double mVsq = std::pow( m_mass, 2 );

  return T2(qsq) + ( T3(qsq) * ( qsq / ( mBsq - mVsq ) ) );

}

double ESS18FFCalculation::A12( const double qsq ) const{

  double mbmk = m_motherMass + m_mass;
 
  double mkmkqsq = m_motherMass*m_motherMass -
               	  m_mass*m_mass - qsq;

  double mbmkmk = 16.0 * m_motherMass * m_mass*m_mass;


  double numer = mbmk*mbmk*mkmkqsq*A1(qsq) -
	  Kinematic::lambda(qsq, m_motherMass, m_mass)*A2(qsq);

  double denom = mbmkmk * mbmk;

	return numer / denom;
}

double ESS18FFCalculation::T23( const double qsq ) const{

  double mbmk = m_motherMass*m_motherMass - m_mass*m_mass;
  
  double mbmkqsq = m_motherMass*m_motherMass +
        	3.0*m_mass*m_mass - qsq;
  
  double denom = 8.0*m_motherMass*m_mass*m_mass*
  	( m_motherMass - m_mass );
  
  double numer = mbmk * mbmkqsq * T2(qsq) -
		Kinematic::lambda(qsq, m_motherMass, m_mass) * T3(qsq);

  return numer / denom ;
}
//---------------------------------------
//YANG10
//---------------------------------------
YANG10FFCalculation::YANG10FFCalculation( ){ } 

YANG10FFCalculation::YANG10FFCalculation( const YANG10FFCalculation& other ) : 
  params_( other.params_ ) {} 

YANG10FFCalculation::YANG10FFCalculation( TVectorD& params ) {
  setParams( params );
  print();
}

YANG10FFCalculation::YANG10FFCalculation( const std::vector<double>& params ) { 
  setParams( params );
  print();
}

YANG10FFCalculation::~YANG10FFCalculation() {}

void YANG10FFCalculation::setParams( const TVectorD& params ) {
    
  assert( params.GetNrows() == 21 );
  
  for ( int i = 0; i < 21; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void YANG10FFCalculation::setParams( const std::vector< double >& params ) {
    
  assert( params.size() == 21 );
  
  for ( int i = 0; i < 21; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void YANG10FFCalculation::print() const {

  std::cout << "Printing Scalar parameters:" << std::endl;
  
  for(int i=0; i < 21; i++){
    std::cout << params_[i] << std::endl;
  }
  
  return;
}

TVectorD YANG10FFCalculation::getParams() const{
  
  TVectorD vec(21);
  
  for(int i=0; i < 21; i++){
    vec[i] = params_[i];
  }
  
  return vec;
}

double YANG10FFCalculation::A1( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass,  params_[0], params_[1], params_[2]);
}

double YANG10FFCalculation::A2( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[3],params_[4],params_[5]);
}

double YANG10FFCalculation::A0( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[6], params_[7], params_[8]);
}

double YANG10FFCalculation::V( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[9],params_[10], params_[11]);
}

double YANG10FFCalculation::T1( const double qsq ) const { 
  return calcFabFF(qsq, m_motherMass, params_[12], params_[13], params_[14]);
}

double YANG10FFCalculation::T2( const double qsq ) const { 
  return calcFabFF(qsq,m_motherMass, params_[15], params_[16], params_[17]);
}

double YANG10FFCalculation::T3( const double qsq ) const { 
  return calcFabFF(qsq,m_motherMass, params_[18], params_[19], params_[20]);
}

double YANG10FFCalculation::T3tilde( const double qsq) const{

  double mBsq = std::pow( m_motherMass, 2 );
  double mVsq = std::pow( m_mass, 2 );

  return T2(qsq) + ( T3(qsq) * ( qsq / ( mBsq - mVsq ) ) );

}

double YANG10FFCalculation::A12( const double qsq ) const{

  double mbmk = m_motherMass + m_mass;
 
  double mkmkqsq = m_motherMass*m_motherMass -
               	  m_mass*m_mass - qsq;

  double mbmkmk = 16.0 * m_motherMass * m_mass*m_mass;


  double numer = mbmk*mbmk*mkmkqsq*A1(qsq) -
	  Kinematic::lambda(qsq, m_motherMass, m_mass)*A2(qsq);

  double denom = mbmkmk * mbmk;

	return numer / denom;
}

double YANG10FFCalculation::T23( const double qsq ) const{

  double mbmk = m_motherMass*m_motherMass - m_mass*m_mass;
  
  double mbmkqsq = m_motherMass*m_motherMass +
        	3.0*m_mass*m_mass - qsq;
  
  double denom = 8.0*m_motherMass*m_mass*m_mass*
  	( m_motherMass - m_mass );
  
  double numer = mbmk * mbmkqsq * T2(qsq) -
		Kinematic::lambda(qsq, m_motherMass, m_mass) * T3(qsq);

  return numer / denom ;
}

//=============================================================================
// Scalar Form-factors calculation
//=============================================================================
//
//---------------------------------------
//AAS07
//---------------------------------------

double AAS07FFCalculation::fplus( const double qsq ) const {
  return calcFabFF(qsq, m_motherMass, params_[0], params_[1], params_[2]);
}

double AAS07FFCalculation::fminus( const double qsq ) const {
  return calcFabFF(qsq, m_motherMass, params_[3], params_[4], params_[5]);
}

double AAS07FFCalculation::fzero( const double qsq ) const { 

  const double fp = fplus( qsq );
  const double fm = fminus( qsq ); 
  
  const double scale = qsq/( m_motherMass*m_motherMass - 
			     m_mass*m_mass );
  
  return fp + scale * fm;
}

double AAS07FFCalculation::fT(const double qsq) const {
  return calcFabFF(qsq, m_motherMass, params_[6], params_[7], params_[8]);
}

//---------------------------------------
//CFW10
//---------------------------------------

CFW10FFCalculation::CFW10FFCalculation( ){ } 

CFW10FFCalculation::CFW10FFCalculation( const CFW10FFCalculation& other ) : 
  params_( other.params_ ) {} 

CFW10FFCalculation::CFW10FFCalculation( TVectorD& params){

  setParams( params );
  print();
}

CFW10FFCalculation::CFW10FFCalculation( const std::vector<double>& params ){
  setParams( params );
  print();
}

CFW10FFCalculation::~CFW10FFCalculation() {}

void CFW10FFCalculation::setParams( const TVectorD& params ) {
    
  assert( params.GetNrows() == 9 );
  
  for ( int i = 0; i < 9; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void CFW10FFCalculation::setParams( const std::vector< double >& params ) {
    
  assert( params.size() == 9 );
  
  for ( int i = 0; i < 9; ++i ){
    params_[i] = params[i];
  }
  
  return;
}

void CFW10FFCalculation::print() const {

  std::cout << "Printing Scalar parameters:" << std::endl;
  
  for(int i=0; i < 9; i++){
    std::cout << params_[i] << std::endl;
  }
  
  return;
}

TVectorD CFW10FFCalculation::getParams() const{
  
  TVectorD vec(9);
  
  for(int i=0; i < 9; i++){
    vec[i] = params_[i];
  }
  
  return vec;
}

double CFW10FFCalculation::fplus( const double qsq ) const {
  return calcFabFF(qsq, m_motherMass, params_[0], params_[1], params_[2]);
}

double CFW10FFCalculation::fzero( const double qsq ) const {
  return calcFabFF(qsq, m_motherMass, params_[3], params_[4], params_[5]);
}

double CFW10FFCalculation::fT( const double qsq ) const {
  return calcFabFF(qsq, m_motherMass, params_[6], params_[7], params_[8]);
}

double CFW10FFCalculation::fminus( const double qsq ) const { 
 
  const double fp = fplus( qsq );
  const double f0 = fzero( qsq ); 
  
  const double scale = ( m_motherMass*m_motherMass - 
			 m_mass*m_mass ) / qsq;

  
  return scale*(f0-fp);
}
