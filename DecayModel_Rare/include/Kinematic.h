#ifndef KINEMATIC_H
#define KINEMATIC_H 1

//Useful kinematic terms used to determine the form factor variation with momenta
namespace Kinematic { 
  
  extern double beta( const double qsq, const double m ); 

  extern double betasq( const double qsq, const double m );
  
  extern double lambda( const double qsq, const double m, const double mB );

  extern double z( const double qsq, const double m, const double mB );

  extern double zseries( const double qsq,
			 const double m, 
			 const double mB,
			 const double a0,
			 const double a1,
			 const double a2 ) ;

  extern double pole( const double qsq, const double mpole );

  extern double ffparam( const double shat, 
			 const double f, 
			 const double a, 
			 const double b );

  extern double shat( const double qsq, 
		      const double mB );

  extern bool isForbidden( const double qsq, 
			   const double mlepton, 
			   const double mhadron, 
			   const double mhone  ,
			   const double mhtwo  ,
			   const double mB );
}


#endif
