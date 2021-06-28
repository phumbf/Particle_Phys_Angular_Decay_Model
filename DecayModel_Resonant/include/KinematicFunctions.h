#ifndef KINEMATICFUNCTIONS_H
#define KINEMATICFUNCTIONS_H

//useful kinematic terms relevant for the maths

namespace KinematicFunctions { 
  
  double meffective( const double mean, const double min, const double max );

  double momentum  ( const double m, const double m1, const double m2 ); 
    
  double barrier   ( const double p0, 
		     const double p , 
		     const double radius, 
		     const unsigned int LX );
  
  double parity    ( const int pz, const unsigned int jz ); 

  double beta      ( const double msq, const double m );
  
}

#endif
