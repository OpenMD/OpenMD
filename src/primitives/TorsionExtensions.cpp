#include "primitives/SRI.hpp"

CubicTorsion::CubicTorsion( Atom &a, Atom &b, Atom &c, Atom &d ){

  set_atoms( a, b, c, d );
  k1 = 0.0;
  k2 = 0.0;
  k3 = 0.0;
  k4 = 0.0;
}

void CubicTorsion::setConstants( double the_k1, double the_k2, double the_k3,
				 double the_k4 ){
  
  k1 = the_k1;
  k2 = the_k2;
  k3 = the_k3;
  k4 = the_k4;
}

double CubicTorsion::torsion_force( double cos_phi ){
  
  double cp, cp2, cp3;
  double force;
  
  cp = cos_phi;
  cp2 = cp * cp;
  cp3 = cp2 * cp;

  c_potential_E = ( k1 * cp3 ) + ( k2 * cp2 ) + ( k3 * cp ) + k4;

  force = -( ( 3.0 * k1 * cp2 ) + ( 2.0 * k2 * cp ) + k3 );
  return force;
}
