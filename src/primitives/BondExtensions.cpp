

#include "SRI.hpp"


ConstrainedBond::ConstrainedBond( Atom &a, Atom &b, double constraint ){

  set_atoms( a, b );
  constrain( constraint );
  d0 = constraint;
  c_potential_E = 0.0;
}


HarmonicBond::HarmonicBond(Atom &a, Atom &b, double theR0, double theK0 ){
  
  set_atoms( a, b );
  d0 = theR0;
  k0 = theK0;
  c_potential_E = 0.0;
}


double HarmonicBond::bond_force( double r_ab ){

  double force;
  double dr, dr2;

  dr = r_ab - d0;
  dr2 = dr * dr;

  c_potential_E = 0.5 * k0 * dr2;
  force = - k0 * dr;
  return force;
 
}
