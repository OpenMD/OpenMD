#include <math.h>

#include "primitives/SRI.hpp"
#include "utils/simError.h"

QuadraticBend::QuadraticBend( Atom &a, Atom &b, Atom &c ){

  set_atoms( a, b, c );
  c1 = 0.0;
  c2 = 0.0;
  c3 = 0.0;
  theta0 = 0.0;
}

void QuadraticBend::setConstants( double the_c1, double the_c2, double the_c3, 
				  double the_Th0 ){
  c1 = the_c1;
  c2 = the_c2;
  c3 = the_c3;
  theta0 = the_Th0;
}


double QuadraticBend::bend_force( double theta ){

  double dt, dt2;
  double force;




  dt = ( theta - theta0 ) * M_PI / 180.0;
  dt2 = dt * dt;

  c_potential_E = ( c1 * dt2 ) + ( c2 * dt ) + c3;
  force = -( ( 2.0 * c1 * dt ) + c2 );
  return force;
}
