#include "primitives/SRI.hpp"
#include "primitives/Atom.hpp"
#include <math.h>
#include <iostream>
#include <stdlib.h>



Bond::Bond(){
  
  c_constraint = NULL;
  c_is_constrained = 0;
}

void Bond::set_atoms( Atom &a, Atom &b ){
  
  c_p_a = &a;
  c_p_b = &b;
}  

void Bond::constrain(double bond_distance){

  double dsqr = bond_distance * bond_distance;
  
  c_is_constrained = 1;
  
  c_constraint = new Constraint();
  c_constraint->set_a( c_p_a->getIndex() );
  c_constraint->set_b( c_p_b->getIndex() );
  c_constraint->set_dsqr( dsqr );
}

Bond::~Bond(){
  delete c_constraint;
  c_constraint = 0;
}

void Bond::calc_forces(){
  
  /* return 0 if the bond is constrained and stop wasting cpu */
  
  if(c_is_constrained){
    
    c_potential_E = 0.0;
    return;
  }
  
  vect r_ab; /*the vector whose origin is a and end is b */
  double force; /* the force scaling factor. */
  double Fab_x; /*the x,y, and z components of the force */
  double Fab_y;
  double Fab_z;

  double aR[3], bR[3];
  double aF[3], bF[3];

  /* initialize the vector */
  
  c_p_a->getPos(aR);
  c_p_b->getPos(bR);

  r_ab.x = bR[0] - aR[0];
  r_ab.y = bR[1] - aR[1];
  r_ab.z = bR[2] - aR[2];

  r_ab.length = sqrt((r_ab.x * r_ab.x + r_ab.y * r_ab.y + r_ab.z * r_ab.z));
  
  /* calculate the force here */

  force = bond_force(r_ab.length);
  
  Fab_x = -force *  r_ab.x / r_ab.length;
  Fab_y = -force *  r_ab.y / r_ab.length;
  Fab_z = -force *  r_ab.z / r_ab.length;

  aF[0] = Fab_x;
  aF[1] = Fab_y;
  aF[2] = Fab_z;

  bF[0] = -Fab_x;
  bF[1] = -Fab_y;
  bF[2] = -Fab_z;

  c_p_a->addFrc(aF);
  c_p_b->addFrc(bF);

  return;
}
