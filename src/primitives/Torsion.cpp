#include "SRI.hpp"
#include "Atom.hpp"
#include <math.h>
#include <iostream>
#include <stdlib.h>

void Torsion::set_atoms( Atom &a, Atom &b, Atom &c, Atom &d){
  c_p_a = &a;
  c_p_b = &b;
  c_p_c = &c;
  c_p_d = &d;
}


void Torsion::calc_forces(){
  
  /**********************************************************************
   * 
   * initialize vectors
   *
   ***********************************************************************/
  
  vect r_ab; /* the vector whose origin is a and end is b */
  vect r_cb; /* the vector whose origin is c and end is b */
  vect r_cd; /* the vector whose origin is c and end is b */
  vect r_cr1; /* the cross product of r_ab and r_cb */
  vect r_cr2; /* the cross product of r_cb and r_cd */

  double r_cr1_x2; /* the components of r_cr1 squared */
  double r_cr1_y2;
  double r_cr1_z2;
  
  double r_cr2_x2; /* the components of r_cr2 squared */
  double r_cr2_y2;
  double r_cr2_z2;

  double r_cr1_sqr; /* the length of r_cr1 squared */
  double r_cr2_sqr; /* the length of r_cr2 squared */
  
  double r_cr1_r_cr2; /* the length of r_cr1 * length of r_cr2 */
  
  double aR[3], bR[3], cR[3], dR[3];
  double aF[3], bF[3], cF[3], dF[3];

  c_p_a->getPos( aR );
  c_p_b->getPos( bR );
  c_p_c->getPos( cR );
  c_p_d->getPos( dR );

  r_ab.x = bR[0] - aR[0];
  r_ab.y = bR[1] - aR[1];
  r_ab.z = bR[2] - aR[2];
  r_ab.length  = sqrt((r_ab.x * r_ab.x + r_ab.y * r_ab.y + r_ab.z * r_ab.z));

  r_cb.x = bR[0] - cR[0];
  r_cb.y = bR[1] - cR[1];
  r_cb.z = bR[2] - cR[2];
  r_cb.length = sqrt((r_cb.x * r_cb.x + r_cb.y * r_cb.y + r_cb.z * r_cb.z));

  r_cd.x = dR[0] - cR[0];
  r_cd.y = dR[1] - cR[1];
  r_cd.z = dR[2] - cR[2];
  r_cd.length = sqrt((r_cd.x * r_cd.x + r_cd.y * r_cd.y + r_cd.z * r_cd.z));

  r_cr1.x = r_ab.y * r_cb.z - r_cb.y * r_ab.z;
  r_cr1.y = r_ab.z * r_cb.x - r_cb.z * r_ab.x;
  r_cr1.z = r_ab.x * r_cb.y - r_cb.x * r_ab.y;
  r_cr1_x2 = r_cr1.x * r_cr1.x;
  r_cr1_y2 = r_cr1.y * r_cr1.y;
  r_cr1_z2 = r_cr1.z * r_cr1.z;
  r_cr1_sqr = r_cr1_x2 + r_cr1_y2 + r_cr1_z2;
  r_cr1.length = sqrt(r_cr1_sqr);

  r_cr2.x = r_cb.y * r_cd.z - r_cd.y * r_cb.z;
  r_cr2.y = r_cb.z * r_cd.x - r_cd.z * r_cb.x;
  r_cr2.z = r_cb.x * r_cd.y - r_cd.x * r_cb.y;
  r_cr2_x2 = r_cr2.x * r_cr2.x;
  r_cr2_y2 = r_cr2.y * r_cr2.y;
  r_cr2_z2 = r_cr2.z * r_cr2.z;
  r_cr2_sqr = r_cr2_x2 + r_cr2_y2 + r_cr2_z2;
  r_cr2.length = sqrt(r_cr2_sqr);

  r_cr1_r_cr2 = r_cr1.length * r_cr2.length;

  /**********************************************************************
   *
   * dot product and angle calculations 
   *
   ***********************************************************************/
  
  double cr1_dot_cr2; /* the dot product of the cr1 and cr2 vectors */
  double cos_phi; /* the cosine of the torsion angle */

  cr1_dot_cr2 = r_cr1.x * r_cr2.x + r_cr1.y * r_cr2.y + r_cr1.z * r_cr2.z;
  
  cos_phi = cr1_dot_cr2 / r_cr1_r_cr2;
  
   /* adjust for the granularity of the numbers for angles near 0 or pi */

  if(cos_phi > 1.0) cos_phi = 1.0;
  if(cos_phi < -1.0) cos_phi = -1.0;


  /********************************************************************
   *
   * This next section calculates derivatives needed for the force
   * calculation
   *
   ********************************************************************/


  /* the derivatives of cos phi with respect to the x, y,
     and z components of vectors cr1 and cr2. */
  double d_cos_dx_cr1;
  double d_cos_dy_cr1;
  double d_cos_dz_cr1;
  double d_cos_dx_cr2;
  double d_cos_dy_cr2;
  double d_cos_dz_cr2;

  d_cos_dx_cr1 = r_cr2.x / r_cr1_r_cr2 - (cos_phi * r_cr1.x) / r_cr1_sqr;
  d_cos_dy_cr1 = r_cr2.y / r_cr1_r_cr2 - (cos_phi * r_cr1.y) / r_cr1_sqr;
  d_cos_dz_cr1 = r_cr2.z / r_cr1_r_cr2 - (cos_phi * r_cr1.z) / r_cr1_sqr;

  d_cos_dx_cr2 = r_cr1.x / r_cr1_r_cr2 - (cos_phi * r_cr2.x) / r_cr2_sqr;
  d_cos_dy_cr2 = r_cr1.y / r_cr1_r_cr2 - (cos_phi * r_cr2.y) / r_cr2_sqr;
  d_cos_dz_cr2 = r_cr1.z / r_cr1_r_cr2 - (cos_phi * r_cr2.z) / r_cr2_sqr;

  /***********************************************************************
   *
   * Calculate the actual forces and place them in the atoms. 
   *
   ***********************************************************************/

  double force; /*the force scaling factor */

  force = torsion_force(cos_phi);

  aF[0] = force * (d_cos_dy_cr1 * r_cb.z - d_cos_dz_cr1 * r_cb.y);
  aF[1] = force * (d_cos_dz_cr1 * r_cb.x - d_cos_dx_cr1 * r_cb.z);
  aF[2] = force * (d_cos_dx_cr1 * r_cb.y - d_cos_dy_cr1 * r_cb.x);

  bF[0] = force * (  d_cos_dy_cr1 * (r_ab.z - r_cb.z)
		   - d_cos_dy_cr2 *  r_cd.z	  
		   + d_cos_dz_cr1 * (r_cb.y - r_ab.y)
		   + d_cos_dz_cr2 *  r_cd.y);
  bF[1] = force * (  d_cos_dx_cr1 * (r_cb.z - r_ab.z)
		   + d_cos_dx_cr2 *  r_cd.z	  
		   + d_cos_dz_cr1 * (r_ab.x - r_cb.x)
		   - d_cos_dz_cr2 *  r_cd.x);
  bF[2] = force * (  d_cos_dx_cr1 * (r_ab.y - r_cb.y)
		   - d_cos_dx_cr2 *  r_cd.y	  
		   + d_cos_dy_cr1 * (r_cb.x - r_ab.x)
		   + d_cos_dy_cr2 *  r_cd.x);

  cF[0] = force * (- d_cos_dy_cr1 *  r_ab.z
		   - d_cos_dy_cr2 * (r_cb.z - r_cd.z)
		   + d_cos_dz_cr1 *  r_ab.y
		   - d_cos_dz_cr2 * (r_cd.y - r_cb.y));
  cF[1] = force * (  d_cos_dx_cr1 *  r_ab.z
		   - d_cos_dx_cr2 * (r_cd.z - r_cb.z)
		   - d_cos_dz_cr1 *  r_ab.x
		   - d_cos_dz_cr2 * (r_cb.x - r_cd.x));
  cF[2] = force * (- d_cos_dx_cr1 *  r_ab.y
		   - d_cos_dx_cr2 * (r_cb.y - r_cd.y)
		   + d_cos_dy_cr1 *  r_ab.x
		   - d_cos_dy_cr2 * (r_cd.x - r_cb.x));

  dF[0] = force * (d_cos_dy_cr2 * r_cb.z - d_cos_dz_cr2 * r_cb.y);
  dF[1] = force * (d_cos_dz_cr2 * r_cb.x - d_cos_dx_cr2 * r_cb.z);
  dF[2] = force * (d_cos_dx_cr2 * r_cb.y - d_cos_dy_cr2 * r_cb.x);


  c_p_a->addFrc(aF);
  c_p_b->addFrc(bF);
  c_p_c->addFrc(cF);
  c_p_d->addFrc(dF);
}
