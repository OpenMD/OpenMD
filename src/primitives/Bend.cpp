#include "primitives/SRI.hpp"
#include "primitives/Atom.hpp"

#include <math.h>
#include <iostream>
#include <stdlib.h>

void Bend::set_atoms( Atom &a, Atom &b, Atom &c){

  c_p_a = &a;
  c_p_b = &b;
  c_p_c = &c;

  c_potential_E = 0.0;
}


void Bend::calc_forces(){
  
  double dx,dy,dz,gx,gy,gz,dx2,dy2,dz2,gx2,gy2,gz2;
  double rij2, rkj2, riji2, rkji2, dot, denom, cosang, angl;
  
  double sina2, sinai;

  double comf2, comf3, comf4;
  double dcsidx, dcsidy, dcsidz, dcskdx, dcskdy, dcskdz;
  // double dcsjdx, dcsjdy, dcsjdz;
  double dadxi, dadyi, dadzi;
  double dadxk, dadyk, dadzk;//, dadxj, dadyj, dadzj;
  double daxi, dayi, dazi, daxk, dayk, dazk, daxj, dayj, dazj;
  
  double aR[3], bR[3], cR[3];
  double aF[3], bF[3], cF[3];

  c_p_a->getPos( aR );
  c_p_b->getPos( bR );
  c_p_c->getPos( cR );
  

  dx = aR[0] - bR[0];
  dy = aR[1] - bR[1];
  dz = aR[2] - bR[2];
 
  gx = cR[0] - bR[0];
  gy = cR[1] - bR[1];
  gz = cR[2] - bR[2];
  
  dx2 = dx * dx;
  dy2 = dy * dy;
  dz2 = dz * dz;

  gx2 = gx * gx;
  gy2 = gy * gy;
  gz2 = gz * gz;
  
  rij2 = dx2 + dy2 + dz2;
  rkj2 = gx2 + gy2 + gz2;
  
  riji2 = 1.0 / rij2;
  rkji2 = 1.0 / rkj2;

  dot = dx * gx + dy * gy + dz * gz;
  denom = sqrt((riji2 * rkji2));
  cosang = dot * denom;

  if(cosang > 1.0)cosang = 1.0;
  if(cosang < -1.0) cosang = -1.0;

  angl = acos(cosang);
  angl = angl * 180.0 / M_PI;

  sina2 = 1.0 - cosang*cosang;
  if(fabs(sina2) < 1.0E-12 ) sina2 = 1.0E-12;
  sinai = 1.0 / sqrt(sina2);

  comf2 = cosang * riji2;
  comf3 = cosang * rkji2;
  comf4 = bend_force(angl);


  dcsidx = gx*denom - comf2*dx;
  dcsidy = gy*denom - comf2*dy;
  dcsidz = gz*denom - comf2*dz;
  
  dcskdx = dx*denom - comf3*gx;
  dcskdy = dy*denom - comf3*gy;
  dcskdz = dz*denom - comf3*gz;
  
//   dcsjdx = -dcsidx - dcskdx;
//   dcsjdy = -dcsidy - dcskdy;
//   dcsjdz = -dcsidz - dcskdz;

  dadxi = -sinai*dcsidx;
  dadyi = -sinai*dcsidy;
  dadzi = -sinai*dcsidz;

  dadxk = -sinai*dcskdx;
  dadyk = -sinai*dcskdy;
  dadzk = -sinai*dcskdz;

//   dadxj = -dadxi - dadxk;
//   dadyj = -dadyi - dadyk;
//   dadzj = -dadzi - dadzk;

  daxi = comf4*dadxi;
  dayi = comf4*dadyi;
  dazi = comf4*dadzi;

  daxk = comf4*dadxk;
  dayk = comf4*dadyk;
  dazk = comf4*dadzk;
  
  daxj = -daxi - daxk;
  dayj = -dayi - dayk;
  dazj = -dazi - dazk;
  
  aF[0] = daxi;
  aF[1] = dayi;
  aF[2] = dazi;

  bF[0] = daxj;
  bF[1] = dayj;
  bF[2] = dazj;

  cF[0] = daxk;
  cF[1] = dayk;
  cF[2] = dazk;

  c_p_a->addFrc(aF);
  c_p_b->addFrc(bF);
  c_p_c->addFrc(cF);

  return;
}
