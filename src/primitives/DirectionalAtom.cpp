#include <math.h>

#include "primitives/Atom.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "utils/simError.h"
#include "math/MatVec3.h"

void DirectionalAtom::zeroForces() {
  if( hasCoords ){

    Atom::zeroForces();
    
    trq[offsetX] = 0.0; 
    trq[offsetY] = 0.0; 
    trq[offsetZ] = 0.0;
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to zero frc and trq for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void DirectionalAtom::setCoords(void){

  if( myConfig->isAllocated() ){

    myConfig->getAtomPointers( index,
		     &pos, 
		     &vel, 
		     &frc, 
		     &trq, 
		     &Amat,
		     &mu,  
		     &ul);
  }
  else{
    sprintf( painCave.errMsg,
	     "Attempted to set Atom %d  coordinates with an unallocated "
	     "SimState object.\n", index );
    painCave.isFatal = 1;
    simError();
  }

  hasCoords = true;

}

void DirectionalAtom::setA( double the_A[3][3] ){

  if( hasCoords ){
    Amat[Axx] = the_A[0][0]; Amat[Axy] = the_A[0][1]; Amat[Axz] = the_A[0][2];
    Amat[Ayx] = the_A[1][0]; Amat[Ayy] = the_A[1][1]; Amat[Ayz] = the_A[1][2];
    Amat[Azx] = the_A[2][0]; Amat[Azy] = the_A[2][1]; Amat[Azz] = the_A[2][2];
    
    this->updateU();  
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to set Amat for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void DirectionalAtom::setI( double the_I[3][3] ){  

  int n_linear_coords, i, j;
  
  Ixx = the_I[0][0]; Ixy = the_I[0][1]; Ixz = the_I[0][2];
  Iyx = the_I[1][0]; Iyy = the_I[1][1]; Iyz = the_I[1][2];
  Izx = the_I[2][0]; Izy = the_I[2][1]; Izz = the_I[2][2];
  
  n_linear_coords = 0;

  for (i = 0; i<3; i++) {
    if (fabs(the_I[i][i]) < momIntTol) {
      is_linear = true;
      n_linear_coords++;
      linear_axis = i;
    }
  }
  
  if (n_linear_coords > 1) {
    sprintf( painCave.errMsg,
             "DirectionalAtom error.\n"
             "\tOOPSE was told to set more than one axis in this\n"
             "\tDirectionalAtom to a vanishing moment of inertia.\n"
             "\tThis should not be a DirectionalAtom.  Use an Atom.\n"
             );
      painCave.isFatal = 1;
      simError();
  }


}

void DirectionalAtom::setQ( double the_q[4] ){

  double q0Sqr, q1Sqr, q2Sqr, q3Sqr;

  if( hasCoords ){
    q0Sqr = the_q[0] * the_q[0];
    q1Sqr = the_q[1] * the_q[1];
    q2Sqr = the_q[2] * the_q[2];
    q3Sqr = the_q[3] * the_q[3];
    
    
    Amat[Axx] = q0Sqr + q1Sqr - q2Sqr - q3Sqr;
    Amat[Axy] = 2.0 * ( the_q[1] * the_q[2] + the_q[0] * the_q[3] );
    Amat[Axz] = 2.0 * ( the_q[1] * the_q[3] - the_q[0] * the_q[2] );
    
    Amat[Ayx] = 2.0 * ( the_q[1] * the_q[2] - the_q[0] * the_q[3] );
    Amat[Ayy] = q0Sqr - q1Sqr + q2Sqr - q3Sqr;
    Amat[Ayz] = 2.0 * ( the_q[2] * the_q[3] + the_q[0] * the_q[1] );
    
    Amat[Azx] = 2.0 * ( the_q[1] * the_q[3] + the_q[0] * the_q[2] );
    Amat[Azy] = 2.0 * ( the_q[2] * the_q[3] - the_q[0] * the_q[1] );
    Amat[Azz] = q0Sqr - q1Sqr -q2Sqr +q3Sqr;
    
    this->updateU();
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to set Q for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }

}

void DirectionalAtom::getA( double the_A[3][3] ){
  
  if( hasCoords ){
    the_A[0][0] = Amat[Axx];
    the_A[0][1] = Amat[Axy];
    the_A[0][2] = Amat[Axz];
    
    the_A[1][0] = Amat[Ayx];
    the_A[1][1] = Amat[Ayy];
    the_A[1][2] = Amat[Ayz];
    
    the_A[2][0] = Amat[Azx];
    the_A[2][1] = Amat[Azy];
    the_A[2][2] = Amat[Azz];
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to get Amat for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }

}

void DirectionalAtom::printAmatIndex( void ){

  if( hasCoords ){
    std::cerr << "Atom[" << index << "] index =>\n" 
	      << "[ " << Axx << ", " << Axy << ", " << Axz << " ]\n"
	      << "[ " << Ayx << ", " << Ayy << ", " << Ayz << " ]\n"
	      << "[ " << Azx << ", " << Azy << ", " << Azz << " ]\n";
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to print Amat indices for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}


void DirectionalAtom::getU( double the_u[3] ){
  
  the_u[0] = sU[2][0];
  the_u[1] = sU[2][1];
  the_u[2] = sU[2][2];
  
  this->body2Lab( the_u );
}

void DirectionalAtom::getQ( double q[4] ){
  
  double t, s;
  double ad1, ad2, ad3;

  if( hasCoords ){
    
    t = Amat[Axx] + Amat[Ayy] + Amat[Azz] + 1.0;
    if( t > 0.0 ){
      
      s = 0.5 / sqrt( t );
      q[0] = 0.25 / s;
      q[1] = (Amat[Ayz] - Amat[Azy]) * s;
      q[2] = (Amat[Azx] - Amat[Axz]) * s;
      q[3] = (Amat[Axy] - Amat[Ayx]) * s;
    }
    else{
      
      ad1 = fabs( Amat[Axx] );
      ad2 = fabs( Amat[Ayy] );
      ad3 = fabs( Amat[Azz] );
      
      if( ad1 >= ad2 && ad1 >= ad3 ){
	
	s = 2.0 * sqrt( 1.0 + Amat[Axx] - Amat[Ayy] - Amat[Azz] );
	q[0] = (Amat[Ayz] + Amat[Azy]) / s;
	q[1] = 0.5 / s;
	q[2] = (Amat[Axy] + Amat[Ayx]) / s;
	q[3] = (Amat[Axz] + Amat[Azx]) / s;
      }
      else if( ad2 >= ad1 && ad2 >= ad3 ){
	
	s = sqrt( 1.0 + Amat[Ayy] - Amat[Axx] - Amat[Azz] ) * 2.0;
	q[0] = (Amat[Axz] + Amat[Azx]) / s;
	q[1] = (Amat[Axy] + Amat[Ayx]) / s;
	q[2] = 0.5 / s;
	q[3] = (Amat[Ayz] + Amat[Azy]) / s;
      }
      else{
	
	s = sqrt( 1.0 + Amat[Azz] - Amat[Axx] - Amat[Ayy] ) * 2.0;
	q[0] = (Amat[Axy] + Amat[Ayx]) / s;
	q[1] = (Amat[Axz] + Amat[Azx]) / s;
	q[2] = (Amat[Ayz] + Amat[Azy]) / s;
	q[3] = 0.5 / s;
      }
    }
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to get Q for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void DirectionalAtom::setUnitFrameFromEuler(double phi, 
                                            double theta, 
                                            double psi) {

  double myA[3][3];
  double uFrame[3][3];
  double len;
  int i, j;
  
  myA[0][0] = (cos(phi) * cos(psi)) - (sin(phi) * cos(theta) * sin(psi));
  myA[0][1] = (sin(phi) * cos(psi)) + (cos(phi) * cos(theta) * sin(psi));
  myA[0][2] = sin(theta) * sin(psi);
  
  myA[1][0] = -(cos(phi) * sin(psi)) - (sin(phi) * cos(theta) * cos(psi));
  myA[1][1] = -(sin(phi) * sin(psi)) + (cos(phi) * cos(theta) * cos(psi));
  myA[1][2] = sin(theta) * cos(psi);
  
  myA[2][0] = sin(phi) * sin(theta);
  myA[2][1] = -cos(phi) * sin(theta);
  myA[2][2] = cos(theta);
  
  // Make the unit Frame:

  for (i=0; i < 3; i++) 
    for (j=0; j < 3; j++)
      uFrame[i][j] = 0.0;

  for (i=0; i < 3; i++)
    uFrame[i][i] = 1.0;

  // rotate by the given rotation matrix:

  matMul3(myA, uFrame, sU);

  // renormalize column vectors:

  for (i=0; i < 3; i++) {
    len = 0.0;
    for (j = 0; j < 3; j++) {
      len += sU[i][j]*sU[i][j];
    }
    len = sqrt(len);
    for (j = 0; j < 3; j++) {
      sU[i][j] /= len;     
    }
  }
   
  // sU now contains the coordinates of the 'special' frame;
    
}

void DirectionalAtom::setEuler( double phi, double theta, double psi ){
  
  if( hasCoords ){
    Amat[Axx] = (cos(phi) * cos(psi)) - (sin(phi) * cos(theta) * sin(psi));
    Amat[Axy] = (sin(phi) * cos(psi)) + (cos(phi) * cos(theta) * sin(psi));
    Amat[Axz] = sin(theta) * sin(psi);
    
    Amat[Ayx] = -(cos(phi) * sin(psi)) - (sin(phi) * cos(theta) * cos(psi));
    Amat[Ayy] = -(sin(phi) * sin(psi)) + (cos(phi) * cos(theta) * cos(psi));
    Amat[Ayz] = sin(theta) * cos(psi);
    
    Amat[Azx] = sin(phi) * sin(theta);
    Amat[Azy] = -cos(phi) * sin(theta);
    Amat[Azz] = cos(theta);
    
    this->updateU();
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to set Euler angles for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}


void DirectionalAtom::lab2Body( double r[3] ){

  double rl[3]; // the lab frame vector 
  
  if( hasCoords ){
    rl[0] = r[0];
    rl[1] = r[1];
    rl[2] = r[2];
    
    r[0] = (Amat[Axx] * rl[0]) + (Amat[Axy] * rl[1]) + (Amat[Axz] * rl[2]);
    r[1] = (Amat[Ayx] * rl[0]) + (Amat[Ayy] * rl[1]) + (Amat[Ayz] * rl[2]);
    r[2] = (Amat[Azx] * rl[0]) + (Amat[Azy] * rl[1]) + (Amat[Azz] * rl[2]);
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to convert lab2body for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }

}

void DirectionalAtom::rotateBy( double by_A[3][3]) {

  // Check this
  
  double r00, r01, r02, r10, r11, r12, r20, r21, r22;

  if( hasCoords ){

    r00 = by_A[0][0]*Amat[Axx] + by_A[0][1]*Amat[Ayx] + by_A[0][2]*Amat[Azx];
    r01 = by_A[0][0]*Amat[Axy] + by_A[0][1]*Amat[Ayy] + by_A[0][2]*Amat[Azy];
    r02 = by_A[0][0]*Amat[Axz] + by_A[0][1]*Amat[Ayz] + by_A[0][2]*Amat[Azz];
    
    r10 = by_A[1][0]*Amat[Axx] + by_A[1][1]*Amat[Ayx] + by_A[1][2]*Amat[Azx];
    r11 = by_A[1][0]*Amat[Axy] + by_A[1][1]*Amat[Ayy] + by_A[1][2]*Amat[Azy];
    r12 = by_A[1][0]*Amat[Axz] + by_A[1][1]*Amat[Ayz] + by_A[1][2]*Amat[Azz];
    
    r20 = by_A[2][0]*Amat[Axx] + by_A[2][1]*Amat[Ayx] + by_A[2][2]*Amat[Azx];
    r21 = by_A[2][0]*Amat[Axy] + by_A[2][1]*Amat[Ayy] + by_A[2][2]*Amat[Azy];
    r22 = by_A[2][0]*Amat[Axz] + by_A[2][1]*Amat[Ayz] + by_A[2][2]*Amat[Azz];
    
    Amat[Axx] = r00; Amat[Axy] = r01; Amat[Axz] = r02;
    Amat[Ayx] = r10; Amat[Ayy] = r11; Amat[Ayz] = r12;
    Amat[Azx] = r20; Amat[Azy] = r21; Amat[Azz] = r22;

  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to rotate frame for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }

}


void DirectionalAtom::body2Lab( double r[3] ){

  double rb[3]; // the body frame vector 
  
  if( hasCoords ){
    rb[0] = r[0];
    rb[1] = r[1];
    rb[2] = r[2];
    
    r[0] = (Amat[Axx] * rb[0]) + (Amat[Ayx] * rb[1]) + (Amat[Azx] * rb[2]);
    r[1] = (Amat[Axy] * rb[0]) + (Amat[Ayy] * rb[1]) + (Amat[Azy] * rb[2]);
    r[2] = (Amat[Axz] * rb[0]) + (Amat[Ayz] * rb[1]) + (Amat[Azz] * rb[2]);
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to convert body2lab for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void DirectionalAtom::updateU( void ){

  if( hasCoords ){
    ul[offsetX] = (Amat[Axx] * sU[2][0]) + 
      (Amat[Ayx] * sU[2][1]) + (Amat[Azx] * sU[2][2]);
    ul[offsetY] = (Amat[Axy] * sU[2][0]) + 
      (Amat[Ayy] * sU[2][1]) + (Amat[Azy] * sU[2][2]);
    ul[offsetZ] = (Amat[Axz] * sU[2][0]) + 
      (Amat[Ayz] * sU[2][1]) + (Amat[Azz] * sU[2][2]);
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to updateU for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void DirectionalAtom::getJ( double theJ[3] ){
  
  theJ[0] = jx;
  theJ[1] = jy;
  theJ[2] = jz;
}

void DirectionalAtom::setJ( double theJ[3] ){
  
  jx = theJ[0];
  jy = theJ[1];
  jz = theJ[2];
}

void DirectionalAtom::getTrq( double theT[3] ){
  
  if( hasCoords ){
    theT[0] = trq[offsetX];
    theT[1] = trq[offsetY];
    theT[2] = trq[offsetZ];
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to get Trq for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void DirectionalAtom::addTrq( double theT[3] ){
  
  if( hasCoords ){
    trq[offsetX] += theT[0];
    trq[offsetY] += theT[1];
    trq[offsetZ] += theT[2];
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to add Trq for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}


void DirectionalAtom::getI( double the_I[3][3] ){
  
  the_I[0][0] = Ixx;
  the_I[0][1] = Ixy;
  the_I[0][2] = Ixz;

  the_I[1][0] = Iyx;
  the_I[1][1] = Iyy;
  the_I[1][2] = Iyz;

  the_I[2][0] = Izx;
  the_I[2][1] = Izy;
  the_I[2][2] = Izz;
}

void DirectionalAtom::getGrad( double grad[6] ) {

  double myEuler[3];
  double phi, theta, psi;
  double cphi, sphi, ctheta, stheta;
  double ephi[3];
  double etheta[3];
  double epsi[3];

  this->getEulerAngles(myEuler);

  phi = myEuler[0];
  theta = myEuler[1];
  psi = myEuler[2];

  cphi = cos(phi);
  sphi = sin(phi);
  ctheta = cos(theta);
  stheta = sin(theta);

  // get unit vectors along the phi, theta and psi rotation axes

  ephi[0] = 0.0;
  ephi[1] = 0.0;
  ephi[2] = 1.0;

  etheta[0] = cphi;
  etheta[1] = sphi;
  etheta[2] = 0.0;
  
  epsi[0] = stheta * cphi;
  epsi[1] = stheta * sphi;
  epsi[2] = ctheta;
  
  for (int j = 0 ; j<3; j++)
    grad[j] = frc[j];

  grad[3] = 0;
  grad[4] = 0;
  grad[5] = 0;

  for (int j = 0; j < 3; j++ ) {
    
    grad[3] += trq[j]*ephi[j];
    grad[4] += trq[j]*etheta[j];
    grad[5] += trq[j]*epsi[j];
    
  }

}

/**
  * getEulerAngles computes a set of Euler angle values consistent
  *  with an input rotation matrix.  They are returned in the following
  * order:
  *  myEuler[0] = phi;
  *  myEuler[1] = theta;
  *  myEuler[2] = psi;
*/
void DirectionalAtom::getEulerAngles(double myEuler[3]) {

  // We use so-called "x-convention", which is the most common definition. 
  // In this convention, the rotation given by Euler angles (phi, theta, psi), where the first 
  // rotation is by an angle phi about the z-axis, the second is by an angle  
  // theta (0 <= theta <= 180)about the x-axis, and thethird is by an angle psi about the
  //z-axis (again). 
  
  
  double phi,theta,psi,eps;
  double ctheta,stheta;

  // set the tolerance for Euler angles and rotation elements
  
  eps = 1.0e-8;

  theta = acos(min(1.0,max(-1.0,Amat[Azz])));
  ctheta = Amat[Azz]; 
  stheta = sqrt(1.0 - ctheta * ctheta);

  // when sin(theta) is close to 0, we need to consider singularity
  // In this case, we can assign an arbitary value to phi (or psi), and then determine 
  // the psi (or phi) or vice-versa. We'll assume that phi always gets the rotation, and psi is 0
  // in cases of singularity.  
  // we use atan2 instead of atan, since atan2 will give us -Pi to Pi. 
  // Since 0 <= theta <= 180, sin(theta) will be always non-negative. Therefore, it never
  // change the sign of both of the parameters passed to atan2.
  
  if (fabs(stheta) <= eps){
    psi = 0.0;
    phi = atan2(-Amat[Ayx], Amat[Axx]);  
  }
  // we only have one unique solution
  else{    
      phi = atan2(Amat[Azx], -Amat[Azy]);
      psi = atan2(Amat[Axz], Amat[Ayz]);
  }

  //wrap phi and psi, make sure they are in the range from 0 to 2*Pi
  //if (phi < 0)
  //  phi += M_PI;

  //if (psi < 0)
  //  psi += M_PI;

  myEuler[0] = phi;
  myEuler[1] = theta;
  myEuler[2] = psi;
  
  return;
}

double DirectionalAtom::getZangle( ){
  
  if( hasCoords ){
    return zAngle;
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to get zAngle for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
    return 0;
  }
}

void DirectionalAtom::setZangle( double zAng ){
  
  if( hasCoords ){
    zAngle = zAng;
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to set zAngle for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void DirectionalAtom::addZangle( double zAng ){
  
  if( hasCoords ){
    zAngle += zAng;
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to add zAngle to atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

double DirectionalAtom::max(double x, double  y) {  
  return (x > y) ? x : y;
}

double DirectionalAtom::min(double x, double  y) {  
  return (x > y) ? y : x;
}
