#include <math.h>
#include "primitives/RigidBody.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "utils/simError.h"
#include "math/MatVec3.h"

RigidBody::RigidBody() : StuntDouble() {
  objType = OT_RIGIDBODY;
  is_linear = false;
  linear_axis =  -1;
  momIntTol = 1e-6;
}

RigidBody::~RigidBody() {
}

void RigidBody::addAtom(Atom* at, AtomStamp* ats) {

  vec3 coords;
  vec3 euler;
  mat3x3 Atmp;

  myAtoms.push_back(at);
 
  if( !ats->havePosition() ){
    sprintf( painCave.errMsg,
             "RigidBody error.\n"
             "\tAtom %s does not have a position specified.\n"
             "\tThis means RigidBody cannot set up reference coordinates.\n",
             ats->getType() );
    painCave.isFatal = 1;
    simError();
  }
  
  coords[0] = ats->getPosX();
  coords[1] = ats->getPosY();
  coords[2] = ats->getPosZ();

  refCoords.push_back(coords);
  
  if (at->isDirectional()) {   

    if( !ats->haveOrientation() ){
      sprintf( painCave.errMsg,
               "RigidBody error.\n"
               "\tAtom %s does not have an orientation specified.\n"
               "\tThis means RigidBody cannot set up reference orientations.\n",
               ats->getType() );
      painCave.isFatal = 1;
      simError();
    }    
    
    euler[0] = ats->getEulerPhi();
    euler[1] = ats->getEulerTheta();
    euler[2] = ats->getEulerPsi();
    
    doEulerToRotMat(euler, Atmp);
    
    refOrients.push_back(Atmp);
    
  }
}

void RigidBody::getPos(double theP[3]){
  for (int i = 0; i < 3 ; i++) 
    theP[i] = pos[i];
}	

void RigidBody::setPos(double theP[3]){
  for (int i = 0; i < 3 ; i++) 
    pos[i] = theP[i];
}	

void RigidBody::getVel(double theV[3]){
  for (int i = 0; i < 3 ; i++) 
    theV[i] = vel[i];
}	

void RigidBody::setVel(double theV[3]){
  for (int i = 0; i < 3 ; i++) 
    vel[i] = theV[i];
}	

void RigidBody::getFrc(double theF[3]){
  for (int i = 0; i < 3 ; i++) 
    theF[i] = frc[i];
}	

void RigidBody::addFrc(double theF[3]){
  for (int i = 0; i < 3 ; i++) 
    frc[i] += theF[i];
}    

void RigidBody::zeroForces() {

  for (int i = 0; i < 3; i++) {
    frc[i] = 0.0;
    trq[i] = 0.0;
  }

}

void RigidBody::setEuler( double phi, double theta, double psi ){
  
    A[0][0] = (cos(phi) * cos(psi)) - (sin(phi) * cos(theta) * sin(psi));
    A[0][1] = (sin(phi) * cos(psi)) + (cos(phi) * cos(theta) * sin(psi));
    A[0][2] = sin(theta) * sin(psi);
    
    A[1][0] = -(cos(phi) * sin(psi)) - (sin(phi) * cos(theta) * cos(psi));
    A[1][1] = -(sin(phi) * sin(psi)) + (cos(phi) * cos(theta) * cos(psi));
    A[1][2] = sin(theta) * cos(psi);
    
    A[2][0] = sin(phi) * sin(theta);
    A[2][1] = -cos(phi) * sin(theta);
    A[2][2] = cos(theta);

}

void RigidBody::getQ( double q[4] ){
  
  double t, s;
  double ad1, ad2, ad3;
    
  t = A[0][0] + A[1][1] + A[2][2] + 1.0;
  if( t > 0.0 ){
    
    s = 0.5 / sqrt( t );
    q[0] = 0.25 / s;
    q[1] = (A[1][2] - A[2][1]) * s;
    q[2] = (A[2][0] - A[0][2]) * s;
    q[3] = (A[0][1] - A[1][0]) * s;
  }
  else{
    
    ad1 = fabs( A[0][0] );
    ad2 = fabs( A[1][1] );
    ad3 = fabs( A[2][2] );
    
    if( ad1 >= ad2 && ad1 >= ad3 ){
      
      s = 2.0 * sqrt( 1.0 + A[0][0] - A[1][1] - A[2][2] );
      q[0] = (A[1][2] + A[2][1]) / s;
      q[1] = 0.5 / s;
      q[2] = (A[0][1] + A[1][0]) / s;
      q[3] = (A[0][2] + A[2][0]) / s;
    }
    else if( ad2 >= ad1 && ad2 >= ad3 ){
      
      s = sqrt( 1.0 + A[1][1] - A[0][0] - A[2][2] ) * 2.0;
      q[0] = (A[0][2] + A[2][0]) / s;
      q[1] = (A[0][1] + A[1][0]) / s;
      q[2] = 0.5 / s;
      q[3] = (A[1][2] + A[2][1]) / s;
    }
    else{
      
      s = sqrt( 1.0 + A[2][2] - A[0][0] - A[1][1] ) * 2.0;
      q[0] = (A[0][1] + A[1][0]) / s;
      q[1] = (A[0][2] + A[2][0]) / s;
      q[2] = (A[1][2] + A[2][1]) / s;
      q[3] = 0.5 / s;
    }
  }
}

void RigidBody::setQ( double the_q[4] ){

  double q0Sqr, q1Sqr, q2Sqr, q3Sqr;
  
  q0Sqr = the_q[0] * the_q[0];
  q1Sqr = the_q[1] * the_q[1];
  q2Sqr = the_q[2] * the_q[2];
  q3Sqr = the_q[3] * the_q[3];
  
  A[0][0] = q0Sqr + q1Sqr - q2Sqr - q3Sqr;
  A[0][1] = 2.0 * ( the_q[1] * the_q[2] + the_q[0] * the_q[3] );
  A[0][2] = 2.0 * ( the_q[1] * the_q[3] - the_q[0] * the_q[2] );
  
  A[1][0] = 2.0 * ( the_q[1] * the_q[2] - the_q[0] * the_q[3] );
  A[1][1] = q0Sqr - q1Sqr + q2Sqr - q3Sqr;
  A[1][2] = 2.0 * ( the_q[2] * the_q[3] + the_q[0] * the_q[1] );
  
  A[2][0] = 2.0 * ( the_q[1] * the_q[3] + the_q[0] * the_q[2] );
  A[2][1] = 2.0 * ( the_q[2] * the_q[3] - the_q[0] * the_q[1] );
  A[2][2] = q0Sqr - q1Sqr -q2Sqr +q3Sqr;   

}

void RigidBody::getA( double the_A[3][3] ){
  
  for (int i = 0; i < 3; i++) 
    for (int j = 0; j < 3; j++) 
      the_A[i][j] = A[i][j];

}

void RigidBody::setA( double the_A[3][3] ){

  for (int i = 0; i < 3; i++) 
    for (int j = 0; j < 3; j++) 
      A[i][j] = the_A[i][j];
   
}

void RigidBody::getJ( double theJ[3] ){
  
  for (int i = 0; i < 3; i++)
    theJ[i] = ji[i];

}

void RigidBody::setJ( double theJ[3] ){
  
  for (int i = 0; i < 3; i++)
    ji[i] = theJ[i];

}

void RigidBody::getTrq(double theT[3]){
  for (int i = 0; i < 3 ; i++) 
    theT[i] = trq[i];
}	

void RigidBody::addTrq(double theT[3]){
  for (int i = 0; i < 3 ; i++) 
    trq[i] += theT[i];
}	

void RigidBody::getI( double the_I[3][3] ){  

    for (int i = 0; i < 3; i++) 
      for (int j = 0; j < 3; j++) 
        the_I[i][j] = I[i][j];

}

void RigidBody::lab2Body( double r[3] ){

  double rl[3]; // the lab frame vector 
  
  rl[0] = r[0];
  rl[1] = r[1];
  rl[2] = r[2];
  
  r[0] = (A[0][0] * rl[0]) + (A[0][1] * rl[1]) + (A[0][2] * rl[2]);
  r[1] = (A[1][0] * rl[0]) + (A[1][1] * rl[1]) + (A[1][2] * rl[2]);
  r[2] = (A[2][0] * rl[0]) + (A[2][1] * rl[1]) + (A[2][2] * rl[2]);

}

void RigidBody::body2Lab( double r[3] ){

  double rb[3]; // the body frame vector 
  
  rb[0] = r[0];
  rb[1] = r[1];
  rb[2] = r[2];
  
  r[0] = (A[0][0] * rb[0]) + (A[1][0] * rb[1]) + (A[2][0] * rb[2]);
  r[1] = (A[0][1] * rb[0]) + (A[1][1] * rb[1]) + (A[2][1] * rb[2]);
  r[2] = (A[0][2] * rb[0]) + (A[1][2] * rb[1]) + (A[2][2] * rb[2]);

}

double RigidBody::getZangle( ){
    return zAngle;
}

void RigidBody::setZangle( double zAng ){
    zAngle = zAng;
}

void RigidBody::addZangle( double zAng ){
    zAngle += zAng;
}

void RigidBody::calcRefCoords( ) {

  int i,j,k, it, n_linear_coords;
  double mtmp;
  vec3 apos;
  double refCOM[3];
  vec3 ptmp;
  double Itmp[3][3];
  double evals[3];
  double evects[3][3];
  double r, r2, len;

  // First, find the center of mass:
  
  mass = 0.0;
  for (j=0; j<3; j++)
    refCOM[j] = 0.0;
  
  for (i = 0; i < myAtoms.size(); i++) {
    mtmp = myAtoms[i]->getMass();
    mass += mtmp;

    apos = refCoords[i];
    
    for(j = 0; j < 3; j++) {
      refCOM[j] += apos[j]*mtmp;     
    }    
  }
  
  for(j = 0; j < 3; j++) 
    refCOM[j] /= mass;

// Next, move the origin of the reference coordinate system to the COM:

  for (i = 0; i < myAtoms.size(); i++) {
    apos = refCoords[i];
    for (j=0; j < 3; j++) {
      apos[j] = apos[j] - refCOM[j];
    }
    refCoords[i] = apos;
  }

// Moment of Inertia calculation

  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++)
      Itmp[i][j] = 0.0;  
  
  for (it = 0; it < myAtoms.size(); it++) {

    mtmp = myAtoms[it]->getMass();
    ptmp = refCoords[it];
    r= norm3(ptmp.vec);
    r2 = r*r;
    
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        
        if (i==j) Itmp[i][j] += mtmp * r2;

        Itmp[i][j] -= mtmp * ptmp.vec[i]*ptmp.vec[j];
      }
    }
  }
  
  diagonalize3x3(Itmp, evals, sU);
  
  // zero out I and then fill the diagonals with the moments of inertia:

  n_linear_coords = 0;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      I[i][j] = 0.0;  
    }
    I[i][i] = evals[i];

    if (fabs(evals[i]) < momIntTol) {
      is_linear = true;
      n_linear_coords++;
      linear_axis = i;
    }
  }

  if (n_linear_coords > 1) {
          sprintf( painCave.errMsg,
               "RigidBody error.\n"
               "\tOOPSE found more than one axis in this rigid body with a vanishing \n"
               "\tmoment of inertia.  This can happen in one of three ways:\n"
               "\t 1) Only one atom was specified, or \n"
               "\t 2) All atoms were specified at the same location, or\n"
               "\t 3) The programmers did something stupid.\n"
               "\tIt is silly to use a rigid body to describe this situation.  Be smarter.\n"
               );
      painCave.isFatal = 1;
      simError();
  }
  
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
}

void RigidBody::doEulerToRotMat(vec3 &euler, mat3x3 &myA ){

  double phi, theta, psi;
  
  phi = euler[0];
  theta = euler[1];
  psi = euler[2];
  
  myA[0][0] = (cos(phi) * cos(psi)) - (sin(phi) * cos(theta) * sin(psi));
  myA[0][1] = (sin(phi) * cos(psi)) + (cos(phi) * cos(theta) * sin(psi));
  myA[0][2] = sin(theta) * sin(psi);
  
  myA[1][0] = -(cos(phi) * sin(psi)) - (sin(phi) * cos(theta) * cos(psi));
  myA[1][1] = -(sin(phi) * sin(psi)) + (cos(phi) * cos(theta) * cos(psi));
  myA[1][2] = sin(theta) * cos(psi);
  
  myA[2][0] = sin(phi) * sin(theta);
  myA[2][1] = -cos(phi) * sin(theta);
  myA[2][2] = cos(theta);

}

void RigidBody::calcForcesAndTorques() {

  // Convert Atomic forces and torques to total forces and torques:
  int i, j;
  double apos[3];
  double afrc[3];
  double atrq[3];
  double rpos[3];

  zeroForces();
  
  for (i = 0; i < myAtoms.size(); i++) {

    myAtoms[i]->getPos(apos);
    myAtoms[i]->getFrc(afrc);

    for (j=0; j<3; j++) {
      rpos[j] = apos[j] - pos[j];
      frc[j] += afrc[j];
    }
    
    trq[0] += rpos[1]*afrc[2] - rpos[2]*afrc[1];
    trq[1] += rpos[2]*afrc[0] - rpos[0]*afrc[2];
    trq[2] += rpos[0]*afrc[1] - rpos[1]*afrc[0];

    // If the atom has a torque associated with it, then we also need to 
    // migrate the torques onto the center of mass:

    if (myAtoms[i]->isDirectional()) {

      myAtoms[i]->getTrq(atrq);
      
      for (j=0; j<3; j++) 
        trq[j] += atrq[j];
    }
  }

  // Convert Torque to Body-fixed coordinates:
  // (Actually, on second thought, don't.  Integrator does this now.)
  // lab2Body(trq);

}

void RigidBody::updateAtoms() {
  int i, j;
  vec3 ref;
  double apos[3];
  DirectionalAtom* dAtom;
  
  for (i = 0; i < myAtoms.size(); i++) {
     
    ref = refCoords[i];

    body2Lab(ref.vec);
    
    for (j = 0; j<3; j++) 
      apos[j] = pos[j] + ref.vec[j];
    
    myAtoms[i]->setPos(apos);
    
    if (myAtoms[i]->isDirectional()) {
      
      dAtom = (DirectionalAtom *) myAtoms[i];
      dAtom->rotateBy( A );
      
    }
  }  
}

void RigidBody::getGrad( double grad[6] ) {

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

  grad[3] = 0.0;
  grad[4] = 0.0;
  grad[5] = 0.0;
  
  for (int j = 0; j < 3; j++ ) {
    
    grad[3] += trq[j]*ephi[j];
    grad[4] += trq[j]*etheta[j];
    grad[5] += trq[j]*epsi[j];
    
  }
  
}

/**
  * getEulerAngles computes a set of Euler angle values consistent
  * with an input rotation matrix.  They are returned in the following
  * order:
  *  myEuler[0] = phi;
  *  myEuler[1] = theta;
  *  myEuler[2] = psi;
*/
void RigidBody::getEulerAngles(double myEuler[3]) {

  // We use so-called "x-convention", which is the most common
  // definition.  In this convention, the rotation given by Euler
  // angles (phi, theta, psi), where the first rotation is by an angle
  // phi about the z-axis, the second is by an angle theta (0 <= theta
  // <= 180) about the x-axis, and the third is by an angle psi about
  // the z-axis (again).
  
  
  double phi,theta,psi,eps;
  double pi;
  double cphi,ctheta,cpsi;
  double sphi,stheta,spsi;
  double b[3];
  int flip[3];
  
  // set the tolerance for Euler angles and rotation elements
  
  eps = 1.0e-8;

  theta = acos(min(1.0,max(-1.0,A[2][2])));
  ctheta = A[2][2]; 
  stheta = sqrt(1.0 - ctheta * ctheta);

  // when sin(theta) is close to 0, we need to consider the
  // possibility of a singularity. In this case, we can assign an
  // arbitary value to phi (or psi), and then determine the psi (or
  // phi) or vice-versa.  We'll assume that phi always gets the
  // rotation, and psi is 0 in cases of singularity.  we use atan2
  // instead of atan, since atan2 will give us -Pi to Pi.  Since 0 <=
  // theta <= 180, sin(theta) will be always non-negative. Therefore,
  // it never changes the sign of both of the parameters passed to
  // atan2.
  
  if (fabs(stheta) <= eps){
    psi = 0.0;
    phi = atan2(-A[1][0], A[0][0]);  
  }
  // we only have one unique solution
  else{    
    phi = atan2(A[2][0], -A[2][1]);
    psi = atan2(A[0][2], A[1][2]);
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

double RigidBody::max(double x, double  y) {  
  return (x > y) ? x : y;
}

double RigidBody::min(double x, double  y) {  
  return (x > y) ? y : x;
}

void RigidBody::findCOM() {
  
  size_t i;
  int j;
  double mtmp;
  double ptmp[3];
  double vtmp[3];
  
  for(j = 0; j < 3; j++) {
    pos[j] = 0.0;
    vel[j] = 0.0;
  }
  mass = 0.0;
  
  for (i = 0; i < myAtoms.size(); i++) {
    
    mtmp = myAtoms[i]->getMass();    
    myAtoms[i]->getPos(ptmp);
    myAtoms[i]->getVel(vtmp);
    
    mass += mtmp;
    
    for(j = 0; j < 3; j++) {
      pos[j] += ptmp[j]*mtmp;
      vel[j] += vtmp[j]*mtmp;
    }
    
  }
  
  for(j = 0; j < 3; j++) {
    pos[j] /= mass;
    vel[j] /= mass;
  }

}

void RigidBody::accept(BaseVisitor* v){
  vector<Atom*>::iterator atomIter;
  v->visit(this);

  //for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
  //  (*atomIter)->accept(v);
}
void RigidBody::getAtomRefCoor(double pos[3], int index){
  vec3 ref;

  ref = refCoords[index];
  pos[0] = ref[0];
  pos[1] = ref[1];
  pos[2] = ref[2];
  
}


void RigidBody::getAtomPos(double theP[3], int index){
  vec3 ref;

  if (index >= myAtoms.size())
    cerr << index << " is an invalid index, current rigid body contains " << myAtoms.size() << "atoms" << endl;

  ref = refCoords[index];
  body2Lab(ref.vec);
  
  theP[0] = pos[0] + ref[0];
  theP[1] = pos[1] + ref[1];
  theP[2] = pos[2] + ref[2];
}


void RigidBody::getAtomVel(double theV[3], int index){
  vec3 ref;
  double velRot[3];
  double skewMat[3][3];
  double aSkewMat[3][3];
  double aSkewTransMat[3][3];
  
  //velRot = $(A\cdot skew(I^{-1}j))^{T}refCoor$

  if (index >= myAtoms.size())
    cerr << index << " is an invalid index, current rigid body contains " << myAtoms.size() << "atoms" << endl;

  ref = refCoords[index];

  skewMat[0][0] =0;
  skewMat[0][1] = ji[2] /I[2][2];
  skewMat[0][2] = -ji[1] /I[1][1];

  skewMat[1][0] = -ji[2] /I[2][2];
  skewMat[1][1] = 0;
  skewMat[1][2] = ji[0]/I[0][0];

  skewMat[2][0] =ji[1] /I[1][1];
  skewMat[2][1] = -ji[0]/I[0][0];
  skewMat[2][2] = 0;
  
  matMul3(A, skewMat, aSkewMat);

  transposeMat3(aSkewMat, aSkewTransMat);

  matVecMul3(aSkewTransMat, ref.vec, velRot);
  theV[0] = vel[0] + velRot[0];
  theV[1] = vel[1] + velRot[1];
  theV[2] = vel[2] + velRot[2];
}


