#ifndef _DIRECTIONALATOM_H_
#define _DIRECTIONALATOM_H_

#include <string.h>
#include <stdlib.h>
#include <iostream>

#include "primitives/StuntDouble.hpp"
#include "primitives/Atom.hpp"

class DirectionalAtom : public Atom {
  
public:
  DirectionalAtom(int theIndex, SimState* theConfig) : Atom(theIndex, 
							    theConfig)
  { 
    objType = OT_DATOM;

    for (int i=0; i < 3; i++) 
      for (int j=0; j < 3; j++)
        sU[i][j] = 0.0;

    is_linear = false;
    linear_axis =  -1;
    momIntTol = 1e-6;
  }
  virtual ~DirectionalAtom() {}

  virtual void setCoords(void);

  void printAmatIndex( void );
  
  void setUnitFrameFromEuler(double phi, double theta, double psi);
  void setEuler( double phi, double theta, double psi );

  void zeroForces();

  void getA( double the_A[3][3] ); // get the full rotation matrix
  void setA( double the_A[3][3] );
  void rotateBy( double by_A[3][3] );  // rotate your frame using this matrix

  void getU( double the_u[3] ); // get the unit vetor
  void updateU( void );

  void getQ( double the_q[4] ); // get the quanternions
  void setQ( double the_q[4] );
 
  void getJ( double theJ[3] );
  void setJ( double theJ[3] );

  void getTrq( double theT[3] );
  void addTrq( double theT[3] );

  void setI( double the_I[3][3] );
  void getI( double the_I[3][3] );

  bool isLinear() {return is_linear;}
  int linearAxis() {return linear_axis;}
  
  void lab2Body( double r[3] );
  void body2Lab( double r[3] );

  double getZangle( );
  void setZangle( double zAng );
  void addZangle( double zAng );

  // Four functions added for derivatives with respect to Euler Angles:
  // (Needed for minimization routines):

  void getGrad(double gradient[6] );
  void getEulerAngles( double myEuler[3] );

  double max(double x, double y);
  double min(double x, double y);

  virtual void accept(BaseVisitor* v) {v->visit(this);}
  
private:
  int dIndex;

  double sU[3][3];       // the standard unit vectors    ( body fixed )
  
  double jx, jy, jz;    // the angular momentum vector ( body fixed )
  
  double Ixx, Ixy, Ixz; // the inertial tensor matrix  ( body fixed )
  double Iyx, Iyy, Iyz;
  double Izx, Izy, Izz;

  bool is_linear;
  int linear_axis;
  double momIntTol;

};

#endif
