#ifndef __RIGIDBODY_HPP__
#define __RIGIDBODY_HPP__

#include <vector>
//#include "primitives/Atom.hpp"
//#include "types/AtomStamp.hpp"
#include "types/RigidBodyStamp.hpp"
#include "primitives/StuntDouble.hpp"
using namespace std;

class Atom;
class AtomStamp;

typedef struct {
  double vec[3];
  double& operator[](int index) {return vec[index];}  
} vec3;

typedef struct {
  double mat[3][3];
  double* operator[](int index) {return mat[index];}  
} mat3x3;

class RigidBody : public StuntDouble {

public:
  
  RigidBody();
  //RigidBody(const RigidBody& rb);
  
  virtual ~RigidBody();
 
  void addAtom(Atom* at, AtomStamp* ats);

  void getPos( double theP[3] );
  void setPos( double theP[3] );

  void getVel( double theV[3] );
  void setVel( double theV[3] );

  void getFrc( double theF[3] );
  void addFrc( double theF[3] );
  void zeroForces();
  
  bool isLinear() {return is_linear;}
  int linearAxis() {return linear_axis;}

  double getMass( void ) { return mass; }

  void printAmatIndex( void );
  void setEuler( double phi, double theta, double psi );
  void getQ( double the_q[4] ); // get the quanternions
  void setQ( double the_q[4] );

  void getA( double the_A[3][3] ); // get the full rotation matrix
  void setA( double the_A[3][3] );

  void getJ( double theJ[3] );
  void setJ( double theJ[3] );

  virtual void setType(char* type) {strcpy(rbName, type);}
  virtual char* getType() { return rbName;}

  void getTrq( double theT[3] );
  void addTrq( double theT[3] );

  void getI( double the_I[3][3] );
  void lab2Body( double r[3] );
  void body2Lab( double r[3] );

  double getZangle( );
  void setZangle( double zAng );
  void addZangle( double zAng );

  void calcRefCoords( void );
  void doEulerToRotMat(vec3 &euler, mat3x3 &myA );
  void calcForcesAndTorques( void );
  void updateAtoms( void );

  //void yourAtomsHaveMoved( void );

  // Four functions added for derivatives with respect to Euler Angles:
  // (Needed for minimization routines):

  void getGrad(double gradient[6] );
  void getEulerAngles( double myEuler[3] ); 
 
  double max(double x, double y);
  double min(double x, double y);


  // utility routines

  void findCOM( void );

  virtual void accept(BaseVisitor* v);

  vector<Atom*> getAtoms() { return myAtoms;}
  int getNumAtoms() {return myAtoms.size();}

  void getAtomPos(double theP[3], int index);
  void getAtomVel(double theV[3], int index);
  void getAtomRefCoor(double pos[3], int index);
protected:

  double mass;     // the total mass
  double pos[3];   // the position array (center of mass)
  double vel[3];   // the velocity array (center of mass)
  double frc[3];   // the force array    (center of mass)
  double trq[3];   // the torque vector  ( space fixed )
  double ji[3];    // the angular momentum vector (body fixed)
  double A[3][3];  // the rotation matrix 
  double I[3][3];  // the inertial tensor (body fixed)
  double sU[3][3]; // the standard unit vectors (body fixed)
  double zAngle;   // the rotation about the z-axis (body fixed)

  bool is_linear;
  int linear_axis;
  double momIntTol;

  vector<Atom*> myAtoms;  // the vector of atoms
  vector<vec3> refCoords;
  vector<mat3x3> refOrients;

  char rbName[100]; //it will eventually be converted into string
};

#endif
