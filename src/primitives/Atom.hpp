#ifndef _ATOM_H_
#define _ATOM_H_

#include <string.h>
#include <stdlib.h>
#include <iostream>

#include "brains/SimState.hpp"
#include "primitives/StuntDouble.hpp"
#include "visitors/BaseVisitor.hpp"

class Atom : public StuntDouble {
public:

  Atom(int theIndex, SimState* theConfig );
  virtual ~Atom() {} 

  virtual void setCoords(void);

  void getPos( double theP[3] );
  void setPos( double theP[3] );

  void getVel( double theV[3] );
  void setVel( double theV[3] );

  void getFrc( double theF[3] );
  void addFrc( double theF[3] );

  virtual void zeroForces();

  double getMass() {return c_mass;}
  void setMass(double mass) {c_mass = mass;}
  
  int getIndex() const {return index;}
  void setIndex(int theIndex);

  char *getType() {return c_name;}
  void setType(char * name) {strcpy(c_name,name);}
  
  int getIdent( void ) { return ident; }
  void setIdent( int info ) { ident = info; }

#ifdef IS_MPI
  int getGlobalIndex( void ) { return myGlobalIndex; }
  void setGlobalIndex( int info ) { myGlobalIndex = info; }
#endif // is_mpi

  void setHasDipole( int value ) { has_dipole = value; }
  int hasDipole( void ) { return has_dipole; }

  void setHasCharge(int value) {has_charge = value;}
  int hasCharge(void) {return has_charge;}


  virtual void accept(BaseVisitor* v) {v->visit(this);}
  
protected:
  
  SimState* myConfig;

  double* pos; // the position array
  double* vel; // the velocity array
  double* frc; // the forc array
  double* trq; // the torque vector  ( space fixed )
  double* Amat; // the rotation matrix
  double* mu;   // the array of dipole moments
  double* ul;   // the lab frame unit directional vector

  double zAngle; // the rotation about the z-axis ( body-fixed )

  double c_mass; /* the mass of the atom in amu */

  int index; /* set the atom's index */
  int offset; // the atom's offset in the storage array
  int offsetX, offsetY, offsetZ;

  int Axx, Axy, Axz; // the rotational matrix indices
  int Ayx, Ayy, Ayz;
  int Azx, Azy, Azz;

  char c_name[100]; /* it's name */
  int ident;  // it's unique numeric identity.
  
  int has_dipole; // dipole boolean
  int has_charge; // charge boolean

  bool hasCoords;

#ifdef IS_MPI
  int myGlobalIndex;
#endif
  
};

#endif
