#ifndef __ATOMSTAMP_H__
#define __ATOMSTAMP_H__

#include "LinkedAssign.hpp"

class AtomStamp{
  
public:
  AtomStamp();
  ~AtomStamp();

  void setPosition( double x, double y, double z );
  void setOrientation( double phi, double theta, double psi );
  char* assignString( char* lhs, char* rhs );
  char* assignDouble( char* lhs, double rhs );
  char* assignInt( char* lhs, int rhs );
  char* checkMe( void );

  char* getType( void ) { return type; }
  short int havePosition( void ) { return have_position; }
  short int haveOrientation( void ) { return have_orientation; }
  double getPosX( void ) { return pos[0]; }
  double getPosY( void ) { return pos[1]; }
  double getPosZ( void ) { return pos[2]; }
  double getEulerPhi( void )   { return ornt[0]; }
  double getEulerTheta( void ) { return ornt[1]; }
  double getEulerPsi( void )   { return ornt[2]; }
  

private:

  double pos[3]; //the position vector
  short int have_position; // boolean for positions
  double ornt[3]; // the Euler angles
  short int have_orientation;
  char type[100]; // the type name of the atom
  short int have_type;
  
  LinkedAssign* unhandled; // the list of unhandled assignments
  short int have_extras;
};

#endif
