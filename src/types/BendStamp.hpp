#ifndef __BENDSTAMP_H__
#define __BENDSTAMP_H__

#include "io/LinkedAssign.hpp"

class BendStamp{

public:
  BendStamp();
  ~BendStamp();
  
  void assignString( char* lhs, char* rhs );
  void assignDouble( char* lhs, double rhs );
  void assignInt( char* lhs, int rhs );
  void members( int the_a, int the_b, int the_c );
  void constrain( double the_constraint );
  char* checkMe( void );

  int getA( void ){ return a; }
  int getB( void ){ return b; }
  int getC( void ){ return c; }

  int haveExtras( void ) { return have_extras; }
  LinkedAssign* getExtras( void ) { return unhandled; }
  
private:

  int a, b, c; //the members
  double constraint;
  short int have_mbrs, have_constraint;

  LinkedAssign* unhandled; // the unhandled assignments
  short int have_extras;

};

#endif
