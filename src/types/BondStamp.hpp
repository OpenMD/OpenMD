#ifndef __BONDSTAMP_H__
#define __BONDSTAMP_H__

#include "io/LinkedAssign.hpp"

class BondStamp{

public:
  BondStamp();
  ~BondStamp();
  
  void assignString( char* lhs, char* rhs );
  void assignDouble( char* lhs, double rhs );
  void assignInt( char* lhs, int rhs );
  void members( int the_a, int the_b );
  void constrain( double the_constraint );
  char* checkMe( void );

  int getA( void ){ return a; }
  int getB( void ){ return b; }

  int haveExtras( void ) { return have_extras; }
  LinkedAssign* getExtras( void ) { return unhandled; }

private:

  int a, b; //the members
  double constraint;
  short int have_mbrs, have_constraint;

  LinkedAssign* unhandled; // the unhandled assignments
  short int have_extras;

};


#endif
