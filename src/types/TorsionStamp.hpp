#ifndef __TORSIONSTAMP_H__
#define __TORSIONSTAMP_H__

#include "LinkedAssign.hpp"

class TorsionStamp{

public:
  TorsionStamp();
  ~TorsionStamp();
  
  void assignString( char* lhs, char* rhs );
  void assignDouble( char* lhs, double rhs );
  void assignInt( char* lhs, int rhs );
  void members( int a, int b, int c, int d );
  void constrain( double constraint );
  char* checkMe( void );

  int getA( void ){ return a; }
  int getB( void ){ return b; }
  int getC( void ){ return c; }
  int getD( void ){ return d; }

  int haveExtras( void ) { return have_extras; }
  LinkedAssign* getExtras( void ) { return unhandled; }

private:

  int a, b, c, d; //the members
  double constraint;
  short int have_mbrs, have_constraint;

  LinkedAssign* unhandled; // the unhandled assignments
  short int have_extras;  
};

#endif
