#ifndef __ZCONSTAMP_H__
#define __ZCONSTAMP_H__

#include "io/LinkedAssign.hpp"

class ZconStamp{

public:
  ZconStamp(int theIndex);
  ~ZconStamp();

  int assignString( char* lhs, char* rhs, char** err );
  int assignDouble( char* lhs, double rhs, char** err );
  int assignInt( char* lhs, int rhs, char** err );

  int haveExtras( void ) { return have_extras; }
  LinkedAssign* getExtras( void ) { return unhandled; }

  char* checkMe( void );

  int    getMolIndex( void ) { return molIndex; }
  double getZpos( void )     { return zPos; }
  double getKratio( void )   { return kRatio; }
  double getCantVel( void )   { return cantVel; }
  
  short int haveZpos( void ) { return have_zPos; }
  short int haveKratio( void ) { return have_kRatio; }
  short int haveCantVel( void ) { return have_cantVel; }
  
private:

  short int have_zPos;
  short int have_molIndex;
  short int have_kRatio;
  short int have_cantVel;

  double zPos;
  double kRatio;
  double cantVel;
  int molIndex;
  
  int index;

  LinkedAssign* unhandled; // the unhandled assignments
  short int have_extras;  
};


#endif // __ZCONSTAMP_H__
