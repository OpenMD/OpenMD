#ifndef __CUTOFFGROUPSTAMP_H__
#define __CUTOFFGROUPSTAMP_H__

#include "LinkedAssign.hpp"
#include "AtomStamp.hpp"


class CutoffGroupStamp{

public:
  CutoffGroupStamp();
  ~CutoffGroupStamp();

  char* assignString( char* lhs, char* rhs );
  char* assignDouble( char* lhs, double rhs );
  char* assignInt( char* lhs, int rhs );
  char* checkMe( void );

  char*      addMember( int atomIndex );
  int        getNMembers( void )    { return n_members; }
  int        getMember( int index ) { return members[index]; }
  
  int haveExtras( void ) { return have_extras; }
  LinkedAssign* getExtras( void ) { return unhandled; }

  static char errMsg[500];
private:

  int n_members;
  int which;
  short int have_members;
  
  int* members;

  LinkedAssign* unhandled; // the unhandled assignments
  short int have_extras;
};

#endif
