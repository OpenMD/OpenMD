#ifndef __RIGIDBODYSTAMP_H__
#define __RIGIDBODYSTAMP_H__

#include "io/LinkedAssign.hpp"
#include "types/AtomStamp.hpp"


class RigidBodyStamp{

public:
  RigidBodyStamp();
  ~RigidBodyStamp();

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
