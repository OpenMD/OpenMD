#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "types/RigidBodyStamp.hpp"

char RigidBodyStamp::errMsg[500];

RigidBodyStamp::RigidBodyStamp(){
  
  unhandled = NULL;
  have_members = 0;
  have_extras = 0;
  n_members = 0;
  which = 0;

}

RigidBodyStamp::~RigidBodyStamp(){
  int i;
  
  if( unhandled != NULL ) delete unhandled;
  
  free(members);
 
}

char* RigidBodyStamp::assignString( char* lhs, char* rhs ){

  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
  return NULL;

}

char* RigidBodyStamp::assignDouble( char* lhs, double rhs ){
  int i;

  if( !strcmp( lhs, "nMembers" ) ){
    n_members = (int)rhs;
    
    if( have_members ){
      sprintf( errMsg,
               "RigidBodyStamp error, nMembers already declared"
               " for this RigidBody.\n");
      return strdup( errMsg );
    }
    have_members = 1;
    members = (int *) calloc(n_members, sizeof(int));    
  } 
  else {
    if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
    else unhandled->add( lhs, rhs );
    have_extras = 1;
  }
  return NULL;
}

char* RigidBodyStamp::assignInt( char* lhs, int rhs ){
  int i;

  if( !strcmp( lhs, "nMembers" ) ){
    n_members = rhs;

    if( have_members ){
      sprintf( errMsg,
               "RigidBodyStamp error, nMembers already declared for"
               " this RigidBody.\n");
      return strdup( errMsg );
    }
    have_members = 1;
    members = (int *) calloc(n_members, sizeof(int));    
  }
  else {  
    if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
    else unhandled->add( lhs, rhs );
    have_extras = 1;
  }
  return NULL;
}

char* RigidBodyStamp::addMember( int atomIndex ){

  if( have_members && which < n_members ) {
    members[which] = atomIndex;
    which++;
  } else {
    if( have_members ){
      sprintf( errMsg, "RigidBodyStamp error, %d out of nMembers range",
               which );
      return strdup( errMsg );
    }
    else return strdup("RigidBodyStamp error, nMembers not given before"
                       " member list declaration." );
  }
  return NULL;
}

char* RigidBodyStamp::checkMe( void ){

  int i;
  short int no_member;
  
  if( !have_members ){
    return strdup( "RigidBodyStamp error. RigidBody contains no members." );
  }
 
  if (which < n_members) {
    sprintf( errMsg,
             "RigidBodyStamp error. Not all of the members were"
             " declared for this RigidBody.");
    return strdup( errMsg );
  }
  
  return NULL;
  
}
