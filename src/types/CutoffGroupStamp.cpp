#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "types/CutoffGroupStamp.hpp"

char CutoffGroupStamp::errMsg[500];

CutoffGroupStamp::CutoffGroupStamp(){
  
  unhandled = NULL;
  have_members = 0;
  have_extras = 0;
  n_members = 0;
  which = 0;

}

CutoffGroupStamp::~CutoffGroupStamp(){
  int i;
  
  if( unhandled != NULL ) delete unhandled;
  
  free(members);
 
}

char* CutoffGroupStamp::assignString( char* lhs, char* rhs ){

  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
  return NULL;

}

char* CutoffGroupStamp::assignDouble( char* lhs, double rhs ){
  int i;

  if( !strcmp( lhs, "nMembers" ) ){
    n_members = (int)rhs;
    
    if( have_members ){
      sprintf( errMsg,
               "CutoffGroupStamp error, nMembers already declared"
               " for this CutoffGroup.\n");
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

char* CutoffGroupStamp::assignInt( char* lhs, int rhs ){
  int i;

  if( !strcmp( lhs, "nMembers" ) ){
    n_members = rhs;

    if( have_members ){
      sprintf( errMsg,
               "CutoffGroupStamp error, nMembers already declared for"
               " this CutoffGroup.\n");
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

char* CutoffGroupStamp::addMember( int atomIndex ){

  if( have_members && which < n_members ) {
    members[which] = atomIndex;
    which++;
  } else {
    if( have_members ){
      sprintf( errMsg, "CutoffGroupStamp error, %d out of nMembers range",
               which );
      return strdup( errMsg );
    }
    else return strdup("CutoffGroupStamp error, nMembers not given before"
                       " member list declaration." );
  }
  return NULL;
}

char* CutoffGroupStamp::checkMe( void ){

  int i;
  short int no_member;
  
  if( !have_members ){
    return strdup( "CutoffGroupStamp error. CutoffGroup contains no members." );
  }
 
  if (which < n_members) {
    sprintf( errMsg,
             "CutoffGroupStamp error. Not all of the members were"
             " declared for this CutoffGroup.");
    return strdup( errMsg );
  }
  
  return NULL;
  
}
