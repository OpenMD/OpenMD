#include <stdlib.h>
#include <string.h>

#include "types/BondStamp.hpp"

BondStamp::BondStamp(){
  
  have_mbrs = 0;
  have_constraint = 0;
  
  unhandled = NULL;
  have_extras = 0;
}

BondStamp::~BondStamp(){
  
  if( unhandled != NULL ) delete unhandled;
}

void BondStamp::members( int the_a, int the_b ){
  
  a = the_a;
  b = the_b;
  have_mbrs = 1;
}

void BondStamp::constrain( double the_constraint ){
  
  constraint = the_constraint;
  have_constraint = 1;
}

void BondStamp::assignString( char* lhs, char* rhs ){

  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

void BondStamp::assignDouble( char* lhs, double rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

void BondStamp::assignInt( char* lhs, int rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

char* BondStamp::checkMe(){

  if( !have_mbrs ){
    return strdup( "BondStamp error. Bond was not given members." );
  }
  
  return NULL;
}
