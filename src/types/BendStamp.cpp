#include <stdlib.h>
#include <string.h>

#include "BendStamp.hpp"

BendStamp::BendStamp(){
  
  have_mbrs = 0;
  have_constraint = 0;
  
  unhandled = NULL;
  have_extras = 0;
}

BendStamp::~BendStamp(){
  
  if( unhandled != NULL ) delete unhandled;
}

void BendStamp::members( int the_a, int the_b, int the_c ){
  
  a = the_a;
  b = the_b;
  c = the_c;
  have_mbrs = 1;
}

void BendStamp::constrain( double the_constraint ){
  
  constraint = the_constraint;
  have_constraint = 1;
}

void BendStamp::assignString( char* lhs, char* rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

void BendStamp::assignDouble( char* lhs, double rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

void BendStamp::assignInt( char* lhs, int rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

char* BendStamp::checkMe(){

  if( !have_mbrs ){
    return strdup( "BendStamp error. Bend was not given members." );
  }
  
  return NULL;
}
