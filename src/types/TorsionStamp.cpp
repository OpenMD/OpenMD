#include <stdlib.h>
#include <string.h>

#include "TorsionStamp.hpp"

TorsionStamp::TorsionStamp(){
  
  have_mbrs = 0;
  have_constraint = 0;
  
  unhandled = NULL;
  have_extras = 0;
}

TorsionStamp::~TorsionStamp(){
  
  if( unhandled != NULL ) delete unhandled;
}

void TorsionStamp::members( int the_a, int the_b, int the_c, int the_d ){
  
  a = the_a;
  b = the_b;
  c = the_c;
  d = the_d;
  have_mbrs = 1;
}

void TorsionStamp::constrain( double the_constraint ){
  
  constraint = the_constraint;
  have_constraint = 1;
}

void TorsionStamp::assignString( char* lhs, char* rhs ){

  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

void TorsionStamp::assignDouble( char* lhs, double rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

void TorsionStamp::assignInt( char* lhs, int rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
}

char* TorsionStamp::checkMe(){

  if( !have_mbrs ){
    return strdup( "TorsionStamp error. Torsion was not given members." );
  }
  
  return NULL;
}
