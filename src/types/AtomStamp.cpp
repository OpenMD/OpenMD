#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "AtomStamp.hpp"

AtomStamp::AtomStamp(){
  
  unhandled = NULL;
  have_position = 0;
  have_orientation = 0;
  have_type = 0;
  have_extras = 0;
}

AtomStamp::~AtomStamp(){
  
  if( unhandled != NULL ) delete unhandled;
}

void AtomStamp::setPosition( double x, double y, double z ){

  pos[0] = x;
  pos[1] = y;
  pos[2] = z;

  have_position = 1;
}

void AtomStamp::setOrientation( double phi, double theta, double psi ){

  ornt[0] = phi;
  ornt[1] = theta;
  ornt[2] = psi;

  have_orientation = 1;
}

char* AtomStamp::assignString( char* lhs, char* rhs ){

  if( !strcmp( lhs, "type" ) ){
    strcpy( type, rhs );
    have_type = 1;
  }
  else{
    if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
    else unhandled->add( lhs, rhs );
    have_extras = 1;
  }
  return NULL;
}

char* AtomStamp::assignDouble( char* lhs, double rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
  return NULL;
}

char* AtomStamp::assignInt( char* lhs, int rhs ){
  
  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
  return NULL;
}

char* AtomStamp::checkMe( void ){
  
  if( !have_type ){
    return strdup( "AtomStamp error. Atom was untyped." );
  }
  return NULL;
}
