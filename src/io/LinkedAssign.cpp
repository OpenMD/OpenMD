#include <stdlib.h>
#include <string.h>

#include "io/LinkedAssign.hpp"

LinkedAssign::LinkedAssign( char* the_text, int the_ival ){
  
  next = NULL;
  strcpy( text, the_text );
  ival = the_ival;
  type = 0;
}

LinkedAssign::LinkedAssign( char* the_text, double the_dval ){
  
  next = NULL;
  strcpy( text, the_text );
  dval = the_dval;
  type = 1;
}

LinkedAssign::LinkedAssign( char* the_text, char* the_sval ){
  
  next = NULL;
  strcpy( text, the_text );
  strcpy( sval, the_sval );
  type = 2;
}

void LinkedAssign::setValues( char* the_text, int the_ival ){
  
  strcpy( text, the_text );
  ival = the_ival;
  type = 0;
}

void LinkedAssign::setValues( char* the_text, double the_dval ){
  
  strcpy( text, the_text );
  dval = the_dval;
  type = 1;
}

void LinkedAssign::setValues( char* the_text, char* the_sval ){
  
  strcpy( text, the_text );
  strcpy( sval, the_sval );
  type = 2;
}

void LinkedAssign::add( char* the_text, int the_ival ){

  if( next != NULL ) next->add( the_text, the_ival );
  
  else next = new LinkedAssign( the_text, the_ival );
}

void LinkedAssign::add( char* the_text, double the_dval ){

  if( next != NULL ) next->add( the_text, the_dval );
  
  else next = new LinkedAssign( the_text, the_dval );
}

void LinkedAssign::add( char* the_text, char* the_sval ){

  if( next != NULL ) next->add( the_text, the_sval );
  
  else next = new LinkedAssign( the_text, the_sval );
}

LinkedAssign* LinkedAssign::find( char* key ){

  if( !strcmp(key, text) ) return this;
  
  if( next != NULL ) return next->find(key);
  
  return NULL;
}

