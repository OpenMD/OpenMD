#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "Component.hpp"

Component::Component(){
  
  start_array = NULL;
  have_type = 0;
  have_nMol = 0;
  have_molFraction = 0;
  have_start_array = 0;
}

Component::~Component(){
  
  if( start_array != NULL ) free( start_array );
}

int Component::assignString( char* lhs, char* rhs, char** err ){
  
  char myErr[1000];


  if( !strcmp( lhs, "type" ) ){
    strcpy( type, rhs );
    have_type = 1;
  }

  else{
    sprintf(myErr, "Invalid assignment type for Component: %s = %s\n", lhs, rhs);
    *err = strdup(myErr);
    return 0;
  }

  return 1;
}

int Component::assignInt( char* lhs, int rhs, char** err ){
  
  char myErr[1000];

  if( !strcmp( lhs, "nMol" ) ){
    nMol = rhs;
    have_nMol = 1;
  }

  else if( !strcmp( lhs, "molFraction" ) ){
    if( rhs > 1 || rhs < 0 ){
      sprintf(myErr,"Component error. %d is an invalid molFraction. It must lie between 0 and 1\n",rhs);
      *err = strdup(myErr);
      return 0;
    }
    molFraction = rhs;
    have_molFraction = 1;
  }

  else{
    sprintf(myErr, "Invalid assignment type for Component: %s = %d\n", lhs, rhs);
    *err = strdup(myErr);
    return 0;
  }

  return 1;
}

int Component::assignDouble( char* lhs, double rhs, char** err ){

  char myErr[1000];

  if( !strcmp( lhs, "molFraction" ) ){
    if( rhs > 1 || rhs < 0 ){
      sprintf(myErr,"Component error. %lf is an invalid molFraction. It must lie between 0 and 1\n",rhs);
      *err = strdup(myErr);
      return 0;
    }
    molFraction = rhs;
    have_molFraction = 1;
  }

  else if( !strcmp( lhs, "nMol" ) ){
    nMol = (int) rhs;
    have_nMol = 1;
  }
  
  else{
    sprintf(myErr, "Invalid assignment type for Component: %s = %lf\n", lhs, rhs);
    *err = strdup(myErr);
    return 0;
  }

  return 1;
}

void Component::startIndex( int* the_start_array, int n_elements ){
  
  start_array = the_start_array;
  n_start = n_elements;
  have_start_array = 1;
}

char* Component::checkMe( void ){
  
  if( !have_type ){
    return strdup( "type was not identified for the Component." );
  }

  return NULL;
}
