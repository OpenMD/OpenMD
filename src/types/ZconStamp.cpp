#include <stdio.h>
#include <string.h>

#include "types/ZconStamp.hpp"


ZconStamp::ZconStamp(int theIndex){
  index = theIndex;
  
  have_zPos = 0;
  have_molIndex = 0;
  have_kRatio = 0;

  unhandled = NULL;
  have_extras = 0;
}

ZconStamp::~ZconStamp(){
  
  if( unhandled != NULL ) delete unhandled;
}

char* ZconStamp::checkMe( void ){
  char myErr[1000];
    
    //   if( !have_zPos ){
    //     sprintf( myErr,
    // 	     "No zPos was given to the Zconstraint[%d].",
    // 	     index );
    //     return strdup( myErr );
    //   }
  
  if( !have_molIndex ){
    sprintf( myErr,
	     "No Index was given to the Zconstraint[%d].",
	     index );
    return strdup( myErr );
  }
  
  return NULL;
}

int ZconStamp::assignString( char* lhs, char* rhs, char** err ){

  if( unhandled == NULL ){
    unhandled = new LinkedAssign( lhs, rhs );
    return 1;
  }
  else {
    unhandled->add( lhs, rhs );
    have_extras = 1;
    return 1;
  }

  return 0;
}

int ZconStamp::assignDouble( char* lhs, double rhs, char** err ){

  if( !strcasecmp( lhs, "zPos" )){
    
    zPos = rhs;
    have_zPos = 1;
    return 1;
  }
  else if( !strcasecmp( lhs, "kRatio" ) ){
    
    kRatio = rhs;
    have_kRatio = 1;
    return 1;
  }
  else if( !strcasecmp( lhs, "cantVel" )){
    
    cantVel = (double)rhs;
    have_cantVel = 1;
    return 1;
  }    
  else if( unhandled == NULL ){
    unhandled = new LinkedAssign( lhs, rhs );
    return 1;
  }
  else {
    unhandled->add( lhs, rhs );
    have_extras = 1;
    return 1;
  }

  return 0;
}

int ZconStamp::assignInt( char* lhs, int rhs, char** err ){

  if( !strcasecmp( lhs, "molIndex" ) ){
    
    molIndex = rhs;
    have_molIndex = 1;
    return 1;
  }
  else if( !strcasecmp( lhs, "kRatio" ) ){
    
    kRatio = (double)rhs;
    have_kRatio = 1;
    return 1;
  }
  else if( !strcasecmp( lhs, "zPos" )){
    
    zPos = (double)rhs;
    have_zPos = 1;
    return 1;
  }  
  else if( !strcasecmp( lhs, "cantVel" )){
    
    cantVel = (double)rhs;
    have_cantVel = 1;
    return 1;
  }  
  else if( unhandled == NULL ){
    unhandled = new LinkedAssign( lhs, rhs );
    return 1;
  }
  else {
    unhandled->add( lhs, rhs );
    have_extras = 1;
    return 1;
  }
  return 0;
}
