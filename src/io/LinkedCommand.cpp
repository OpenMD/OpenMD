#include <stdlib.h>
#include <string.h>

#include "LinkedCommand.hpp"

int LinkedCommand::match( char* the_text ){
  
  if( !strcmp( the_text, text ) ) return token; //search successful
  else if( next != NULL ) return next->match( the_text ); // search the next
  
  return 0; // search failed
}

void LinkedCommand::setValues( char* the_text, int the_token ){

  strcpy( text, the_text );
  token = the_token;
}
