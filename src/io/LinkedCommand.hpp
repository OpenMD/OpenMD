#ifndef __LINKEDCOMMAND_H__
#define __LINKEDCOMMAND_H__

#include <stdlib.h>

class LinkedCommand{

public:
  
  LinkedCommand(){ next = NULL; token = 0; }
  ~LinkedCommand(){ if( next != NULL ) delete next; }
  
  int match( char* the_text );
  void setValues( char* the_text, int the_token );
  void setNext( LinkedCommand* the_next ){ next = the_next; }
  LinkedCommand* getNext( void ){ return next; }
  
private:
  
  LinkedCommand* next;
  char text[100];
  int token;
};

#endif
