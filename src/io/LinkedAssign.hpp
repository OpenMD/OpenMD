#ifndef __LINKEDASSIGN_H__
#define __LINKEDASSIGN_H__

class LinkedAssign{

public:
  LinkedAssign() { next = NULL; }
  LinkedAssign( char* the_text, int the_ival );
  LinkedAssign( char* the_text, double the_dval );
  LinkedAssign( char* the_text, char* the_sval );
  ~LinkedAssign(){ if( next != NULL ) delete next; }


  void setValues( char* the_text, int the_ival );
  void setValues( char* the_text, double the_dval );
  void setValues( char* the_text, char* the_sval );
  void add( char* the_text, int the_ival );
  void add( char* the_text, double the_dval );
  void add( char* the_text, char* the_sval );
  
  char* getlhs( void ) { return text; }
  short int getType( void ) { return type; }
  int getInt( void ) { return ival; }
  double getDouble( void ) { return dval; }
  char* getString( void ) { return sval; }
  void setNext( LinkedAssign* the_next ) { next = the_next; }
  LinkedAssign* getNext( void ) { return next; }

  LinkedAssign* find( char* key );

private:
  char text[100];
  int ival;
  double dval;
  char sval[100];
  short int type; // 0 = int, 1 = double, 2 = string
  LinkedAssign* next;
};

#endif
