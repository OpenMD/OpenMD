#ifndef __COMPONENT_H__
#define __COMPONENT_H__

class Component{
  
public:
  Component();
  ~Component();
  
  int assignString( char* lhs, char* rhs, char** err );
  int assignDouble( char* lhs, double rhs, char** err );
  int assignInt( char* lhs, int rhs, char** err );
  void startIndex( int* the_start_array, int n_elements );
  char* checkMe( void );

  int   getNMol( void ) { return nMol; }
  char* getType( void ) { return type; }
  
  short int haveNMol( void ) { return have_nMol; }

private:

  
  char type[100];
  int nMol;
  double molFraction;
  int* start_array;
  int n_start;

  short int have_type, have_nMol, have_molFraction, have_start_array;
};

#endif
