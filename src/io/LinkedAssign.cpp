 /*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */
 
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

