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
