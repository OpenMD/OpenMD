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
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "types/CutoffGroupStamp.hpp"

char CutoffGroupStamp::errMsg[500];

CutoffGroupStamp::CutoffGroupStamp(){
  
  unhandled = NULL;
  have_members = 0;
  have_extras = 0;
  n_members = 0;
  which = 0;

}

CutoffGroupStamp::~CutoffGroupStamp(){
  if( unhandled != NULL ) delete unhandled;
  
  free(members);
 
}

char* CutoffGroupStamp::assignString( char* lhs, char* rhs ){

  if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
  else unhandled->add( lhs, rhs );
  have_extras = 1;
  return NULL;

}

char* CutoffGroupStamp::assignDouble( char* lhs, double rhs ){
  if( !strcmp( lhs, "nMembers" ) ){
    n_members = (int)rhs;
    
    if( have_members ){
      sprintf( errMsg,
               "CutoffGroupStamp error, nMembers already declared"
               " for this CutoffGroup.\n");
      return strdup( errMsg );
    }
    have_members = 1;
    members = (int *) calloc(n_members, sizeof(int));    
  } 
  else {
    if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
    else unhandled->add( lhs, rhs );
    have_extras = 1;
  }
  return NULL;
}

char* CutoffGroupStamp::assignInt( char* lhs, int rhs ){
  if( !strcmp( lhs, "nMembers" ) ){
    n_members = rhs;

    if( have_members ){
      sprintf( errMsg,
               "CutoffGroupStamp error, nMembers already declared for"
               " this CutoffGroup.\n");
      return strdup( errMsg );
    }
    have_members = 1;
    members = (int *) calloc(n_members, sizeof(int));    
  }
  else {  
    if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
    else unhandled->add( lhs, rhs );
    have_extras = 1;
  }
  return NULL;
}

char* CutoffGroupStamp::addMember( int atomIndex ){

  if( have_members && which < n_members ) {
    members[which] = atomIndex;
    which++;
  } else {
    if( have_members ){
      sprintf( errMsg, "CutoffGroupStamp error, %d out of nMembers range",
               which );
      return strdup( errMsg );
    }
    else return strdup("CutoffGroupStamp error, nMembers not given before"
                       " member list declaration." );
  }
  return NULL;
}

char* CutoffGroupStamp::checkMe( void ){

  short int no_member;
  
  if( !have_members ){
    return strdup( "CutoffGroupStamp error. CutoffGroup contains no members." );
  }
 
  if (which < n_members) {
    sprintf( errMsg,
             "CutoffGroupStamp error. Not all of the members were"
             " declared for this CutoffGroup.");
    return strdup( errMsg );
  }
  
  return NULL;
  
}
