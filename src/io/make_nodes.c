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
#include <stdlib.h>
#include <string.h>

#include "io/node_list.h"
#include "io/make_nodes.h"
#include "utils/simError.h"

/*
  walks to to the top of a stmt_list. Needed when assigning the
  statment list to a block of code.
*/

struct node_tag* walk_to_top( struct node_tag* walk_me ){

  while( walk_me->prev_stmt != NULL ){
    walk_me = walk_me->prev_stmt;
  }
  return walk_me;
}

/*
  the next three functions are for creating and initializing the
  assignment statements.
*/

struct node_tag* assign_i( char* lhs, int rhs ){

  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = ASSIGNMENT_STMT;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = NULL;
  
  the_node->the_data.asmt.type = INT_ASSN;
  the_node->the_data.asmt.identifier = lhs;
  the_node->the_data.asmt.rhs.i_val = rhs;

  return the_node;
}

struct node_tag* assign_d( char* lhs, double rhs ){

  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = ASSIGNMENT_STMT;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = NULL;
  
  the_node->the_data.asmt.type = DOUBLE_ASSN;
  the_node->the_data.asmt.identifier = lhs;
  the_node->the_data.asmt.rhs.d_val = rhs;

  return the_node;
}

struct node_tag* assign_s( char* lhs, char* rhs ){

  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = ASSIGNMENT_STMT;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = NULL;
  
  the_node->the_data.asmt.type = STR_ASSN;
  the_node->the_data.asmt.identifier = lhs;
  the_node->the_data.asmt.rhs.str_ptr = rhs;

  return the_node;
}

/*
  The next function deals with creating and initializing of the
  members statements used by the bonds, bends, torsions, rigid bodies, 
  and cutoff groups.
*/

struct node_tag* members( char* list_str ){

  int nTokens, i;
  char *foo;

  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = MEMBERS_STMT;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = NULL;
  
  nTokens = count_tokens(list_str, " ,;\t");

  if (nTokens < 1) {
    yyerror("can't find any members in this members statement");
  }

  the_node->the_data.mbrs.nMembers = nTokens;

  the_node->the_data.mbrs.memberList = (int *) calloc(nTokens, sizeof(int));

  foo = strtok(list_str, " ,;\t");
  the_node->the_data.mbrs.memberList[0] = atoi( foo );

  for (i = 1; i < nTokens; i++) {
    foo = strtok(NULL, " ,;\t");
    the_node->the_data.mbrs.memberList[i] = atoi( foo );
  }

  free(list_str);
  return the_node;
}

/*
  This function does the C & I of the constraint statements
*/

struct node_tag* constraint( char* list_str ){

  int nTokens;
  char *foo;

  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = CONSTRAINT_STMT;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = NULL;
  
  nTokens = count_tokens(list_str, " ,;\t");
  if (nTokens != 1) {
    yyerror("wrong number of arguments given in constraint statement");
  }
      
  foo = strtok(list_str, " ,;\t");
  the_node->the_data.cnstr.constraint_val = atof( foo );

  free(list_str);
  return the_node;
}

/*
  This function does the C & I for the orientation statements
*/

struct node_tag* orientation( char* list_str ){

  int nTokens;
  char *foo;

  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = ORIENTATION_STMT;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = NULL;
 
  nTokens = count_tokens(list_str, " ,;\t");
  if (nTokens != 3) {
    yyerror("wrong number of arguments given in orientation");
  }
      
  foo = strtok(list_str, " ,;\t");
  the_node->the_data.ort.phi = atof( foo );

  foo = strtok(NULL, " ,;\t");
  the_node->the_data.ort.theta = atof( foo );

  foo = strtok(NULL, " ,;\t");
  the_node->the_data.ort.psi = atof( foo );

  free(list_str);
  return the_node;
}

/*
  This function does the C & I for the position statements
*/

struct node_tag* position( char* list_str ){

  int nTokens;
  char *foo;

  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = POSITION_STMT;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = NULL;

  nTokens = count_tokens(list_str, " ,;\t");
  if (nTokens != 3) {
    yyerror("wrong number of arguments given in position");
  }
      
  foo = strtok(list_str, " ,;\t");
  the_node->the_data.pos.x = atof( foo );

  foo = strtok(NULL, " ,;\t");
  the_node->the_data.pos.y = atof( foo );

  foo = strtok(NULL, " ,;\t");
  the_node->the_data.pos.z = atof( foo );

  free(list_str);
  return the_node;
}

/*
  The following six functions initialize the statement nodes
  coresponding to the different code block types.
*/

struct node_tag* molecule_blk( struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = MOLECULE_HEAD;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
    
  return the_node;
}

struct node_tag* rigidbody_blk( int index, struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = RIGIDBODY_HEAD;
  the_node->index = index;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
    
  return the_node;
}

struct node_tag* cutoffgroup_blk( int index, struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  // The guillotine statement:
  the_node->type = CUTOFFGROUP_HEAD;
  the_node->index = index;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
    
  return the_node;
}

struct node_tag* atom_blk( int index, struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = ATOM_HEAD;
  the_node->index = index;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
    
  return the_node;
}

struct node_tag* bond_blk( int index, struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = BOND_HEAD;
  the_node->index = index;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
    
  return the_node;
}
 
struct node_tag* bend_blk( int index, struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = BEND_HEAD;
  the_node->index = index;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
    
  return the_node;
}

struct node_tag* torsion_blk( int index, struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = TORSION_HEAD;
  the_node->index = index;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
    
  return the_node;
}

struct node_tag* zconstraint_blk( int index, struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = ZCONSTRAINT_HEAD;
  the_node->index = index;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
    
  return the_node;
} 
 
struct node_tag* component_blk( struct node_tag* stmt_list ){
  
  struct node_tag* the_node;
  the_node = ( struct node_tag* )malloc( sizeof( node ) );
  
  the_node->type = COMPONENT_HEAD;
  the_node->index = 0;
  the_node->next_stmt = NULL;
  the_node->prev_stmt = NULL;
  the_node->stmt_list = walk_to_top( stmt_list );
  
  return the_node;
}


int count_tokens(char *line, char *delimiters) {
  /* PURPOSE: RETURN A COUNT OF THE NUMBER OF TOKENS ON THE LINE. */

  char *working_line;   /* WORKING COPY OF LINE. */
  int ntokens;          /* NUMBER OF TOKENS FOUND IN LINE. */
  char *strtok_ptr;     /* POINTER FOR STRTOK. */
  
  strtok_ptr= working_line= strdup(line);
  
  ntokens=0;
  while (strtok(strtok_ptr,delimiters)!=NULL)
    {
      ntokens++;
      strtok_ptr=NULL;
    }
  
  free(working_line);
  return(ntokens);
}
