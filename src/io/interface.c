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

#include "io/parse_interface.h"
#include "io/BASS_interface.h"
#include "utils/simError.h"
#ifdef IS_MPI
#include "io/mpiBASS.h"
#endif

void interface_error( event* the_event );

void init_component( int comp_index ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = COMPONENT;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = comp_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif

  free( the_event );  
}

void init_molecule( int mol_index ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = MOLECULE;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = mol_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );  
}

void init_rigidbody( int rigidbody_index ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = RIGIDBODY;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = rigidbody_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );  
}

void init_cutoffgroup( int cutoffgroup_index ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = CUTOFFGROUP;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = cutoffgroup_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );  
}

void init_atom( int atom_index ){
  event* the_event;
  
  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = ATOM;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = atom_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}

void init_bond( int bond_index ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = BOND;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = bond_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}

void init_bend( int bend_index ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = BEND;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = bend_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}

void init_torsion( int torsion_index ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = TORSION;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = torsion_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}

void init_zconstraint( int zconstraint_index ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = ZCONSTRAINT;
  the_event->err_msg = NULL;
  the_event->evt.blk_index = zconstraint_index;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}


/*
 * the next few deal with the statement initializations 
 */

void init_members( struct node_tag* the_node, 
		   struct namespc the_namespc ){
  event* the_event;
  int i;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = MEMBERS;
  the_event->err_msg = NULL;

  the_event->evt.mbrs.nMembers = the_node->the_data.mbrs.nMembers;

  the_event->evt.mbrs.memberList = (int *) calloc(the_node->the_data.mbrs.nMembers, 
                                                  sizeof(int));

  for (i = 0; i < the_node->the_data.mbrs.nMembers; i++) {
    the_event->evt.mbrs.memberList[i] = the_node->the_data.mbrs.memberList[i];
  }
  
  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event->evt.mbrs.memberList );
  free( the_event );
}

void init_constraint( struct node_tag* the_node, 
		      struct namespc the_namespc ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = CONSTRAINT;
  the_event->err_msg = NULL;
  the_event->evt.cnstr = the_node->the_data.cnstr.constraint_val;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}

void init_assignment( struct node_tag* the_node, 
		      struct namespc the_namespc ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = ASSIGNMENT;
  the_event->err_msg = NULL;

  strcpy( the_event->evt.asmt.lhs, the_node->the_data.asmt.identifier );
  switch( the_node->the_data.asmt.type ){
    
  case STR_ASSN:
    the_event->evt.asmt.asmt_type = STRING;
    strcpy( the_event->evt.asmt.rhs.sval,
	    the_node->the_data.asmt.rhs.str_ptr );
    break;

  case INT_ASSN:
    the_event->evt.asmt.asmt_type = INT;
    the_event->evt.asmt.rhs.ival = the_node->the_data.asmt.rhs.i_val;
    break;

  case DOUBLE_ASSN:
    the_event->evt.asmt.asmt_type = DOUBLE;
    the_event->evt.asmt.rhs.dval = the_node->the_data.asmt.rhs.d_val;
    break;
  }

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}

void init_position( struct node_tag* the_node, 
		    struct namespc the_namespc ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = POSITION;
  the_event->err_msg = NULL;
  the_event->evt.pos.x = the_node->the_data.pos.x;
  the_event->evt.pos.y = the_node->the_data.pos.y;
  the_event->evt.pos.z = the_node->the_data.pos.z;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}

void init_orientation( struct node_tag* the_node, 
		       struct namespc the_namespc ){
  event* the_event;

  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = ORIENTATION;
  the_event->err_msg = NULL;
  the_event->evt.ornt.phi   = the_node->the_data.ort.phi;
  the_event->evt.ornt.theta = the_node->the_data.ort.theta;
  the_event->evt.ornt.psi   = the_node->the_data.ort.psi;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}


void end_of_block( void ){
  event* the_event;
  
  the_event = (event* )malloc( sizeof( event ) );
  
  the_event->event_type = BLOCK_END;
  the_event->err_msg = NULL;

  if( !event_handler( the_event ) ) interface_error( the_event );
#ifdef IS_MPI
  throwMPIEvent(the_event);
#endif
  
  free( the_event );
}

void interface_error( event* the_event ){

  sprintf( painCave.errMsg,
	   "Error in parsing meta-data file!\n"
	   "\t%s\n",
	   the_event->err_msg );
  painCave.severity = OOPSE_ERROR;
  painCave.isFatal = 1;
  simError();
#ifdef IS_MPI
  mpiInterfaceExit();
#endif //IS_MPI
}
