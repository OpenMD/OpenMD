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
 
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "io/Globals.hpp"
#include "io/parse_me.h"
#include "utils/simError.h"

#ifdef IS_MPI
#include "io/mpiBASS.h"
#endif


// Globals ************************************************

typedef enum { GLOBAL_BLK, MOLECULE_BLK, ATOM_BLK, BOND_BLK, BEND_BLK,
	       TORSION_BLK, COMPONENT_BLK, ZCONSTRAINT_BLK, 
               RIGIDBODY_BLK, CUTOFFGROUP_BLK } block_type;

block_type current_block = GLOBAL_BLK;
#define MAX_NEST 20 // the max number of nested blocks
block_type block_stack[MAX_NEST];
int block_stack_ptr = 0;

void incr_block( block_type new_block );
void decr_block();

MakeStamps *the_stamps;
Globals *the_simParams;

// Functions **********************************************

/*
 * This is the main event handler
 *
 * returns 1 if the event was handled. If the event isn't handled, the
 * event's err_msg is set, and 0 is returned
 */

int event_handler( event* the_event ){

  int handled = 0;

  switch( current_block ){
    
  case GLOBAL_BLK:
    switch( the_event->event_type ){
      
    case MOLECULE:
      incr_block( MOLECULE_BLK );
      handled = the_stamps->newMolecule( the_event );
      break;

    case ZCONSTRAINT:
      incr_block( ZCONSTRAINT_BLK );
      handled = the_simParams->newZconstraint( the_event );
      break;
      
    case COMPONENT:
      incr_block( COMPONENT_BLK );
      handled = the_simParams->newComponent( the_event );
      break;

    case ASSIGNMENT:
      handled = the_simParams->globalAssign( the_event );
      break;

    case BLOCK_END:
      handled = the_simParams->globalEnd( the_event );
      break;
      
    default:
      the_event->err_msg = 
	strdup("Not a valid global event\n" );
      return 0;
    }
    break;

  case MOLECULE_BLK:
    
    switch( the_event->event_type ){
      
    case ATOM:
      incr_block( ATOM_BLK );
      handled = the_stamps->newAtom( the_event );
      break;
      
    case BOND:
      incr_block( BOND_BLK );
      handled = the_stamps->newBond( the_event );
      break;
      
    case BEND:
      incr_block( BEND_BLK );
      handled = the_stamps->newBend( the_event );
      break;

    case TORSION:
      incr_block( TORSION_BLK );
      handled = the_stamps->newTorsion( the_event );
      break;

    case RIGIDBODY:
      incr_block( RIGIDBODY_BLK );
      handled = the_stamps->newRigidBody( the_event );
      break;

    case CUTOFFGROUP:
      incr_block( CUTOFFGROUP_BLK );
      handled = the_stamps->newCutoffGroup( the_event );
      break;

    case ASSIGNMENT:
      handled = the_stamps->moleculeAssign( the_event );
      break;

    case BLOCK_END:
      decr_block();
      handled = the_stamps->moleculeEnd( the_event );
      break;

    default:
      the_event->err_msg = 
	strdup( "Not a valid molecule event\n" );
      return 0;
    }
    break;


  case RIGIDBODY_BLK:
    
    switch( the_event->event_type ){
      
    case ASSIGNMENT:
      handled = the_stamps->rigidBodyAssign( the_event );
      break;

    case MEMBERS:
      handled = the_stamps->rigidBodyMembers( the_event );
      break;
      
    case BLOCK_END:
      decr_block();
      handled = the_stamps->rigidBodyEnd( the_event );
      break;
      
    default:
      the_event->err_msg = 
	strdup( "Not a valid RigidBody event\n" );
      return 0;
    }
    break;

  case CUTOFFGROUP_BLK:
    
    switch( the_event->event_type ){
      
    case ASSIGNMENT:
      handled = the_stamps->cutoffGroupAssign( the_event );
      break;

    case MEMBERS:
      handled = the_stamps->cutoffGroupMembers( the_event );
      break;
      
    case BLOCK_END:
      decr_block();
      handled = the_stamps->cutoffGroupEnd( the_event );
      break;
      
    default:
      the_event->err_msg = 
	strdup( "Not a valid CutoffGroup event\n" );
      return 0;
    }
    break;
    

  case ATOM_BLK:
    
    switch( the_event->event_type ){
      
    case POSITION:
      handled = the_stamps->atomPosition( the_event );
      break;
      
    case ORIENTATION:
      handled = the_stamps->atomOrientation( the_event );
      break;

    case ASSIGNMENT:
      handled = the_stamps->atomAssign( the_event );
      break;
      
    case BLOCK_END:
      decr_block();
      handled = the_stamps->atomEnd( the_event );
      break;
      
    default:
      the_event->err_msg = 
	strdup( "Not a valid atom event\n" );
      return 0;
    }
    break;

  case BOND_BLK:

    switch( the_event->event_type ){
      
    case ASSIGNMENT:
      handled = the_stamps->bondAssign( the_event );
      break;
      
    case MEMBERS:
      handled = the_stamps->bondMembers( the_event );
      break;
      
    case CONSTRAINT:
      handled = the_stamps->bondConstraint( the_event );
      break;

    case BLOCK_END:
      decr_block();
      handled = the_stamps->bondEnd(the_event );
      break;

    default:
      the_event->err_msg = 
	strdup( "not a valid bond event\n" );
      return 0;
    }
    break;
    
  case BEND_BLK:
    
    switch( the_event->event_type ){
      
    case ASSIGNMENT:
      handled = the_stamps->bendAssign( the_event );
      break;
      
    case MEMBERS:
      handled = the_stamps->bendMembers( the_event );
      break;
      
    case CONSTRAINT:
      handled = the_stamps->bendConstraint( the_event );
      break;

    case BLOCK_END:
      decr_block();
      handled = the_stamps->bendEnd(the_event );
      break;

    default:
      the_event->err_msg = 
	strdup( "not a valid bend event\n" );
      return 0;
    }
    break;

  case TORSION_BLK:
    
    switch( the_event->event_type ){
      
    case ASSIGNMENT:
      handled = the_stamps->torsionAssign( the_event );
      break;
      
    case MEMBERS:
      handled = the_stamps->torsionMembers( the_event );
      break;
      
    case CONSTRAINT:
      handled = the_stamps->torsionConstraint( the_event );
      break;

    case BLOCK_END:
      decr_block();
      handled = the_stamps->torsionEnd(the_event );
      break;

    default:
      the_event->err_msg = 
	strdup( "not a valid torsion event\n" );
      return 0;
    }
    break;

  case ZCONSTRAINT_BLK:
 
    switch( the_event->event_type ){
      
    case ASSIGNMENT:
      handled = the_simParams->zConstraintAssign( the_event );
      break;
      
    case BLOCK_END:
      decr_block();
      handled = the_simParams->zConstraintEnd( the_event );
      break;

    default:
      the_event->err_msg = 
	strdup( "not a valid zConstraint event\n" );
      return 0;
    }
    break;
    
  case COMPONENT_BLK:
    
    switch( the_event->event_type ){
      
    case ASSIGNMENT:
      handled = the_simParams->componentAssign( the_event );
      break;
      
    case BLOCK_END:
      decr_block();
      handled = the_simParams->componentEnd(the_event );
      break;

    default:
      the_event->err_msg = 
	strdup( "not a valid component event\n" );
      return 0;
    }
   break;
   
  default:
    the_event->err_msg = 
      strdup( "event is not in a valid block type\n" );
    return 0;
  }
  
  return handled;
}

void incr_block( block_type new_block ){

  block_stack[block_stack_ptr] = current_block;
  block_stack_ptr++;
  
  if( block_stack_ptr >= MAX_NEST ){
    sprintf( painCave.errMsg,
	     "Event blocks nested too deeply\n" );
    painCave.isFatal = 1;
    simError();

#ifdef IS_MPI
    if( worldRank == 0 ) mpiInterfaceExit();
#endif //is_mpi
  }

  else current_block = new_block; 
}


void decr_block( void ){

  block_stack_ptr--;

  if( block_stack_ptr < 0 ){
    
    sprintf( painCave.errMsg,
	     "Too many event blocks closed\n" );
    painCave.isFatal = 1;
    simError();

#ifdef IS_MPI
    if( worldRank == 0 ) mpiInterfaceExit();
#endif //is_mpi

  }
  
  else current_block = block_stack[block_stack_ptr];
}

void set_interface_stamps( MakeStamps* ms, Globals* g ){

  the_stamps = ms;
  the_simParams = g;
  
}
