
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "Globals.hpp"
#include "BASS_interface.h"
#include "simError.h"

#ifdef IS_MPI
#include "mpiBASS.h"
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
Globals *the_globals;

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
      handled = the_globals->newZconstraint( the_event );
      break;
      
    case COMPONENT:
      incr_block( COMPONENT_BLK );
      handled = the_globals->newComponent( the_event );
      break;

    case ASSIGNMENT:
      handled = the_globals->globalAssign( the_event );
      break;

    case BLOCK_END:
      handled = the_globals->globalEnd( the_event );
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
      handled = the_globals->zConstraintAssign( the_event );
      break;
      
    case BLOCK_END:
      decr_block();
      handled = the_globals->zConstraintEnd( the_event );
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
      handled = the_globals->componentAssign( the_event );
      break;
      
    case BLOCK_END:
      decr_block();
      handled = the_globals->componentEnd(the_event );
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
  the_globals = g;
  
}
