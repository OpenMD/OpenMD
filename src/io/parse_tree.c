#include <stdlib.h>
#include <stdio.h>

#include "parse_tree.h"
#include "simError.h"

#ifdef IS_MPI
#define __is_lex__
#include "mpiBASS.h"
#endif

void walk_down( struct node_tag* the_node, struct namespc the_namespc );
int mol_index; // keeps track of the number of molecules
int comp_index; // keeps track of the number of components.

/*
 * This is the parse tree function that is called by the yacc
 * routine. It passes the global node and namespace to the recursive
 * walk_down routine to explore the node tree.
 */

void pt_me( struct node_tag* head_node ){

  struct namespc global_namespc;
  
  if( head_node->type != GLOBAL_HEAD ){
    sprintf( painCave.errMsg, 
	     "Parse tree error: The head node was not the global node.\n" );
    painCave.isFatal = 1;
    simError();
#ifdef IS_MPI
    mpiInterfaceExit();
#endif //is_mpi
  }

  global_namespc.index = 0;
  global_namespc.type = GLOBAL_HEAD;

  mol_index = 0;
  comp_index = 0;
  walk_down( head_node->next_stmt, global_namespc ); 
  // closed global namespace and exit

}

/*
 * This is the main logic workhorse routine for the node tree
 * parser. It recursively walks down the node list and calls the
 * appropriate interface functions acording to the node types. It will
 * also set the appropriate namespace for the nodes.  
 *
 *       Note: the nodes only know about the namespace of their
 *       current level, the namespace above themselves is hidden.
 */

void walk_down( struct node_tag* the_node, struct namespc the_namespc ){
  
  struct namespc current_namespc;

  if( the_node != NULL ){
    
    if( the_node->type == GLOBAL_HEAD ){
      print_tree_error( the_node, "Too many global regions" );
    }

    if( the_node->stmt_list != NULL ){

      // the statement is a block node of some sort
      
      switch( the_node->type ){

      case COMPONENT_HEAD:
	if( the_namespc.type != GLOBAL_HEAD ){
	  print_tree_error( the_node, 
			    "component block is not in the global namespace" );
	}
	else{
	  init_component( comp_index );
	  current_namespc.index = comp_index;
	  current_namespc.type = COMPONENT_HEAD;
	  walk_down( the_node->stmt_list, current_namespc );
	  comp_index++;
	}
	break;

      case MOLECULE_HEAD:
	if( the_namespc.type != GLOBAL_HEAD ){
	  print_tree_error( the_node, 
			    "Molecule block is not in the global namespace" );
	}
	else{
	  init_molecule( mol_index );
	  current_namespc.index = mol_index;
	  current_namespc.type = MOLECULE_HEAD;
	  walk_down( the_node->stmt_list, current_namespc );
	  mol_index++;
	}
	break;

      case ATOM_HEAD:
	if( the_namespc.type != MOLECULE_HEAD ){
	     print_tree_error( the_node,
			    "The atom block is not in a molecule namespace" );
	}
	else{
	  init_atom( the_node->index );
	  current_namespc.index = the_node->index;
	  current_namespc.type = the_node->type;
	  walk_down( the_node->stmt_list, current_namespc );
	}
	break;

      case RIGIDBODY_HEAD:
	if( the_namespc.type != MOLECULE_HEAD ){
	  print_tree_error( the_node,
			    "The RigidBody block is not in a Molecule namespace" );
	}
	else{
	  init_rigidbody( the_node->index );
	  current_namespc.index = the_node->index;
	  current_namespc.type = the_node->type;
	  walk_down( the_node->stmt_list, current_namespc );
	}
	break;

      case CUTOFFGROUP_HEAD:
	if( the_namespc.type != MOLECULE_HEAD ){
	  print_tree_error( the_node,
			    "The CutoffGroup block is not in a Molecule namespace" );
	}
	else{
	  init_cutoffgroup( the_node->index );
	  current_namespc.index = the_node->index;
	  current_namespc.type = the_node->type;
	  walk_down( the_node->stmt_list, current_namespc );
	}
	break;
	
      case BOND_HEAD:
	if( the_namespc.type != MOLECULE_HEAD ){
	  print_tree_error( the_node,
			    "The bond block is not in a molecule namespace" );
	}
	else{
	  init_bond( the_node->index );
	  current_namespc.index = the_node->index;
	  current_namespc.type = the_node->type;
	  walk_down( the_node->stmt_list, current_namespc );
	}
	break;

      case BEND_HEAD:
	if( the_namespc.type != MOLECULE_HEAD ){
	  print_tree_error( the_node,
			    "The bend block is not in a molecule namespace" );
	}
	else{
	  init_bend( the_node->index );
	  current_namespc.index = the_node->index;
	  current_namespc.type = the_node->type;
	  walk_down( the_node->stmt_list, current_namespc );
	}
	break;

      case TORSION_HEAD:
	if( the_namespc.type != MOLECULE_HEAD ){
	  print_tree_error( the_node,
			    "The torsion block is not in "
			    "a molecule namespace" );
	}
	else{
	  init_torsion( the_node->index );
	  current_namespc.index = the_node->index;
	  current_namespc.type = the_node->type;
	  walk_down( the_node->stmt_list, current_namespc );
	}
	break;
      
      case ZCONSTRAINT_HEAD:
	if( the_namespc.type != GLOBAL_HEAD ){
	  print_tree_error( the_node,
			    "The Zconstraint block is not in "
			    "the global namespace" );
	}
	else{
	  init_zconstraint( the_node->index );
	  current_namespc.index = the_node->index;
	  current_namespc.type = the_node->type;
	  walk_down( the_node->stmt_list, current_namespc );
	}
	break;
	
      default:
	print_tree_error( the_node, "Not a valid code block" );
      }
    }
    
    else{
      
      // the node is a statement 

      switch( the_node->type ){

      case MEMBERS_STMT:
	switch( the_namespc.type ){
	case BOND_HEAD: // fall through
	case BEND_HEAD: // fall through
	case TORSION_HEAD: 
        case RIGIDBODY_HEAD:
        case CUTOFFGROUP_HEAD: // same for the first four
	  init_members( the_node, the_namespc );
	  break;

	default:
	  print_tree_error( the_node, 
			    "Members statement not in a bond, bend, "
			    "torsion, RigidBody, or CutoffGroup" );
	  break;
	}
	break;

      case CONSTRAINT_STMT:
	switch( the_namespc.type ){
	case BOND_HEAD: // fall through
	case BEND_HEAD: // fall through
	case TORSION_HEAD: // same for the first three
	  init_constraint( the_node, the_namespc );
	  break;

	default:
	  print_tree_error( the_node, 
			    "Constraint statement not in a bond, bend, "
			    "or torsion" );
	  break;
	}
	break;

      case ASSIGNMENT_STMT:
	init_assignment( the_node, the_namespc );
	break;

      case POSITION_STMT:
	if( the_namespc.type != ATOM_HEAD ){
	  print_tree_error( the_node,
			    "position statement is not located in an "
			    "atom block" );
	}
	
	init_position( the_node, the_namespc );
	break;

      case ORIENTATION_STMT:
	if( the_namespc.type != ATOM_HEAD ){
	  print_tree_error( the_node,
			    "orientation statement is not located in an "
			    "atom block" );
	}
	
	init_orientation( the_node, the_namespc );
	break;

      default:
	print_tree_error( the_node, "unrecognized statement" );
	break;
      }
    }
    
    // recurse down to the next node

    walk_down( the_node->next_stmt, the_namespc );
  }

  // send an end of block signal
  else end_of_block();
  
  // we're done
}



/*
 * This is a routine utilized by the node parsers to make printing
 * error messages easy.
 */

void print_tree_error( struct node_tag* err_node, char* err_msg ){

  switch( err_node->type ){

  case GLOBAL_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: global head node error -> %s\n",
	     err_msg );
    break;
    
  case COMPONENT_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: component head node error -> %s\n",
	     err_msg );
    break;

  case MOLECULE_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: molecule head node error -> %s\n",
	     err_msg );
    break;

  case RIGIDBODY_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: rigidBody head node error -> %s\n",
	     err_msg );
    break;
    
  case CUTOFFGROUP_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: CutoffGroup head node error -> %s\n",
	     err_msg );
    break;

  case ATOM_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: atom head node error [%d] -> %s\n",
	     err_node->index,
	     err_msg );
    break;
    
  case BOND_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: bond head node error [%d] -> %s\n",
	     err_node->index,
	     err_msg );
    break;
    
  case BEND_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: bend head node error [%d] -> %s\n",
	     err_node->index,
	     err_msg );
    break;
      
  case ZCONSTRAINT_HEAD:
    sprintf( painCave.errMsg,
	     "Parse tree error: Zconstraint head node error [%d] -> %s\n",
	     err_node->index,
	     err_msg );
    break;

  case MEMBERS_STMT:
    sprintf( painCave.errMsg,
	     "Parse tree error: members node error (nMembers = %d)\n"
	     "                  -> %s\n",
             err_node->the_data.mbrs.nMembers,
	     err_msg );
    break;

  case CONSTRAINT_STMT:
    sprintf( painCave.errMsg,
	    "Parse tree error: constraint node error => ( %lf )\n"
	    "                 -> %s\n",
	    err_node->the_data.cnstr.constraint_val,
	    err_msg );
    break;
    
  case ASSIGNMENT_STMT:
    sprintf( painCave.errMsg,
	     "Parse tree error: assignment node error\n"
	     "                  => %s = ",
	     err_node->the_data.asmt.identifier );
    
    switch( err_node->the_data.asmt.type ){
      
    case STR_ASSN:
      sprintf( painCave.errMsg,
	       "%s",
	       err_node->the_data.asmt.rhs.str_ptr );
      break;
      
    case INT_ASSN:
      sprintf( painCave.errMsg,
	       "%d",
	       err_node->the_data.asmt.rhs.i_val );
      break;
      
    case DOUBLE_ASSN:
      sprintf( painCave.errMsg,
	       "%lf",
	       err_node->the_data.asmt.rhs.d_val );
      break;
    }
    
    sprintf( painCave.errMsg,
	     "\n"
	     "                  -> %s\n",
	     err_msg );
    break;

  case POSITION_STMT:
    sprintf( painCave.errMsg,
	     "Parse tree error: position node error => ( %lf, %lf, %lf )\n"
	     "                  -> %s\n",
	     err_node->the_data.pos.x, 
	     err_node->the_data.pos.y, 
	     err_node->the_data.pos.z, 
	     err_msg );
    break;

  case ORIENTATION_STMT:
    sprintf( painCave.errMsg,
	     "Parse tree error: orientation node error => ( %lf, %lf, %lf )\n"
	     "                  -> %s\n",
	     err_node->the_data.ort.phi, 
	     err_node->the_data.ort.theta, 
	     err_node->the_data.ort.psi, 
	     err_msg );
    break;

  default:
    sprintf( painCave.errMsg,
	     "Parse tree error: unknown node type -> %s\n",
	     err_msg );
  }

  painCave.isFatal = 1;
  simError();
#ifdef IS_MPI
  mpiInterfaceExit();
#endif //is_mpi

}


/*
 * recursive walkdown and kill of the node tree
 * note: looks mighty similar to the walkdown routine.
 */ 
  
void kill_tree( struct node_tag* the_node ){
  

  if( the_node != NULL ){
    
    if( the_node->stmt_list != NULL ){

      // the statement is a block node of some sort
      
      kill_tree( the_node->stmt_list );
    }
    
    else{
      
      // the node is a statement 

      switch( the_node->type ){

      case ASSIGNMENT_STMT:
	
	if( the_node->the_data.asmt.type == STR_ASSN )
	  free( the_node->the_data.asmt.rhs.str_ptr );
	
	free( the_node->the_data.asmt.identifier );
	break;

      default:
	// nothing to do here, everyone else can be freed normally.
	break;
      }
    }
    
    // recurse down to the next node

    kill_tree( the_node->next_stmt );
    free( the_node );    
  }

  // we're done
}
