#ifndef __INTERFACE_H__
#define __INTERFACE_H__

#include "node_list.h"

/*
 * the following is used to clarify the namespace of a given statement
 * or block node.
 */

struct namespc{
  int index;
  node_type type;
};

/*
 * The following are the function calls to initialize the different
 * code blocks.
 */

extern void init_component( int component_index );
extern void init_molecule( int mol_index );
extern void init_atom( int atom_index );
extern void init_bond( int bond_index );
extern void init_bend( int bend_index );
extern void init_torsion( int torsion_index );
extern void init_zconstraint( int zconstraint_index );
extern void init_rigidbody( int rigidbody_index );
extern void init_cutoffgroup( int cutoffgroup_index );


/*
 * the next few deal with the statement initializations 
 */

extern void init_members( struct node_tag* the_node, 
			  struct namespc the_namespc );
extern void init_constraint( struct node_tag* the_node, 
			     struct namespc the_namespc );
extern void init_assignment( struct node_tag* the_node, 
			     struct namespc the_namespc );
extern void init_position( struct node_tag* the_node, 
			   struct namespc the_namespc );
extern void init_orientation( struct node_tag* the_node, 
			      struct namespc the_namespc );
extern void init_start_index( struct node_tag* the_node, 
			      struct namespc the_namespc );
/*
 * This let's us know when we encounter the end of a block.
 */

extern void end_of_block( void );


#endif
