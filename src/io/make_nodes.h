#ifndef __MAKE_NODES_H__
#define __MAKE_NODES_H__

#include "io/node_list.h"

/* walks to the top node of the current list */

extern struct node_tag* walk_to_top( struct node_tag* walk_me );

/* handles the assignment functions */

extern struct node_tag* assign_i( char * lhs, int rhs );
extern struct node_tag* assign_d( char * lhs, double rhs );
extern struct node_tag* assign_s( char * lhs, char* rhs );

/* handles the members functions */

extern struct node_tag* members( char * list_str );

/* handles the constraint funtion */

extern struct node_tag* constraint( char * list_str );

/* handles the orientation function */

extern struct node_tag* orientation( char * list_str );

/* handles the position function */

extern struct node_tag* position( char * list_str );

/* handles the various block modes */

extern struct node_tag* molecule_blk( struct node_tag* stmt_list );
extern struct node_tag* zconstraint_blk( int index, struct node_tag* stmt_list );
extern struct node_tag* rigidbody_blk( int index, struct node_tag* stmt_list );
extern struct node_tag* cutoffgroup_blk( int index, struct node_tag* stmt_list );
extern struct node_tag* atom_blk( int index, struct node_tag* stmt_list );
extern struct node_tag* bond_blk( int index, struct node_tag* stmt_list );
extern struct node_tag* bend_blk( int index, struct node_tag* stmt_list );
extern struct node_tag* torsion_blk( int index, struct node_tag* stmt_list );
extern struct node_tag* component_blk( struct node_tag* stmt_list );

int count_tokens(char *line, char *delimiters);
#endif
