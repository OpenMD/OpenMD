#ifndef __PARSE_TREE_H__
#define __PARSE_TREE_H__

#include "node_list.h"
#include "parse_interface.h"

/* 
   This is the function that takes the statement node tree and gives
   it to the interface to our program. 
*/

extern void pt_me( struct node_tag* head_node );

/* use the following to kill the node tree when done. */

extern void kill_tree( struct node_tag* the_node );

// function to print out a node causing errors 

extern void print_tree_error( struct node_tag* err_node, char* err_msg );


#endif
