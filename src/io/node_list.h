#ifndef __NODE_LIST_H__
#define __NODE_LIST_H__

/* enumerates the different types of nodes the statments can be */

typedef enum { GLOBAL_HEAD, COMPONENT_HEAD,
	       MOLECULE_HEAD, ATOM_HEAD, BOND_HEAD, BEND_HEAD, TORSION_HEAD,
	       MEMBERS_STMT, CONSTRAINT_STMT, ASSIGNMENT_STMT, POSITION_STMT,
	       ORIENTATION_STMT, ZCONSTRAINT_HEAD, 
               RIGIDBODY_HEAD, CUTOFFGROUP_HEAD } node_type;

/* a structure to hold the index of the members of a bond, bend, or torsion */

typedef struct members_data_tag{
  int nMembers;
  int* memberList;
} members_data;

/* a structure to hold constraint information */

typedef struct constraint_data_tag{
  double constraint_val;
} constraint_data;

/* a structure to hold assignment info */

typedef enum{ STR_ASSN, INT_ASSN, DOUBLE_ASSN } assign_type;

typedef union{
  int i_val;
  double d_val;
  char* str_ptr;
} assignment_value;

typedef struct assignment_data_tag{
  assign_type type;
  char* identifier; // left hand side
  assignment_value rhs; // right hand side
} assignment_data;

 /* a structure to hold the position information */

typedef struct position_data_tag{
  double x;
  double y;
  double z;
} position_data;

 /* a structure to hold the orientation information */

typedef struct orientation_data_tag{
  double phi;
  double theta;
  double psi;
} orientation_data;

/* here's the master node type. This is the node that will be strung
   together by the yacc parser. Each statement will an individual
   node. Block statements will act as branch nodes, pointing to their
   own set of statement lists.*/

typedef struct node_tag{
  node_type type;
  int index; // needed for atoms, bonds, etc.
  struct node_tag* next_stmt; // the next statement
  struct node_tag* prev_stmt; // the previous statement
  struct node_tag* stmt_list; // the statment list if this is a block.
  
  /* this is a union to hold the statement data */

  union{ 
    members_data mbrs;
    constraint_data cnstr;
    assignment_data asmt;
    position_data pos;
    orientation_data ort;
  } the_data;
    
} node;



#endif
