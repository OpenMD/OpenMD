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
 
#ifndef IO_NODELIST_H
#define IO_NODELIST_H

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
