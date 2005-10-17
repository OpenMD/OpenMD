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
 
#ifndef IO_BASS_INTERFACE_H
#define IO_BASS_INTERFACE_H



typedef enum { MOLECULE, ATOM, BOND, BEND, TORSION, COMPONENT, 
	       POSITION, ASSIGNMENT, MEMBERS, CONSTRAINT, ORIENTATION,
	       ZCONSTRAINT, RIGIDBODY, CUTOFFGROUP, BLOCK_END } event_enum;


typedef struct{
  double x;
  double y;
  double z;
} position_event;

typedef struct{
  double phi;
  double theta;
  double psi;
} orientation_event;

typedef enum { STRING, INT, DOUBLE } interface_assign_type;

typedef struct{
  interface_assign_type asmt_type;
  char lhs[80];
  union{
    int ival;
    double dval;
    char sval[1024];
  }rhs;
} assignment_event;

typedef struct{
  int nMembers;
  int *memberList;
} members_event;

typedef struct{
  event_enum event_type;
  char* err_msg;

  union{
    int               blk_index; /* block index*/
    position_event    pos;
    orientation_event ornt; /* use the same structure for orientation*/
    assignment_event  asmt;
    members_event     mbrs;
    double            cnstr; /* the constraint value*/
  } evt;
} event;

#ifdef __cplusplus
extern "C" {
#endif

  int event_handler( event* the_event );

#ifdef __cplusplus
}
#endif


#endif /* ifndef __BASS_INTERFACE_H__*/
