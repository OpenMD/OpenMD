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
 
#ifndef TYPES_MAKESTAMPS_HPP
#define TYPES_MAKESTAMPS_HPP

#include <stdlib.h>
#include <string.h>

#include "io/BASS_interface.h"
#include "types/MoleculeStamp.hpp"
#include "types/AtomStamp.hpp"
#include "types/BondStamp.hpp"
#include "types/BendStamp.hpp"
#include "types/TorsionStamp.hpp"
#include "types/RigidBodyStamp.hpp"
#include "types/CutoffGroupStamp.hpp"

class LinkedMolStamp{

public:
  LinkedMolStamp(){ mol_stamp = NULL; next = NULL; prev = NULL; }
  ~LinkedMolStamp();
  
  MoleculeStamp* match( char* id );
  LinkedMolStamp* extract( char* id );
  void setStamp( MoleculeStamp* the_stamp ){ mol_stamp = the_stamp; }
  MoleculeStamp* getStamp(){ return mol_stamp; }
  void add( LinkedMolStamp* newbie );
  void setPrev( LinkedMolStamp* thePrev ){ prev = thePrev; }
  void setNext( LinkedMolStamp* theNext ){ next = theNext; }
  LinkedMolStamp* getNext() { return next; }
  
private:
  MoleculeStamp* mol_stamp;
  LinkedMolStamp* next;
  LinkedMolStamp* prev;
};

class MakeStamps{

public:
  MakeStamps();
  ~MakeStamps();

  int newMolecule( event* the_event );
  int moleculeAssign( event* the_event );
  int moleculeEnd( event* the_event );

  int newAtom( event* the_event );
  int atomPosition( event* the_event );
  int atomOrientation( event* the_event );
  int atomAssign( event* the_event );
  int atomEnd( event* the_event );

  int newRigidBody( event* the_event );
  int rigidBodyAssign( event* the_event );
  int rigidBodyMembers( event* the_event );
  int rigidBodyEnd( event* the_event );

  int newCutoffGroup( event* the_event );
  int cutoffGroupAssign( event* the_event );
  int cutoffGroupMembers( event* the_event );
  int cutoffGroupEnd( event* the_event );

  int newBond( event* the_event );
  int bondAssign( event* the_event );
  int bondMembers( event* the_event );
  int bondConstraint( event* the_event );
  int bondEnd( event* the_event );
  
  int newBend( event* the_event );
  int bendAssign( event* the_event );
  int bendMembers( event* the_event );
  int bendConstraint( event* the_event );
  int bendEnd( event* the_event );

  int newTorsion( event* the_event );
  int torsionAssign( event* the_event );
  int torsionMembers( event* the_event );
  int torsionConstraint( event* the_event );
  int torsionEnd( event* the_event );

  LinkedMolStamp* extractMolStamp( char* the_id );

private:

  int hash_size;
  int hash_shift;
  int hash( char* text );
  LinkedMolStamp** my_mols;
  void addMolStamp( MoleculeStamp* the_stamp );

  MoleculeStamp* current_mol;
  AtomStamp* current_atom;
  BondStamp* current_bond;
  BendStamp* current_bend;
  TorsionStamp* current_torsion;
  RigidBodyStamp* current_rigidbody;
  CutoffGroupStamp* current_cutoffgroup;

};

#endif
