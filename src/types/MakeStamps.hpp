#ifndef __MAKESTAMPS_H__
#define __MAKESTAMPS_H__

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
