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
 
#ifndef TYPES_MOLECULESTAMP_HPP
#define TYPES_MOLECULESTAMP_HPP
#include <vector>
#include <utility>
#include "types/AtomStamp.hpp"
#include "types/BondStamp.hpp"
#include "types/BendStamp.hpp"
#include "types/TorsionStamp.hpp"
#include "types/RigidBodyStamp.hpp"
#include "types/CutoffGroupStamp.hpp"
#include "io/LinkedAssign.hpp"


class MoleculeStamp{

public:
  MoleculeStamp();
  ~MoleculeStamp();

  char* assignString( char* lhs, char* rhs );
  char* assignDouble( char* lhs, double rhs );
  char* assignInt( char* lhs, int rhs );
  char* checkMe( void );

  char* addAtom( AtomStamp* the_atom, int atomIndex );
  char* addRigidBody( RigidBodyStamp* the_rigidbody, int rigidBodyIndex );
  char* addCutoffGroup( CutoffGroupStamp* the_cutoffgroup, int cutoffGroupIndex );
  char* addBond( BondStamp* the_bond, int bondIndex );
  char* addBend( BendStamp* the_bend, int bendIndex );
  char* addTorsion( TorsionStamp* the_torsion, int torsionIndex );  

  char* getID( void )         { return name; }
  int   getNAtoms( void )     { return n_atoms; }
  int   getNBonds( void )     { return n_bonds; }
  int   getNBends( void )     { return n_bends; }
  int   getNTorsions( void )  { return n_torsions; }
  int   getNRigidBodies(void) { return n_rigidbodies; }
  int   getNCutoffGroups(void){ return n_cutoffgroups; }  
  int   getNIntegrable(void)  { return n_integrable; }

  AtomStamp* getAtom( int index ) { return atoms[index]; }
  BondStamp* getBond( int index ) { return bonds[index]; }
  BendStamp* getBend( int index ) { return bends[index]; }
  TorsionStamp* getTorsion( int index ) { return torsions[index]; }
  RigidBodyStamp* getRigidBody( int index ) { return rigidBodies[index]; }
  CutoffGroupStamp* getCutoffGroup( int index ) { return cutoffGroups[index]; }


  bool isBondInSameRigidBody(BondStamp*bond);
  bool isAtomInRigidBody(int atomIndex);  
  bool isAtomInRigidBody(int atomIndex, int& whichRigidBody, int& consAtomIndex);  
   std::vector<std::pair<int, int> > getJointAtoms(int rb1, int rb2);
  
  int haveExtras( void ) { return have_extras; }
  LinkedAssign* getUnhandled( void ) { return unhandled; }
  
  static char errMsg[500];
private:
  
  
  char name[100];
  int n_atoms;
  int n_bonds;
  int n_bends;
  int n_torsions;
  int n_rigidbodies;
  int n_cutoffgroups;
  int n_integrable;
  
  int have_name, have_atoms, have_bonds, have_bends, have_torsions;
  int have_rigidbodies, have_cutoffgroups;

  AtomStamp** atoms;
  BondStamp** bonds;
  BendStamp** bends;
  TorsionStamp** torsions;  
  RigidBodyStamp** rigidBodies;  
  CutoffGroupStamp** cutoffGroups;  

  LinkedAssign* unhandled; // the unhandled assignments
  short int have_extras;
};

#endif
