#ifndef __MOLECULESTAMP_H__
#define __MOLECULESTAMP_H__
#include <vector>
#include <utility>
#include "AtomStamp.hpp"
#include "BondStamp.hpp"
#include "BendStamp.hpp"
#include "TorsionStamp.hpp"
#include "RigidBodyStamp.hpp"
#include "CutoffGroupStamp.hpp"
#include "LinkedAssign.hpp"

using namespace std;
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
  vector<pair<int, int> > getJointAtoms(int rb1, int rb2);
  
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
