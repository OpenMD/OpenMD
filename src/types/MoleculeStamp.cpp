#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "types/MoleculeStamp.hpp"

char MoleculeStamp::errMsg[500];

MoleculeStamp::MoleculeStamp(){
  
  n_atoms = 0;
  n_bonds = 0;
  n_bends = 0;
  n_torsions = 0;
  n_rigidbodies = 0;
  n_cutoffgroups = 0;
  n_integrable = 0;

  unhandled = NULL;
  atoms = NULL;
  bonds = NULL;
  bends = NULL;
  torsions = NULL;
  rigidBodies = NULL;
  cutoffGroups = NULL;

  have_name = 0;
  have_atoms = 0;
  have_bonds = 0;
  have_bends = 0;
  have_torsions = 0;
  have_rigidbodies = 0;
  have_cutoffgroups = 0;

}

MoleculeStamp::~MoleculeStamp(){
  int i;
  
  if( unhandled != NULL) delete unhandled;
  
  if( rigidBodies != NULL ) {
    for( i=0; i<n_rigidbodies; i++ ) delete rigidBodies[i];
  }
  delete[] rigidBodies;

  if( cutoffGroups != NULL ) {
    for( i=0; i<n_cutoffgroups; i++ ) delete cutoffGroups[i];
  }
  delete[] cutoffGroups;
  
  if( atoms != NULL ){
    for( i=0; i<n_atoms; i++ ) delete atoms[i];
  }
  delete[] atoms;
  
  if( bonds != NULL ){
    for( i=0; i<n_bonds; i++ ) delete bonds[i];
  }
  delete[] bonds;
  
  if( bends != NULL ){
    for( i=0; i<n_bends; i++ ) delete bends[i];
  }
  delete[] bends;
  
  if( torsions != NULL ){
    for( i=0; i<n_torsions; i++ ) delete torsions[i];
  }
  delete[] torsions;
  



}

char* MoleculeStamp::assignString( char* lhs, char* rhs ){
  
  if( !strcmp( lhs, "name" ) ){
    strcpy( name, rhs );
    have_name = 1;
  }
  else{
    if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
    else unhandled->add( lhs, rhs );
    have_extras = 1;
  }
  return NULL;
}

char* MoleculeStamp::assignDouble( char* lhs, double rhs ){
  int i;

  if( !strcmp( lhs, "nAtoms" ) ){
    n_atoms = (int)rhs;
    
    if( have_atoms ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_atoms already declared"
	       " for molecule: %s\n",
	       name);
      return strdup( errMsg );
    }
    have_atoms = 1;
    atoms = new AtomStamp*[n_atoms];
    for( i=0; i<n_atoms; i++ ) atoms[i] = NULL;
  }
  
  else if( !strcmp( lhs, "nBonds" ) ){
    n_bonds = (int)rhs;

    if( have_bonds ){
      sprintf( errMsg,  
	       "MoleculeStamp error, n_bonds already declared for"
	       " molecule: %s\n",
	       name);
      return strdup( errMsg );
    }
    have_bonds = 1;
    bonds = new BondStamp*[n_bonds];
    for( i=0; i<n_bonds; i++ ) bonds[i] = NULL;
  }
  
  else if( !strcmp( lhs, "nBends" ) ){
    n_bends = (int)rhs;

    if( have_bends ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_bends already declared for"
	       " molecule: %s\n",
	       name);
      return strdup( errMsg );
    }
    have_bends = 1;
    bends = new BendStamp*[n_bends];
    for( i=0; i<n_bends; i++ ) bends[i] = NULL;
  }
  
  else if( !strcmp( lhs, "nTorsions" ) ){
    n_torsions = (int)rhs;

    if( have_torsions ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_torsions already declared for"
	       " molecule: %s\n",
	       name );
      return strdup( errMsg );
    }
    have_torsions = 1;
    torsions = new TorsionStamp*[n_torsions];
    for( i=0; i<n_torsions; i++ ) torsions[i] = NULL;
  }

  else if( !strcmp( lhs, "nRigidBodies" ) ){
    n_rigidbodies = (int)rhs;

    if( have_rigidbodies ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_rigidbodies already declared for"
	       " molecule: %s\n",
	       name );
      return strdup( errMsg );
    }
    have_rigidbodies = 1;
    rigidBodies = new RigidBodyStamp*[n_rigidbodies];
    for( i=0; i<n_rigidbodies; i++ ) rigidBodies[i] = NULL;
  }

  else if( !strcmp( lhs, "nCutoffGroups" ) ){
    n_cutoffgroups = (int)rhs;

    if( have_cutoffgroups ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_cutoffgroups already declared for"
	       " molecule: %s\n",
	       name );
      return strdup( errMsg );
    }
    have_cutoffgroups = 1;
    cutoffGroups = new CutoffGroupStamp*[n_cutoffgroups];
    for( i=0; i<n_cutoffgroups; i++ ) cutoffGroups[i] = NULL;
  }
  
  else{
    if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
    else unhandled->add( lhs, rhs );
    have_extras = 1;
  }
  return NULL;
}

char*  MoleculeStamp::assignInt( char* lhs, int rhs ){
  int i;
  
  if( !strcmp( lhs, "nAtoms" ) ){
    n_atoms = rhs;
    
    if( have_atoms ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_atoms already declared for"
	       " molecule: %s\n",
	       name);
      return strdup( errMsg );
    }
    have_atoms = 1;
    atoms = new AtomStamp*[n_atoms];
    for( i=0; i<n_atoms; i++ ) atoms[i] = NULL;
  }
  
  else if( !strcmp( lhs, "nBonds" ) ){
    n_bonds = rhs;

    if( have_bonds ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_bonds already declared for"
	       " molecule: %s\n",
	       name);
      return strdup( errMsg );
    }
    have_bonds = 1;
    bonds = new BondStamp*[n_bonds];
    for( i=0; i<n_bonds; i++ ) bonds[i] = NULL;
  }
  
  else if( !strcmp( lhs, "nBends" ) ){
    n_bends = rhs;

    if( have_bends ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_bends already declared for"
	       " molecule: %s\n",
	       name );
      return strdup( errMsg );
    }
    have_bends = 1;
    bends = new BendStamp*[n_bends];
    for( i=0; i<n_bends; i++ ) bends[i] = NULL;
  }
  
  else if( !strcmp( lhs, "nTorsions" ) ){
    n_torsions = rhs;

    if( have_torsions ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_torsions already declared for"
	       " molecule: %s\n",
	       name);
      return strdup( errMsg );
    }
    have_torsions = 1;
    torsions = new TorsionStamp*[n_torsions];
    for( i=0; i<n_torsions; i++ ) torsions[i] = NULL;
  }

  else if( !strcmp( lhs, "nRigidBodies" ) ){
    n_rigidbodies = rhs;

    if( have_rigidbodies ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_rigidbodies already declared for"
	       " molecule: %s\n",
	       name);
      return strdup( errMsg );
    }
    have_rigidbodies = 1;
    rigidBodies = new RigidBodyStamp*[n_rigidbodies];
    for( i=0; i<n_rigidbodies; i++ ) rigidBodies[i] = NULL;
  }
  else if( !strcmp( lhs, "nCutoffGroups" ) ){
    n_cutoffgroups = rhs;

    if( have_cutoffgroups ){
      sprintf( errMsg,
	       "MoleculeStamp error, n_cutoffgroups already declared for"
	       " molecule: %s\n",
	       name);
      return strdup( errMsg );
    }
    have_cutoffgroups = 1;
    cutoffGroups = new CutoffGroupStamp*[n_cutoffgroups];
    for( i=0; i<n_cutoffgroups; i++ ) cutoffGroups[i] = NULL;
  }
  else{
    if( unhandled == NULL ) unhandled = new LinkedAssign( lhs, rhs );
    else unhandled->add( lhs, rhs );
    have_extras = 1;
  }
  return NULL;
}

char* MoleculeStamp::addAtom( AtomStamp* the_atom, int atomIndex ){
  
  if( have_atoms && atomIndex < n_atoms ) atoms[atomIndex] = the_atom;
  else {
    if( have_atoms ){
      sprintf( errMsg, "MoleculeStamp error, %d out of nAtoms range", 
	       atomIndex );
      return strdup( errMsg );
    }
    else return strdup("MoleculeStamp error, nAtoms not given before"
		       " first atom declaration." );
  }

  return NULL;
}

char* MoleculeStamp::addRigidBody( RigidBodyStamp* the_rigidbody, 
                                   int rigidBodyIndex ){
  
  if( have_rigidbodies && rigidBodyIndex < n_rigidbodies ) 
    rigidBodies[rigidBodyIndex] = the_rigidbody;
  else {
    if( have_rigidbodies ){
      sprintf( errMsg, "MoleculeStamp error, %d out of nRigidBodies range", 
	       rigidBodyIndex );
      return strdup( errMsg );
    }
    else return strdup("MoleculeStamp error, nRigidBodies not given before"
		       " first rigidBody declaration." );
  }
  
  return NULL;
}

char* MoleculeStamp::addCutoffGroup( CutoffGroupStamp* the_cutoffgroup, 
                                     int cutoffGroupIndex ){
  
  if( have_cutoffgroups && cutoffGroupIndex < n_cutoffgroups ) 
    cutoffGroups[cutoffGroupIndex] = the_cutoffgroup;
  else {
    if( have_cutoffgroups ){
      sprintf( errMsg, "MoleculeStamp error, %d out of nCutoffGroups range", 
	       cutoffGroupIndex );
      return strdup( errMsg );
    }
    else return strdup("MoleculeStamp error, nCutoffGroups not given before"
		       " first CutoffGroup declaration." );
  }
  
  return NULL;
}

char* MoleculeStamp::addBond( BondStamp* the_bond, int bondIndex ){
  
  
  if( have_bonds && bondIndex < n_bonds ) bonds[bondIndex] = the_bond;
  else{
    if( have_bonds ){
      sprintf( errMsg, "MoleculeStamp error, %d out of nBonds range", 
	       bondIndex );
      return strdup( errMsg );
    }
    else return strdup("MoleculeStamp error, nBonds not given before"
		       "first bond declaration." );
  }

  return NULL;
}

char* MoleculeStamp::addBend( BendStamp* the_bend, int bendIndex ){
  
  
  if( have_bends && bendIndex < n_bends ) bends[bendIndex] = the_bend;
  else{
    if( have_bends ){
      sprintf( errMsg, "MoleculeStamp error, %d out of nBends range", 
	       bendIndex );
      return strdup( errMsg );
    }
    else return strdup("MoleculeStamp error, nBends not given before"
		       "first bend declaration." );
  }

  return NULL;
}

char* MoleculeStamp::addTorsion( TorsionStamp* the_torsion, int torsionIndex ){
  
  
  if( have_torsions && torsionIndex < n_torsions ) 
    torsions[torsionIndex] = the_torsion;
  else{
    if( have_torsions ){
      sprintf( errMsg, "MoleculeStamp error, %d out of nTorsions range", 
	       torsionIndex );
      return strdup( errMsg );
    }
    else return strdup("MoleculeStamp error, nTorsions not given before"
		       "first torsion declaration." );
  }

  return NULL;
}

char* MoleculeStamp::checkMe( void ){
  
  int i;
  short int no_atom, no_rigidbody, no_cutoffgroup;

  if( !have_name ) return strdup( "MoleculeStamp error. Molecule's name"
                                  " was not given.\n" );
  
  if( !have_atoms ){
    return strdup( "MoleculeStamp error. Molecule contains no atoms." );
  }
  
  no_rigidbody = 0;
  for( i=0; i<n_rigidbodies; i++ ){
    if( rigidBodies[i] == NULL ) no_rigidbody = 1;
  }

  if( no_rigidbody ){
    sprintf( errMsg, 
	     "MoleculeStamp error. Not all of the RigidBodies were"
	     " declared in molecule \"%s\".\n", name );
    return strdup( errMsg );
  }

  no_cutoffgroup = 0;
  for( i=0; i<n_cutoffgroups; i++ ){
    if( cutoffGroups[i] == NULL ) no_cutoffgroup = 1;
  }

  if( no_cutoffgroup ){
    sprintf( errMsg, 
	     "MoleculeStamp error. Not all of the CutoffGroups were"
	     " declared in molecule \"%s\".\n", name );
    return strdup( errMsg );
  }
  
  no_atom = 0;
  for( i=0; i<n_atoms; i++ ){
    if( atoms[i] == NULL ) no_atom = 1;
  }

  if( no_atom ){
    sprintf( errMsg, 
	     "MoleculeStamp error. Not all of the atoms were"
	     " declared in molecule \"%s\".\n", name );
    return strdup( errMsg );
  }

  n_integrable = n_atoms;
  for (i = 0; i < n_rigidbodies; i++) 
    n_integrable = n_integrable - rigidBodies[i]->getNMembers() + 1; //rigidbody is an integrable object
  
  if (n_integrable <= 0 || n_integrable > n_atoms) {
    sprintf( errMsg, 
	     "MoleculeStamp error. n_integrable is either <= 0 or"
	     " greater than n_atoms in molecule \"%s\".\n", name );
    return strdup( errMsg );
  }

  return NULL;
}  


//Function Name: isBondInSameRigidBody
//Return true is both atoms of the bond belong to the same rigid body, otherwise return false
bool MoleculeStamp::isBondInSameRigidBody(BondStamp* bond){
  int rbA;
  int rbB;
  int consAtomA;
  int consAtomB;

  if (!isAtomInRigidBody(bond->getA(),rbA, consAtomA))
    return false;

  if(!isAtomInRigidBody(bond->getB(),rbB, consAtomB) )
    return false;

  if(rbB == rbA)
    return true;
  else
    return false;
}

// Function Name: isAtomInRigidBody 
//return false if atom does not belong to a rigid body, otherwise return true 
bool MoleculeStamp::isAtomInRigidBody(int atomIndex){
  int whichRigidBody;
  int consAtomIndex;

  return isAtomInRigidBody(atomIndex, whichRigidBody, consAtomIndex);
   
}

// Function Name: isAtomInRigidBody 
//return false if atom does not belong to a rigid body otherwise return true and set whichRigidBody 
//and consAtomIndex
//atomIndex : the index of atom in component
//whichRigidBody: the index of rigidbody in component
//consAtomIndex:  the position of joint atom apears in  rigidbody's definition
bool MoleculeStamp::isAtomInRigidBody(int atomIndex, int& whichRigidBody, int& consAtomIndex){
  RigidBodyStamp* rbStamp;
  int numRb;
  int numAtom;

  whichRigidBody = -1;
  consAtomIndex = -1;

  numRb = this->getNRigidBodies();
  
  for(int i = 0 ; i < numRb; i++){
    rbStamp = this->getRigidBody(i);
    numAtom = rbStamp->getNMembers();
    for(int j = 0; j < numAtom; j++)
      if (rbStamp->getMember(j) == atomIndex){
        whichRigidBody = i;
        consAtomIndex = j;
        return true;
      }
  }

  return false;
   
}

//return the position of joint atom apears in  rigidbody's definition
//for the time being, we will use the most inefficient algorithm, the complexity is O(N2)
//actually we could improve the complexity to O(NlgN) by sorting the atom index in rigid body first
vector<pair<int, int> > MoleculeStamp::getJointAtoms(int rb1, int rb2){
  RigidBodyStamp* rbStamp1;
  RigidBodyStamp* rbStamp2;
  int natomInRb1;
  int natomInRb2;
  int atomIndex1;
  int atomIndex2;
  vector<pair<int, int> > jointAtomIndexPair;
  
  rbStamp1 = this->getRigidBody(rb1);
  natomInRb1 =rbStamp1->getNMembers();

  rbStamp2 = this->getRigidBody(rb2);
  natomInRb2 =rbStamp2->getNMembers();

  for(int i = 0; i < natomInRb1; i++){
    atomIndex1 = rbStamp1->getMember(i);
      
    for(int j= 0; j < natomInRb1; j++){
      atomIndex2 = rbStamp2->getMember(j);

      if(atomIndex1 == atomIndex2){
        jointAtomIndexPair.push_back(make_pair(i, j));
        break;
      }
      
    }//end for(j =0)

  }//end for (i = 0)

  return jointAtomIndexPair;
}
