#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#include <set>
#include <vector>

#include "primitives/Atom.hpp"
#include "primitives/SRI.hpp"
#include "types/MoleculeStamp.hpp"
#include "primitives/RigidBody.hpp"
#include "primitives/CutoffGroup.hpp"

using namespace std;

typedef struct{
  
  int stampID;   // the ID in the BASS component stamp array
  int nAtoms;    // the number of atoms in the molecule
  int nBonds;    // ... .. ..  . .bonds .. .. . . . . 
  int nBends;    // . . . . .. . .bends . . . . .. . 
  int nTorsions; // .. . . .. . . torsions . . .. . . 
  int nRigidBodies; // .. .. .. . rigid bodies ... .. 
  int nOriented; // .. . . . .. . oriented atoms . . . 
  
  Atom** myAtoms;      // the array of atoms
  Bond** myBonds;      // arrays of all the short range interactions
  Bend** myBends;
  Torsion** myTorsions;
  vector<RigidBody*>   myRigidBodies;
  vector<StuntDouble*> myIntegrableObjects;
  vector<CutoffGroup*> myCutoffGroups;
} molInit;

class Molecule{

public:
  
  Molecule( void );
  ~Molecule( void );

  void initialize( molInit &theInit );

  void setMyIndex( int theIndex ){ myIndex = theIndex;}
  int getMyIndex( void ) { return myIndex; }

  int getGlobalIndex( void ) { return globalIndex; }
  void setGlobalIndex( int theIndex ) { globalIndex = theIndex; }

  int getNAtoms   ( void )    {return nAtoms;}
  int getNBonds   ( void )    {return nBonds;}
  int getNBends   ( void )    {return nBends;}
  int getNTorsions( void )    {return nTorsions;}
  int getNRigidBodies( void ) {return myRigidBodies.size();}
  int getNOriented( void )    {return nOriented;}
  int getNMembers ( void )    {return nMembers;}
  int getStampID  ( void )    {return stampID;}

  Atom**      getMyAtoms   ( void )    {return myAtoms;}
  Bond**      getMyBonds   ( void )    {return myBonds;}
  Bend**      getMyBends   ( void )    {return myBends;}
  Torsion**   getMyTorsions( void )    {return myTorsions;}
  vector<RigidBody*> getMyRigidBodies( void ) {return myRigidBodies;}
  vector<StuntDouble*>& getIntegrableObjects(void) {return myIntegrableObjects;}

  //beginCutoffGroup return the first group and initialize the iterator
  CutoffGroup* beginCutoffGroup(vector<CutoffGroup*>::iterator& i){
    i = myCutoffGroups.begin();
    return i != myCutoffGroups.end()? *i : NULL;
  }

  //nextCutoffGroup return next cutoff group based on the iterator
  CutoffGroup* nextCutoffGroup(vector<CutoffGroup*>::iterator& i){
    i++;
    return i != myCutoffGroups.end()? *i : NULL;
  }

  int getNCutoffGroups() {return nCutoffGroups;}

  void setStampID( int info ) {stampID = info;} 

  void calcForces( void );
  void atoms2rigidBodies( void );
  double getPotential( void );
  
  void printMe( void );

  void getCOM( double COM[3] );
  void moveCOM( double delta[3] );
  double getCOMvel( double COMvel[3] );
  
  double getTotalMass();

private:

  int stampID;   // the ID in the BASS component stamp array
  int nAtoms;    // the number of atoms in the molecule
  int nBonds;    // ... .. ..  . .bonds .. .. . . . . 
  int nBends;    // . . . . .. . .bends . . . . .. . 
  int nTorsions; // .. . . .. . . torsions . . .. . . 
  int nRigidBodies; // .. . . .. .rigid bodies . . .. . . 
  int nOriented; // .. . . . .. . oriented atoms . . . 
  int nMembers;  // .. . . . . . .atoms (legacy code) . . . 
  int nCutoffGroups;
  
  int myIndex; // mostly just for debug (and for making pressure calcs work)
  int globalIndex;

  Atom** myAtoms;     // the array of atoms
  Bond** myBonds;     // arrays of all the short range interactions
  Bend** myBends;
  Torsion** myTorsions;
  vector<RigidBody*>   myRigidBodies;
  vector<StuntDouble*> myIntegrableObjects;
  vector<CutoffGroup*> myCutoffGroups;
};

#endif
