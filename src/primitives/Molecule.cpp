#include <stdlib.h>


#include "primitives/Molecule.hpp"
#include "utils/simError.h"



Molecule::Molecule( void ){

  myAtoms = NULL;
  myBonds = NULL;
  myBends = NULL;
  myTorsions = NULL;
}

Molecule::~Molecule( void ){
  int i;
  CutoffGroup* cg;
  vector<CutoffGroup*>::iterator iter;
  
  if( myAtoms != NULL ){
    for(i=0; i<nAtoms; i++) if(myAtoms[i] != NULL ) delete myAtoms[i];
    delete[] myAtoms;
  }

  if( myBonds != NULL ){
    for(i=0; i<nBonds; i++) if(myBonds[i] != NULL ) delete myBonds[i];
    delete[] myBonds;
  }

  if( myBends != NULL ){
    for(i=0; i<nBends; i++) if(myBends[i] != NULL ) delete myBends[i];
    delete[] myBends;
  }

  if( myTorsions != NULL ){
    for(i=0; i<nTorsions; i++) if(myTorsions[i] != NULL ) delete myTorsions[i];
    delete[] myTorsions;
  }

  for(cg = beginCutoffGroup(iter);  cg != NULL; cg = nextCutoffGroup(iter))
    delete cg;
  myCutoffGroups.clear();
  
}


void Molecule::initialize( molInit &theInit ){

  CutoffGroup* curCutoffGroup;
  vector<CutoffGroup*>::iterator iterCutoff;
  Atom* cutoffAtom;
  vector<Atom*>::iterator iterAtom;
  int atomIndex;
  
  nAtoms = theInit.nAtoms;
  nMembers = nAtoms;
  nBonds = theInit.nBonds;
  nBends = theInit.nBends;
  nTorsions = theInit.nTorsions;
  nRigidBodies = theInit.nRigidBodies;
  nOriented = theInit.nOriented;

  myAtoms = theInit.myAtoms;
  myBonds = theInit.myBonds;
  myBends = theInit.myBends;
  myTorsions = theInit.myTorsions;
  myRigidBodies = theInit.myRigidBodies;

  myIntegrableObjects = theInit.myIntegrableObjects;

  for (int i = 0; i < myRigidBodies.size(); i++) 
      myRigidBodies[i]->calcRefCoords();

  myCutoffGroups = theInit.myCutoffGroups;
  nCutoffGroups = myCutoffGroups.size();

}

void Molecule::calcForces( void ){
  
  int i;
  double com[3];

  for(i=0; i<myRigidBodies.size(); i++) {
    myRigidBodies[i]->updateAtoms();
  }

  for(i=0; i<nBonds; i++){
    myBonds[i]->calc_forces();
  }

  for(i=0; i<nBends; i++){
    myBends[i]->calc_forces();
  }

  for(i=0; i<nTorsions; i++){
    myTorsions[i]->calc_forces();
  }

  // Rigid Body forces and torques are done after the fortran force loop

}


double Molecule::getPotential( void ){
  
  int i;
  double myPot = 0.0;

  for(i=0; i<myRigidBodies.size(); i++) {
    myRigidBodies[i]->updateAtoms();
  }
   
  for(i=0; i<nBonds; i++){
    myPot += myBonds[i]->get_potential();
  }

  for(i=0; i<nBends; i++){
    myPot += myBends[i]->get_potential();
  }

  for(i=0; i<nTorsions; i++){
    myPot += myTorsions[i]->get_potential();
  }

  return myPot;
}

void Molecule::printMe( void ){
  
  int i;

  for(i=0; i<nBonds; i++){
    myBonds[i]->printMe();
  }

  for(i=0; i<nBends; i++){
    myBends[i]->printMe();
  }

  for(i=0; i<nTorsions; i++){
    myTorsions[i]->printMe();
  }

}

void Molecule::moveCOM(double delta[3]){
  double aPos[3];
  int i, j;

  for(i=0; i<myIntegrableObjects.size(); i++) {
    if(myIntegrableObjects[i] != NULL ) {
      
      myIntegrableObjects[i]->getPos( aPos );
      
      for (j=0; j< 3; j++) 
        aPos[j] += delta[j];

      myIntegrableObjects[i]->setPos( aPos );
    }
  }

  for(i=0; i<myRigidBodies.size(); i++) {

      myRigidBodies[i]->getPos( aPos );

      for (j=0; j< 3; j++) 
        aPos[j] += delta[j];
      
      myRigidBodies[i]->setPos( aPos );
    }
}

void Molecule::atoms2rigidBodies( void ) {
  int i;
  for (i = 0; i < myRigidBodies.size(); i++) {
    myRigidBodies[i]->calcForcesAndTorques();   
  }
}

void Molecule::getCOM( double COM[3] ) {

  double mass, mtot;
  double aPos[3];
  int i, j;

  for (j=0; j<3; j++) 
    COM[j] = 0.0;

  mtot   = 0.0;

  for (i=0; i < myIntegrableObjects.size(); i++) {
    if (myIntegrableObjects[i] != NULL) {

      mass = myIntegrableObjects[i]->getMass();
      mtot   += mass;
      
      myIntegrableObjects[i]->getPos( aPos );

      for( j = 0; j < 3; j++) 
        COM[j] += aPos[j] * mass;

    }
  }

  for (j = 0; j < 3; j++) 
    COM[j] /= mtot; 
}

double Molecule::getCOMvel( double COMvel[3] ) {

  double mass, mtot;
  double aVel[3];
  int i, j;


  for (j=0; j<3; j++) 
    COMvel[j] = 0.0;

  mtot   = 0.0;

  for (i=0; i < myIntegrableObjects.size(); i++) {
    if (myIntegrableObjects[i] != NULL) {

      mass = myIntegrableObjects[i]->getMass();
      mtot   += mass;

      myIntegrableObjects[i]->getVel(aVel);

      for (j=0; j<3; j++) 
        COMvel[j] += aVel[j]*mass;

    }
  }

  for (j=0; j<3; j++) 
    COMvel[j] /= mtot;
 
  return mtot;

}

double Molecule::getTotalMass()
{

  double totalMass;
  
  totalMass = 0;
  for(int i =0; i < myIntegrableObjects.size(); i++){
    totalMass += myIntegrableObjects[i]->getMass();
  }

  return totalMass;
}
