#include <iostream>

#include <cstdlib>
#include <cmath>

#include "simError.h"
#include "MoLocator.hpp"
#include "MatVec3.h"

MoLocator::MoLocator( MoleculeStamp* theStamp, ForceFields* theFF){

  myStamp = theStamp;
  myFF = theFF;
  nIntegrableObjects = myStamp->getNIntegrable();
  calcRefCoords();
}

void MoLocator::placeMol( const Vector3d& offset, const Vector3d& ort, Molecule* mol){
  double newCoor[3];
  double curRefCoor[3];
  double zeroVector[3];
  vector<StuntDouble*> myIntegrableObjects;
  double rotMat[3][3];
  
  zeroVector[0] = 0.0;
  zeroVector[1] = 0.0;
  zeroVector[2] = 0.0;
  
  latVec2RotMat(ort, rotMat);
  
  myIntegrableObjects = mol->getIntegrableObjects();

  if(myIntegrableObjects.size() != nIntegrableObjects){
      sprintf( painCave.errMsg,
	       "MoLocator error.\n"
	       "  The number of integrable objects of MoleculeStamp is not the same as  that of Molecule\n");
	  painCave.isFatal = 1;
      simError();

  }
    
  for(int i=0; i<nIntegrableObjects; i++) {

    //calculate the reference coordinate for integrable objects after rotation
    curRefCoor[0] = refCoords[i][0];
    curRefCoor[1] = refCoords[i][1];
    curRefCoor[2] = refCoords[i][2];
    
    matVecMul3(rotMat, curRefCoor, newCoor);    

    newCoor[0] +=  offset[0];
    newCoor[1] +=  offset[1];
    newCoor[2] +=  offset[2];

    myIntegrableObjects[i]->setPos( newCoor);
    myIntegrableObjects[i]->setVel(zeroVector);

    if(myIntegrableObjects[i]->isDirectional()){
      myIntegrableObjects[i]->setA(rotMat);
      myIntegrableObjects[i]->setJ(zeroVector);  
    }
  }

}
 
void MoLocator::calcRefCoords( void ){
  AtomStamp* currAtomStamp;
  int nAtoms; 
  int nRigidBodies;
  vector<double> mass;
  Vector3d coor;
  Vector3d refMolCom;  
  int nAtomsInRb;
  double totMassInRb;
  double currAtomMass;
  double molMass;
  
  nAtoms= myStamp->getNAtoms();
  nRigidBodies = myStamp->getNRigidBodies();

  for(size_t i=0; i<nAtoms; i++){

    currAtomStamp = myStamp->getAtom(i);

    if( !currAtomStamp->havePosition() ){
      sprintf( painCave.errMsg,
                  "MoLocator error.\n"
                  "  Component %s, atom %s does not have a position specified.\n"
                  "  This means MoLocator cannot initalize it's position.\n",
                  myStamp->getID(),
                  currAtomStamp->getType() );

      painCave.isFatal = 1;
      simError();
    }

    //if atom belongs to rigidbody, just skip it
    if(myStamp->isAtomInRigidBody(i))
      continue;
    //get mass and the reference coordinate 
    else{
      currAtomMass = myFF->getAtomTypeMass(currAtomStamp->getType());
      mass.push_back(currAtomMass);
      coor.x = currAtomStamp->getPosX();
      coor.y = currAtomStamp->getPosY();
      coor.z = currAtomStamp->getPosZ();
      refCoords.push_back(coor);

    }
  }

  for(int i = 0; i < nRigidBodies; i++){
    coor.x = 0;
    coor.y = 0;
    coor.z = 0;
    totMassInRb = 0;

    for(int j = 0; j < nAtomsInRb; j++){

      currAtomMass = myFF->getAtomTypeMass(currAtomStamp->getType());
      totMassInRb +=  currAtomMass;
      
      coor.x += currAtomStamp->getPosX() * currAtomMass;
      coor.y += currAtomStamp->getPosY() * currAtomMass;
      coor.z += currAtomStamp->getPosZ() * currAtomMass;
    }

    mass.push_back(totMassInRb);
    coor /= totMassInRb;
    refCoords.push_back(coor);
  }


  //calculate the reference center of mass
  molMass = 0;
  refMolCom.x = 0;
  refMolCom.y = 0;
  refMolCom.z = 0;
  
  for(int i = 0; i < nIntegrableObjects; i++){
    refMolCom += refCoords[i] * mass[i];
   molMass += mass[i];
  }
  
  refMolCom /= molMass;

  //move the reference center of mass to (0,0,0) and adjust the reference coordinate 
  //of the integrabel objects
  for(int i = 0; i < nIntegrableObjects; i++)
    refCoords[i] -= refMolCom;
}


void latVec2RotMat(const Vector3d& lv, double rotMat[3][3]){

  double theta, phi, psi;
  
  theta =acos(lv.z);
  phi = atan2(lv.y, lv.x);
  psi = 0;

  rotMat[0][0] = (cos(phi) * cos(psi)) - (sin(phi) * cos(theta) * sin(psi));
  rotMat[0][1] = (sin(phi) * cos(psi)) + (cos(phi) * cos(theta) * sin(psi));
  rotMat[0][2] = sin(theta) * sin(psi);
  
  rotMat[1][0] = -(cos(phi) * sin(psi)) - (sin(phi) * cos(theta) * cos(psi));
  rotMat[1][1] = -(sin(phi) * sin(psi)) + (cos(phi) * cos(theta) * cos(psi));
  rotMat[1][2] = sin(theta) * cos(psi);
  
  rotMat[2][0] = sin(phi) * sin(theta);
  rotMat[2][1] = -cos(phi) * sin(theta);
  rotMat[2][2] = cos(theta);
}

