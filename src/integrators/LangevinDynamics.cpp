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
 
/**
 * @file LangevinDynamics.cpp
 * @author tlin
 * @date 11/08/2004
 * @time 15:13am
 * @version 1.0
 */

#include "integrators/LangevinDynamics.hpp"
#include "primitives/Molecule.hpp"
#include "utils/OOPSEConstant.hpp"
#include "integrators/LDForceManager.hpp"
namespace oopse {

  
  LangevinDynamics::LangevinDynamics(SimInfo* info) : VelocityVerletIntegrator(info){
    setForceManager(new LDForceManager(info));
  }
  
  void LangevinDynamics::moveA(){
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;
    Vector3d vel;
    Vector3d pos;
    Vector3d frc;
    Vector3d Tb;
    Vector3d ji;
    RealType mass;
    
    for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j)) {

	vel =integrableObject->getVel();
	pos = integrableObject->getPos();
	frc = integrableObject->getFrc();
	mass = integrableObject->getMass();
                
	// velocity half step
	vel += (dt2 /mass * OOPSEConstant::energyConvert) * frc;

	// position whole step
	pos += dt * vel;

	integrableObject->setVel(vel);
	integrableObject->setPos(pos);

	if (integrableObject->isDirectional()){

	  // get and convert the torque to body frame

	  Tb = integrableObject->lab2Body(integrableObject->getTrq());

	  // get the angular momentum, and propagate a half step

	  ji = integrableObject->getJ();

	  ji += (dt2  * OOPSEConstant::energyConvert) * Tb;

	  rotAlgo->rotate(integrableObject, ji, dt);

	  integrableObject->setJ(ji);
	}

            
      }
    } //end for(mol = info_->beginMolecule(i))
    
    rattle->constraintA();
    
  }    

  void LangevinDynamics::moveB(){
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;
    Vector3d vel;
    Vector3d frc;
    Vector3d Tb;
    Vector3d ji;
    RealType mass;
    
    for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j)) {

	vel =integrableObject->getVel();
	frc = integrableObject->getFrc();
	mass = integrableObject->getMass();
                
	// velocity half step
	vel += (dt2 /mass * OOPSEConstant::energyConvert) * frc;
                
	integrableObject->setVel(vel);

	if (integrableObject->isDirectional()){

	  // get and convert the torque to body frame

	  Tb = integrableObject->lab2Body(integrableObject->getTrq());

	  // get the angular momentum, and propagate a half step

	  ji = integrableObject->getJ();

	  ji += (dt2  * OOPSEConstant::energyConvert) * Tb;

	  integrableObject->setJ(ji);
	}

            
      }
    } //end for(mol = info_->beginMolecule(i))
  

    rattle->constraintB();

  }


  RealType LangevinDynamics::calcConservedQuantity() {
    return 0.0;
  }

} //end namespace oopse

