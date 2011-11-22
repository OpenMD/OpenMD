/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "integrators/Velocitizer.hpp"
#include "math/SquareMatrix3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"

#ifndef IS_MPI
#include "math/SeqRandNumGen.hpp"
#else
#include "math/ParallelRandNumGen.hpp"
#endif

/* Remove me after testing*/
/*
#include <cstdio>
#include <iostream>
*/
/*End remove me*/

namespace OpenMD {
  
  Velocitizer::Velocitizer(SimInfo* info) : info_(info) {
    
    int seedValue;
    Globals * simParams = info->getSimParams();
    
#ifndef IS_MPI
    if (simParams->haveSeed()) {
      seedValue = simParams->getSeed();
      randNumGen_ = new SeqRandNumGen(seedValue);
    }else {
      randNumGen_ = new SeqRandNumGen();
    }    
#else
    if (simParams->haveSeed()) {
      seedValue = simParams->getSeed();
      randNumGen_ = new ParallelRandNumGen(seedValue);
    }else {
      randNumGen_ = new ParallelRandNumGen();
    }    
#endif 
  }
  
  Velocitizer::~Velocitizer() {
    delete randNumGen_;
  }
  
  void Velocitizer::velocitize(RealType temperature) {
    Vector3d aVel;
    Vector3d aJ;
    Mat3x3d I;
    int l;
    int m;
    int n; 
    Vector3d vdrift;
    RealType vbar;
    /**@todo refactory kb */
    const RealType kb = 8.31451e-7; // kb in amu, angstroms, fs, etc.
    RealType av2;
    RealType kebar;
    
    Globals * simParams = info_->getSimParams();
    
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule * mol;
    StuntDouble * integrableObject;

    kebar = kb * temperature * info_->getNdfRaw() / (2.0 * info_->getNdf());
    for( mol = info_->beginMolecule(i); mol != NULL;
	 mol = info_->nextMolecule(i) ) {
      for( integrableObject = mol->beginIntegrableObject(j);
	   integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j) ) {
	
	// uses equipartition theory to solve for vbar in angstrom/fs
	
	av2 = 2.0 * kebar / integrableObject->getMass();
	vbar = sqrt(av2);
	
	// picks random velocities from a gaussian distribution
	// centered on vbar
	
	for( int k = 0; k < 3; k++ ) {
	  aVel[k] = vbar * randNumGen_->randNorm(0.0, 1.0);
	}
	integrableObject->setVel(aVel);
	
	if (integrableObject->isDirectional()) {
	  I = integrableObject->getI();
	  
	  if (integrableObject->isLinear()) {
	    l = integrableObject->linearAxis();
	    m = (l + 1) % 3;
	    n = (l + 2) % 3;
	    
	    aJ[l] = 0.0;
	    vbar = sqrt(2.0 * kebar * I(m, m));
	    aJ[m] = vbar * randNumGen_->randNorm(0.0, 1.0);
	    vbar = sqrt(2.0 * kebar * I(n, n));
	    aJ[n] = vbar * randNumGen_->randNorm(0.0, 1.0);
	  } else {
	    for( int k = 0; k < 3; k++ ) {
	      vbar = sqrt(2.0 * kebar * I(k, k));
	      aJ[k] = vbar *randNumGen_->randNorm(0.0, 1.0);
	    }
	  } // else isLinear
	  
	  integrableObject->setJ(aJ);
	}     //isDirectional 
      }
    }             //end for (mol = beginMolecule(i); ...)
    
    
    
    removeComDrift();
    // Remove angular drift if we are not using periodic boundary conditions.
    if(!simParams->getUsePeriodicBoundaryConditions()) removeAngularDrift();
    
  }
  
  
  
  void Velocitizer::removeComDrift() {
    // Get the Center of Mass drift velocity.
    Vector3d vdrift = info_->getComVel();
    
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule * mol;
    StuntDouble * integrableObject;
    
    //  Corrects for the center of mass drift.
    // sums all the momentum and divides by total mass.
    for( mol = info_->beginMolecule(i); mol != NULL;
	 mol = info_->nextMolecule(i) ) {
      for( integrableObject = mol->beginIntegrableObject(j);
	   integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j) ) {
	integrableObject->setVel(integrableObject->getVel() - vdrift);
      }
    }
    
  }
  
   
  void Velocitizer::removeAngularDrift() {
    // Get the Center of Mass drift velocity.
      
    Vector3d vdrift;
    Vector3d com; 
      
    info_->getComAll(com,vdrift);
         
    Mat3x3d inertiaTensor;
    Vector3d angularMomentum;
    Vector3d omega;
      
      
      
    info_->getInertiaTensor(inertiaTensor,angularMomentum);
    // We now need the inverse of the inertia tensor.
    /*
    std::cerr << "Angular Momentum before is "
              << angularMomentum <<  std::endl;
    std::cerr << "Inertia Tensor before is "
              << inertiaTensor <<  std::endl;
    */
    inertiaTensor =inertiaTensor.inverse();
    /*
    std::cerr << "Inertia Tensor after inverse is "
              << inertiaTensor <<  std::endl;
    */
    omega = inertiaTensor*angularMomentum;
      
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule * mol;
    StuntDouble * integrableObject;
    Vector3d tempComPos;
      
    //  Corrects for the center of mass angular drift.
    // sums all the angular momentum and divides by total mass.
    for( mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i) ) {
      for( integrableObject = mol->beginIntegrableObject(j);
           integrableObject != NULL;
           integrableObject = mol->nextIntegrableObject(j) ) {
        tempComPos = integrableObject->getPos()-com;
        integrableObject->setVel((integrableObject->getVel() - vdrift)-cross(omega,tempComPos));
      }
    }
      
    angularMomentum = info_->getAngularMomentum();
    /*
    std::cerr << "Angular Momentum after is "
              << angularMomentum <<  std::endl;
    */ 
  }
   
   
   
   
}
