/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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

#include "optimization/PotentialEnergyObjectiveFunction.hpp"

namespace OpenMD{

  PotentialEnergyObjectiveFunction::PotentialEnergyObjectiveFunction(SimInfo* info, ForceManager* forceMan)
    : info_(info), forceMan_(forceMan), thermo(info) {   
    shake_ = new Shake(info_);    
  }
  

  
  RealType PotentialEnergyObjectiveFunction::value(const DynamicVector<RealType>& x) {
    setCoor(x);
    shake_->constraintR();
    forceMan_->calcForces();
    shake_->constraintF();
    return thermo.getPotential();
  }
  
  void PotentialEnergyObjectiveFunction::gradient(DynamicVector<RealType>& grad, const DynamicVector<RealType>& x) {
    
    setCoor(x);       
    shake_->constraintR();
    forceMan_->calcForces();
    shake_->constraintF();
    getGrad(grad);
  }
  
  RealType PotentialEnergyObjectiveFunction::valueAndGradient(DynamicVector<RealType>& grad,
                                                              const DynamicVector<RealType>& x) {
    setCoor(x);
    shake_->constraintR();
    forceMan_->calcForces();
    shake_->constraintF();
    getGrad(grad); 
    return thermo.getPotential();
  }
  
  void PotentialEnergyObjectiveFunction::setCoor(const DynamicVector<RealType> &x) const {
    Vector3d position;
    Vector3d eulerAngle;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;    
    int index = 0;
    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); 
           integrableObject != NULL;
           integrableObject = mol->nextIntegrableObject(j)) {
        
        position[0] = x[index++];
        position[1] = x[index++];
        position[2] = x[index++];
        
        integrableObject->setPos(position);
        
        if (integrableObject->isDirectional()) {
          eulerAngle[0] = x[index++];
          eulerAngle[1] = x[index++];
          eulerAngle[2] = x[index++];
          
          integrableObject->setEuler(eulerAngle);

          if (integrableObject->isRigidBody()) {
            RigidBody* rb = static_cast<RigidBody*>(integrableObject);
            rb->updateAtoms();
          }        
        }
      }
    }    
  }
  
  void PotentialEnergyObjectiveFunction::getGrad(DynamicVector<RealType> &grad) {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;    
    std::vector<RealType> myGrad;
    
    int index = 0;
    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); 
           integrableObject != NULL;
           integrableObject = mol->nextIntegrableObject(j)) {        
        myGrad = integrableObject->getGrad();

        for (size_t k = 0; k < myGrad.size(); ++k) {   
          grad[index++] = myGrad[k];
        }
      }            
    }         
  }

  DynamicVector<RealType> PotentialEnergyObjectiveFunction::setInitialCoords() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;    

    Vector3d pos;
    Vector3d eulerAngle;
    
    DynamicVector<RealType> xinit(info_->getNdfLocal(), 0.0);

    int index = 0;
    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); 
           integrableObject != NULL;
           integrableObject = mol->nextIntegrableObject(j)) {        

        pos = integrableObject->getPos();

        xinit[index++] = pos[0];
        xinit[index++] = pos[1];
        xinit[index++] = pos[2];
                
        if (integrableObject->isDirectional()) {
          eulerAngle = integrableObject->getEuler();
          xinit[index++] = eulerAngle[0];
          xinit[index++] = eulerAngle[1];
          xinit[index++] = eulerAngle[2];          
        }
      }
    }
    return xinit;
  }
}
