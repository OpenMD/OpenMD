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

#include "flucq/FluctuatingChargeObjectiveFunction.hpp"

namespace OpenMD{
  
  FluctuatingChargeObjectiveFunction::FluctuatingChargeObjectiveFunction(SimInfo* info, ForceManager* forceMan, FluctuatingChargeConstraints* fqConstraints)
    : info_(info), forceMan_(forceMan), fqConstraints_(fqConstraints), 
      thermo(info) {       
  }
  
  RealType FluctuatingChargeObjectiveFunction::value(const DynamicVector<RealType>& x) {
    setCoor(x);
    forceMan_->calcForces();

    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    potVec pot = curSnapshot->getLongRangePotentials();
    potVec exPot = curSnapshot->getExcludedPotentials();  
    cerr << "val p= " <<   pot[ELECTROSTATIC_FAMILY] << "\n";
    cerr << "val e= " << exPot[ELECTROSTATIC_FAMILY] << "\n";
  
    return pot[ELECTROSTATIC_FAMILY] + exPot[ELECTROSTATIC_FAMILY];
  }
  
  void FluctuatingChargeObjectiveFunction::gradient(DynamicVector<RealType>& grad, const DynamicVector<RealType>& x) {
    int shakeStatus;
    
    setCoor(x);     
    
    forceMan_->calcForces(); 
    fqConstraints_->applyConstraints();
    cerr << "grad\n";
    getGrad(grad);      
  }
  
  RealType FluctuatingChargeObjectiveFunction::valueAndGradient(DynamicVector<RealType>& grad,
                                                              const DynamicVector<RealType>& x) {

    setCoor(x);     

    forceMan_->calcForces();
    fqConstraints_->applyConstraints();
    
    getGrad(grad); 

    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    potVec pot = curSnapshot->getLongRangePotentials();
    potVec exPot = curSnapshot->getExcludedPotentials();    
    cerr << "vang p= " <<   pot[ELECTROSTATIC_FAMILY] << "\n";
    cerr << "vang e= " << exPot[ELECTROSTATIC_FAMILY] << "\n";

    return pot[ELECTROSTATIC_FAMILY] + exPot[ELECTROSTATIC_FAMILY];
  }
  
  void FluctuatingChargeObjectiveFunction::setCoor(const DynamicVector<RealType> &x) const {
    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator  j;
    Molecule* mol;
    Atom* atom;
    int index = 0;
    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {
      
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
       
        atom->setFlucQPos(x[index++]);
        cerr << "setting charge = " << x[index -1] << "\n";
      }
    }    
  }
  
  void FluctuatingChargeObjectiveFunction::getGrad(DynamicVector<RealType> &grad) {
    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator  j;
    Molecule* mol;
    Atom* atom;
    
    int index = 0;
    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {

        grad[index++] = atom->getFlucQFrc();
        cerr << "getting grad = " << grad[index -1] << "\n";
      }
    }
  }

  DynamicVector<RealType> FluctuatingChargeObjectiveFunction::setInitialCoords() {
    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator  j;
    Molecule* mol;
    Atom* atom;
    
    int index = 0;    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        
        index++;
      }
    }

    DynamicVector<RealType> initCoords(index);
    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        
        initCoords[index++] = atom->getFlucQPos();
      }
    }   
    return initCoords;
  }
}
