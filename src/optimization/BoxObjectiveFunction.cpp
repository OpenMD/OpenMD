/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "optimization/BoxObjectiveFunction.hpp"
#include "math/CholeskyDecomposition.hpp"

namespace OpenMD{  

  BoxObjectiveFunction::BoxObjectiveFunction(SimInfo* info,
                                             ForceManager* forceMan)
    : info_(info), forceMan_(forceMan), thermo(info) {   
    shake_ = new Shake(info_);
    
    if (info_->usesFluctuatingCharges()) {
      if (info_->getNFluctuatingCharges() > 0) {
        hasFlucQ_ = true;
        fqConstraints_ = new FluctuatingChargeConstraints(info_);
        bool cr = info_->getSimParams()->getFluctuatingChargeParameters()->getConstrainRegions();    
        fqConstraints_->setConstrainRegions(cr);
      }
    }
  }
  
  RealType BoxObjectiveFunction::value(const DynamicVector<RealType>& x) {
    info_->getSnapshotManager()->advance();
    
    if (setCoor(x) == 0) {
      shake_->constraintR();
      forceMan_->calcForces();
      if (hasFlucQ_) fqConstraints_->applyConstraints();
      shake_->constraintF();
      return thermo.getPotential();
    } else {
      // The deformation was too large, so return an infinite potential
      return std::numeric_limits<RealType>::infinity();
    }
  }
  
  void BoxObjectiveFunction::gradient(DynamicVector<RealType>& grad,
                                      const DynamicVector<RealType>& x) {    

    info_->getSnapshotManager()->advance();            
    if (setCoor(x) == 0) {
      shake_->constraintR();
      forceMan_->calcForces();
      if (hasFlucQ_) fqConstraints_->applyConstraints();
      shake_->constraintF();
      getGrad(grad);
    } else {
      // The deformation was too large, so return an infinite gradient:
      for (int j = 0; j < 6; j++) 
        grad[j] = std::numeric_limits<RealType>::infinity();
    }

  }
  
  RealType BoxObjectiveFunction::valueAndGradient(DynamicVector<RealType>& grad,
                                                  const DynamicVector<RealType>& x) {
    
    info_->getSnapshotManager()->advance();            
    if (setCoor(x) == 0) {
      shake_->constraintR();
      forceMan_->calcForces();
      if (hasFlucQ_) fqConstraints_->applyConstraints();
      shake_->constraintF();
      getGrad(grad);
      return thermo.getPotential();;
    } else {
      // The deformation was too large, so return infinite
      // potential and gradient
      for (int j = 0; j < 6; j++) 
        grad[j] = std::numeric_limits<RealType>::infinity();      
      return std::numeric_limits<RealType>::infinity();
    }

  }
  
  int BoxObjectiveFunction::setCoor(const DynamicVector<RealType> &x) {
    Vector3d posO;
    Vector3d posN;
    Vector3d delta;
    SimInfo::MoleculeIterator miter;
    Molecule* mol;
    Mat3x3d eta(0.0);
    Mat3x3d eps(0.0);
    Mat3x3d y(0.0);
    Mat3x3d test(0.0);
    RealType norm;

    // η is the Lagrangian strain tensor:
    eta.setupVoigtTensor(x[0], x[1], x[2], x[3]/2., x[4]/2., x[5]/2.);
    
    // Make sure the deformation isn't too large:
    if (eta.frobeniusNorm() > 0.7) {
      // Deformation is too large, return an error condition:
      return -1;
    }

    // Find the physical strain tensor, ε, from the Lagrangian strain, η:
    // η = ε + 0.5 * ε^2
    norm = 1.0;
    eps = eta;
    while (norm > 1.0e-10) {
      y = eta - eps*eps / 2.0;
      test = y - eps;
      norm = test.frobeniusNorm();       
      eps = y;
    }
    deformation_ = SquareMatrix3<RealType>::identity() + eps;

    int index = 0;
    
    for (mol = info_->beginMolecule(miter); mol != NULL; 
         mol = info_->nextMolecule(miter)) {
      posO = refPos_[index++];
      posN = mol->getCom();
      delta = deformation_ * posO - posN;
      mol->moveCom(delta);
    }

    Mat3x3d Hmat = deformation_ * refHmat_;
    info_->getSnapshotManager()->getCurrentSnapshot()->setHmat(Hmat);
    return 0;
  }
  
  void BoxObjectiveFunction::getGrad(DynamicVector<RealType> &grad) {
    Mat3x3d pressureTensor;
    Vector<RealType, 6> lstress(0.0);

    // Find the Lagragian stress tensor, τ, from the physical
    // stress tensor, σ, that was computed from the pressureTensor
    // in this code.       
    // τ = det(1+ε) (1+ε)^−1 · σ · (1+ε)^−1
    // (Note that 1+ε is the deformation tensor computed above.)
    
    Mat3x3d idm = deformation_.inverse();
    RealType ddm = deformation_.determinant();
    
    pressureTensor = thermo.getPressureTensor();        
    pressureTensor.negate();
    pressureTensor *= Constants::elasticConvert;
    
    Mat3x3d tao = idm * (pressureTensor * idm);
    tao *= ddm;                    
    
    lstress = tao.toVoigtTensor();
    RealType V = thermo.getVolume();
    
    for (int j = 0; j < 6; j++) {
      grad[j] = V * lstress[j];
    }
  }

  DynamicVector<RealType> BoxObjectiveFunction::setInitialCoords() {    
    DynamicVector<RealType> xinit(6, 0.0);
    SimInfo::MoleculeIterator miter;
    Molecule* mol;

    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    refHmat_ = snap->getHmat();
    V0_ = snap->getVolume();

    refPos_.clear();
    for (mol = info_->beginMolecule(miter); mol != NULL; 
         mol = info_->nextMolecule(miter)) {
      refPos_.push_back( mol->getCom() );
    }
    
    return xinit;
  }
}
