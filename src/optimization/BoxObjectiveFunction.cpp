/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "optimization/BoxObjectiveFunction.hpp"

#include "math/CholeskyDecomposition.hpp"

namespace OpenMD {

  BoxObjectiveFunction::BoxObjectiveFunction(SimInfo* info,
                                             ForceManager* forceMan) :
      info_(info),
      forceMan_(forceMan), thermo(info) {
    shake_ = new Shake(info_);

    if (info_->usesFluctuatingCharges()) {
      if (info_->getNFluctuatingCharges() > 0) {
        hasFlucQ_      = true;
        fqConstraints_ = new FluctuatingChargeConstraints(info_);
        bool cr        = info_->getSimParams()
                      ->getFluctuatingChargeParameters()
                      ->getConstrainRegions();
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

  RealType BoxObjectiveFunction::valueAndGradient(
      DynamicVector<RealType>& grad, const DynamicVector<RealType>& x) {
    info_->getSnapshotManager()->advance();
    if (setCoor(x) == 0) {
      shake_->constraintR();
      forceMan_->calcForces();
      if (hasFlucQ_) fqConstraints_->applyConstraints();
      shake_->constraintF();
      getGrad(grad);
      return thermo.getPotential();
    } else {
      // The deformation was too large, so return infinite
      // potential and gradient
      for (int j = 0; j < 6; j++)
        grad[j] = std::numeric_limits<RealType>::infinity();
      return std::numeric_limits<RealType>::infinity();
    }
  }

  int BoxObjectiveFunction::setCoor(const DynamicVector<RealType>& x) {
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
    eta.setupVoigtTensor(x[0], x[1], x[2], x[3] / 2., x[4] / 2., x[5] / 2.);

    // Make sure the deformation isn't too large:
    if (eta.frobeniusNorm() > 0.7) {
      // Deformation is too large, return an error condition:
      return -1;
    }

    // Find the physical strain tensor, ε, from the Lagrangian strain, η:
    // η = ε + 0.5 * ε^2
    norm = 1.0;
    eps  = eta;
    while (norm > 1.0e-10) {
      y    = eta - eps * eps / 2.0;
      test = y - eps;
      norm = test.frobeniusNorm();
      eps  = y;
    }
    deformation_ = SquareMatrix3<RealType>::identity() + eps;

    int index = 0;

    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {
      posO  = refPos_[index++];
      posN  = mol->getCom();
      delta = deformation_ * posO - posN;
      mol->moveCom(delta);
    }

    Mat3x3d Hmat = deformation_ * refHmat_;
    info_->getSnapshotManager()->getCurrentSnapshot()->setHmat(Hmat);
    return 0;
  }

  void BoxObjectiveFunction::getGrad(DynamicVector<RealType>& grad) {
    Mat3x3d pressureTensor;
    Vector<RealType, 6> lstress(0.0);

    // Find the Lagragian stress tensor, τ, from the physical
    // stress tensor, σ, that was computed from the pressureTensor
    // in this code.
    // τ = det(1+ε) (1+ε)^−1 · σ · (1+ε)^−1
    // (Note that 1+ε is the deformation tensor computed above.)

    Mat3x3d idm  = deformation_.inverse();
    RealType ddm = deformation_.determinant();

    pressureTensor = thermo.getPressureTensor();
    pressureTensor.negate();
    pressureTensor *= Constants::elasticConvert;

    Mat3x3d tao = idm * (pressureTensor * idm);
    tao *= ddm;

    lstress    = tao.toVoigtTensor();
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
    refHmat_       = snap->getHmat();
    V0_            = snap->getVolume();

    refPos_.clear();
    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {
      refPos_.push_back(mol->getCom());
    }

    return xinit;
  }
}  // namespace OpenMD
