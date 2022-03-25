/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#include "primitives/DirectionalAtom.hpp"

#include "types/DirectionalAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  DirectionalAtom::DirectionalAtom(AtomType* dAtomType) : Atom(dAtomType) {
    objType_ = otDAtom;

    DirectionalAdapter da = DirectionalAdapter(dAtomType);
    I_                    = da.getI();

    MultipoleAdapter ma = MultipoleAdapter(dAtomType);
    if (ma.isDipole()) { dipole_ = ma.getDipole(); }
    if (ma.isQuadrupole()) { quadrupole_ = ma.getQuadrupole(); }

    // Check if one of the diagonal inertia tensor of this directional
    // atom is zero:
    int nLinearAxis       = 0;
    Mat3x3d inertiaTensor = getI();
    for (int i = 0; i < 3; i++) {
      if (fabs(inertiaTensor(i, i)) < OpenMD::epsilon) {
        linear_     = true;
        linearAxis_ = i;
        ++nLinearAxis;
      }
    }

    if (nLinearAxis > 1) {
      sprintf(
          painCave.errMsg,
          "Directional Atom warning.\n"
          "\tOpenMD found more than one axis in this directional atom with a "
          "vanishing \n"
          "\tmoment of inertia.");
      painCave.isFatal = 0;
      simError();
    }
  }

  Mat3x3d DirectionalAtom::getI() { return I_; }

  void DirectionalAtom::setPrevA(const RotMat3x3d& a) {
    ((snapshotMan_->getPrevSnapshot())->*storage_).aMat[localIndex_] = a;

    if (atomType_->isMultipole()) {
      RotMat3x3d atrans = a.transpose();

      if (atomType_->isDipole()) {
        ((snapshotMan_->getPrevSnapshot())->*storage_).dipole[localIndex_] =
            atrans * dipole_;
      }

      if (atomType_->isQuadrupole()) {
        ((snapshotMan_->getPrevSnapshot())->*storage_).quadrupole[localIndex_] =
            atrans * quadrupole_ * a;
      }
    }
  }

  void DirectionalAtom::setA(const RotMat3x3d& a) {
    ((snapshotMan_->getCurrentSnapshot())->*storage_).aMat[localIndex_] = a;

    if (atomType_->isMultipole()) {
      RotMat3x3d atrans = a.transpose();

      if (atomType_->isDipole()) {
        ((snapshotMan_->getCurrentSnapshot())->*storage_).dipole[localIndex_] =
            atrans * dipole_;
      }

      if (atomType_->isQuadrupole()) {
        ((snapshotMan_->getCurrentSnapshot())->*storage_)
            .quadrupole[localIndex_] = atrans * quadrupole_ * a;
      }
    }
  }

  void DirectionalAtom::setA(const RotMat3x3d& a, int snapshotNo) {
    ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).aMat[localIndex_] = a;

    if (atomType_->isMultipole()) {
      RotMat3x3d atrans = a.transpose();

      if (atomType_->isDipole()) {
        ((snapshotMan_->getSnapshot(snapshotNo))->*storage_)
            .dipole[localIndex_] = atrans * dipole_;
      }

      if (atomType_->isQuadrupole()) {
        ((snapshotMan_->getSnapshot(snapshotNo))->*storage_)
            .quadrupole[localIndex_] = atrans * quadrupole_ * a;
      }
    }
  }

  void DirectionalAtom::rotateBy(const RotMat3x3d& m) { setA(m * getA()); }

  std::vector<RealType> DirectionalAtom::getGrad() {
    std::vector<RealType> grad(6, 0.0);
    Vector3d force;
    Vector3d torque;
    Vector3d myEuler;
    RealType phi, theta;
    RealType cphi, sphi, ctheta, stheta;
    Vector3d ephi;
    Vector3d etheta;
    Vector3d epsi;

    force  = getFrc();
    torque = getTrq();

    myEuler = getA().toEulerAngles();

    phi   = myEuler[0];
    theta = myEuler[1];

    cphi   = cos(phi);
    sphi   = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);

    if (fabs(stheta) < 1.0E-9) { stheta = 1.0E-9; }

    ephi[0] = -sphi * ctheta / stheta;
    ephi[1] = cphi * ctheta / stheta;
    ephi[2] = 1.0;

    etheta[0] = cphi;
    etheta[1] = sphi;
    etheta[2] = 0.0;

    epsi[0] = sphi / stheta;
    epsi[1] = -cphi / stheta;
    epsi[2] = 0.0;

    // gradient is equal to -force
    for (int j = 0; j < 3; j++)
      grad[j] = -force[j];

    for (int j = 0; j < 3; j++) {
      grad[3] -= torque[j] * ephi[j];
      grad[4] -= torque[j] * etheta[j];
      grad[5] -= torque[j] * epsi[j];
    }

    return grad;
  }

  void DirectionalAtom::accept(BaseVisitor* v) { v->visit(this); }
}  // namespace OpenMD
