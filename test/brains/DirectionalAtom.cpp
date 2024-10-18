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

#include "primitives/DirectionalAtom.hpp"

namespace OpenMD {

  DirectionalAtom::DirectionalAtom(DirectionalAtomType* dAtomType) :
      Atom(dAtomType) {
    objType_ = otDAtom;
  }

  Mat3x3d DirectionalAtom::getI() {
    return static_cast<DirectionalAtomType*>(getAtomType())->getI();
  }

  void DirectionalAtom::setPrevA(const RotMat3x3d& a) {
    ((snapshotMan_->getPrevSnapshot())->*storage_).aMat[localIndex_] = a;
    ((snapshotMan_->getPrevSnapshot())->*storage_).unitVector[localIndex_] =
        a.inverse() * sU_.getColumn(2);
  }

  void DirectionalAtom::setA(const RotMat3x3d& a) {
    ((snapshotMan_->getCurrentSnapshot())->*storage_).aMat[localIndex_] = a;
    ((snapshotMan_->getCurrentSnapshot())->*storage_).unitVector[localIndex_] =
        a.inverse() * sU_.getColumn(2);
  }

  void DirectionalAtom::setA(const RotMat3x3d& a, int snapshotNo) {
    ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).aMat[localIndex_] = a;
    ((snapshotMan_->getSnapshot(snapshotNo))->*storage_)
        .unitVector[localIndex_] = a.inverse() * sU_.getColumn(2);
  }

  void DirectionalAtom::rotateBy(const RotMat3x3d& m) { setA(m * getA()); }

  void DirectionalAtom::setUnitFrameFromEuler(double phi, double theta,
                                              double psi) {
    sU_.setupRotMat(phi, theta, psi);
  }

  std::vector<double> DirectionalAtom::getGrad() {
    vector<double> grad(6, 0.0);
    Vector3d force;
    Vector3d torque;
    Vector3d myEuler;
    double phi, theta, psi;
    double cphi, sphi, ctheta, stheta;
    Vector3d ephi;
    Vector3d etheta;
    Vector3d epsi;

    force   = getFrc();
    torque  = getTrq();
    myEuler = getA().toEulerAngles();

    phi   = myEuler[0];
    theta = myEuler[1];
    psi   = myEuler[2];

    cphi   = cos(phi);
    sphi   = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);

    // get unit vectors along the phi, theta and psi rotation axes

    ephi[0] = 0.0;
    ephi[1] = 0.0;
    ephi[2] = 1.0;

    etheta[0] = cphi;
    etheta[1] = sphi;
    etheta[2] = 0.0;

    epsi[0] = stheta * cphi;
    epsi[1] = stheta * sphi;
    epsi[2] = ctheta;

    // gradient is equal to -force
    for (int j = 0; j < 3; j++)
      grad[j] = -force[j];

    for (int j = 0; j < 3; j++) {
      grad[3] += torque[j] * ephi[j];
      grad[4] += torque[j] * etheta[j];
      grad[5] += torque[j] * epsi[j];
    }

    return grad;
  }

  void DirectionalAtom::accept(BaseVisitor* v) { v->visit(this); }

}  // namespace OpenMD
