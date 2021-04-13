/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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

#include "primitives/DirectionalAtom.hpp"

namespace OpenMD {

DirectionalAtom::DirectionalAtom(DirectionalAtomType* dAtomType)
    : Atom(dAtomType) {
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
  ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).unitVector[localIndex_] =
      a.inverse() * sU_.getColumn(2);
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

  force = getFrc();
  torque = getTrq();
  myEuler = getA().toEulerAngles();

  phi = myEuler[0];
  theta = myEuler[1];
  psi = myEuler[2];

  cphi = cos(phi);
  sphi = sin(phi);
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
  for (int j = 0; j < 3; j++) grad[j] = -force[j];

  for (int j = 0; j < 3; j++) {
    grad[3] += torque[j] * ephi[j];
    grad[4] += torque[j] * etheta[j];
    grad[5] += torque[j] * epsi[j];
  }

  return grad;
}

void DirectionalAtom::accept(BaseVisitor* v) { v->visit(this); }

}  // namespace OpenMD
