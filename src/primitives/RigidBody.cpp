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
#include "primitives/RigidBody.hpp"

#include <algorithm>
#include <cmath>

#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  RigidBody::RigidBody() :
      StuntDouble(otRigidBody, &Snapshot::rigidbodyData), inertiaTensor_(0.0) {}

  void RigidBody::setPrevA(const RotMat3x3d& a) {
    ((snapshotMan_->getPrevSnapshot())->*storage_).aMat[localIndex_] = a;

    for (unsigned int i = 0; i < atoms_.size(); ++i) {
      if (atoms_[i]->isDirectional()) {
        atoms_[i]->setPrevA(refOrients_[i].transpose() * a);
      }
    }
  }

  void RigidBody::setA(const RotMat3x3d& a) {
    ((snapshotMan_->getCurrentSnapshot())->*storage_).aMat[localIndex_] = a;

    for (unsigned int i = 0; i < atoms_.size(); ++i) {
      if (atoms_[i]->isDirectional()) {
        atoms_[i]->setA(refOrients_[i].transpose() * a);
      }
    }
  }

  void RigidBody::setA(const RotMat3x3d& a, int snapshotNo) {
    ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).aMat[localIndex_] = a;

    for (unsigned int i = 0; i < atoms_.size(); ++i) {
      if (atoms_[i]->isDirectional()) {
        atoms_[i]->setA(refOrients_[i].transpose() * a, snapshotNo);
      }
    }
  }

  Mat3x3d RigidBody::getI() { return inertiaTensor_; }

  std::vector<RealType> RigidBody::getGrad() {
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

  void RigidBody::accept(BaseVisitor* v) { v->visit(this); }

  /**@todo need modification */
  void RigidBody::calcRefCoords() {
    RealType mtmp;
    Vector3d refCOM(0.0);
    mass_ = 0.0;
    for (std::size_t i = 0; i < atoms_.size(); ++i) {
      mtmp = atoms_[i]->getMass();
      mass_ += mtmp;
      refCOM += refCoords_[i] * mtmp;
    }
    refCOM_ = refCOM / mass_;

    // Next, move the origin of the reference coordinate system to the COM:
    for (std::size_t i = 0; i < atoms_.size(); ++i) {
      refCoords_[i] -= refCOM_;
    }

    // Moment of Inertia calculation
    Mat3x3d Itmp(0.0);
    for (std::size_t i = 0; i < atoms_.size(); i++) {
      Mat3x3d IAtom(0.0);
      mtmp = atoms_[i]->getMass();
      IAtom -= outProduct(refCoords_[i], refCoords_[i]) * mtmp;
      RealType r2 = refCoords_[i].lengthSquare();
      IAtom(0, 0) += mtmp * r2;
      IAtom(1, 1) += mtmp * r2;
      IAtom(2, 2) += mtmp * r2;
      Itmp += IAtom;

      // project the inertial moment of directional atoms into this rigid body
      if (atoms_[i]->isDirectional()) {
        Itmp += refOrients_[i].transpose() * atoms_[i]->getI() * refOrients_[i];
      }
    }

    // diagonalize
    Vector3d evals;
    Mat3x3d::diagonalize(Itmp, evals, sU_);

    // zero out I and then fill the diagonals with the moments of inertia:
    inertiaTensor_(0, 0) = evals[0];
    inertiaTensor_(1, 1) = evals[1];
    inertiaTensor_(2, 2) = evals[2];

    int nLinearAxis = 0;
    for (int i = 0; i < 3; i++) {
      if (fabs(evals[i]) < OpenMD::epsilon) {
        linear_     = true;
        linearAxis_ = i;
        ++nLinearAxis;
      }
    }

    if (nLinearAxis > 1) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "RigidBody error.\n"
          "\tOpenMD found more than one axis in this rigid body with a "
          "vanishing \n"
          "\tmoment of inertia.  This can happen in one of three ways:\n"
          "\t 1) Only one atom was specified, or \n"
          "\t 2) All atoms were specified at the same location, or\n"
          "\t 3) The programmers did something stupid.\n"
          "\tIt is silly to use a rigid body to describe this situation.  Be "
          "smarter.\n");
      painCave.isFatal = 1;
      simError();
    }
  }

  void RigidBody::calcForcesAndTorques() {
    Vector3d afrc;
    Vector3d atrq;
    Vector3d apos;
    Vector3d rpos;
    Vector3d frc(0.0);
    Vector3d trq(0.0);
    Vector3d ef(0.0);
    Vector3d pos = this->getPos();
    AtomType* atype;
    int eCount = 0;

    int sl =
        ((snapshotMan_->getCurrentSnapshot())->*storage_).getStorageLayout();

    for (unsigned int i = 0; i < atoms_.size(); i++) {
      atype = atoms_[i]->getAtomType();

      afrc = atoms_[i]->getFrc();
      apos = atoms_[i]->getPos();
      rpos = apos - pos;

      frc += afrc;

      trq[0] += rpos[1] * afrc[2] - rpos[2] * afrc[1];
      trq[1] += rpos[2] * afrc[0] - rpos[0] * afrc[2];
      trq[2] += rpos[0] * afrc[1] - rpos[1] * afrc[0];

      // If the atom has a torque associated with it, then we also need to
      // migrate the torques onto the center of mass:

      if (atoms_[i]->isDirectional()) {
        atrq = atoms_[i]->getTrq();
        trq += atrq;
      }

      if ((sl & DataStorage::dslElectricField) && (atype->isElectrostatic())) {
        ef += atoms_[i]->getElectricField();
        eCount++;
      }
    }
    addFrc(frc);
    addTrq(trq);

    if (sl & DataStorage::dslElectricField) {
      ef /= eCount;
      setElectricField(ef);
    }
  }

  Mat3x3d RigidBody::calcForcesAndTorquesAndVirial() {
    Vector3d afrc;
    Vector3d atrq;
    Vector3d apos;
    Vector3d rpos;
    Vector3d dfrc;
    Vector3d frc(0.0);
    Vector3d trq(0.0);
    Vector3d ef(0.0);
    AtomType* atype;
    int eCount = 0;

    Vector3d pos = this->getPos();
    Mat3x3d tau_(0.0);

    int sl =
        ((snapshotMan_->getCurrentSnapshot())->*storage_).getStorageLayout();

    for (unsigned int i = 0; i < atoms_.size(); i++) {
      atype = atoms_[i]->getAtomType();

      afrc = atoms_[i]->getFrc();
      apos = atoms_[i]->getPos();
      rpos = apos - pos;

      frc += afrc;

      trq[0] += rpos[1] * afrc[2] - rpos[2] * afrc[1];
      trq[1] += rpos[2] * afrc[0] - rpos[0] * afrc[2];
      trq[2] += rpos[0] * afrc[1] - rpos[1] * afrc[0];

      // If the atom has a torque associated with it, then we also need to
      // migrate the torques onto the center of mass:

      if (atoms_[i]->isDirectional()) {
        atrq = atoms_[i]->getTrq();
        trq += atrq;
      }

      if ((sl & DataStorage::dslElectricField) && (atype->isElectrostatic())) {
        ef += atoms_[i]->getElectricField();
        eCount++;
      }

      tau_(0, 0) -= rpos[0] * afrc[0];
      tau_(0, 1) -= rpos[0] * afrc[1];
      tau_(0, 2) -= rpos[0] * afrc[2];
      tau_(1, 0) -= rpos[1] * afrc[0];
      tau_(1, 1) -= rpos[1] * afrc[1];
      tau_(1, 2) -= rpos[1] * afrc[2];
      tau_(2, 0) -= rpos[2] * afrc[0];
      tau_(2, 1) -= rpos[2] * afrc[1];
      tau_(2, 2) -= rpos[2] * afrc[2];
    }
    addFrc(frc);
    addTrq(trq);

    if (sl & DataStorage::dslElectricField) {
      ef /= eCount;
      setElectricField(ef);
    }

    return tau_;
  }

  void RigidBody::updateAtoms() {
    unsigned int i;
    Vector3d ref;
    Vector3d apos;
    DirectionalAtom* dAtom;
    Vector3d pos = getPos();
    RotMat3x3d a = getA();

    for (i = 0; i < atoms_.size(); i++) {
      ref = body2Lab(refCoords_[i]);

      apos = pos + ref;

      atoms_[i]->setPos(apos);

      if (atoms_[i]->isDirectional()) {
        dAtom = dynamic_cast<DirectionalAtom*>(atoms_[i]);
        dAtom->setA(refOrients_[i].transpose() * a);
      }
    }
  }

  void RigidBody::updateAtoms(int frame) {
    unsigned int i;
    Vector3d ref;
    Vector3d apos;
    DirectionalAtom* dAtom;
    Vector3d pos = getPos(frame);
    RotMat3x3d a = getA(frame);

    for (i = 0; i < atoms_.size(); i++) {
      ref = body2Lab(refCoords_[i], frame);

      apos = pos + ref;

      atoms_[i]->setPos(apos, frame);

      if (atoms_[i]->isDirectional()) {
        dAtom = dynamic_cast<DirectionalAtom*>(atoms_[i]);
        dAtom->setA(refOrients_[i].transpose() * a, frame);
      }
    }
  }

  void RigidBody::updateAtomVel() {
    Mat3x3d skewMat;

    Vector3d ji = getJ();
    Mat3x3d I   = getI();

    skewMat(0, 0) = 0;
    skewMat(0, 1) = ji[2] / I(2, 2);
    skewMat(0, 2) = -ji[1] / I(1, 1);

    skewMat(1, 0) = -ji[2] / I(2, 2);
    skewMat(1, 1) = 0;
    skewMat(1, 2) = ji[0] / I(0, 0);

    skewMat(2, 0) = ji[1] / I(1, 1);
    skewMat(2, 1) = -ji[0] / I(0, 0);
    skewMat(2, 2) = 0;

    Mat3x3d mat    = (getA() * skewMat).transpose();
    Vector3d rbVel = getVel();

    Vector3d velRot;
    for (unsigned int i = 0; i < refCoords_.size(); ++i) {
      atoms_[i]->setVel(rbVel + mat * refCoords_[i]);
    }
  }

  void RigidBody::updateAtomVel(int frame) {
    Mat3x3d skewMat;
    ;

    Vector3d ji = getJ(frame);
    Mat3x3d I   = getI();

    skewMat(0, 0) = 0;
    skewMat(0, 1) = ji[2] / I(2, 2);
    skewMat(0, 2) = -ji[1] / I(1, 1);

    skewMat(1, 0) = -ji[2] / I(2, 2);
    skewMat(1, 1) = 0;
    skewMat(1, 2) = ji[0] / I(0, 0);

    skewMat(2, 0) = ji[1] / I(1, 1);
    skewMat(2, 1) = -ji[0] / I(0, 0);
    skewMat(2, 2) = 0;

    Mat3x3d mat    = (getA(frame) * skewMat).transpose();
    Vector3d rbVel = getVel(frame);

    Vector3d velRot;
    for (unsigned int i = 0; i < refCoords_.size(); ++i) {
      atoms_[i]->setVel(rbVel + mat * refCoords_[i], frame);
    }
  }

  bool RigidBody::getAtomPos(Vector3d& pos, unsigned int index) {
    if (index < atoms_.size()) {
      Vector3d ref = body2Lab(refCoords_[index]);
      pos          = getPos() + ref;
      return true;
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "%d is an invalid index. The current rigid body contans %zu atoms.\n",
          index, atoms_.size());
      painCave.isFatal = 0;
      simError();
      return false;
    }
  }

  bool RigidBody::getAtomPos(Vector3d& pos, Atom* atom) {
    std::vector<Atom*>::iterator i;
    i = std::find(atoms_.begin(), atoms_.end(), atom);
    if (i != atoms_.end()) {
      // RigidBody class makes sure refCoords_ and atoms_ match each other
      Vector3d ref = body2Lab(refCoords_[i - atoms_.begin()]);
      pos          = getPos() + ref;
      return true;
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Atom %d does not belong to rigid body %d.\n",
               atom->getGlobalIndex(), getGlobalIndex());
      painCave.isFatal = 0;
      simError();
      return false;
    }
  }
  bool RigidBody::getAtomVel(Vector3d& vel, unsigned int index) {
    // velRot = $(A\cdot skew(I^{-1}j))^{T}refCoor$

    if (index < atoms_.size()) {
      Vector3d velRot;
      Mat3x3d skewMat;
      ;
      Vector3d ref = refCoords_[index];
      Vector3d ji  = getJ();
      Mat3x3d I    = getI();

      skewMat(0, 0) = 0;
      skewMat(0, 1) = ji[2] / I(2, 2);
      skewMat(0, 2) = -ji[1] / I(1, 1);

      skewMat(1, 0) = -ji[2] / I(2, 2);
      skewMat(1, 1) = 0;
      skewMat(1, 2) = ji[0] / I(0, 0);

      skewMat(2, 0) = ji[1] / I(1, 1);
      skewMat(2, 1) = -ji[0] / I(0, 0);
      skewMat(2, 2) = 0;

      velRot = (getA() * skewMat).transpose() * ref;

      vel = getVel() + velRot;
      return true;

    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "Index %d is an invalid index, the current rigid body contains %zu "
          "atoms.\n",
          index, atoms_.size());
      painCave.isFatal = 0;
      simError();
      return false;
    }
  }

  bool RigidBody::getAtomVel(Vector3d& vel, Atom* atom) {
    std::vector<Atom*>::iterator i;
    i = std::find(atoms_.begin(), atoms_.end(), atom);
    if (i != atoms_.end()) {
      return getAtomVel(vel, i - atoms_.begin());
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Atom %d does not belong to rigid body %d.\n",
               atom->getGlobalIndex(), getGlobalIndex());
      painCave.isFatal = 0;
      simError();
      return false;
    }
  }

  bool RigidBody::getAtomRefCoor(Vector3d& coor, unsigned int index) {
    if (index < atoms_.size()) {
      coor = refCoords_[index];
      return true;
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "Index %d is an invalid index, the current rigid body contains %zu "
          "atoms.\n",
          index, atoms_.size());
      painCave.isFatal = 0;
      simError();
      return false;
    }
  }

  bool RigidBody::getAtomRefCoor(Vector3d& coor, Atom* atom) {
    std::vector<Atom*>::iterator i;
    i = std::find(atoms_.begin(), atoms_.end(), atom);
    if (i != atoms_.end()) {
      // RigidBody class makes sure refCoords_ and atoms_ match each other
      coor = refCoords_[i - atoms_.begin()];
      return true;
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Atom %d does not belong to rigid body %d.\n",
               atom->getGlobalIndex(), getGlobalIndex());
      painCave.isFatal = 0;
      simError();
      return false;
    }
  }

  void RigidBody::addAtom(Atom* at, AtomStamp* ats) {
    Vector3d coords;
    Vector3d euler;

    atoms_.push_back(at);

    if (!ats->havePosition()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RigidBody error.\n"
               "\tAtom %s does not have a position specified.\n"
               "\tThis means RigidBody cannot set up reference coordinates.\n",
               ats->getType().c_str());
      painCave.isFatal = 1;
      simError();
    }

    coords[0] = ats->getPosX();
    coords[1] = ats->getPosY();
    coords[2] = ats->getPosZ();

    refCoords_.push_back(coords);

    RotMat3x3d identMat = RotMat3x3d::identity();

    if (at->isDirectional()) {
      if (!ats->haveOrientation()) {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "RigidBody error.\n"
            "\tAtom %s does not have an orientation specified.\n"
            "\tThis means RigidBody cannot set up reference orientations.\n",
            ats->getType().c_str());
        painCave.isFatal = 1;
        simError();
      }

      euler[0] = ats->getEulerPhi() * Constants::PI / 180.0;
      euler[1] = ats->getEulerTheta() * Constants::PI / 180.0;
      euler[2] = ats->getEulerPsi() * Constants::PI / 180.0;

      RotMat3x3d Atmp(euler);
      refOrients_.push_back(Atmp);

    } else {
      refOrients_.push_back(identMat);
    }
  }

}  // namespace OpenMD
