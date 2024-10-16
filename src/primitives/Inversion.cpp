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

#include "primitives/Inversion.hpp"

#include <config.h>

#include <cmath>

#include "utils/Constants.hpp"

namespace OpenMD {

  Inversion::Inversion(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4,
                       InversionType* it) :
      ShortRangeInteraction(),
      inversionType_(it) {
    atoms_.resize(4);
    atoms_[0] = atom1;
    atoms_[1] = atom2;
    atoms_[2] = atom3;
    atoms_[3] = atom4;

    inversionKey_ = inversionType_->getKey();
  }

  void Inversion::calcForce(RealType& angle, bool doParticlePot) {
    // In OpenMD's version of an inversion, the central atom
    // comes first.  However, to get the planarity in a typical cosine
    // version of this potential (i.e. Amber-style), the central atom
    // is treated as atom *3* in a standard torsion form:

    Vector3d pos1 = atoms_[1]->getPos();
    Vector3d pos2 = atoms_[2]->getPos();
    Vector3d pos3 = atoms_[0]->getPos();
    Vector3d pos4 = atoms_[3]->getPos();

    Vector3d r31 = pos1 - pos3;
    snapshotMan_->getCurrentSnapshot()->wrapVector(r31);
    Vector3d r23 = pos3 - pos2;
    snapshotMan_->getCurrentSnapshot()->wrapVector(r23);
    Vector3d r43 = pos3 - pos4;
    snapshotMan_->getCurrentSnapshot()->wrapVector(r43);

    //  Calculate the cross products and distances
    Vector3d A  = cross(r31, r43);
    RealType rA = A.length();
    Vector3d B  = cross(r43, r23);
    RealType rB = B.length();

    A.normalize();
    B.normalize();

    //  Calculate the sin and cos
    RealType cos_phi = dot(A, B);
    if (cos_phi > 1.0) cos_phi = 1.0;
    if (cos_phi < -1.0) cos_phi = -1.0;

    RealType dVdcosPhi;
    switch (inversionKey_) {
    case itCosAngle:
      inversionType_->calcForce(cos_phi, potential_, dVdcosPhi);
      break;
    case itAngle:
      RealType phi = acos(cos_phi);
      RealType dVdPhi;
      inversionType_->calcForce(phi, potential_, dVdPhi);
      RealType sin_phi = sqrt(1.0 - cos_phi * cos_phi);
      if (fabs(sin_phi) < 1.0E-6) { sin_phi = 1.0E-6; }
      dVdcosPhi = -dVdPhi / sin_phi;

      break;
    }

    Vector3d f1;
    Vector3d f2;
    Vector3d f3;

    Vector3d dcosdA = (cos_phi * A - B) / rA;
    Vector3d dcosdB = (cos_phi * B - A) / rB;

    f1 = dVdcosPhi * cross(r43, dcosdA);
    f2 = dVdcosPhi * (cross(r23, dcosdB) - cross(r31, dcosdA));
    f3 = dVdcosPhi * cross(dcosdB, r43);

    // In OpenMD's version of an improper torsion, the central atom
    // comes first.  However, to get the planarity in a typical cosine
    // version of this potential (i.e. Amber-style), the central atom
    // is treated as atom *3* in a standard torsion form:

    //  AMBER:   I - J - K - L   (e.g. K is sp2 hybridized carbon)
    //  OpenMD:  I - (J - K - L)  (e.g. I is sp2 hybridized carbon)

    // Confusing enough?  Good.

    atoms_[1]->addFrc(f1);
    atoms_[0]->addFrc(f2 - f1 + f3);
    atoms_[3]->addFrc(-f2);
    atoms_[2]->addFrc(-f3);

    if (doParticlePot) {
      atoms_[0]->addParticlePot(potential_);
      atoms_[1]->addParticlePot(potential_);
      atoms_[2]->addParticlePot(potential_);
      atoms_[3]->addParticlePot(potential_);
    }

    angle = acos(cos_phi) / Constants::PI * 180.0;
  }

}  // namespace OpenMD
