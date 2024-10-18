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

#include "primitives/GhostTorsion.hpp"

#include <config.h>

#include <cmath>

#include "utils/Constants.hpp"

namespace OpenMD {

  GhostTorsion::GhostTorsion(Atom* atom1, Atom* atom2,
                             DirectionalAtom* ghostAtom, TorsionType* tt) :
      Torsion(atom1, atom2, ghostAtom, ghostAtom, tt) {}

  void GhostTorsion::calcForce(RealType& angle, bool doParticlePot) {
    DirectionalAtom* ghostAtom = static_cast<DirectionalAtom*>(atoms_[2]);

    Vector3d pos1 = atoms_[0]->getPos();
    Vector3d pos2 = atoms_[1]->getPos();
    Vector3d pos3 = ghostAtom->getPos();

    Vector3d r21 = pos1 - pos2;
    snapshotMan_->getCurrentSnapshot()->wrapVector(r21);
    Vector3d r32 = pos2 - pos3;
    snapshotMan_->getCurrentSnapshot()->wrapVector(r32);
    Vector3d r43 = ghostAtom->getA().transpose().getColumn(2);

    //  Calculate the cross products and distances
    Vector3d A  = cross(r21, r32);
    RealType rA = A.length();
    Vector3d B  = cross(r32, r43);
    RealType rB = B.length();

    /*
     If either of the two cross product vectors is tiny, that means
     the three atoms involved are colinear, and the torsion angle is
     going to be undefined.  The easiest check for this problem is
     to use the product of the two lengths.
  */
    if (rA * rB < OpenMD::epsilon) return;

    A.normalize();
    B.normalize();

    //  Calculate the sin and cos
    RealType cos_phi = dot(A, B);

    RealType dVdcosPhi;
    torsionType_->calcForce(cos_phi, potential_, dVdcosPhi);

    Vector3d dcosdA = (cos_phi * A - B) / rA;
    Vector3d dcosdB = (cos_phi * B - A) / rB;

    Vector3d f1 = dVdcosPhi * cross(r32, dcosdA);
    Vector3d f2 = dVdcosPhi * (cross(r43, dcosdB) - cross(r21, dcosdA));
    Vector3d f3 = dVdcosPhi * cross(dcosdB, r32);

    atoms_[0]->addFrc(f1);
    atoms_[1]->addFrc(f2 - f1);

    ghostAtom->addFrc(-f2);

    f3.negate();
    ghostAtom->addTrq(cross(r43, f3));

    if (doParticlePot) {
      atoms_[0]->addParticlePot(potential_);
      atoms_[1]->addParticlePot(potential_);
      ghostAtom->addParticlePot(potential_);
    }

    angle = acos(cos_phi) / Constants::PI * 180.0;
  }
}  // namespace OpenMD
