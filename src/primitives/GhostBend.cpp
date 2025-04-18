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

#include "primitives/GhostBend.hpp"

#include <config.h>

#include <cmath>

#include "primitives/DirectionalAtom.hpp"
#include "utils/Constants.hpp"

namespace OpenMD {

  /**@todo still a lot left to improve*/
  void GhostBend::calcForce(RealType& angle, bool doParticlePot) {
    DirectionalAtom* ghostAtom = static_cast<DirectionalAtom*>(atoms_[1]);

    Vector3d pos1 = atoms_[0]->getPos();
    Vector3d pos2 = ghostAtom->getPos();

    Vector3d r21 = pos1 - pos2;
    snapshotMan_->getCurrentSnapshot()->wrapVector(r21);
    RealType d21 = r21.length();

    RealType d21inv = 1.0 / d21;

    // we need the transpose of A to get the lab fixed vector:
    Vector3d r23 = ghostAtom->getA().transpose().getColumn(2);
    RealType d23 = r23.length();

    RealType d23inv = 1.0 / d23;

    RealType cosTheta = dot(r21, r23) / (d21 * d23);

    // check roundoff
    if (cosTheta > 1.0) {
      cosTheta = 1.0;
    } else if (cosTheta < -1.0) {
      cosTheta = -1.0;
    }

    RealType theta = acos(cosTheta);

    RealType dVdTheta;

    bendType_->calcForce(theta, potential_, dVdTheta);

    RealType sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    if (fabs(sinTheta) < 1.0E-6) { sinTheta = 1.0E-6; }

    RealType commonFactor1 = dVdTheta / sinTheta * d21inv;
    RealType commonFactor2 = dVdTheta / sinTheta * d23inv;

    Vector3d force1 = commonFactor1 * (r23 * d23inv - r21 * d21inv * cosTheta);
    Vector3d force3 = commonFactor2 * (r21 * d21inv - r23 * d23inv * cosTheta);

    // Total force in current bend is zero

    atoms_[0]->addFrc(force1);
    ghostAtom->addFrc(-force1);

    ghostAtom->addTrq(cross(r23, force3));
    if (doParticlePot) {
      atoms_[0]->addParticlePot(potential_);
      ghostAtom->addParticlePot(potential_);
    }

    angle = theta / Constants::PI * 180.0;
  }
}  // namespace OpenMD
