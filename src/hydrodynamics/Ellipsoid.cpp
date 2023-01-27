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

#include "hydrodynamics/Ellipsoid.hpp"

#include "math/LU.hpp"
#include "utils/Constants.hpp"

namespace OpenMD {

  Ellipsoid::Ellipsoid(Vector3d origin, RealType rAxial, RealType rEquatorial,
                       Mat3x3d rotMat) :
      origin_(origin),
      rAxial_(rAxial), rEquatorial_(rEquatorial), rotMat_(rotMat) {
    if (rAxial_ > rEquatorial_) {
      rMajor_ = rAxial_;
      rMinor_ = rEquatorial_;
    } else {
      rMajor_ = rEquatorial_;
      rMinor_ = rAxial_;
    }
  }

  bool Ellipsoid::isInterior(Vector3d pos) {
    Vector3d r     = pos - origin_;
    Vector3d rbody = rotMat_ * r;

    RealType xoverb = rbody[0] / rEquatorial_;
    RealType yoverb = rbody[1] / rEquatorial_;
    RealType zovera = rbody[2] / rAxial_;

    bool result;
    if (xoverb * xoverb + yoverb * yoverb + zovera * zovera < 1)
      result = true;
    else
      result = false;

    return result;
  }

  std::pair<Vector3d, Vector3d> Ellipsoid::getBoundingBox() {
    std::pair<Vector3d, Vector3d> boundary;
    // make a cubic box
    RealType rad = rAxial_ > rEquatorial_ ? rAxial_ : rEquatorial_;
    Vector3d r(rad, rad, rad);
    boundary.first  = origin_ - r;
    boundary.second = origin_ + r;
    return boundary;
  }

  HydroProp* Ellipsoid::getHydroProp(RealType viscosity) {
    RealType a  = rAxial_;
    RealType b  = rEquatorial_;
    RealType a2 = a * a;
    RealType b2 = b * b;

    RealType p = a / b;
    RealType S;
    if (p > 1.0) {
      // Ellipsoid is prolate:
      S = 2.0 / sqrt(a2 - b2) * log((a + sqrt(a2 - b2)) / b);
    } else {
      // Ellipsoid is oblate:
      S = 2.0 / sqrt(b2 - a2) * atan(sqrt(b2 - a2) / a);
    }

    RealType pi = Constants::PI;
    RealType XittA =
        16.0 * pi * viscosity * (a2 - b2) / ((2.0 * a2 - b2) * S - 2.0 * a);
    RealType XittB = 32.0 * pi * viscosity * (a2 - b2) /
                     ((2.0 * a2 - 3.0 * b2) * S + 2.0 * a);
    RealType XirrA =
        32.0 / 3.0 * pi * viscosity * (a2 - b2) * b2 / (2.0 * a - b2 * S);
    RealType XirrB = 32.0 / 3.0 * pi * viscosity * (a2 * a2 - b2 * b2) /
                     ((2.0 * a2 - b2) * S - 2.0 * a);

    Mat6x6d Xi;

    Xi(0, 0) = XittB;
    Xi(1, 1) = XittB;
    Xi(2, 2) = XittA;
    Xi(3, 3) = XirrB;
    Xi(4, 4) = XirrB;
    Xi(5, 5) = XirrA;

    Xi *= Constants::viscoConvert;

    HydroProp* hprop = new HydroProp(V3Zero, Xi);
    hprop->setName(getName());

    return hprop;
  }
}  // namespace OpenMD
