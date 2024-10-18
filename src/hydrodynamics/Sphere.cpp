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

#include "hydrodynamics/Sphere.hpp"

#include "math/LU.hpp"
#include "utils/Constants.hpp"

namespace OpenMD {

  Sphere::Sphere(Vector3d origin, RealType radius) :
      origin_(origin), radius_(radius) {}

  bool Sphere::isInterior(Vector3d pos) {
    Vector3d r = pos - origin_;

    bool result;
    if (r.length() < radius_)
      result = true;
    else
      result = false;

    return result;
  }

  std::pair<Vector3d, Vector3d> Sphere::getBoundingBox() {
    std::pair<Vector3d, Vector3d> boundary;
    Vector3d r(radius_, radius_, radius_);
    boundary.first  = origin_ - r;
    boundary.second = origin_ + r;
    return boundary;
  }

  HydroProp* Sphere::getHydroProp(RealType viscosity) {
    RealType Xitt = 6.0 * Constants::PI * viscosity * radius_;
    RealType Xirr = 8.0 * Constants::PI * viscosity * pow(radius_, 3);

    Mat6x6d Xi;

    Xi(0, 0) = Xitt;
    Xi(1, 1) = Xitt;
    Xi(2, 2) = Xitt;
    Xi(3, 3) = Xirr;
    Xi(4, 4) = Xirr;
    Xi(5, 5) = Xirr;

    Xi *= Constants::viscoConvert;

    HydroProp* hprop = new HydroProp(V3Zero, Xi);
    hprop->setName(getName());
    return hprop;
  }

}  // namespace OpenMD
