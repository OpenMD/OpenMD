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
    hprop->setName( getName() );
    return hprop;
  }

}  // namespace OpenMD
