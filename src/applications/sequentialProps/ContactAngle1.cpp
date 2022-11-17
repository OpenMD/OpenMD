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

#include "applications/sequentialProps/ContactAngle1.hpp"

#include <algorithm>
#include <functional>
#include <sstream>

#include "io/DumpReader.hpp"
#include "math/Polynomial.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  ContactAngle1::ContactAngle1(SimInfo* info, const std::string& filename,
                               const std::string& sele1,
                               const std::string& sele2, RealType solidZ,
                               RealType dropletRadius) :
      SequentialAnalyzer(info, filename, sele1, sele2),
      solidZ_(solidZ), dropletRadius_(dropletRadius) {
    setOutputName(getPrefix(filename) + ".ca1");

    std::stringstream params;
    params << " solid Z = " << solidZ_
           << ", droplet radius = " << dropletRadius_;

    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  void ContactAngle1::doFrame(int) {
    StuntDouble* sd;
    int i;

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    RealType mtot = 0.0;
    Vector3d com(V3Zero);
    RealType mass;

    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      mass = sd->getMass();
      mtot += mass;
      com += sd->getPos() * mass;
    }

    com /= mtot;

    RealType dz = com.z() - solidZ_;

    if (dz < 0.0) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ContactAngle1: Z-center of mass of selection, %lf, was\n"
               "\tlocated below the solid reference plane, %lf\n",
               com.z(), solidZ_);
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (dz > dropletRadius_) {
      values_.push_back(180.0);
    } else {
      RealType k = pow(2.0, -4.0 / 3.0) * dropletRadius_;

      RealType z2 = dz * dz;
      RealType z3 = z2 * dz;
      RealType k2 = k * k;
      RealType k3 = k2 * k;

      Polynomial<RealType> poly;
      poly.setCoefficient(4, z3 + k3);
      poly.setCoefficient(3, 8.0 * z3 + 8.0 * k3);
      poly.setCoefficient(2, 24.0 * z3 + 18.0 * k3);
      poly.setCoefficient(1, 32.0 * z3);
      poly.setCoefficient(0, 16.0 * z3 - 27.0 * k3);
      vector<RealType> realRoots = poly.FindRealRoots();

      RealType ct;

      vector<RealType>::iterator ri;

      RealType maxct = -1.0;
      for (ri = realRoots.begin(); ri != realRoots.end(); ++ri) {
        ct = *ri;
        if (ct > 1.0) ct = 1.0;
        if (ct < -1.0) ct = -1.0;

        // use the largest magnitude of ct that it finds:
        if (ct > maxct) { maxct = ct; }
      }

      values_.push_back(acos(maxct) * (180.0 / Constants::PI));
    }
  }
}  // namespace OpenMD
