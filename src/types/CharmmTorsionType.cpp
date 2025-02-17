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

#include "types/CharmmTorsionType.hpp"

#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>

#include "math/ChebyshevT.hpp"
#include "math/ChebyshevU.hpp"

namespace OpenMD {

  CharmmTorsionType::CharmmTorsionType(
      std::vector<CharmmTorsionParameter>& parameters) :
      TorsionType(),
      C_(0.0) {
    std::vector<CharmmTorsionParameter>::iterator i;
    i = std::max_element(parameters.begin(), parameters.end(),
                         LessThanPeriodicityFunctor());
    if (i != parameters.end()) {
      int maxPower = i->n;
      ChebyshevT T(maxPower);
      ChebyshevU U(maxPower);

      // convert parameters of charmm type torsion into
      // Polynomial parameters

      for (i = parameters.begin(); i != parameters.end(); ++i) {
        DoublePolynomial cosTerm = T.getChebyshevPolynomial(i->n);
        cosTerm *= (cos(i->delta) * i->kchi);

        // should check that i->n is >= 1
        DoublePolynomial sinTerm = U.getChebyshevPolynomial(i->n - 1);
        sinTerm *= -(sin(i->delta) * i->kchi);

        T_ += cosTerm;
        U_ += sinTerm;
        C_ += i->kchi;
      }
    }
  }

  void CharmmTorsionType::calcForce(RealType cosPhi, RealType& V,
                                    RealType& dVdCosPhi) {
    // check roundoff
    if (cosPhi > 1.0) {
      cosPhi = 1.0;
    } else if (cosPhi < -1.0) {
      cosPhi = -1.0;
    }

    RealType sinPhi = sqrt(1.0 - cosPhi * cosPhi);

    // trick to avoid divergence in angles near 0 and pi:

    if (fabs(sinPhi) < 1.0E-6) { sinPhi = copysign(1.0E-6, sinPhi); }

    V         = C_ + T_.evaluate(cosPhi) + U_.evaluate(cosPhi) * sinPhi;
    dVdCosPhi = T_.evaluateDerivative(cosPhi);
    // Chain rule for U * sinPhi term:
    dVdCosPhi += U_.evaluateDerivative(cosPhi) * sinPhi;
    dVdCosPhi += U_.evaluate(cosPhi) / (2.0 * sinPhi);
  }
}  // namespace OpenMD
