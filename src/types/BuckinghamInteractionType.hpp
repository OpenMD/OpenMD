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

#ifndef TYPES_BUCKINGHAMINTERACTIONTYPE_HPP
#define TYPES_BUCKINGHAMINTERACTIONTYPE_HPP

#include "types/NonBondedInteractionType.hpp"

namespace OpenMD {
  /**
   * @class BuckinghamInteractionType
   *
   * BuckinghamInteractionType is one of the basic non-bonded
   * interactions for representing the non-Coulombic contribution in
   * ionic and networked solids like silica and alumina. The
   * Traditional form has an exponential-6 form, which diverges to
   * negative infinity at short distances. The modified form corrects
   * the divergence with a short range 30-6 potential.
   *
   * Traditional:
   *
   * \f[ V = A \exp( -B r) - \frac{C}{r^6} \f]
   *
   * Modified:
   *
   * \f[ V = A \exp( -B r) - \frac{C}{r^6} + 4 \epsilon \left( \left(
   * \frac{\sigma}{r} \right)^{30} - \left( \frac{\sigma}{r} \right)^6 \right)
   * \f]
   */

  enum BuckinghamType { btTraditional, btModified, btUnknown };

  class BuckinghamInteractionType : public NonBondedInteractionType {
  public:
    BuckinghamInteractionType(RealType myA, RealType myB, RealType myC,
                              BuckinghamType myType) {
      A               = myA;
      B               = myB;
      C               = myC;
      interactionType = myType;
      setBuckingham();
    }

    BuckinghamInteractionType(RealType myA, RealType myB, RealType myC,
                              RealType mySigma, RealType myEpsilon,
                              BuckinghamType myType) {
      A               = myA;
      B               = myB;
      C               = myC;
      sigma           = mySigma;
      epsilon         = myEpsilon;
      interactionType = myType;
      setBuckingham();
    }

    RealType getA() { return A; }

    RealType getB() { return B; }

    RealType getC() { return C; }

    RealType getSigma() { return sigma; }

    RealType getEpsilon() { return epsilon; }

    BuckinghamType getInteractionType() { return interactionType; }

  private:
    RealType A;
    RealType B;
    RealType C;
    RealType sigma;
    RealType epsilon;
    BuckinghamType interactionType;
  };
}  // namespace OpenMD

#endif
