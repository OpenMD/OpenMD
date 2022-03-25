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

/**
 * @file PolynomialBondType.hpp
 * @author    teng lin
 * @date  11/16/2004
 * @version 1.0
 */

#ifndef TYPES_POLYNOMIALBONDTYPE_HPP
#define TYPES_POLYNOMIALBONDTYPE_HPP

#include "math/Polynomial.hpp"
#include "types/BondType.hpp"

namespace OpenMD {

  /**
   * @class PolynomialBondType PolynomialBondType.hpp
   * "types/PolynomialBondType.hpp"
   * @todo documentation
   */
  class PolynomialBondType : public BondType {
  public:
    PolynomialBondType(RealType r0) : BondType(r0) {}

    void setCoefficient(int power, RealType coefficient) {
      polynomial_.setCoefficient(power, coefficient);
    }

    RealType getCoefficient(int power) {
      return polynomial_.getCoefficient(power);
    }

    void calcForce(RealType r, RealType& V, RealType& dVdr) {
      RealType delta = r - r0;
      V              = polynomial_.evaluate(delta);
      dVdr           = polynomial_.evaluateDerivative(delta);
    }

    friend std::ostream& operator<<(std::ostream& os, PolynomialBondType& pbt);

  private:
    DoublePolynomial polynomial_;
  };

  std::ostream& operator<<(std::ostream& os, PolynomialBondType& pbt) {
    DoublePolynomial::const_iterator i;

    i = pbt.polynomial_.begin();

    if (i == pbt.polynomial_.end()) {
      os << "This PolynomialBondType contains nothing" << std::endl;
      return os;
    }

    os << "This PolynomialBondType contains below terms:" << std::endl;

    while (true) {
      os << i->second << "*"
         << "(r - " << pbt.getEquilibriumBondLength() << ")"
         << "^" << i->first;

      if (++i == pbt.polynomial_.end()) {
        // If we reach the end of the polynomial pair, write out a
        // newline and then escape the loop
        os << std::endl;
        break;
      } else {
        // otherwise, write out a "+"
        os << " + ";
      }
    }

    return os;
  }

}  // namespace OpenMD

#endif  // TYPES_POLYNOMIALBONDTYPE_HPP
