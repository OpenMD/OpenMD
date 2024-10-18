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

/**
 * @file ChebyshevU.hpp
 * @author    teng lin
 * @date  11/16/2004
 * @version 1.0
 */

#ifndef MATH_CHEBYSHEVU_HPP
#define MATH_CHEBYSHEVU_HPP

#include <cassert>
#include <vector>

#include "math/Polynomial.hpp"

namespace OpenMD {

  /**
   * @class ChebyshevU
   * A collection of Chebyshev Polynomials.
   * @todo document
   */
  class ChebyshevU {
  public:
    ChebyshevU(int maxPower);
    virtual ~ChebyshevU() {}
    /**
     * Calculates the value of the nth Chebyshev Polynomial evaluated at the
     * given x value.
     * @return The value of the nth Chebyshev Polynomial evaluates at the given
     * x value
     * @param n
     * @param x the value of the independent variable for the nth Chebyshev
     * Polynomial  function
     */

    RealType evaluate(int n, RealType x) {
      assert(n <= maxPower_ && n >= 0);
      return polyList_[n].evaluate(x);
    }

    /**
     * Returns the first derivative of the nth Chebyshev Polynomial.
     * @return the first derivative of the nth Chebyshev Polynomial
     * @param n
     * @param x the value of the independent variable for the nth Chebyshev
     * Polynomial  function
     */
    RealType evaluateDerivative(int n, RealType x) {
      assert(n <= maxPower_ && n >= 0);
      return polyList_[n].evaluateDerivative(x);
    }

    /**
     * Returns the nth Chebyshev Polynomial
     * @return the nth Chebyshev Polynomial
     * @param n
     */
    const DoublePolynomial& getChebyshevPolynomial(int n) const {
      assert(n <= maxPower_ && n >= 0);
      return polyList_[n];
    }

  protected:
    std::vector<DoublePolynomial> polyList_;
    void GeneratePolynomials(int maxPower);

  private:
    void GenerateFirstTwoTerms();

    int maxPower_;
  };
}  // namespace OpenMD

#endif  // MATH_CHEBYSHEVU_HPP
