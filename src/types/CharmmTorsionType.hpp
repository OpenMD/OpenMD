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

#ifndef TYPES_CHARMMTORSIONTYPE_HPP
#define TYPES_CHARMMTORSIONTYPE_HPP

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "math/Polynomial.hpp"
#include "types/TorsionType.hpp"

namespace OpenMD {

  struct CharmmTorsionParameter {
    RealType kchi;
    int n;
    RealType delta;
  };

  class LessThanPeriodicityFunctor {
  public:
    bool operator()(const CharmmTorsionParameter& p1,
                    const CharmmTorsionParameter& p2) {
      return p1.n < p2.n;
    }
  };

  /**
   * @class CharmmTorsionType CharmmTorsionType.hpp
   * "types/CharmmTorsionType.hpp" These torsion types are defined identically
   * with functional form given in the following paper:
   *
   * "A. D. MacKerell, Jr. et al., CHARMM: The energy function and its
   * parameterization with an overview of the program," in The
   * Encyclopedia of Computational Chemistry, edited by
   * P. v. R. Schleyer, et al., volume 1, pages 271–277, John Wiley &
   * Sons, New York, 1998.
   *
   * This torsion potential has the form:
   *
   * \f[ V_{\mathrm{torsion}}(\phi) = \sum_n K_n \left( 1 + \cos(n \phi -
   * \delta_n) \right) \f]
   *
   * Notes:
   *
   * 1. OpenMD converts internally to Chebyshev polynomials for
   *    computational efficiency.
   * 2. Coefficients \f$ K_n \f$ are assumed to be in kcal / mol.
   * 3. Phase angles \f$ \delta_n \f$ are assumed to be in degrees.
   * 4. Periodicity values \f$ n \f$ are positive integers.
   *
   * Internally convert CHARMM torsion functions to two polynomials
   * based on Chebyshev polynomials in cos(phi):
   *
   * \f[ V_{\mathrm{torsion}}(\phi) = \sum_n K_n + \sum_n K_n \cos(\delta_n)
   * T_n(\cos(\phi)) - \sum_n K_n \sin(\delta_n) U_{n-1}((\cos \phi)) \sin(\phi)
   * \f]
   *
   * This conversion has used the cosine addition formula, and two
   * identities of Chebyshev polynomials:
   *
   * \f[ T_n (\cos \phi) = \cos(n \phi) \f]
   *
   * for Chebyshev polynomials of the first type, and:
   *
   * \f[ U_{n-1} (\cos \phi) \sin(\phi) = \sin( n \phi ) \f]
   *
   * for Chebyshev polynomials of the second type. We're left with a
   * simpler equation for the torsion potential in terms of only
   * polynomials of the cosine and an additional sine of the angle:
   *
   * \f[ V_{\mathrm{torsion}}(\phi) = C + T(\cos(\phi)) + U(\cos(\phi)) \sin(\phi)
   * \f] where: \f[ C = \sum_n K_n \f] \f[ T(\cos(\phi)) = \sum_n K_n
   * \cos(\delta_n) T_n(\cos(\phi)) \f] \f[ U(\cos(\phi)) = \sum_n -K_n
   * \sin(\delta_n) U_{n-1}(\cos(\phi)) \f]
   */
  class CharmmTorsionType : public TorsionType {
  public:
    CharmmTorsionType(std::vector<CharmmTorsionParameter>& parameters);
    virtual void calcForce(RealType cosPhi, RealType& V, RealType& dVdCosPhi);

  private:
    DoublePolynomial T_;
    DoublePolynomial U_;
    RealType C_;
  };
}  // namespace OpenMD

#endif  // TYPES_CHARMMTORSIONTYPE_HPP
