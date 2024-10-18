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

#ifndef MATH_VECTOR2_HPP
#define MATH_VECTOR2_HPP

#include <cassert>
#include <cmath>

#include "Vector.hpp"

namespace OpenMD {

  template<typename Real>
  class Vector2 : public Vector<Real, 2> {
  public:
    typedef Real ElemType;
    typedef Real* ElemPoinerType;
    Vector2() : Vector<Real, 2>() {}

    /** Constructs and initializes a Vector2 from x and y coordinates */
    inline Vector2(Real x, Real y) {
      this->data_[0] = x;
      this->data_[1] = y;
    }

    /** Constructs and initializes from an array*/
    inline Vector2(Real* array) : Vector<Real, 2>(array) {}

    inline Vector2(const Vector<Real, 2>& v) : Vector<Real, 2>(v) {}

    inline Vector2<Real>& operator=(const Vector<Real, 2>& v) {
      if (this == &v) { return *this; }
      Vector<Real, 2>::operator=(v);
      return *this;
    }

    /**
     * Returns reference of the first element of Vector2.
     * @return reference of the first element of Vector2
     */
    inline Real& x() { return this->data_[0]; }

    /**
     * Returns the first element of Vector2.
     * @return  the first element of Vector2
     */
    inline Real x() const { return this->data_[0]; }

    /**
     * Returns reference of the second element of Vector2.
     * @return reference of the second element of Vector2
     */
    inline Real& y() { return this->data_[1]; }

    /**
     * Returns  the second element of Vector2.
     * @return c the second element of Vector2
     */
    inline Real y() const { return this->data_[1]; }
  };

  /**
   * Returns the linear indexing for size_t vectors. Compare to
   * Rapaport's VLinear
   *
   * @param p first vector
   * @param s second vector
   */
  inline std::size_t Vlinear(const Vector2<std::size_t>& p,
                             const Vector2<std::size_t>& s) {
    return std::size_t(p.y() * s.x() + p.x());
  }

  typedef Vector2<int> Vector2i;

  typedef Vector2<RealType> Vector2d;

  const Vector2d V2Zero(0.0, 0.0);
  const Vector2d V2X(1.0, 0.0);
  const Vector2d V2Y(0.0, 1.0);

}  // namespace OpenMD

#endif
