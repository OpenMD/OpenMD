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

/**
 * @file Vector3.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */

#ifndef MATH_VECTOR3_HPP
#define MATH_VECTOR3_HPP

#include <cassert>
#include <cmath>

#include "Vector.hpp"

namespace OpenMD {

  /**
   * @class Vector3 Vector3.hpp "math/Vector3.hpp"
   * @brief
   */
  template<typename Real>
  class Vector3 : public Vector<Real, 3> {
  public:
    using ElemType       = Real;
    using ElemPoinerType = Real*;

    Vector3() : Vector<Real, 3>() {}

    /** Constructs and initializes a Vector3 from x, y, z coordinates */
    inline Vector3(Real x, Real y, Real z) {
      this->data_[0] = x;
      this->data_[1] = y;
      this->data_[2] = z;
    }

    /** Constructs and initializes from an array*/
    inline Vector3(Real* array) : Vector<Real, 3>(array) {}

    inline Vector3(const Vector<Real, 3>& v) : Vector<Real, 3>(v) {}

    inline Vector3<Real>& operator=(const Vector<Real, 3>& v) {
      if (this == &v) { return *this; }
      Vector<Real, 3>::operator=(v);
      return *this;
    }

    /**
     * Returns reference of the first element of Vector3.
     * @return reference of the first element of Vector3
     */
    inline Real& x() { return this->data_[0]; }

    /**
     * Returns the first element of Vector3.
     * @return  the first element of Vector3
     */
    inline Real x() const { return this->data_[0]; }

    /**
     * Returns reference of the second element of Vector3.
     * @return reference of the second element of Vector3
     */
    inline Real& y() { return this->data_[1]; }

    /**
     * Returns  the second element of Vector3.
     * @return c the second element of Vector3
     */
    inline Real y() const { return this->data_[1]; }

    /**
     * Returns reference of the third element of Vector3.
     * @return reference of the third element of Vector3
     */
    inline Real& z() { return this->data_[2]; }

    /**
     * Returns  the third element of Vector3.
     * @return f the third element of Vector3
     */
    inline Real z() const { return this->data_[2]; }
  };

  /**
   * Returns the cross product of two Vectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the cross product  of v1 and v2
   */
  template<typename Real>
  inline Vector3<Real> cross(const Vector3<Real>& v1, const Vector3<Real>& v2) {
    Vector3<Real> result;

    result.x() = v1.y() * v2.z() - v1.z() * v2.y();
    result.y() = v1.z() * v2.x() - v1.x() * v2.z();
    result.z() = v1.x() * v2.y() - v1.y() * v2.x();

    return result;
  }

  /**
   * Returns the linear indexing for integer vectors. Compare to
   * Rapaport's VLinear
   *
   * @param p first vector
   * @param s second vector
   */
  template<typename Real>
  inline Real Vlinear(const Vector3<Real>& p, const Vector3<Real>& s) {
    return (p.z() * s.y() + p.y()) * s.x() + p.x();
  }

  using Vector3i = Vector3<int>;

  using Vector3d = Vector3<RealType>;

  const Vector3d V3Zero(0.0, 0.0, 0.0);
  const Vector3d V3X(1.0, 0.0, 0.0);
  const Vector3d V3Y(0.0, 1.0, 0.0);
  const Vector3d V3Z(0.0, 0.0, 1.0);
}  // namespace OpenMD

#endif
