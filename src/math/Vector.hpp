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
 * @file Vector.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */

#ifndef MATH_VECTOR_HPP
#define MATH_VECTOR_HPP

#include <config.h>

#include <cassert>
#include <cmath>
#include <iostream>

namespace OpenMD {

  static constexpr RealType epsilon = 0.000001;

  template<typename T>
  inline bool equal(T e1, T e2) {
    if constexpr (std::is_same_v<T, RealType>)
      return std::fabs(e1 - e2) < epsilon;
    else
      return e1 == e2;
  }

  /**
   * @class Vector Vector.hpp "math/Vector.hpp"
   * @brief Fix length vector class
   */
  template<typename Real, unsigned int Dim>
  class Vector {
  public:
    using ElemType       = Real;
    using ElemPoinerType = Real*;

    /** default constructor */
    inline Vector() {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = 0;
    }

    /** Constructs and initializes a Vector from a vector */
    inline Vector(const Vector<Real, Dim>& v) { *this = v; }

    /** copy assignment operator */
    inline Vector<Real, Dim>& operator=(const Vector<Real, Dim>& v) {
      if (this == &v) return *this;

      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = v[i];

      return *this;
    }

    /** array assignment operator */
    inline Vector<Real, Dim>& operator=(const Real* v) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = v[i];

      return *this;
    }

    // template<typename T>
    // inline Vector(const T& s){
    inline Vector(const Real& s) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = s;
    }

    /** Constructs and initializes a Vector from an array */
    inline Vector(Real* v) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = v[i];
    }

    /**
     * Returns reference of ith element.
     * @return reference of ith element
     * @param i index
     */
    inline Real& operator[](unsigned int i) {
      assert(i < Dim);
      return this->data_[i];
    }

    /**
     * Returns reference of ith element.
     * @return reference of ith element
     * @param i index
     */
    inline Real& operator()(unsigned int i) {
      assert(i < Dim);
      return this->data_[i];
    }

    /**
     * Returns constant reference of ith element.
     * @return reference of ith element
     * @param i index
     */
    inline const Real& operator[](unsigned int i) const {
      assert(i < Dim);
      return this->data_[i];
    }

    /**
     * Returns constant reference of ith element.
     * @return reference of ith element
     * @param i index
     */
    inline const Real& operator()(unsigned int i) const {
      assert(i < Dim);
      return this->data_[i];
    }

    /** Copy the internal data to an array*/
    void getArray(Real* array) {
      for (unsigned int i = 0; i < Dim; i++) {
        array[i] = this->data_[i];
      }
    }

    /** Returns the pointer of internal array */
    Real* getArrayPointer() { return this->data_; }

    /**
     * Tests if this vetor is equal to other vector
     * @return true if equal, otherwise return false
     * @param v vector to be compared
     */
    inline bool operator==(const Vector<Real, Dim>& v) {
      for (unsigned int i = 0; i < Dim; i++) {
        if (!equal(this->data_[i], v[i])) { return false; }
      }

      return true;
    }

    /**
     * Tests if this vetor is not equal to other vector
     * @return true if equal, otherwise return false
     * @param v vector to be compared
     */
    inline bool operator!=(const Vector<Real, Dim>& v) { return !(*this == v); }

    /** Zeros out the values in this vector in place */
    inline void zero() {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = 0;
    }

    /** Negates the value of this vector in place. */
    inline void negate() {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = -this->data_[i];
    }

    /**
     * Sets the value of this vector to the negation of vector v1.
     * @param v1 the source vector
     */
    inline void negate(const Vector<Real, Dim>& v1) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = -v1.data_[i];
    }

    /**
     * Sets the value of this vector to the sum of itself and v1 (*this += v1).
     * @param v1 the other vector
     */
    inline void add(const Vector<Real, Dim>& v1) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] += v1.data_[i];
    }

    /**
     * Sets the value of this vector to the sum of v1 and v2 (*this = v1 + v2).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    inline void add(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = v1.data_[i] + v2.data_[i];
    }

    /**
     * Sets the value of this vector to the difference  of itself and v1 (*this
     * -= v1).
     * @param v1 the other vector
     */
    inline void sub(const Vector<Real, Dim>& v1) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] -= v1.data_[i];
    }

    /**
     * Sets the value of this vector to the difference of vector v1 and v2
     * (*this = v1 - v2).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    inline void sub(const Vector<Real, Dim>& v1, const Vector& v2) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = v1.data_[i] - v2.data_[i];
    }

    /**
     * Sets the value of this vector to the scalar multiplication of itself
     * (*this *= s).
     * @param s the scalar value
     */
    inline void mul(Real s) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] *= s;
    }

    /**
     * Sets the value of this vector to the scalar multiplication of vector v1
     * (*this = s * v1).
     * @param v1 the vector
     * @param s the scalar value
     */
    inline void mul(const Vector<Real, Dim>& v1, Real s) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = s * v1.data_[i];
    }

    /**
     * Sets the elements of this vector to the multiplication of
     * elements of two other vectors.  Not to be confused with scalar
     * multiplication (mul) or dot products.
     *
     * (*this.data_[i] =  v1.data_[i] * v2.data_[i]).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    inline void Vmul(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = v1.data_[i] * v2.data_[i];
    }

    /* replaces the elements with the absolute values of those elements */
    inline Vector<Real, Dim>& abs() {
      for (unsigned int i = 0; i < Dim; i++) {
        this->data_[i] = std::abs(this->data_[i]);
      }
      return *this;
    }

    /* returns the maximum value in this vector */
    inline Real max() {
      Real val = this->data_[0];
      for (unsigned int i = 0; i < Dim; i++) {
        if (this->data_[i] > val) val = this->data_[i];
      }
      return val;
    }

    /**
     * Sets the value of this vector to the scalar division of itself  (*this /=
     * s ).
     * @param s the scalar value
     */
    inline void div(Real s) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] /= s;
    }

    /**
     * Sets the value of this vector to the scalar division of vector v1  (*this
     * = v1 / s ).
     * @param v1 the source vector
     * @param s the scalar value
     */
    inline void div(const Vector<Real, Dim>& v1, Real s) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = v1.data_[i] / s;
    }

    /**
     * Sets the elements of this vector to the division of
     * elements of two other vectors.  Not to be confused with scalar
     * division (div)
     *
     * (*this.data_[i] =  v1.data_[i] / v2.data_[i]).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    inline void Vdiv(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
      for (unsigned int i = 0; i < Dim; i++)
        this->data_[i] = v1.data_[i] / v2.data_[i];
    }

    /** @see #add */
    inline Vector<Real, Dim>& operator+=(const Vector<Real, Dim>& v1) {
      add(v1);
      return *this;
    }

    /** @see #sub */
    inline Vector<Real, Dim>& operator-=(const Vector<Real, Dim>& v1) {
      sub(v1);
      return *this;
    }

    /** @see #mul */
    inline Vector<Real, Dim>& operator*=(Real s) {
      mul(s);
      return *this;
    }

    /** @see #div */
    inline Vector<Real, Dim>& operator/=(Real s) {
      div(s);
      return *this;
    }

    /**
     * Returns the sum of all elements of this vector.
     * @return the sum of all elements of this vector
     */
    inline Real sum() {
      Real tmp;
      tmp = 0;
      for (unsigned int i = 0; i < Dim; i++)
        tmp += this->data_[i];
      return tmp;
    }

    /**
     * Returns the product of all elements of this vector.
     * @return the product of all elements of this vector
     */
    inline Real componentProduct() {
      Real tmp;
      tmp = 1;
      for (unsigned int i = 0; i < Dim; i++)
        tmp *= this->data_[i];
      return tmp;
    }

    /**
     * Returns the length of this vector.
     * @return the length of this vector
     */
    inline Real length() { return sqrt(lengthSquare()); }

    /**
     * Returns the squared length of this vector.
     * @return the squared length of this vector
     */
    inline Real lengthSquare() { return dot(*this, *this); }

    /** Normalizes this vector in place */
    inline void normalize() {
      Real len;

      len = length();

      // if (len < OpenMD::Constants::epsilon)
      //  throw();

      *this /= len;
    }

    /**
     * Tests if this vector is normalized
     * @return true if this vector is normalized, otherwise return false
     */
    inline bool isNormalized() { return equal(lengthSquare(), (RealType)1); }

    unsigned int size() const { return Dim; }

  protected:
    Real data_[Dim] {};
  };

  /** unary minus*/
  template<typename Real, unsigned int Dim>
  inline Vector<Real, Dim> operator-(const Vector<Real, Dim>& v1) {
    Vector<Real, Dim> tmp(v1);
    tmp.negate();
    return tmp;
  }

  /**
   * Return the sum of two vectors  (v1 - v2).
   * @return the sum of two vectors
   * @param v1 the first vector
   * @param v2 the second vector
   */
  template<typename Real, unsigned int Dim>
  inline Vector<Real, Dim> operator+(const Vector<Real, Dim>& v1,
                                     const Vector<Real, Dim>& v2) {
    Vector<Real, Dim> result;

    result.add(v1, v2);
    return result;
  }

  /**
   * Return the difference of two vectors  (v1 - v2).
   * @return the difference of two vectors
   * @param v1 the first vector
   * @param v2 the second vector
   */
  template<typename Real, unsigned int Dim>
  Vector<Real, Dim> operator-(const Vector<Real, Dim>& v1,
                              const Vector<Real, Dim>& v2) {
    Vector<Real, Dim> result;
    result.sub(v1, v2);
    return result;
  }

  /**
   * Returns the vaule of scalar multiplication of this vector v1 (v1 * r).
   * @return  the vaule of scalar multiplication of this vector
   * @param v1 the source vector
   * @param s the scalar value
   */
  template<typename Real, unsigned int Dim>
  Vector<Real, Dim> operator*(const Vector<Real, Dim>& v1, Real s) {
    Vector<Real, Dim> result;
    result.mul(v1, s);
    return result;
  }

  /**
   * Returns the vaule of scalar multiplication of this vector v1 (v1 * r).
   * @return  the vaule of scalar multiplication of this vector
   * @param s the scalar value
   * @param v1 the source vector
   */
  template<typename Real, unsigned int Dim>
  Vector<Real, Dim> operator*(Real s, const Vector<Real, Dim>& v1) {
    Vector<Real, Dim> result;
    result.mul(v1, s);
    return result;
  }

  /**
   * Returns the  value of division of a vector by a scalar.
   * @return  the vaule of scalar division of this vector
   * @param v1 the source vector
   * @param s the scalar value
   */
  template<typename Real, unsigned int Dim>
  Vector<Real, Dim> operator/(const Vector<Real, Dim>& v1, Real s) {
    Vector<Real, Dim> result;
    result.div(v1, s);
    return result;
  }

  /**
   * Returns the dot product of two Vectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the dot product of v1 and v2
   */
  template<typename Real, unsigned int Dim>
  inline Real dot(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
    Real tmp;
    tmp = 0;

    for (unsigned int i = 0; i < Dim; i++)
      tmp += v1[i] * v2[i];

    return tmp;
  }

  /**
   * Returns the wide dot product of three Vectors.  Compare with
   * Rapaport's VWDot function.
   *
   * @param v1 first vector
   * @param v2 second vector
   * @param v3 third vector
   * @return the wide dot product of v1, v2, and v3.
   */
  template<typename Real, unsigned int Dim>
  inline Real dot(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2,
                  const Vector<Real, Dim>& v3) {
    Real tmp;
    tmp = 0;

    for (unsigned int i = 0; i < Dim; i++)
      tmp += v1[i] * v2[i] * v3[i];

    return tmp;
  }

  /**
   * Returns the distance between  two Vectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the distance between v1 and v2
   */
  template<typename Real, unsigned int Dim>
  inline Real distance(const Vector<Real, Dim>& v1,
                       const Vector<Real, Dim>& v2) {
    Vector<Real, Dim> tempVector = v1 - v2;
    return tempVector.length();
  }

  /**
   * Returns the squared distance between  two Vectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the squared distance between v1 and v2
   */
  template<typename Real, unsigned int Dim>
  inline Real distanceSquare(const Vector<Real, Dim>& v1,
                             const Vector<Real, Dim>& v2) {
    Vector<Real, Dim> tempVector = v1 - v2;
    return tempVector.lengthSquare();
  }

  /**
   * Write to an output stream
   */
  template<typename Real, unsigned int Dim>
  std::ostream& operator<<(std::ostream& o, const Vector<Real, Dim>& v) {
    o << "[ ";

    for (unsigned int i = 0; i < Dim; i++) {
      o << v[i];

      if (i != Dim - 1) { o << ", "; }
    }

    o << " ]";
    return o;
  }

}  // namespace OpenMD

#endif
