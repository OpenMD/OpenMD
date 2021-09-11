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

/**
 * @file DynamicVector.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */

#ifndef MATH_DYNAMICVECTOR_HPP
#define MATH_DYNAMICVECTOR_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <vector>

#include "math/Vector.hpp"

namespace OpenMD {

  /**
   * @class DynamicVector DynamicVector.hpp "math/DynamicVector.hpp"
   * @brief Dynamically-sized vector class
   */
  template<typename Real, typename Alloc = std::allocator<Real>>
  class DynamicVector {
  public:
    using value_type             = Real;
    using allocator_type         = Alloc;
    using VectorType             = std::vector<Real, Alloc>;
    using size_type              = typename VectorType::size_type;
    using difference_type        = typename VectorType::difference_type;
    using reference              = typename VectorType::reference;
    using const_reference        = typename VectorType::const_reference;
    using pointer                = typename VectorType::pointer;
    using const_pointer          = typename VectorType::const_pointer;
    using iterator               = typename VectorType::iterator;
    using const_iterator         = typename VectorType::const_iterator;
    using reverse_iterator       = typename VectorType::reverse_iterator;
    using const_reverse_iterator = typename VectorType::const_reverse_iterator;

    /**
     * @brief  Default constructor creates no elements.
     * @param  alloc  The allocator_type to use
     */
    explicit DynamicVector(const allocator_type& alloc = allocator_type()) :
        data_(alloc) {}

    /**
     * @brief  Create a %DynamicVector with copies of an exemplar element.
     * @param  n  The number of elements to initially create.
     * @param  value  An element to copy.
     * @param  alloc  The allocator_type to use
     *
     * This constructor fills the %DynamicVector with @a n copies of @a value.
     */
    DynamicVector(size_type n, const value_type& value,
                  const allocator_type& alloc = allocator_type()) :
        data_(n, value, alloc) {}

    /**
     * @brief  Create a %DynamicVector with default elements.
     * @param  n  The number of elements to initially create.
     * @param  alloc  The allocator_type to use
     *
     *  This constructor fills the %DynamicVector with @a n copies of a
     *  default-constructed element.
     */
    explicit DynamicVector(size_type n,
                           const allocator_type& alloc = allocator_type()) :
        data_(n, alloc) {}

    /**
     * @brief  Create a %DynamicVector using an iterator range
     * @param  first  The beginning of the range to copy the elements from
     * @param  last   The end of the range to copy the elements from
     * @param  alloc  The allocator_type to use
     */
    template<typename InputIterator>
    DynamicVector(InputIterator first, InputIterator last,
                  const allocator_type& alloc = allocator_type()) :
        data_(first, last, alloc) {}

    /**
     * @brief  Create a %DynamicVector with the contents of an initializer_list
     * @param  init   Initializer list to initialize the elements with
     * @param  alloc  The allocator_type to use
     */
    DynamicVector(std::initializer_list<value_type> init,
                  const allocator_type& alloc = allocator_type()) :
        data_(init, alloc) {}

    // Element access functions
    reference operator[](size_type i) { return data_[i]; }
    const_reference operator[](size_type i) const { return data_[i]; }

    reference operator()(size_type i) { return data_[i]; }
    const_reference operator()(size_type i) const { return data_[i]; }

    // Iterator functions
    iterator begin() noexcept { return data_.begin(); }
    const_iterator begin() const noexcept { return data_.begin(); }
    const_iterator cbegin() const noexcept { return data_.cbegin(); }

    iterator end() noexcept { return data_.end(); }
    const_iterator end() const noexcept { return data_.end(); }
    const_iterator cend() const noexcept { return data_.cend(); }

    // Capacity functions
    bool empty() const noexcept { return data_.empty(); }
    size_type size() const noexcept { return data_.size(); }

    // Modifier functions
    void resize(size_type n) { return data_.resize(n); }
    void resize(size_type n, const value_type& value) {
      data_.resize(n, value);
    }

    /**
     * Tests if this vetor is equal to other vector
     * @return true if equal, otherwise return false
     * @param v vector to be compared
     */
    bool operator==(const DynamicVector<Real>& v) {
      if (this->size() != v.size()) return false;

      return std::equal(
          this->begin(), this->end(), v.begin(),
          [](Real val1, Real val2) { return OpenMD::equal(val1, val2); });
    }

    /**
     * Tests if this vetor is not equal to other vector
     * @return true if equal, otherwise return false
     * @param v vector to be compared
     */
    bool operator!=(const DynamicVector<Real>& v) { return !(*this == v); }

    /** Negates the value of this vector in place. */
    void negate() {
      std::transform(this->begin(), this->end(), this->begin(),
                     [](Real val) { return -val; });
    }

    /**
     * Sets the value of this vector to the negation of vector v1.
     * @param v1 the source vector
     */
    void negate(const DynamicVector<Real>& v1) {
      std::transform(v1.begin(), v1.end(), this->begin(),
                     [](Real val) { return -val; });
    }

    /**
     * Sets the value of this vector to the sum of itself and v1 (*this += v1).
     * @param v1 the other vector
     */
    void add(const DynamicVector<Real>& v1) {
      std::transform(this->begin(), this->end(), v1.begin(), this->begin(),
                     [](Real val1, Real val2) { return val1 + val2; });
    }

    /**
     * Sets the value of this vector to the sum of v1 and v2 (*this = v1 + v2).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    void add(const DynamicVector<Real>& v1, const DynamicVector<Real>& v2) {
      std::transform(v1.begin(), v1.end(), v2.begin(), this->begin(),
                     [](Real val1, Real val2) { return val1 + val2; });
    }

    /**
     * Sets the value of this vector to the difference  of itself and v1 (*this
     * -= v1).
     * @param v1 the other vector
     */
    void sub(const DynamicVector<Real>& v1) {
      std::transform(this->begin(), this->end(), v1.begin(), this->begin(),
                     [](Real val1, Real val2) { return val1 - val2; });
    }

    /**
     * Sets the value of this vector to the difference of vector v1 and v2
     * (*this = v1 - v2).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    void sub(const DynamicVector<Real>& v1, const DynamicVector<Real>& v2) {
      std::transform(v1.begin(), v1.end(), v2.begin(), this->begin(),
                     [](Real val1, Real val2) { return val1 - val2; });
    }

    /**
     * Sets the value of this vector to the scalar multiplication of itself
     * (*this *= s).
     * @param s the scalar value
     */
    void mul(Real s) {
      std::transform(this->begin(), this->end(), this->begin(),
                     [s](Real val) { return val * s; });
    }

    /**
     * Sets the value of this vector to the scalar multiplication of vector v1
     * (*this = s * v1).
     * @param v1 the vector
     * @param s the scalar value
     */
    void mul(const DynamicVector<Real>& v1, Real s) {
      if (this->size() != v1.size()) this->resize(v1.size());

      std::transform(v1.begin(), v1.end(), this->begin(),
                     [s](Real val) { return val * s; });
    }

    /**
     * Sets the value of this vector to the scalar division of itself
     * (*this /= s).
     * @param s the scalar value
     */
    void div(Real s) {
      std::transform(this->begin(), this->end(), this->begin(),
                     [s](Real val) { return val / s; });
    }

    /**
     * Sets the value of this vector to the scalar division of vector v1
     * (*this = v1 / s).
     * @param v1 the source vector
     * @param s the scalar value
     */
    void div(const DynamicVector<Real>& v1, Real s) {
      if (this->size() != v1.size()) this->resize(v1.size());

      std::transform(v1.begin(), v1.end(), this->begin(),
                     [s](Real val) { return val / s; });
    }

    /** @see #add */
    DynamicVector<Real>& operator+=(const DynamicVector<Real>& v1) {
      add(v1);
      return *this;
    }

    /** @see #sub */
    DynamicVector<Real>& operator-=(const DynamicVector<Real>& v1) {
      sub(v1);
      return *this;
    }

    /** @see #mul */
    DynamicVector<Real>& operator*=(Real s) {
      mul(s);
      return *this;
    }

    /** @see #div */
    DynamicVector<Real>& operator/=(Real s) {
      div(s);
      return *this;
    }

    /** zero out the vector */
    void setZero() { std::fill(this->begin(), this->end(), 0); }

    /**
     * Returns the length of this vector.
     * @return the length of this vector
     */
    Real length() { return std::sqrt(lengthSquare()); }

    /**
     * Returns the squared length of this vector.
     * @return the squared length of this vector
     */
    Real lengthSquare() { return dot(*this, *this); }

    /** Normalizes this vector in place */
    void normalize() {
      Real len = length();

      // if (len < OpenMD::Constants::epsilon)
      //  throw();

      *this /= len;
    }

    /**
     * Tests if this vector is normalized
     * @return true if this vector is normalized, otherwise return false
     */
    bool isNormalized() { return OpenMD::equal(lengthSquare(), 1.0); }

    template<class VectorType>
    void getSubVector(size_type beginning, VectorType& v) {
      assert(beginning + v.size() - 1 <= this->size());

      for (size_type i {}; i < v.size(); ++i)
        v(i) = (*this)[beginning + i];
    }

  private:
    std::vector<Real, Alloc> data_;
  };

  /** unary minus*/
  template<typename Real>
  inline DynamicVector<Real> operator-(const DynamicVector<Real>& v1) {
    DynamicVector<Real> tmp(v1);
    tmp.negate();
    return tmp;
  }

  /**
   * Return the sum of two vectors  (v1 - v2).
   * @return the sum of two vectors
   * @param v1 the first vector
   * @param v2 the second vector
   */
  template<typename Real>
  inline DynamicVector<Real> operator+(const DynamicVector<Real>& v1,
                                       const DynamicVector<Real>& v2) {
    assert(v1.size() == v2.size());
    DynamicVector<Real> result(v1.size());
    result.add(v1, v2);
    return result;
  }

  /**
   * Return the difference of two vectors  (v1 - v2).
   * @return the difference of two vectors
   * @param v1 the first vector
   * @param v2 the second vector
   */
  template<typename Real>
  DynamicVector<Real> operator-(const DynamicVector<Real>& v1,
                                const DynamicVector<Real>& v2) {
    assert(v1.size() == v2.size());
    DynamicVector<Real> result(v1.size());
    result.sub(v1, v2);
    return result;
  }

  /**
   * Returns the vaule of scalar multiplication of this vector v1 (v1 * r).
   * @return  the vaule of scalar multiplication of this vector
   * @param v1 the source vector
   * @param s the scalar value
   */
  template<typename Real>
  DynamicVector<Real> operator*(const DynamicVector<Real>& v1, Real s) {
    DynamicVector<Real> result(v1.size());
    result.mul(v1, s);
    return result;
  }

  /**
   * Returns the vaule of scalar multiplication of this vector v1 (v1 * r).
   * @return  the vaule of scalar multiplication of this vector
   * @param s the scalar value
   * @param v1 the source vector
   */
  template<typename Real>
  DynamicVector<Real> operator*(Real s, const DynamicVector<Real>& v1) {
    DynamicVector<Real> result(v1.size());
    result.mul(v1, s);
    return result;
  }

  /**
   * Returns the  value of division of a vector by a scalar.
   * @return  the vaule of scalar division of this vector
   * @param v1 the source vector
   * @param s the scalar value
   */
  template<typename Real>
  DynamicVector<Real> operator/(const DynamicVector<Real>& v1, Real s) {
    DynamicVector<Real> result(v1.size());
    result.div(v1, s);
    return result;
  }

  /**
   * Returns the dot product of two DynamicVectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the dot product of v1 and v2
   */
  template<typename Real>
  inline Real dot(const DynamicVector<Real>& v1,
                  const DynamicVector<Real>& v2) {
    Real tmp;
    tmp = 0;
    assert(v1.size() == v2.size());
    for (typename DynamicVector<Real>::size_type i {}; i < v1.size(); i++)
      tmp += v1[i] * v2[i];

    return tmp;
  }

  /**
   * Returns the distance between  two DynamicVectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the distance between v1 and v2
   */
  template<typename Real>
  inline Real distance(const DynamicVector<Real>& v1,
                       const DynamicVector<Real>& v2) {
    DynamicVector<Real> tempDynamicVector = v1 - v2;
    return tempDynamicVector.length();
  }

  /**
   * Returns the squared distance between  two DynamicVectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the squared distance between v1 and v2
   */
  template<typename Real>
  inline Real distanceSquare(const DynamicVector<Real>& v1,
                             const DynamicVector<Real>& v2) {
    DynamicVector<Real> tempDynamicVector = v1 - v2;
    return tempDynamicVector.lengthSquare();
  }

  /**
   * Write to an output stream
   */
  template<typename Real>
  std::ostream& operator<<(std::ostream& strm, const DynamicVector<Real>& v) {
    strm << "[ ";

    std::for_each(v.begin(), v.end() - 1,
                  [&strm](auto elem) { strm << elem << ", "; });

    strm << *(v.end() - 1) << " ]";

    return strm;
  }
}  // namespace OpenMD

#endif
