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
 * @file DynamicRectMatrix.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 2.0
 */

#ifndef MATH_DYNAMICRECTMATRIX_HPP
#define MATH_DYNAMICRECTMATRIX_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "math/DynamicVector.hpp"
#include "math/RectMatrix.hpp"

namespace OpenMD {

  /**
   * @class DynamicRectMatrix DynamicRectMatrix.hpp "math/DynamicRectMatrix.hpp"
   * @brief Rectangular matrix class with contiguous flat storage.
   *
   * Elements are stored in row-major order in a single contiguous
   * std::vector<Real>.  Element (i,j) is at data_[i*ncol_ + j].
   * The contiguous layout allows direct interop with BLAS/LAPACK
   * (via getArrayPointer()) without an intermediate copy.
   */
  template<typename Real>
  class DynamicRectMatrix {
  public:
    using ElemType       = Real;
    using ElemPoinerType = Real*;
    using SelfType       = DynamicRectMatrix<Real>;

    /** default constructor */
    DynamicRectMatrix() : nrow_(0), ncol_(0) {}

    DynamicRectMatrix(unsigned int nrow, unsigned int ncol) :
      nrow_(nrow), ncol_(ncol), data_(nrow * ncol, Real(0)) {}

    /** Constructs and initializes every element of this matrix to a scalar */
    DynamicRectMatrix(unsigned int nrow, unsigned int ncol, Real s) :
      nrow_(nrow), ncol_(ncol), data_(nrow * ncol, s) {}

    DynamicRectMatrix(unsigned int nrow, unsigned int ncol, Real* array) :
      nrow_(nrow), ncol_(ncol), data_(array, array + nrow * ncol) {}

    /** copy constructor */
    DynamicRectMatrix(const SelfType& m) = default;

    /** destructor */
    ~DynamicRectMatrix() = default;

    /** copy assignment operator */
    DynamicRectMatrix<Real>& operator=(const DynamicRectMatrix<Real>& m) = default;

    /** move constructor */
    DynamicRectMatrix(SelfType&& m) noexcept = default;

    /** move assignment operator */
    DynamicRectMatrix<Real>& operator=(DynamicRectMatrix<Real>&& m) noexcept = default;

    /**
     * Returns the reference of a single element of this matrix.
     * @param i row index
     * @param j column index
     */
    Real& operator()(unsigned int i, unsigned int j) {
      return data_[i * ncol_ + j];
    }

    /**
     * Returns the value of a single element of this matrix.
     * @param i row index
     * @param j column index
     */
    Real operator()(unsigned int i, unsigned int j) const {
      return data_[i * ncol_ + j];
    }

    /**
     * Copies the internal data to an array (row-major order).
     * @param array the pointer of destination array
     */
    void getArray(Real* array) const {
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        array[i] = data_[i];
    }

    /** Returns a pointer to the contiguous internal storage (row-major). */
    Real* getArrayPointer() { return data_.data(); }
    const Real* getArrayPointer() const { return data_.data(); }

    /**
     * Returns a row of this matrix as a vector.
     * @param row the row index
     */
    DynamicVector<Real> getRow(unsigned int row) {
      DynamicVector<Real> v(ncol_);
      for (unsigned int i = 0; i < ncol_; i++)
        v[i] = data_[row * ncol_ + i];
      return v;
    }

    /**
     * Sets a row of this matrix
     * @param row the row index
     * @param v the vector to be set
     */
    void setRow(unsigned int row, const DynamicVector<Real>& v) {
      assert(v.size() == ncol_);
      for (unsigned int i = 0; i < ncol_; i++)
        data_[row * ncol_ + i] = v[i];
    }

    /**
     * Returns a column of this matrix as a vector.
     * @param col the column index
     */
    DynamicVector<Real> getColumn(unsigned int col) {
      DynamicVector<Real> v(nrow_);
      for (unsigned int j = 0; j < nrow_; j++)
        v[j] = data_[j * ncol_ + col];
      return v;
    }

    /**
     * Sets a column of this matrix
     * @param col the column index
     * @param v the vector to be set
     */
    void setColumn(unsigned int col, const DynamicVector<Real>& v) {
      for (unsigned int j = 0; j < nrow_; j++)
        data_[j * ncol_ + col] = v[j];
    }

    /**
     * swap two rows of this matrix
     * @param i the first row
     * @param j the second row
     */
    void swapRow(unsigned int i, unsigned int j) {
      assert(i < nrow_ && j < nrow_);
      for (unsigned int k = 0; k < ncol_; k++)
        std::swap(data_[i * ncol_ + k], data_[j * ncol_ + k]);
    }

    /**
     * swap two columns of this matrix
     * @param i the first column
     * @param j the second column
     */
    void swapColumn(unsigned int i, unsigned int j) {
      assert(i < ncol_ && j < ncol_);
      for (unsigned int k = 0; k < nrow_; k++)
        std::swap(data_[k * ncol_ + i], data_[k * ncol_ + j]);
    }

    /**
     * Tests if this matrix is identical to matrix m
     * @param m matrix to be compared
     */
    bool operator==(const DynamicRectMatrix<Real>& m) const {
      assert(nrow_ == m.nrow_ && ncol_ == m.ncol_);
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        if (!equal(data_[i], m.data_[i])) return false;
      return true;
    }

    /**
     * Tests if this matrix is not equal to matrix m
     * @param m matrix to be compared
     */
    bool operator!=(const DynamicRectMatrix<Real>& m) const {
      return !(*this == m);
    }

    /** Negates the value of this matrix in place. */
    inline void negate() {
      for (auto& v : data_) v = -v;
    }

    /**
     * Sets the value of this matrix to the negation of matrix m.
     * @param m the source matrix
     */
    inline void negate(const DynamicRectMatrix<Real>& m) {
      assert(nrow_ == m.nrow_ && ncol_ == m.ncol_);
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        data_[i] = -m.data_[i];
    }

    /**
     * Sets the value of this matrix to the sum of itself and m (*this += m).
     * @param m the other matrix
     */
    inline void add(const DynamicRectMatrix<Real>& m) {
      assert(nrow_ == m.nrow_ && ncol_ == m.ncol_);
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        data_[i] += m.data_[i];
    }

    /**
     * Sets the value of this matrix to the sum of m1 and m2 (*this = m1 + m2).
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    inline void add(const DynamicRectMatrix<Real>& m1,
                    const DynamicRectMatrix<Real>& m2) {
      assert(m1.nrow_ == m2.nrow_ && m1.ncol_ == m2.ncol_);
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        data_[i] = m1.data_[i] + m2.data_[i];
    }

    /**
     * Sets the value of this matrix to the difference of itself and m
     * (*this -= m).
     * @param m the other matrix
     */
    inline void sub(const DynamicRectMatrix<Real>& m) {
      assert(nrow_ == m.nrow_ && ncol_ == m.ncol_);
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        data_[i] -= m.data_[i];
    }

    /**
     * Sets the value of this matrix to the difference of matrix m1 and m2
     * (*this = m1 - m2).
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    inline void sub(const DynamicRectMatrix<Real>& m1,
                    const DynamicRectMatrix<Real>& m2) {
      assert(m1.nrow_ == m2.nrow_ && m1.ncol_ == m2.ncol_);
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        data_[i] = m1.data_[i] - m2.data_[i];
    }

    /**
     * Sets the value of this matrix to the scalar multiplication of itself
     * (*this *= s).
     * @param s the scalar value
     */
    inline void mul(Real s) {
      for (auto& v : data_) v *= s;
    }

    /**
     * Sets the value of this matrix to the scalar multiplication of matrix m
     * (*this = s * m).
     * @param s the scalar value
     * @param m the matrix
     */
    inline void mul(Real s, const DynamicRectMatrix<Real>& m) {
      assert(nrow_ == m.nrow_ && ncol_ == m.ncol_);
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        data_[i] = s * m.data_[i];
    }

    /**
     * Sets the value of this matrix to the scalar division of itself
     * (*this /= s).
     * @param s the scalar value
     */
    inline void div(Real s) {
      for (auto& v : data_) v /= s;
    }

    /**
     * Sets the value of this matrix to the scalar division of matrix m
     * (*this = m / s).
     * @param s the scalar value
     * @param m the matrix
     */
    inline void div(Real s, const DynamicRectMatrix<Real>& m) {
      assert(nrow_ == m.nrow_ && ncol_ == m.ncol_);
      for (unsigned int i = 0; i < nrow_ * ncol_; i++)
        data_[i] = m.data_[i] / s;
    }

    /**
     * Multiples a scalar onto every element of this matrix.
     * @param s the scalar value
     */
    DynamicRectMatrix<Real>& operator*=(const Real s) {
      this->mul(s);
      return *this;
    }

    /**
     * Divides every element of this matrix by a scalar.
     * @param s the scalar value
     */
    DynamicRectMatrix<Real>& operator/=(const Real s) {
      this->div(s);
      return *this;
    }

    /**
     * Sets the value of this matrix to the sum of the other matrix and itself
     * (*this += m).
     * @param m the other matrix
     */
    DynamicRectMatrix<Real>& operator+=(const DynamicRectMatrix<Real>& m) {
      add(m);
      return *this;
    }

    /**
     * Sets the value of this matrix to the difference of itself and the other
     * matrix (*this -= m).
     * @param m the other matrix
     */
    DynamicRectMatrix<Real>& operator-=(const DynamicRectMatrix<Real>& m) {
      sub(m);
      return *this;
    }

    /** Return the transpose of this matrix */
    DynamicRectMatrix<Real> transpose() const {
      DynamicRectMatrix<Real> result(ncol_, nrow_);
      for (unsigned int i = 0; i < nrow_; i++)
        for (unsigned int j = 0; j < ncol_; j++)
          result(j, i) = data_[i * ncol_ + j];
      return result;
    }

    unsigned int getNRow() const { return nrow_; }
    unsigned int getNCol() const { return ncol_; }

    template<class MatrixType>
    void setSubMatrix(unsigned int beginRow, unsigned int beginCol,
                      const MatrixType& m) {
      assert(beginRow + m.getNRow() - 1 <= nrow_);
      assert(beginCol + m.getNCol() - 1 <= ncol_);
      for (unsigned int i = 0; i < m.getNRow(); ++i)
        for (unsigned int j = 0; j < m.getNCol(); ++j)
          data_[(beginRow + i) * ncol_ + (beginCol + j)] = m(i, j);
    }

    template<class MatrixType>
    void getSubMatrix(unsigned int beginRow, unsigned int beginCol,
                      MatrixType& m) {
      assert(beginRow + m.getNRow() - 1 <= nrow_);
      assert(beginCol + m.getNCol() - 1 <= ncol_);
      for (unsigned int i = 0; i < m.getNRow(); ++i)
        for (unsigned int j = 0; j < m.getNCol(); ++j)
          m(i, j) = data_[(beginRow + i) * ncol_ + (beginCol + j)];
    }

  protected:
    unsigned int nrow_;
    unsigned int ncol_;
    std::vector<Real> data_;
  };

  /** Negate the value of every element of this matrix. */
  template<typename Real>
  inline DynamicRectMatrix<Real> operator-(const DynamicRectMatrix<Real>& m) {
    DynamicRectMatrix<Real> result(m);
    result.negate();
    return result;
  }

  /**
   * Return the sum of two matrices  (m1 + m2).
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  template<typename Real>
  inline DynamicRectMatrix<Real> operator+(const DynamicRectMatrix<Real>& m1,
                                           const DynamicRectMatrix<Real>& m2) {
    DynamicRectMatrix<Real> result(m1.getNRow(), m1.getNCol());
    result.add(m1, m2);
    return result;
  }

  /**
   * Return the difference of two matrices  (m1 - m2).
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  template<typename Real>
  inline DynamicRectMatrix<Real> operator-(const DynamicRectMatrix<Real>& m1,
                                           const DynamicRectMatrix<Real>& m2) {
    DynamicRectMatrix<Real> result(m1.getNRow(), m1.getNCol());
    result.sub(m1, m2);
    return result;
  }

  /**
   * Return the multiplication of scalar and matrix  (m * s).
   * @param m the matrix
   * @param s the scalar
   */
  template<typename Real>
  inline DynamicRectMatrix<Real> operator*(const DynamicRectMatrix<Real>& m,
                                           Real s) {
    DynamicRectMatrix<Real> result(m.getNRow(), m.getNCol());
    result.mul(s, m);
    return result;
  }

  /**
   * Return the multiplication of a scalar and a matrix  (s * m).
   * @param s the scalar
   * @param m the matrix
   */
  template<typename Real>
  inline DynamicRectMatrix<Real> operator*(Real s,
                                           const DynamicRectMatrix<Real>& m) {
    DynamicRectMatrix<Real> result(m.getNRow(), m.getNCol());
    result.mul(s, m);
    return result;
  }

  /**
   * Return the multiplication of two matrices  (m1 * m2).
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  template<typename Real>
  inline DynamicRectMatrix<Real> operator*(const DynamicRectMatrix<Real>& m1,
                                           const DynamicRectMatrix<Real>& m2) {
    assert(m1.getNCol() == m2.getNRow());
    unsigned int sameDim = m1.getNCol();
    unsigned int nrow    = m1.getNRow();
    unsigned int ncol    = m2.getNCol();
    DynamicRectMatrix<Real> result(nrow, ncol);
    for (unsigned int i = 0; i < nrow; i++)
      for (unsigned int j = 0; j < ncol; j++)
        for (unsigned int k = 0; k < sameDim; k++)
          result(i, j) += m1(i, k) * m2(k, j);
    return result;
  }

  /**
   * Return the multiplication of a matrix and a vector  (m * v).
   * @param m the matrix
   * @param v the vector
   */
  template<typename Real>
  inline DynamicVector<Real> operator*(const DynamicRectMatrix<Real>& m,
                                       const DynamicVector<Real>& v) {
    unsigned int nrow = m.getNRow();
    unsigned int ncol = m.getNCol();
    assert(ncol == v.size());
    DynamicVector<Real> result(nrow);
    for (unsigned int i = 0; i < nrow; i++)
      for (unsigned int j = 0; j < ncol; j++)
        result[i] += m(i, j) * v[j];
    return result;
  }

  /**
   * Return the scalar division of matrix   (m / s).
   * @param m the matrix
   * @param s the scalar
   */
  template<typename Real>
  inline DynamicRectMatrix<Real> operator/(const DynamicRectMatrix<Real>& m,
                                           Real s) {
    DynamicRectMatrix<Real> result(m.getNRow(), m.getNCol());
    result.div(s, m);
    return result;
  }

  /**
   * Write to an output stream
   */
  template<typename Real>
  std::ostream& operator<<(std::ostream& o, const DynamicRectMatrix<Real>& m) {
    for (unsigned int i = 0; i < m.getNRow(); i++) {
      o << "(";
      for (unsigned int j = 0; j < m.getNCol(); j++) {
        o << m(i, j);
        if (j != m.getNCol() - 1) o << "\t";
      }
      o << ")" << std::endl;
    }
    return o;
  }
}  // namespace OpenMD

#endif  // MATH_DYNAMICRECTMATRIX_HPP
