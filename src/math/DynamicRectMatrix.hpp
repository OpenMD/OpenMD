/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file DynamicRectMatrix.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */

#ifndef MATH_DYNAMICRECTMATRIX_HPP
#define MATH_DYNAMICRECTMATRIX_HPP
#include <math.h>
#include <cmath>
#include "math/DynamicVector.hpp"

namespace OpenMD {

  /**
   * @class DynamicRectMatrix DynamicRectMatrix.hpp "math/DynamicRectMatrix.hpp"
   * @brief rectangular matrix class
   */
  template<typename Real> 
  class DynamicRectMatrix {
  public:
    typedef Real ElemType;
    typedef Real* ElemPoinerType;
    typedef DynamicRectMatrix<Real> SelfType;
            
    /** default constructor */
    DynamicRectMatrix(){
      nrow_ = 0;
      ncol_ = 0;
      data_ = NULL;
    }

    DynamicRectMatrix(int nrow, int ncol) {
      allocate(nrow, ncol);

      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)
	  this->data_[i][j] = 0.0;
    }

    /** Constructs and initializes every element of this matrix to a scalar */ 
    DynamicRectMatrix(int nrow, int ncol, Real s) {
      allocate(nrow, ncol);
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)
	  this->data_[i][j] = s;
    }

    DynamicRectMatrix(int nrow, int ncol, Real* array) {
      allocate(nrow, ncol);
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)
	  this->data_[i][j] = array[i * nrow_ + j];
    }

    /** copy constructor */
    DynamicRectMatrix(const SelfType& m) {
      allocate(m.getNRow(), m.getNCol());
      
      for (unsigned int i = 0; i < nrow_; i++)
        for (unsigned int j = 0; j < ncol_; j++)
	    this->data_[i][j] = m.data_[i][j];
    }
            
    /** destructor*/
    ~DynamicRectMatrix() { deallocate();}

    /** copy assignment operator */
    DynamicRectMatrix<Real> operator =(const DynamicRectMatrix<Real> &m) {
      if (this == &m)
	  return *this;
      if (nrow_ != m.getNRow() || ncol_ != m.getNCol()) {
        deallocate();
        allocate(m.getNRow(), m.getNCol());
      }
      
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)
	  this->data_[i][j] = m.data_[i][j];
      return *this;
    }
            
    /**
     * Return the reference of a single element of this matrix.
     * @return the reference of a single element of this matrix 
     * @param i row index
     * @param j Column index
     */
    Real& operator()(unsigned int i, unsigned int j) {
      return this->data_[i][j];
    }

    /**
     * Return the value of a single element of this matrix.
     * @return the value of a single element of this matrix 
     * @param i row index
     * @param j Column index
     */        
    Real operator()(unsigned int i, unsigned int j) const  {
                
      return this->data_[i][j];  
    }

    /** 
     * Copy the internal data to an array
     * @param array the pointer of destination array
     */
    void getArray(Real* array) {
      for (unsigned int i = 0; i < nrow_; i++) {
	for (unsigned int j = 0; j < ncol_; j++) {
	  array[i * nrow_ + j] = this->data_[i][j];
	}
      }
    }

    /**
     * Returns a row of  this matrix as a vector.
     * @return a row of  this matrix as a vector 
     * @param row the row index
     */                
    DynamicVector<Real> getRow(unsigned int row) {
      DynamicVector<Real> v;

      for (unsigned int i = 0; i < ncol_; i++)
	v[i] = this->data_[row][i];

      return v;
    }

    /**
     * Sets a row of  this matrix
     * @param row the row index
     * @param v the vector to be set
     */                
    void setRow(unsigned int row, const DynamicVector<Real>& v) {
      assert(v.size() == nrow_);
      for (unsigned int i = 0; i < ncol_; i++)
	this->data_[row][i] = v[i];
    }

    /**
     * Returns a column of  this matrix as a vector.
     * @return a column of  this matrix as a vector 
     * @param col the column index
     */                
    DynamicVector<Real> getColumn(unsigned int col) {
      DynamicVector<Real> v(ncol_);

      for (unsigned int j = 0; j < nrow_; j++)
	v[j] = this->data_[j][col];

      return v;
    }

    /**
     * Sets a column of  this matrix
     * @param col the column index
     * @param v the vector to be set
     */                
    void setColumn(unsigned int col, const DynamicVector<Real>& v){

      for (unsigned int j = 0; j < nrow_; j++)
	this->data_[j][col] = v[j];
    }         

    /**
     * swap two rows of this matrix
     * @param i the first row
     * @param j the second row
     */
    void swapRow(unsigned int i, unsigned int j){
      assert(i < nrow_ && j < nrow_);

      for (unsigned int k = 0; k < ncol_; k++)
	std::swap(this->data_[i][k], this->data_[j][k]);
    }

    /**
     * swap two Columns of this matrix
     * @param i the first Column
     * @param j the second Column
     */
    void swapColumn(unsigned int i, unsigned int j){
      assert(i < ncol_ && j < ncol_);
                    
      for (unsigned int k = 0; k < nrow_; k++)
	std::swap(this->data_[k][i], this->data_[k][j]);
    }

    /**
     * Tests if this matrix is identical to matrix m
     * @return true if this matrix is equal to the matrix m, return false otherwise
     * @param m matrix to be compared
     *
     * @todo replace operator == by template function equal
     */
    bool operator ==(const DynamicRectMatrix<Real> &m) {
      assert(nrow_ == m.getNRow() && ncol_ == m.getNCol());
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)
	  if (!equal(this->data_[i][j], m.data_[i][j]))
	    return false;

      return true;
    }

    /**
     * Tests if this matrix is not equal to matrix m
     * @return true if this matrix is not equal to the matrix m, return false otherwise
     * @param m matrix to be compared
     */
    bool operator !=(const DynamicRectMatrix<Real> &m) {
      return !(*this == m);
    }

    /** Negates the value of this matrix in place. */           
    inline void negate() {
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)
	  this->data_[i][j] = -this->data_[i][j];
    }
            
    /**
     * Sets the value of this matrix to the negation of matrix m.
     * @param m the source matrix
     */
    inline void negate(const DynamicRectMatrix<Real> &m) {
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)
	  this->data_[i][j] = -m.data_[i][j];        
    }
            
    /**
     * Sets the value of this matrix to the sum of itself and m (*this += m).
     * @param m the other matrix
     */
    inline void add( const DynamicRectMatrix<Real> &m ) {
      assert(nrow_ == m.getNRow() && ncol_ == m.getNCol());
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)        
	  this->data_[i][j] += m.data_[i][j];
    }
            
    /**
     * Sets the value of this matrix to the sum of m1 and m2 (*this = m1 + m2).
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    inline void add( const DynamicRectMatrix<Real> &m1, const DynamicRectMatrix<Real> &m2 ) {
      assert(m1.getNRow() == m2.getNRow() && m1.getNCol() == m2.getNCol());
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)        
	  this->data_[i][j] = m1.data_[i][j] + m2.data_[i][j];
    }
            
    /**
     * Sets the value of this matrix to the difference  of itself and m (*this -= m).
     * @param m the other matrix
     */
    inline void sub( const DynamicRectMatrix<Real> &m ) {
      assert(nrow_ == m.getNRow() && ncol_ == m.getNCol());
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)        
	  this->data_[i][j] -= m.data_[i][j];
    }
            
    /**
     * Sets the value of this matrix to the difference of matrix m1 and m2 (*this = m1 - m2).
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    inline void sub( const DynamicRectMatrix<Real> &m1, const DynamicRectMatrix<Real> &m2){
      assert(m1.getNRow() == m2.getNRow() && m1.getNCol() == m2.getNCol());
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)        
	  this->data_[i][j] = m1.data_[i][j] - m2.data_[i][j];
    }

    /**
     * Sets the value of this matrix to the scalar multiplication of itself (*this *= s).
     * @param s the scalar value
     */
    inline void mul( Real s ) {
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)  
	  this->data_[i][j] *= s;
    }

    /**
     * Sets the value of this matrix to the scalar multiplication of matrix m  (*this = s * m).
     * @param s the scalar value
     * @param m the matrix
     */
    inline void mul( Real s, const DynamicRectMatrix<Real> &m ) {
      assert(nrow_ == m.getNRow() && ncol_ == m.getNCol());    
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)  
	  this->data_[i][j] = s * m.data_[i][j];
    }

    /**
     * Sets the value of this matrix to the scalar division of itself  (*this /= s ).
     * @param s the scalar value
     */             
    inline void div( Real s) {
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)  
	  this->data_[i][j] /= s;
    }

    /**
     * Sets the value of this matrix to the scalar division of matrix m  (*this = m /s).
     * @param s the scalar value
     * @param m the matrix
     */
    inline void div( Real s, const DynamicRectMatrix<Real> &m ) {
      assert(nrow_ == m.getNRow() && ncol_ == m.getNCol());
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)  
	  this->data_[i][j] = m.data_[i][j] / s;
    }

    /**
     *  Multiples a scalar into every element of this matrix.
     * @param s the scalar value
     */
    DynamicRectMatrix<Real> operator *=(const Real s) {
      this->mul(s);
      return *this;
    }

    /**
     *  Divides every element of this matrix by a scalar.
     * @param s the scalar value
     */
    DynamicRectMatrix<Real> operator /=(const Real s) {
      this->div(s);
      return *this;
    }

    /**
     * Sets the value of this matrix to the sum of the other matrix and itself (*this += m).
     * @param m the other matrix
     */
    DynamicRectMatrix<Real> operator += (const DynamicRectMatrix<Real> m) {
      assert(nrow_ == m.getNRow() && ncol_ == m.getNCol());
      add(m);
      return *this;
    }

    /**
     * Sets the value of this matrix to the differerence of itself and the other matrix (*this -= m) 
     * @param m the other matrix
     */
    DynamicRectMatrix<Real> operator -= (const DynamicRectMatrix<Real> m){
      assert(nrow_ == m.getNRow() && ncol_ == m.getNCol());
      sub(m);
      return *this;
    }

    /** Return the transpose of this matrix */
    DynamicRectMatrix<Real> transpose() const{
      DynamicRectMatrix<Real> result(ncol_,nrow_);
                
      for (unsigned int i = 0; i < nrow_; i++)
	for (unsigned int j = 0; j < ncol_; j++)              
	  result(j, i) = this->data_[i][j];

      return result;
    }

    unsigned int getNRow() const {return nrow_;}
    unsigned int getNCol() const {return ncol_;}

    template<class MatrixType>
    void setSubMatrix(unsigned int beginRow, unsigned int beginCol, const MatrixType& m) {
        assert(beginRow + m.getNRow() -1 <= nrow_);
        assert(beginCol + m.getNCol() -1 <= ncol_);

        for (unsigned int i = 0; i < m.getNRow(); ++i)
            for (unsigned int j = 0; j < m.getNCol(); ++j)
                this->data_[beginRow+i][beginCol+j] = m(i, j);
    }

    template<class MatrixType>
    void getSubMatrix(unsigned int beginRow, unsigned int beginCol, MatrixType& m) {
        assert(beginRow + m.getNRow() -1 <= nrow_);
        assert(beginCol + m.getNCol() - 1 <= ncol_);

        for (unsigned int i = 0; i < m.getNRow(); ++i)
            for (unsigned int j = 0; j < m.getNCol(); ++j)
                m(i, j) = this->data_[beginRow+i][beginCol+j];
    }
    
  protected:
    Real** data_;
    unsigned int nrow_;
    unsigned int ncol_;
  private:
    void allocate( int nrow,  int ncol ) {
      nrow_ = (unsigned int) nrow;
      ncol_ = (unsigned int) ncol;
      data_ = new Real*[nrow_];
      for (unsigned int i = 0; i < nrow_; ++i)
        data_[i] = new Real[ncol_];
    }
    
    void deallocate() {
      for (unsigned int i = 0; i < nrow_; ++i)
        delete data_[i];
      delete []data_;
      
      nrow_ = 0;
      ncol_ = 0;
      data_ = NULL;
    }
    
  };

  /** Negate the value of every element of this matrix. */
  template<typename Real> 
  inline DynamicRectMatrix<Real> operator -(const DynamicRectMatrix<Real> &m) {
    DynamicRectMatrix<Real> result(m);

    result.negate();

    return result;
  }
    
  /**
   * Return the sum of two matrixes  (m1 + m2). 
   * @return the sum of two matrixes
   * @param m1 the first matrix
   * @param m2 the second matrix
   */ 
  template<typename Real> 
  inline DynamicRectMatrix<Real> operator + (const DynamicRectMatrix<Real> &m1, const DynamicRectMatrix<Real> &m2) {
    
    DynamicRectMatrix<Real> result(m1.getNRow(), m1.getNCol());

    result.add(m1, m2);

    return result;
  }
    
  /**
   * Return the difference of two matrixes  (m1 - m2). 
   * @return the sum of two matrixes
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  template<typename Real> 
  inline DynamicRectMatrix<Real> operator - (const DynamicRectMatrix<Real> &m1, const DynamicRectMatrix<Real> &m2) {
    DynamicRectMatrix<Real> result(m1.getNRow(), m1.getNCol());

    result.sub(m1, m2);

    return result;
  }

  /**
   * Return the multiplication of scalra and  matrix  (m * s). 
   * @return the multiplication of a scalra and  a matrix 
   * @param m the matrix
   * @param s the scalar
   */
  template<typename Real> 
  inline DynamicRectMatrix<Real> operator *(const DynamicRectMatrix<Real> &m, Real s) {
    DynamicRectMatrix<Real> result(m.getNRow(), m.getNCol());

    result.mul(s, m);

    return result;
  }

  /**
   * Return the multiplication of a scalra and  a matrix  (s * m). 
   * @return the multiplication of a scalra and  a matrix 
   * @param s the scalar
   * @param m the matrix
   */
  template<typename Real> 
  inline DynamicRectMatrix<Real> operator *(Real s, const DynamicRectMatrix<Real> &m) {
    DynamicRectMatrix<Real> result(m.getNRow(), m.getNCol());

    result.mul(s, m);

    return result;
  }
    
  /**
   * Return the multiplication of two matrixes  (m1 * m2). 
   * @return the multiplication of two matrixes
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  template<typename Real> 
  inline DynamicRectMatrix<Real> operator *(const DynamicRectMatrix<Real>& m1, const DynamicRectMatrix<Real>& m2) {
    assert(m1.getNCol() == m2.getNRow());
    unsigned int sameDim = m1.getNCol();
    int nrow = m1.getNRow();
    int ncol = m2.getNCol();
    DynamicRectMatrix<Real> result(nrow, ncol );
    for (unsigned int i = 0; i < nrow; i++)
      for (unsigned int j = 0; j < ncol; j++)
	for (unsigned int k = 0; k < sameDim; k++)
	  result(i, j)  += m1(i, k) * m2(k, j);                

    return result;
  }
    
  /**
   * Return the multiplication of  a matrix and a vector  (m * v). 
   * @return the multiplication of a matrix and a vector
   * @param m the matrix
   * @param v the vector
   */
  template<typename Real>
  inline DynamicVector<Real> operator *(const DynamicRectMatrix<Real> &m, const DynamicVector<Real> &v) {
    int nrow = m.getNRow();
    int ncol = m.getNCol();
    assert(ncol == v.size());
    DynamicVector<Real> result(nrow);
    
    for (unsigned int i = 0; i < nrow ; i++)
      for (unsigned int j = 0; j < ncol ; j++)            
	result[i] += m(i, j) * v[j];
            
    return result;                                                                 
  }

  /**
   * Return the scalar division of matrix   (m / s). 
   * @return the scalar division of matrix  
   * @param m the matrix
   * @param s the scalar
   */
  template<typename Real> 
  inline DynamicRectMatrix<Real> operator /(const DynamicRectMatrix<Real> &m, Real s) {
    DynamicRectMatrix<Real> result(m.getNRow(), m.getNCol());

    result.div(s, m);

    return result;
  }    

  /**
   * Write to an output stream
   */
  template<typename Real>
  std::ostream &operator<< ( std::ostream& o, const DynamicRectMatrix<Real> &m) {
    for (unsigned int i = 0; i < m.getNRow() ; i++) {
      o << "(";
      for (unsigned int j = 0; j < m.getNCol() ; j++) {
	o << m(i, j);
	if (j != m.getNCol() -1)
	  o << "\t";
      }
      o << ")" << std::endl;
    }
    return o;        
  }    
}
#endif //MATH_RECTMATRIX_HPP

