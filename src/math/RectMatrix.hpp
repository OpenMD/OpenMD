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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file RectMatrix.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */

#ifndef MATH_RECTMATRIX_HPP
#define MATH_RECTMATRIX_HPP
#include <math.h>
#include <cmath>
#include "Vector.hpp"

namespace OpenMD {

  /**
   * @class RectMatrix RectMatrix.hpp "math/RectMatrix.hpp"
   * @brief rectangular matrix class
   */
  template<typename Real, unsigned int Row, unsigned int Col> 
  class RectMatrix {
  public:
    typedef Real ElemType;
    typedef Real* ElemPoinerType;
            
    /** default constructor */
    RectMatrix() {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)
	  this->data_[i][j] = 0.0;
    }

    /** Constructs and initializes every element of this matrix to a scalar */ 
    RectMatrix(Real s) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)
	  this->data_[i][j] = s;
    }

    RectMatrix(Real* array) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)
	  this->data_[i][j] = array[i * Row + j];
    }

    /** copy constructor */
    RectMatrix(const RectMatrix<Real, Row, Col>& m) {
      *this = m;
    }
            
    /** destructor*/
    ~RectMatrix() {}

    /** copy assignment operator */
    RectMatrix<Real, Row, Col>& operator =(const RectMatrix<Real, Row, Col>& m) {
      if (this == &m)
	return *this;
                
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)
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
      //assert( i < Row && j < Col);
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
      for (unsigned int i = 0; i < Row; i++) {
	for (unsigned int j = 0; j < Col; j++) {
	  array[i * Row + j] = this->data_[i][j];
	}
      }
    }


    /** Returns the pointer of internal array */
    Real* getArrayPointer() {
      return &this->data_[0][0];
    }

    /**
     * Returns a row of  this matrix as a vector.
     * @return a row of  this matrix as a vector 
     * @param row the row index
     */                
    Vector<Real, Row> getRow(unsigned int row) {
      Vector<Real, Row> v;

      for (unsigned int i = 0; i < Col; i++)
	v[i] = this->data_[row][i];

      return v;
    }

    /**
     * Sets a row of  this matrix
     * @param row the row index
     * @param v the vector to be set
     */                
    void setRow(unsigned int row, const Vector<Real, Row>& v) {

      for (unsigned int i = 0; i < Col; i++)
	this->data_[row][i] = v[i];
    }

    /**
     * Returns a column of  this matrix as a vector.
     * @return a column of  this matrix as a vector 
     * @param col the column index
     */                
    Vector<Real, Col> getColumn(unsigned int col) {
      Vector<Real, Col> v;

      for (unsigned int j = 0; j < Row; j++)
	v[j] = this->data_[j][col];

      return v;
    }

    /**
     * Sets a column of  this matrix
     * @param col the column index
     * @param v the vector to be set
     */                
    void setColumn(unsigned int col, const Vector<Real, Col>& v){

      for (unsigned int j = 0; j < Row; j++)
	this->data_[j][col] = v[j];
    }         

    /**
     * swap two rows of this matrix
     * @param i the first row
     * @param j the second row
     */
    void swapRow(unsigned int i, unsigned int j){
      assert(i < Row && j < Row);

      for (unsigned int k = 0; k < Col; k++)
	std::swap(this->data_[i][k], this->data_[j][k]);
    }

    /**
     * swap two Columns of this matrix
     * @param i the first Column
     * @param j the second Column
     */
    void swapColumn(unsigned int i, unsigned int j){
      assert(i < Col && j < Col);
                    
      for (unsigned int k = 0; k < Row; k++)
	std::swap(this->data_[k][i], this->data_[k][j]);
    }

    /**
     * Tests if this matrix is identical to matrix m
     * @return true if this matrix is equal to the matrix m, return false otherwise
     * @m matrix to be compared
     *
     * @todo replace operator == by template function equal
     */
    bool operator ==(const RectMatrix<Real, Row, Col>& m) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)
	  if (!equal(this->data_[i][j], m.data_[i][j]))
	    return false;

      return true;
    }

    /**
     * Tests if this matrix is not equal to matrix m
     * @return true if this matrix is not equal to the matrix m, return false otherwise
     * @m matrix to be compared
     */
    bool operator !=(const RectMatrix<Real, Row, Col>& m) {
      return !(*this == m);
    }

    /** Negates the value of this matrix in place. */           
    inline void negate() {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)
	  this->data_[i][j] = -this->data_[i][j];
    }
            
    /**
     * Sets the value of this matrix to the negation of matrix m.
     * @param m the source matrix
     */
    inline void negate(const RectMatrix<Real, Row, Col>& m) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)
	  this->data_[i][j] = -m.data_[i][j];        
    }
            
    /**
     * Sets the value of this matrix to the sum of itself and m (*this += m).
     * @param m the other matrix
     */
    inline void add( const RectMatrix<Real, Row, Col>& m ) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)        
	  this->data_[i][j] += m.data_[i][j];
    }
            
    /**
     * Sets the value of this matrix to the sum of m1 and m2 (*this = m1 + m2).
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    inline void add( const RectMatrix<Real, Row, Col>& m1, const RectMatrix<Real, Row, Col>& m2 ) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)        
	  this->data_[i][j] = m1.data_[i][j] + m2.data_[i][j];
    }
            
    /**
     * Sets the value of this matrix to the difference  of itself and m (*this -= m).
     * @param m the other matrix
     */
    inline void sub( const RectMatrix<Real, Row, Col>& m ) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)        
	  this->data_[i][j] -= m.data_[i][j];
    }
            
    /**
     * Sets the value of this matrix to the difference of matrix m1 and m2 (*this = m1 - m2).
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    inline void sub( const RectMatrix<Real, Row, Col>& m1, const RectMatrix<Real, Row, Col>& m2){
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)        
	  this->data_[i][j] = m1.data_[i][j] - m2.data_[i][j];
    }

    /**
     * Sets the value of this matrix to the scalar multiplication of itself (*this *= s).
     * @param s the scalar value
     */
    inline void mul( Real s ) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)  
	  this->data_[i][j] *= s;
    }

    /**
     * Sets the value of this matrix to the scalar multiplication of matrix m  (*this = s * m).
     * @param s the scalar value
     * @param m the matrix
     */
    inline void mul( Real s, const RectMatrix<Real, Row, Col>& m ) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)  
	  this->data_[i][j] = s * m.data_[i][j];
    }

    /**
     * Sets the value of this matrix to the scalar division of itself  (*this /= s ).
     * @param s the scalar value
     */             
    inline void div( Real s) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)  
	  this->data_[i][j] /= s;
    }

    /**
     * Sets the value of this matrix to the scalar division of matrix m  (*this = m /s).
     * @param s the scalar value
     * @param m the matrix
     */
    inline void div( Real s, const RectMatrix<Real, Row, Col>& m ) {
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)  
	  this->data_[i][j] = m.data_[i][j] / s;
    }

    /**
     *  Multiples a scalar into every element of this matrix.
     * @param s the scalar value
     */
    RectMatrix<Real, Row, Col>& operator *=(const Real s) {
      this->mul(s);
      return *this;
    }

    /**
     *  Divides every element of this matrix by a scalar.
     * @param s the scalar value
     */
    RectMatrix<Real, Row, Col>& operator /=(const Real s) {
      this->div(s);
      return *this;
    }

    /**
     * Sets the value of this matrix to the sum of the other matrix and itself (*this += m).
     * @param m the other matrix
     */
    RectMatrix<Real, Row, Col>& operator += (const RectMatrix<Real, Row, Col>& m) {
      add(m);
      return *this;
    }

    /**
     * Sets the value of this matrix to the differerence of itself and the other matrix (*this -= m) 
     * @param m the other matrix
     */
    RectMatrix<Real, Row, Col>& operator -= (const RectMatrix<Real, Row, Col>& m){
      sub(m);
      return *this;
    }

    /** Return the transpose of this matrix */
    RectMatrix<Real,  Col, Row> transpose() const{
      RectMatrix<Real,  Col, Row> result;
                
      for (unsigned int i = 0; i < Row; i++)
	for (unsigned int j = 0; j < Col; j++)              
	  result(j, i) = this->data_[i][j];

      return result;
    }

    template<class MatrixType>
    void setSubMatrix(unsigned int beginRow, unsigned int beginCol, const MatrixType& m) {
        assert(beginRow + m.getNRow() -1 <= getNRow());
        assert(beginCol + m.getNCol() -1 <= getNCol());

        for (unsigned int i = 0; i < m.getNRow(); ++i)
            for (unsigned int j = 0; j < m.getNCol(); ++j)
                this->data_[beginRow+i][beginCol+j] = m(i, j);
    }

    template<class MatrixType>
    void getSubMatrix(unsigned int beginRow, unsigned int beginCol, MatrixType& m) {
        assert(beginRow + m.getNRow() -1 <= getNRow());
        assert(beginCol + m.getNCol() - 1 <= getNCol());

        for (unsigned int i = 0; i < m.getNRow(); ++i)
            for (unsigned int j = 0; j < m.getNCol(); ++j)
                m(i, j) = this->data_[beginRow+i][beginCol+j];
    }
    
    unsigned int getNRow() const {return Row;}
    unsigned int getNCol() const {return Col;}        

  protected:
    Real data_[Row][Col];
  };

  /** Negate the value of every element of this matrix. */
  template<typename Real, unsigned int Row, unsigned int Col> 
  inline RectMatrix<Real, Row, Col> operator -(const RectMatrix<Real, Row, Col>& m) {
    RectMatrix<Real, Row, Col> result(m);

    result.negate();

    return result;
  }
    
  /**
   * Return the sum of two matrixes  (m1 + m2). 
   * @return the sum of two matrixes
   * @param m1 the first matrix
   * @param m2 the second matrix
   */ 
  template<typename Real, unsigned int Row, unsigned int Col> 
  inline RectMatrix<Real, Row, Col> operator + (const RectMatrix<Real, Row, Col>& m1,const RectMatrix<Real, Row, Col>& m2) {
    RectMatrix<Real, Row, Col> result;

    result.add(m1, m2);

    return result;
  }
    
  /**
   * Return the difference of two matrixes  (m1 - m2). 
   * @return the sum of two matrixes
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  template<typename Real, unsigned int Row, unsigned int Col> 
  inline RectMatrix<Real, Row, Col> operator - (const RectMatrix<Real, Row, Col>& m1, const RectMatrix<Real, Row, Col>& m2) {
    RectMatrix<Real, Row, Col> result;

    result.sub(m1, m2);

    return result;
  }

  /**
   * Return the multiplication of scalra and  matrix  (m * s). 
   * @return the multiplication of a scalra and  a matrix 
   * @param m the matrix
   * @param s the scalar
   */
  template<typename Real, unsigned int Row, unsigned int Col> 
  inline RectMatrix<Real, Row, Col> operator *(const RectMatrix<Real, Row, Col>& m, Real s) {
    RectMatrix<Real, Row, Col> result;

    result.mul(s, m);

    return result;
  }

  /**
   * Return the multiplication of a scalra and  a matrix  (s * m). 
   * @return the multiplication of a scalra and  a matrix 
   * @param s the scalar
   * @param m the matrix
   */
  template<typename Real, unsigned int Row, unsigned int Col> 
  inline RectMatrix<Real, Row, Col> operator *(Real s, const RectMatrix<Real, Row, Col>& m) {
    RectMatrix<Real, Row, Col> result;

    result.mul(s, m);

    return result;
  }
    
  /**
   * Return the multiplication of two matrixes  (m1 * m2). 
   * @return the multiplication of two matrixes
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  template<typename Real, unsigned int Row, unsigned int Col, unsigned int SameDim> 
  inline RectMatrix<Real, Row, Col> operator *(const RectMatrix<Real, Row, SameDim>& m1, const RectMatrix<Real, SameDim, Col>& m2) {
    RectMatrix<Real, Row, Col> result;

    for (unsigned int i = 0; i < Row; i++)
      for (unsigned int j = 0; j < Col; j++)
	for (unsigned int k = 0; k < SameDim; k++)
	  result(i, j)  += m1(i, k) * m2(k, j);                

    return result;
  }
    
  /**
   * Return the multiplication of  a matrix and a vector  (m * v). 
   * @return the multiplication of a matrix and a vector
   * @param m the matrix
   * @param v the vector
   */
  template<typename Real, unsigned int Row, unsigned int Col>
  inline Vector<Real, Row> operator *(const RectMatrix<Real, Row, Col>& m, const Vector<Real, Col>& v) {
    Vector<Real, Row> result;

    for (unsigned int i = 0; i < Row ; i++)
      for (unsigned int j = 0; j < Col ; j++)            
	result[i] += m(i, j) * v[j];
            
    return result;                                                                 
  }

  /**
   * Return the scalar division of matrix   (m / s). 
   * @return the scalar division of matrix  
   * @param m the matrix
   * @param s the scalar
   */
  template<typename Real, unsigned int Row, unsigned int Col> 
  inline RectMatrix<Real, Row, Col> operator /(const RectMatrix<Real, Row, Col>& m, Real s) {
    RectMatrix<Real, Row, Col> result;

    result.div(s, m);

    return result;
  }    

  /**
   * Write to an output stream
   */
  template<typename Real,  unsigned int Row, unsigned int Col>
  std::ostream &operator<< ( std::ostream& o, const RectMatrix<Real, Row, Col>& m) {
    for (unsigned int i = 0; i < Row ; i++) {
      o << "(";
      for (unsigned int j = 0; j < Col ; j++) {
	o << m(i, j);
	if (j != Col -1)
	  o << "\t";
      }
      o << ")" << std::endl;
    }
    return o;        
  }    
}
#endif //MATH_RECTMATRIX_HPP
