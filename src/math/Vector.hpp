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
 * @file Vector.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */
 
#ifndef MATH_VECTOR_HPP
#define MATH_VECTOR_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <math.h>
#include "config.h"
namespace OpenMD {

  static const RealType epsilon = 0.000001;

  template<typename T>
  inline bool equal(T e1, T e2) {
    return e1 == e2;
  }

  //template<>
  //inline bool equal(float e1, float e2) {
  //  return fabs(e1 - e2) < epsilon;
  //}

  template<>
  inline bool equal(RealType e1, RealType e2) {
    return fabs(e1 - e2) < epsilon;
  }
    
  /**
   * @class Vector Vector.hpp "math/Vector.hpp"
   * @brief Fix length vector class
   */
  template<typename Real, unsigned int Dim>
  class Vector{
  public:

    typedef Real ElemType;
    typedef Real* ElemPoinerType;

    /** default constructor */
    inline Vector(){
      for (unsigned int i = 0; i < Dim; i++)
	this->data_[i] = 0;
    }

    /** Constructs and initializes a Vector from a vector */
    inline Vector(const Vector<Real, Dim>& v) {
      *this  = v;
    }

    /** copy assignment operator */
    inline Vector<Real, Dim>& operator=(const Vector<Real, Dim>& v) {
      if (this == &v)
	return *this;
                
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
    inline Vector( Real* v) {
      for (unsigned int i = 0; i < Dim; i++)
	this->data_[i] = v[i];
    }

    /** 
     * Returns reference of ith element.
     * @return reference of ith element
     * @param i index
     */
    inline Real& operator[](unsigned int  i) {
      assert( i < Dim);
      return this->data_[i];
    }

    /** 
     * Returns reference of ith element.
     * @return reference of ith element
     * @param i index
     */
    inline Real& operator()(unsigned int  i) {
      assert( i < Dim);
      return this->data_[i];
    }

    /** 
     * Returns constant reference of ith element.
     * @return reference of ith element
     * @param i index
     */
    inline  const Real& operator[](unsigned int i) const {
      assert( i < Dim);
      return this->data_[i];
    }

    /** 
     * Returns constant reference of ith element.
     * @return reference of ith element
     * @param i index
     */
    inline  const Real& operator()(unsigned int i) const {
      assert( i < Dim);
      return this->data_[i];
    }

    /** Copy the internal data to an array*/
    void getArray(Real* array) {
      for (unsigned int i = 0; i < Dim; i ++) {
	array[i] = this->data_[i];
      }                
    }

    /** Returns the pointer of internal array */
    Real* getArrayPointer() {
      return this->data_;
    }
            
    /**
     * Tests if this vetor is equal to other vector
     * @return true if equal, otherwise return false
     * @param v vector to be compared
     */
    inline bool operator ==(const Vector<Real, Dim>& v) {

      for (unsigned int i = 0; i < Dim; i ++) {
	if (!equal(this->data_[i], v[i])) {
	  return false;
	}
      }
                
      return true;
    }

    /**
     * Tests if this vetor is not equal to other vector
     * @return true if equal, otherwise return false
     * @param v vector to be compared
     */
    inline bool operator !=(const Vector<Real, Dim>& v) {
      return !(*this == v);
    }

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
    inline void add( const Vector<Real, Dim>& v1 ) {
      for (unsigned int i = 0; i < Dim; i++)
	this->data_[i] += v1.data_[i];
    }

    /**
     * Sets the value of this vector to the sum of v1 and v2 (*this = v1 + v2).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    inline void add( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {
      for (unsigned int i = 0; i < Dim; i++)
	this->data_[i] = v1.data_[i] + v2.data_[i];
    }

    /**
     * Sets the value of this vector to the difference  of itself and v1 (*this -= v1).
     * @param v1 the other vector
     */
    inline void sub( const Vector<Real, Dim>& v1 ) {
      for (unsigned int i = 0; i < Dim; i++)
	this->data_[i] -= v1.data_[i];
    }

    /**
     * Sets the value of this vector to the difference of vector v1 and v2 (*this = v1 - v2).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    inline void sub( const Vector<Real, Dim>& v1, const Vector  &v2 ){
      for (unsigned int i = 0; i < Dim; i++)
	this->data_[i] = v1.data_[i] - v2.data_[i];
    }

    /**
     * Sets the value of this vector to the scalar multiplication of itself (*this *= s).
     * @param s the scalar value
     */
    inline void mul( Real s ) {
      for (unsigned int i = 0; i < Dim; i++)
	this->data_[i] *= s;
    }

    /**
     * Sets the value of this vector to the scalar multiplication of vector v1  
     * (*this = s * v1).
     * @param v1 the vector            
     * @param s the scalar value
     */
    inline void mul( const Vector<Real, Dim>& v1, Real s) {
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
    inline void Vmul( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
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
     * Sets the value of this vector to the scalar division of itself  (*this /= s ).
     * @param s the scalar value
     */             
    inline void div( Real s) {
      for (unsigned int i = 0; i < Dim; i++)            
	this->data_[i] /= s;
    }

    /**
     * Sets the value of this vector to the scalar division of vector v1  (*this = v1 / s ).
     * @param v1 the source vector
     * @param s the scalar value
     */                         
    inline void div( const Vector<Real, Dim>& v1, Real s ) {
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
    inline void Vdiv( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
      for (unsigned int i = 0; i < Dim; i++)
	this->data_[i] = v1.data_[i] / v2.data_[i];
    }


    /** @see #add */
    inline Vector<Real, Dim>& operator +=( const Vector<Real, Dim>& v1 ) {
      add(v1);
      return *this;
    }

    /** @see #sub */
    inline Vector<Real, Dim>& operator -=( const Vector<Real, Dim>& v1 ) {
      sub(v1);
      return *this;
    }

    /** @see #mul */
    inline Vector<Real, Dim>& operator *=( Real s) {
      mul(s);
      return *this;
    }

    /** @see #div */
    inline Vector<Real, Dim>& operator /=( Real s ) {
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
    inline Real length() {
      return sqrt(lengthSquare());  
    }
            
    /**
     * Returns the squared length of this vector.
     * @return the squared length of this vector
     */
    inline Real lengthSquare() {
      return dot(*this, *this);
    }
            
    /** Normalizes this vector in place */
    inline void normalize() {
      Real len;

      len = length();
                
      //if (len < OpenMD::NumericConstant::epsilon)
      //  throw();
                
      *this /= len;
    }

    /**
     * Tests if this vector is normalized
     * @return true if this vector is normalized, otherwise return false
     */
    inline bool isNormalized() {
      return equal(lengthSquare(), (RealType)1);
    }           

    unsigned int size() {return Dim;}
  protected:
    Real data_[Dim];
        
  };

  /** unary minus*/
  template<typename Real, unsigned int Dim>    
  inline Vector<Real, Dim> operator -(const Vector<Real, Dim>& v1){
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
  inline Vector<Real, Dim> operator +(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
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
  Vector<Real, Dim> operator -(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
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
  Vector<Real, Dim> operator * ( const Vector<Real, Dim>& v1, Real s) {       
    Vector<Real, Dim> result;
    result.mul(v1,s);
    return result;           
  }
    
  /**
   * Returns the vaule of scalar multiplication of this vector v1 (v1 * r). 
   * @return  the vaule of scalar multiplication of this vector
   * @param s the scalar value
   * @param v1 the source vector
   */  
  template<typename Real, unsigned int Dim>
  Vector<Real, Dim> operator * ( Real s, const Vector<Real, Dim>& v1 ) {
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
  Vector<Real, Dim> operator / ( const Vector<Real, Dim>& v1, Real s) {       
    Vector<Real, Dim> result;
    result.div( v1,s);
    return result;           
  }
    
  /**
   * Returns the dot product of two Vectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the dot product of v1 and v2
   */
  template<typename Real, unsigned int Dim>    
  inline Real dot( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {
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
  inline Real dot( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2, const Vector<Real, Dim>& v3 ) {
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
  inline Real distance( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {
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
  inline Real distanceSquare( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {
    Vector<Real, Dim> tempVector = v1 - v2;
    return tempVector.lengthSquare();
  }

  /**
   * Write to an output stream
   */
  template<typename Real, unsigned int Dim>
  std::ostream &operator<< ( std::ostream& o, const Vector<Real, Dim>& v) {

    o << "[ ";
        
    for (unsigned int i = 0 ; i< Dim; i++) {
      o << v[i];

      if (i  != Dim -1) {
	o<< ", ";
      }
    }

    o << " ]";
    return o;        
  }
    
}
#endif
