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
 * @file DynamicVector.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */
 
#ifndef MATH_DYNAMICVECTOR_HPP
#define MATH_DYNAMICVECTOR_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>

namespace OpenMD {
    
  /**
   * @class DynamicVector DynamicVector.hpp "math/DynamicVector.hpp"
   * @brief Fix length vector class
   */
  template<typename Real, typename Alloc = std::allocator<Real> >
  class DynamicVector : public std::vector<Real, Alloc> {
  
  public:
      typedef Real value_type;
      typedef std::vector<Real, Alloc> VectorType;
      typedef typename VectorType::pointer pointer;
      typedef typename VectorType::const_pointer const_pointer;
      typedef typename VectorType::reference reference;
      typedef typename VectorType::const_reference const_reference;
      typedef typename VectorType::iterator iterator;
      typedef typename VectorType::const_iterator const_iterator;
      typedef typename VectorType::const_reverse_iterator  const_reverse_iterator;
      typedef typename VectorType::reverse_iterator reverse_iterator;
      typedef typename VectorType::size_type size_type;
      typedef typename VectorType::difference_type difference_type;
      typedef typename VectorType::allocator_type allocator_type;


      // [23.2.4.1] construct/copy/destroy
      // (assign() and get_allocator() are also listed in this section)
      /**
       *  @brief  Default constructor creates no elements.
       */      explicit
      DynamicVector(const allocator_type& alloc = allocator_type())
      : std::vector<Real, Alloc>(alloc) { }
  
      /**
       *  @brief  Create a %DynamicVector with copies of an exemplar element.
       *  @param  n  The number of elements to initially create.
       *  @param  value  An element to copy.
       *
       *  This constructor fills the %DynamicVector with @a n copies of @a value.
       */
      DynamicVector(size_type n, const value_type& value,
             const allocator_type& alloc = allocator_type())
      : std::vector<Real, Alloc>(n, value, alloc){ }

      /**
       *  @brief  Create a %DynamicVector with default elements.
       *  @param  n  The number of elements to initially create.
       *
       *  This constructor fills the %DynamicVector with @a n copies of a
       *  default-constructed element.
       */
      explicit
      DynamicVector(size_type n) : std::vector<Real, Alloc>(n) { }

      /**
       *  @brief  %Vector copy constructor.
       *  @param  x  A %DynamicVector of identical element and allocator types.
       *
       *  The newly-created %DynamicVector uses a copy of the allocation
       *  object used by @a x.  All the elements of @a x are copied,
       *  but any extra memory in
       *  @a x (for fast expansion) will not be copied.
       */
      DynamicVector(const DynamicVector& x)
      : std::vector<Real, Alloc>(x) {}

      template<typename _InputIterator>
        DynamicVector(_InputIterator first, _InputIterator last,
               const allocator_type& alloc = allocator_type())
        : std::vector<Real, Alloc>(first, last, alloc) {}
       
    inline Real operator()(unsigned int i) const{
      return (*this)[i];
    }
    
    inline Real& operator()(unsigned int i){
      return (*this)[i];
    }     
    /**
     * Tests if this vetor is equal to other vector
     * @return true if equal, otherwise return false
     * @param v vector to be compared
     */
    inline bool operator ==(const DynamicVector<Real>& v) {

      for (unsigned int i = 0; i < this->size(); i ++) {
	if (!equal((*this)[i], v[i])) {
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
    inline bool operator !=(const DynamicVector<Real>& v) {
      return !(*this == v);
    }
             
    /** Negates the value of this vector in place. */           
    inline void negate() {
      for (unsigned int i = 0; i < this->size(); i++)
	(*this)[i] = -(*this)[i];
    }

    /**
     * Sets the value of this vector to the negation of vector v1.
     * @param v1 the source vector
     */
    inline void negate(const DynamicVector<Real>& v1) {
      for (unsigned int i = 0; i < this->size(); i++)
	(*this)[i] = -v1[i];

    }
            
    /**
     * Sets the value of this vector to the sum of itself and v1 (*this += v1).
     * @param v1 the other vector
     */
    inline void add( const DynamicVector<Real>& v1 ) {
      std::transform(this->begin(), this->end(), v1.begin(),this->begin(),std::plus<Real>());
    }

    /**
     * Sets the value of this vector to the sum of v1 and v2 (*this = v1 + v2).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    inline void add( const DynamicVector<Real>& v1, const DynamicVector<Real>& v2 ) {
      std::transform(v1.begin(), v1.end(), v2.begin(),this->begin(),std::plus<Real>());
    }

    /**
     * Sets the value of this vector to the difference  of itself and v1 (*this -= v1).
     * @param v1 the other vector
     */
    inline void sub( const DynamicVector<Real>& v1 ) {
        std::transform(this->begin(), this->end(), v1.begin(),this->begin(),std::minus<Real>());
    }

    /**
     * Sets the value of this vector to the difference of vector v1 and v2 (*this = v1 - v2).
     * @param v1 the first vector
     * @param v2 the second vector
     */
    inline void sub( const DynamicVector<Real>& v1, const DynamicVector  &v2 ){
      std::transform(v1.begin(), v1.end(), v2.begin(),this->begin(),std::minus<Real>());    
    }

    /**
     * Sets the value of this vector to the scalar multiplication of itself (*this *= s).
     * @param s the scalar value
     */
    inline void mul( Real s ) {
      for (unsigned int i = 0; i < this->size(); i++)
	(*this)[i] *= s;
    }

    /**
     * Sets the value of this vector to the scalar multiplication of vector v1  
     * (*this = s * v1).
     * @param v1 the vector            
     * @param s the scalar value
     */
    inline void mul( const DynamicVector<Real>& v1, Real s) {
      this->resize(v1.size());
      for (unsigned int i = 0; i < this->size(); i++)
	(*this)[i] = s * v1[i];
    }

    /**
     * Sets the value of this vector to the scalar division of itself  (*this /= s ).
     * @param s the scalar value
     */             
    inline void div( Real s) {
      for (unsigned int i = 0; i < this->size(); i++)            
	(*this)[i] /= s;
    }

    /**
     * Sets the value of this vector to the scalar division of vector v1  (*this = v1 / s ).
     * @param v1 the source vector
     * @param s the scalar value
     */                         
    inline void div( const DynamicVector<Real>& v1, Real s ) {
      for (unsigned int i = 0; i < this->size(); i++)
	(*this)[i] = v1[i] / s;
    }

    /** @see #add */
    inline DynamicVector<Real>& operator +=( const DynamicVector<Real>& v1 ) {
      add(v1);
      return *this;
    }

    /** @see #sub */
    inline DynamicVector<Real>& operator -=( const DynamicVector<Real>& v1 ) {
      sub(v1);
      return *this;
    }

    /** @see #mul */
    inline DynamicVector<Real>& operator *=( Real s) {
      mul(s);
      return *this;
    }

    /** @see #div */
    inline DynamicVector<Real>& operator /=( Real s ) {
      div(s);
      return *this;
    }

    /** zero out the vector */
    inline void setZero( ) {
      for (unsigned int i = 0; i < this->size(); i++)
	(*this)[i] = 0;
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
      return equal(lengthSquare(), 1.0);
    }           

    template<class VectorType>
    void getSubVector(unsigned int beginning, VectorType& v) {
        assert(beginning + v.size() -1 <= this->size());

        for (unsigned int i = 0; i < v.size(); ++i)
             v(i) = (*this)[beginning+i];
    }

    
  };

  /** unary minus*/
  template<typename Real>    
  inline DynamicVector<Real> operator -(const DynamicVector<Real>& v1){
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
  inline DynamicVector<Real> operator +(const DynamicVector<Real>& v1, const DynamicVector<Real>& v2) {
    assert(v1.size() == v2.size());
    DynamicVector<Real>result(v1.size());
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
  DynamicVector<Real> operator -(const DynamicVector<Real>& v1, const DynamicVector<Real>& v2) {
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
  DynamicVector<Real> operator *( const DynamicVector<Real>& v1, Real s) {      
    DynamicVector<Real> result(v1.size());
    result.mul(v1,s);
    return result;           
  }
    
  /**
   * Returns the vaule of scalar multiplication of this vector v1 (v1 * r). 
   * @return  the vaule of scalar multiplication of this vector
   * @param s the scalar value
   * @param v1 the source vector
   */  
  template<typename Real>
  DynamicVector<Real> operator *( Real s, const DynamicVector<Real>& v1 ) {
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
  DynamicVector<Real> operator / ( const DynamicVector<Real>& v1, Real s) {     
    DynamicVector<Real> result(v1.size());
    result.div( v1,s);
    return result;           
  }
    
  /**
   * Returns the dot product of two DynamicVectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the dot product of v1 and v2
   */
  template<typename Real>    
  inline Real dot( const DynamicVector<Real>& v1, const DynamicVector<Real>& v2 ) {
    Real tmp;
    tmp = 0;
    assert(v1.size() == v2.size());
    for (unsigned int i = 0; i < v1.size(); i++)
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
  inline Real distance( const DynamicVector<Real>& v1, const DynamicVector<Real>& v2 ) {
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
  inline Real distanceSquare( const DynamicVector<Real>& v1, const DynamicVector<Real>& v2 ) {
    DynamicVector<Real> tempDynamicVector = v1 - v2;
    return tempDynamicVector.lengthSquare();
  }

  /**
   * Write to an output stream
   */
  template<typename Real>
  std::ostream &operator<< ( std::ostream& o, const DynamicVector<Real>& v) {

    o << "[ ";
        
    for (unsigned int i = 0 ; i< v.size(); i++) {
      o << v[i];

      if (i  != v.size() -1) {
	o<< ", ";
      }
    }

    o << " ]";
    return o;        
  }
    
}
#endif

