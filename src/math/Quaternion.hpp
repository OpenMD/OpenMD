/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
/**
 * @file Quaternion.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */

#ifndef MATH_QUATERNION_HPP
#define MATH_QUATERNION_HPP

#include "math/Vector.hpp"
#include "math/SquareMatrix.hpp"

namespace oopse{

  /**
   * @class Quaternion Quaternion.hpp "math/Quaternion.hpp"
   * Quaternion is a sort of a higher-level complex number.
   * It is defined as Q = w + x*i + y*j + z*k,
   * where w, x, y, and z are numbers of type T (e.g. RealType), and
   * i*i = -1; j*j = -1; k*k = -1;
   * i*j = k; j*k = i; k*i = j;
   */
  template<typename Real>
  class Quaternion : public Vector<Real, 4> {
  public:
    Quaternion() : Vector<Real, 4>() {}

    /** Constructs and initializes a Quaternion from w, x, y, z values */     
    Quaternion(Real w, Real x, Real y, Real z) {
      this->data_[0] = w;
      this->data_[1] = x;
      this->data_[2] = y;
      this->data_[3] = z;                
    }
            
    /** Constructs and initializes a Quaternion from a  Vector<Real,4> */     
    Quaternion(const Vector<Real,4>& v) 
      : Vector<Real, 4>(v){
      }

    /** copy assignment */
    Quaternion& operator =(const Vector<Real, 4>& v){
      if (this == & v)
	return *this;

      Vector<Real, 4>::operator=(v);
                
      return *this;
    }

    /**
     * Returns the value of the first element of this quaternion.
     * @return the value of the first element of this quaternion
     */
    Real w() const {
      return this->data_[0];
    }

    /**
     * Returns the reference of the first element of this quaternion.
     * @return the reference of the first element of this quaternion
     */
    Real& w() {
      return this->data_[0];    
    }

    /**
     * Returns the value of the first element of this quaternion.
     * @return the value of the first element of this quaternion
     */
    Real x() const {
      return this->data_[1];
    }

    /**
     * Returns the reference of the second element of this quaternion.
     * @return the reference of the second element of this quaternion
     */
    Real& x() {
      return this->data_[1];    
    }

    /**
     * Returns the value of the thirf element of this quaternion.
     * @return the value of the third element of this quaternion
     */
    Real y() const {
      return this->data_[2];
    }

    /**
     * Returns the reference of the third element of this quaternion.
     * @return the reference of the third element of this quaternion
     */           
    Real& y() {
      return this->data_[2];    
    }

    /**
     * Returns the value of the fourth element of this quaternion.
     * @return the value of the fourth element of this quaternion
     */
    Real z() const {
      return this->data_[3];
    }
    /**
     * Returns the reference of the fourth element of this quaternion.
     * @return the reference of the fourth element of this quaternion
     */
    Real& z() {
      return this->data_[3];    
    }

    /**
     * Tests if this quaternion is equal to other quaternion
     * @return true if equal, otherwise return false
     * @param q quaternion to be compared
     */
    inline bool operator ==(const Quaternion<Real>& q) {

      for (unsigned int i = 0; i < 4; i ++) {
	if (!equal(this->data_[i], q[i])) {
	  return false;
	}
      }
                
      return true;
    }
            
    /**
     * Returns the inverse of this quaternion
     * @return inverse
     * @note since quaternion is a complex number, the inverse of quaternion
     * q = w + xi + yj+ zk is inv_q = (w -xi - yj - zk)/(|q|^2)
     */
    Quaternion<Real> inverse() {
      Quaternion<Real> q;
      Real d = this->lengthSquare();
                
      q.w() = w() / d;
      q.x() = -x() / d;
      q.y() = -y() / d;
      q.z() = -z() / d;
                
      return q;
    }

    /**
     * Sets the value to the multiplication of itself and another quaternion
     * @param q the other quaternion
     */
    void mul(const Quaternion<Real>& q) {
      Quaternion<Real> tmp(*this);

      this->data_[0] = (tmp[0]*q[0]) -(tmp[1]*q[1]) - (tmp[2]*q[2]) - (tmp[3]*q[3]);
      this->data_[1] = (tmp[0]*q[1]) + (tmp[1]*q[0]) + (tmp[2]*q[3]) - (tmp[3]*q[2]);
      this->data_[2] = (tmp[0]*q[2]) + (tmp[2]*q[0]) + (tmp[3]*q[1]) - (tmp[1]*q[3]);
      this->data_[3] = (tmp[0]*q[3]) + (tmp[3]*q[0]) + (tmp[1]*q[2]) - (tmp[2]*q[1]);                
    }

    void mul(const Real& s) {
      this->data_[0] *= s;
      this->data_[1] *= s;
      this->data_[2] *= s;
      this->data_[3] *= s;
    }

    /** Set the value of this quaternion to the division of itself by another quaternion */
    void div(Quaternion<Real>& q) {
      mul(q.inverse());
    }

    void div(const Real& s) {
      this->data_[0] /= s;
      this->data_[1] /= s;
      this->data_[2] /= s;
      this->data_[3] /= s;
    }
            
    Quaternion<Real>& operator *=(const Quaternion<Real>& q) {
      mul(q);
      return *this;
    }

    Quaternion<Real>& operator *=(const Real& s) {
      mul(s);
      return *this;
    }
            
    Quaternion<Real>& operator /=(Quaternion<Real>& q) {                
      *this *= q.inverse();
      return *this;
    }

    Quaternion<Real>& operator /=(const Real& s) {
      div(s);
      return *this;
    }            
    /**
     * Returns the conjugate quaternion of this quaternion
     * @return the conjugate quaternion of this quaternion
     */
    Quaternion<Real> conjugate() {
      return Quaternion<Real>(w(), -x(), -y(), -z());            
    }

    /**
     * Returns the corresponding rotation matrix (3x3)
     * @return a 3x3 rotation matrix
     */
    SquareMatrix<Real, 3> toRotationMatrix3() {
      SquareMatrix<Real, 3> rotMat3;

      Real w2;
      Real x2;
      Real y2;
      Real z2;

      if (!this->isNormalized())
	this->normalize();
                
      w2 = w() * w();
      x2 = x() * x();
      y2 = y() * y();
      z2 = z() * z();

      rotMat3(0, 0) = w2 + x2 - y2 - z2;
      rotMat3(0, 1) = 2.0 * ( x() * y() + w() * z() );
      rotMat3(0, 2) = 2.0 * ( x() * z() - w() * y() );

      rotMat3(1, 0) = 2.0 * ( x() * y() - w() * z() );
      rotMat3(1, 1) = w2 - x2 + y2 - z2;
      rotMat3(1, 2) = 2.0 * ( y() * z() + w() * x() );

      rotMat3(2, 0) = 2.0 * ( x() * z() + w() * y() );
      rotMat3(2, 1) = 2.0 * ( y() * z() - w() * x() );
      rotMat3(2, 2) = w2 - x2 -y2 +z2;

      return rotMat3;
    }

  };//end Quaternion


    /**
     * Returns the vaule of scalar multiplication of this quaterion q (q * s). 
     * @return  the vaule of scalar multiplication of this vector
     * @param q the source quaternion
     * @param s the scalar value
     */
  template<typename Real, unsigned int Dim>                 
  Quaternion<Real> operator * ( const Quaternion<Real>& q, Real s) {       
    Quaternion<Real> result(q);
    result.mul(s);
    return result;           
  }
    
  /**
   * Returns the vaule of scalar multiplication of this quaterion q (q * s). 
   * @return  the vaule of scalar multiplication of this vector
   * @param s the scalar value
   * @param q the source quaternion
   */  
  template<typename Real, unsigned int Dim>
  Quaternion<Real> operator * ( const Real& s, const Quaternion<Real>& q ) {
    Quaternion<Real> result(q);
    result.mul(s);
    return result;           
  }    

  /**
   * Returns the multiplication of two quaternion
   * @return the multiplication of two quaternion
   * @param q1 the first quaternion
   * @param q2 the second quaternion
   */
  template<typename Real>
  inline Quaternion<Real> operator *(const Quaternion<Real>& q1, const Quaternion<Real>& q2) {
    Quaternion<Real> result(q1);
    result *= q2;
    return result;
  }

  /**
   * Returns the division of two quaternion
   * @param q1 divisor
   * @param q2 dividen
   */

  template<typename Real>
  inline Quaternion<Real> operator /( Quaternion<Real>& q1,  Quaternion<Real>& q2) {
    return q1 * q2.inverse();
  }

  /**
   * Returns the value of the division of a scalar by a quaternion
   * @return the value of the division of a scalar by a quaternion
   * @param s scalar
   * @param q quaternion
   * @note for a quaternion q, 1/q = q.inverse()
   */
  template<typename Real>
  Quaternion<Real> operator /(const Real& s,  Quaternion<Real>& q) {

    Quaternion<Real> x;
    x = q.inverse();
    x *= s;
    return x;
  }
    
  template <class T>
  inline bool operator==(const Quaternion<T>& lhs, const Quaternion<T>& rhs) {
    return equal(lhs[0] ,rhs[0]) && equal(lhs[1] , rhs[1]) && equal(lhs[2], rhs[2]) && equal(lhs[3], rhs[3]);
  }
    
  typedef Quaternion<RealType> Quat4d;
}
#endif //MATH_QUATERNION_HPP 
