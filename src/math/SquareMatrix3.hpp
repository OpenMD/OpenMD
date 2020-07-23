/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * @file SquareMatrix3.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */

#ifndef MATH_SQUAREMATRIX3_HPP
#define  MATH_SQUAREMATRIX3_HPP
#include "config.h"
#include <cmath>
#include <vector>
#include <limits>
#include "Quaternion.hpp"
#include "SquareMatrix.hpp"
#include "Vector3.hpp"
namespace OpenMD {

  template<typename Real>
  class SquareMatrix3 : public SquareMatrix<Real, 3> {
  public:

    typedef Real ElemType;
    typedef Real* ElemPoinerType;
            
    /** default constructor */
    SquareMatrix3() : SquareMatrix<Real, 3>() {
    }

    /** Constructs and initializes every element of this matrix to a scalar */ 
    SquareMatrix3(Real s) : SquareMatrix<Real,3>(s){
    }

    /** Constructs and initializes from an array */ 
    SquareMatrix3(Real* array) : SquareMatrix<Real,3>(array){
    }


    /** copy  constructor */
    SquareMatrix3(const SquareMatrix<Real, 3>& m)  : SquareMatrix<Real, 3>(m) {
    }
            
    SquareMatrix3( const Vector3<Real>& eulerAngles) {
      setupRotMat(eulerAngles);
    }
            
    SquareMatrix3(Real phi, Real theta, Real psi) {
      setupRotMat(phi, theta, psi);
    }

    SquareMatrix3(const Quaternion<Real>& q) {
      setupRotMat(q);

    }

    SquareMatrix3(Real w, Real x, Real y, Real z) {
      setupRotMat(w, x, y, z);
    }
            
    /** copy assignment operator */
    SquareMatrix3<Real>& operator =(const SquareMatrix<Real, 3>& m) {
      if (this == &m)
	return *this;
      SquareMatrix<Real, 3>::operator=(m);
      return *this;
    }


    SquareMatrix3<Real>& operator =(const Quaternion<Real>& q) {
      this->setupRotMat(q);
      return *this;
    }

    
    /**
     * Sets this matrix to a rotation matrix by three euler angles
     * @ param euler
     */
    void setupRotMat(const Vector3<Real>& eulerAngles) {
      setupRotMat(eulerAngles[0], eulerAngles[1], eulerAngles[2]);
    }

    /**
     * Sets this matrix to a rotation matrix by three euler angles
     * @param phi
     * @param theta
     * @param psi
     */
    void setupRotMat(Real phi, Real theta, Real psi) {
      Real sphi, stheta, spsi;
      Real cphi, ctheta, cpsi;

      sphi = sin(phi);
      stheta = sin(theta);
      spsi = sin(psi);
      cphi = cos(phi);
      ctheta = cos(theta);
      cpsi = cos(psi);

      this->data_[0][0] = cpsi * cphi - ctheta * sphi * spsi;
      this->data_[0][1] = cpsi * sphi + ctheta * cphi * spsi;
      this->data_[0][2] = spsi * stheta;
                
      this->data_[1][0] = -spsi * cphi - ctheta * sphi * cpsi;
      this->data_[1][1] = -spsi * sphi + ctheta * cphi * cpsi;
      this->data_[1][2] = cpsi * stheta;

      this->data_[2][0] = stheta * sphi;
      this->data_[2][1] = -stheta * cphi;
      this->data_[2][2] = ctheta;
    }


    /**
     * Sets this matrix to a rotation matrix by quaternion
     * @param quat
     */
    void setupRotMat(const Quaternion<Real>& quat) {
      setupRotMat(quat.w(), quat.x(), quat.y(), quat.z());
    }

    /**
     * Sets this matrix to a rotation matrix by quaternion
     * @param w the first element 
     * @param x the second element
     * @param y the third element
     * @param z the fourth element
     */
    void setupRotMat(Real w, Real x, Real y, Real z) {
      Quaternion<Real> q(w, x, y, z);
      *this = q.toRotationMatrix3();
    }

    void setupSkewMat(Vector3<Real> v) {
      setupSkewMat(v[0], v[1], v[2]);
    }

    void setupSkewMat(Real v1, Real v2, Real v3) {
      this->data_[0][0] = 0;
      this->data_[0][1] = -v3;
      this->data_[0][2] = v2;
      this->data_[1][0] = v3;
      this->data_[1][1] = 0;
      this->data_[1][2] = -v1;
      this->data_[2][0] = -v2;
      this->data_[2][1] = v1;
      this->data_[2][2] = 0;               
    }

    /**
     * Sets this matrix to a symmetric tensor using Voigt Notation
     * @param vt
     */
    void setupVoigtTensor(Vector<Real, 6> vt) {
      setupVoigtTensor(vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
    }
    
    void setupVoigtTensor(Real v1, Real v2, Real v3, Real v4, Real v5,
                          Real v6) {
      this->data_[0][0] = v1;
      this->data_[1][1] = v2;
      this->data_[2][2] = v3;
      this->data_[1][2] = v4;
      this->data_[2][1] = v4;
      this->data_[0][2] = v5;
      this->data_[2][0] = v5;
      this->data_[0][1] = v6;
      this->data_[1][0] = v6;
    }
    
    /**
     * Sets this matrix to an upper-triangular (asymmetric) tensor
     * using Voigt Notation
     * @param vt
     */
    void setupUpperTriangularVoigtTensor(Vector<Real, 6> vt) {
      setupUpperTriangularVoigtTensor(vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
    }
    
    void setupUpperTriangularVoigtTensor(Real v1, Real v2, Real v3, Real v4,
                                         Real v5, Real v6) {
      this->data_[0][0] = v1;
      this->data_[1][1] = v2;
      this->data_[2][2] = v3;
      this->data_[1][2] = v4;
      this->data_[0][2] = v5;
      this->data_[0][1] = v6;
    }
         
    /**
     * Uses Rodrigues' rotation formula for a rotation matrix.
     * @param axis the axis to rotate around
     * @param angle the angle to rotate (in radians)
     */
    void axisAngle(Vector3d axis, RealType angle) {
      axis.normalize();
      RealType ct = cos(angle);
      RealType st = sin(angle);

      *this = ct * SquareMatrix3<Real>::identity();
      *this += st * SquareMatrix3<Real>::setupSkewMat(axis);
      *this += (1-ct) * SquareMatrix3<Real>::outProduct(axis, axis);      
    }
    
    /**
     * Returns the quaternion from this rotation matrix
     * @return the quaternion from this rotation matrix
     * @exception invalid rotation matrix
     */            
    Quaternion<Real> toQuaternion() {
      Quaternion<Real> q;
      Real t, s;
      Real ad1, ad2, ad3;    
      t = this->data_[0][0] + this->data_[1][1] + this->data_[2][2] + 1.0;

      if( t > std::numeric_limits<RealType>::epsilon() ){

	s = 0.5 / sqrt( t );
	q[0] = 0.25 / s;
	q[1] = (this->data_[1][2] - this->data_[2][1]) * s;
	q[2] = (this->data_[2][0] - this->data_[0][2]) * s;
	q[3] = (this->data_[0][1] - this->data_[1][0]) * s;
      } else {

	ad1 = this->data_[0][0];
	ad2 = this->data_[1][1];
	ad3 = this->data_[2][2];

	if( ad1 >= ad2 && ad1 >= ad3 ){

	  s = 0.5 / sqrt( 1.0 + this->data_[0][0] - this->data_[1][1] -
                          this->data_[2][2] );
	  q[0] = (this->data_[1][2] - this->data_[2][1]) * s;
	  q[1] = 0.25 / s;
	  q[2] = (this->data_[0][1] + this->data_[1][0]) * s;
	  q[3] = (this->data_[0][2] + this->data_[2][0]) * s;
	} else if ( ad2 >= ad1 && ad2 >= ad3 ) {
	  s = 0.5 / sqrt( 1.0 + this->data_[1][1] - this->data_[0][0] -
                          this->data_[2][2] );
	  q[0] = (this->data_[2][0] - this->data_[0][2] ) * s;
	  q[1] = (this->data_[0][1] + this->data_[1][0]) * s;
	  q[2] = 0.25 / s;
	  q[3] = (this->data_[1][2] + this->data_[2][1]) * s;
	} else {

	  s = 0.5 / sqrt( 1.0 + this->data_[2][2] - this->data_[0][0] -
                          this->data_[1][1] );
	  q[0] = (this->data_[0][1] - this->data_[1][0]) * s;
	  q[1] = (this->data_[0][2] + this->data_[2][0]) * s;
	  q[2] = (this->data_[1][2] + this->data_[2][1]) * s;
	  q[3] = 0.25 / s;
	}
      }             

      return q;
                
    }

    /**
     * Returns the euler angles from this rotation matrix
     * @return the euler angles in a vector 
     * @exception invalid rotation matrix
     * We use so-called "x-convention", which is the most common definition. 
     * In this convention, the rotation given by Euler angles (phi, theta, 
     * psi), where the first rotation is by an angle phi about the z-axis, 
     * the second is by an angle theta (0 <= theta <= 180) about the x-axis, 
     * and the third is by an angle psi about the z-axis (again). 
     */            
    Vector3<Real> toEulerAngles() {
      Vector3<Real> myEuler;
      Real phi;
      Real theta;
      Real psi;
      Real ctheta;
      Real stheta;
                
      // set the tolerance for Euler angles and rotation elements

      theta = acos(std::min((RealType)1.0,
                            std::max((RealType)-1.0, this->data_[2][2])));
      ctheta = this->data_[2][2]; 
      stheta = sqrt(1.0 - ctheta * ctheta);

      // when sin(theta) is close to 0, we need to consider
      // singularity In this case, we can assign an arbitary value to
      // phi (or psi), and then determine the psi (or phi) or
      // vice-versa. We'll assume that phi always gets the rotation,
      // and psi is 0 in cases of singularity.
      // we use atan2 instead of atan, since atan2 will give us -Pi to Pi. 
      // Since 0 <= theta <= 180, sin(theta) will be always
      // non-negative. Therefore, it will never change the sign of both of
      // the parameters passed to atan2.

      if (fabs(stheta) < 1e-6){
	psi = 0.0;
	phi = atan2(-this->data_[1][0], this->data_[0][0]);  
      }
      // we only have one unique solution
      else{    
	phi = atan2(this->data_[2][0], -this->data_[2][1]);
	psi = atan2(this->data_[0][2], this->data_[1][2]);
      }

      //wrap phi and psi, make sure they are in the range from 0 to 2*Pi
      if (phi < 0)
	phi += 2.0 * Constants::PI;

      if (psi < 0)
	psi += 2.0 * Constants::PI;

      myEuler[0] = phi;
      myEuler[1] = theta;
      myEuler[2] = psi;

      return myEuler;
    }

    
    Vector<Real, 6> toVoigtTensor() {
      Vector<Real, 6> voigt;
      voigt[0] = this->data_[0][0];
      voigt[1] = this->data_[1][1];
      voigt[2] = this->data_[2][2];
      voigt[3] = 0.5 * (this->data_[1][2] + this->data_[2][1]);
      voigt[4] = 0.5 * (this->data_[0][2] + this->data_[2][0]);
      voigt[5] = 0.5 * (this->data_[0][1] + this->data_[1][0]);
      return voigt;
    }
            
    /** Returns the determinant of this matrix. */
    Real determinant() const {
      Real x,y,z;

      x = this->data_[0][0] * (this->data_[1][1] * this->data_[2][2] -
                               this->data_[1][2] * this->data_[2][1]);
      y = this->data_[0][1] * (this->data_[1][2] * this->data_[2][0] -
                               this->data_[1][0] * this->data_[2][2]);
      z = this->data_[0][2] * (this->data_[1][0] * this->data_[2][1] -
                               this->data_[1][1] * this->data_[2][0]);
      return(x + y + z);
    }            

    /** Returns the trace of this matrix. */
    Real trace() const {
      return this->data_[0][0] + this->data_[1][1] + this->data_[2][2];
    }
            
    /**
     * Sets the value of this matrix to the inverse of itself.
     * @note since this simple algorithm can be applied to invert a 3 by 3
     * matrix, we hide the implementation of inverse in SquareMatrix
     * class
     */
    SquareMatrix3<Real>  inverse() const {
      SquareMatrix3<Real> m;
      RealType det = determinant();
      m(0, 0) = this->data_[1][1]*this->data_[2][2]
        - this->data_[1][2]*this->data_[2][1];
      m(1, 0) = this->data_[1][2]*this->data_[2][0]
        - this->data_[1][0]*this->data_[2][2];
      m(2, 0) = this->data_[1][0]*this->data_[2][1]
        - this->data_[1][1]*this->data_[2][0];
      m(0, 1) = this->data_[2][1]*this->data_[0][2]
        - this->data_[2][2]*this->data_[0][1];
      m(1, 1) = this->data_[2][2]*this->data_[0][0]
        - this->data_[2][0]*this->data_[0][2];
      m(2, 1) = this->data_[2][0]*this->data_[0][1]
        - this->data_[2][1]*this->data_[0][0];
      m(0, 2) = this->data_[0][1]*this->data_[1][2]
        - this->data_[0][2]*this->data_[1][1];
      m(1, 2) = this->data_[0][2]*this->data_[1][0]
        - this->data_[0][0]*this->data_[1][2];
      m(2, 2) = this->data_[0][0]*this->data_[1][1]
        - this->data_[0][1]*this->data_[1][0];
      m /= det;
      return m;
    }

    SquareMatrix3<Real> transpose() const{
      SquareMatrix3<Real> result;
                
      for (unsigned int i = 0; i < 3; i++)
	for (unsigned int j = 0; j < 3; j++)              
	  result(j, i) = this->data_[i][j];

      return result;
    }
    /**
     * Extract the eigenvalues and eigenvectors from a 3x3 matrix.
     * The eigenvectors (the columns of V) will be normalized. 
     * The eigenvectors are aligned optimally with the x, y, and z
     * axes respectively.
     * @param a symmetric matrix whose eigenvectors are to be computed. On return, the matrix is overwritten             
     * @param w will contain the eigenvalues of the matrix On return of this function
     * @param v the columns of this matrix will contain the eigenvectors. The eigenvectors are normalized and mutually orthogonal.              
     * @warning a will be overwritten
     */
    static void diagonalize(SquareMatrix3<Real>& a, Vector3<Real>& w,
                            SquareMatrix3<Real>& v); 
  };
  /*=========================================================================

    Program:   Visualization Toolkit
    Module:    $RCSfile: SquareMatrix3.hpp,v $

    Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
    All rights reserved.
    See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

    This software is distributed WITHOUT ANY WARRANTY; without even
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
    PURPOSE.  See the above copyright notice for more information.

    =========================================================================*/
  template<typename Real>
  void SquareMatrix3<Real>::diagonalize(SquareMatrix3<Real>& a,
                                        Vector3<Real>& w, 
					SquareMatrix3<Real>& v) {
    int i,j,k,maxI;
    Real tmp, maxVal;
    Vector3<Real> v_maxI, v_k, v_j;

    // diagonalize using Jacobi
    SquareMatrix3<Real>::jacobi(a, w, v);
    // if all the eigenvalues are the same, return identity matrix
    if (w[0] == w[1] && w[0] == w[2] ) {
      v = SquareMatrix3<Real>::identity();
      return;
    }

    // transpose temporarily, it makes it easier to sort the eigenvectors
    v = v.transpose(); 
        
    // if two eigenvalues are the same, re-orthogonalize to optimally line
    // up the eigenvectors with the x, y, and z axes
    for (i = 0; i < 3; i++) {
      if (w((i+1)%3) == w((i+2)%3)) {// two eigenvalues are the same
	// find maximum element of the independant eigenvector
	maxVal = fabs(v(i, 0));
	maxI = 0;
	for (j = 1; j < 3; j++) {
	  if (maxVal < (tmp = fabs(v(i, j)))){
	    maxVal = tmp;
	    maxI = j;
	  }
	}
            
	// swap the eigenvector into its proper position
	if (maxI != i) {
	  tmp = w(maxI);
	  w(maxI) = w(i);
	  w(i) = tmp;

	  v.swapRow(i, maxI);
	}
	// maximum element of eigenvector should be positive
	if (v(maxI, maxI) < 0) {
	  v(maxI, 0) = -v(maxI, 0);
	  v(maxI, 1) = -v(maxI, 1);
	  v(maxI, 2) = -v(maxI, 2);
	}

	// re-orthogonalize the other two eigenvectors
	j = (maxI+1)%3;
	k = (maxI+2)%3;

	v(j, 0) = 0.0; 
	v(j, 1) = 0.0; 
	v(j, 2) = 0.0;
	v(j, j) = 1.0;

	/** @todo */
	v_maxI = v.getRow(maxI);
	v_j = v.getRow(j);
	v_k = cross(v_maxI, v_j);
	v_k.normalize();
	v_j = cross(v_k, v_maxI);
	v.setRow(j, v_j);
	v.setRow(k, v_k);


	// transpose vectors back to columns
	v = v.transpose();
	return;
      }
    }

    // the three eigenvalues are different, just sort the eigenvectors
    // to align them with the x, y, and z axes

    // find the vector with the largest x element, make that vector
    // the first vector
    maxVal = fabs(v(0, 0));
    maxI = 0;
    for (i = 1; i < 3; i++) {
      if (maxVal < (tmp = fabs(v(i, 0)))) {
	maxVal = tmp;
	maxI = i;
      }
    }

    // swap eigenvalue and eigenvector
    if (maxI != 0) {
      tmp = w(maxI);
      w(maxI) = w(0);
      w(0) = tmp;
      v.swapRow(maxI, 0);
    }
    // do the same for the y element
    if (fabs(v(1, 1)) < fabs(v(2, 1))) {
      tmp = w(2);
      w(2) = w(1);
      w(1) = tmp;
      v.swapRow(2, 1);
    }

    // ensure that the sign of the eigenvectors is correct
    for (i = 0; i < 2; i++) {
      if (v(i, i) < 0) {
	v(i, 0) = -v(i, 0);
	v(i, 1) = -v(i, 1);
	v(i, 2) = -v(i, 2);
      }
    }

    // set sign of final eigenvector to ensure that determinant is positive
    if (v.determinant() < 0) {
      v(2, 0) = -v(2, 0);
      v(2, 1) = -v(2, 1);
      v(2, 2) = -v(2, 2);
    }

    // transpose the eigenvectors back again
    v = v.transpose();
    return ;
  }

  /**
   * Return the multiplication of two matrixes  (m1 * m2). 
   * @return the multiplication of two matrixes
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  template<typename Real> 
  inline SquareMatrix3<Real> operator *(const SquareMatrix3<Real>& m1,
                                        const SquareMatrix3<Real>& m2) {
    SquareMatrix3<Real> result;

    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
	for (unsigned int k = 0; k < 3; k++)
	  result(i, j)  += m1(i, k) * m2(k, j);                

    return result;
  }

  template<typename Real> 
  inline SquareMatrix3<Real> outProduct(const Vector3<Real>& v1,
                                        const Vector3<Real>& v2) {
    SquareMatrix3<Real> result;

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
	result(i, j)  = v1[i] * v2[j];                
      }
    }
            
    return result;        
  }

    
  typedef SquareMatrix3<RealType> Mat3x3d;
  typedef SquareMatrix3<RealType> RotMat3x3d;

  const Mat3x3d M3Zero(0.0);


} //namespace OpenMD
#endif // MATH_SQUAREMATRIX3_HPP

