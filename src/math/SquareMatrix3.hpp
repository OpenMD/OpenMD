/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

/**
 * @file SquareMatrix3.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */
 #ifndef MATH_SQUAREMATRIX3_HPP
#define  MATH_SQUAREMATRIX3_HPP

#include "Quaternion.hpp"
#include "SquareMatrix.hpp"
#include "Vector3.hpp"

namespace oopse {

    template<typename Real>
    class SquareMatrix3 : public SquareMatrix<Real, 3> {
        public:
            
            /** default constructor */
            SquareMatrix3() : SquareMatrix<Real, 3>() {
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
             * @psi theta
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

                data_[0][0] = cpsi * cphi - ctheta * sphi * spsi;
                data_[0][1] = cpsi * sphi + ctheta * cphi * spsi;
                data_[0][2] = spsi * stheta;
                
                data_[1][0] = -spsi * ctheta - ctheta * sphi * cpsi;
                data_[1][1] = -spsi * stheta + ctheta * cphi * cpsi;
                data_[1][2] = cpsi * stheta;

                data_[2][0] = stheta * sphi;
                data_[2][1] = -stheta * cphi;
                data_[2][2] = ctheta;
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

            /**
             * Returns the quaternion from this rotation matrix
             * @return the quaternion from this rotation matrix
             * @exception invalid rotation matrix
            */            
            Quaternion<Real> toQuaternion() {
                Quaternion<Real> q;
                Real t, s;
                Real ad1, ad2, ad3;    
                t = data_[0][0] + data_[1][1] + data_[2][2] + 1.0;

                if( t > 0.0 ){

                    s = 0.5 / sqrt( t );
                    q[0] = 0.25 / s;
                    q[1] = (data_[1][2] - data_[2][1]) * s;
                    q[2] = (data_[2][0] - data_[0][2]) * s;
                    q[3] = (data_[0][1] - data_[1][0]) * s;
                } else {

                    ad1 = fabs( data_[0][0] );
                    ad2 = fabs( data_[1][1] );
                    ad3 = fabs( data_[2][2] );

                    if( ad1 >= ad2 && ad1 >= ad3 ){

                        s = 2.0 * sqrt( 1.0 + data_[0][0] - data_[1][1] - data_[2][2] );
                        q[0] = (data_[1][2] + data_[2][1]) / s;
                        q[1] = 0.5 / s;
                        q[2] = (data_[0][1] + data_[1][0]) / s;
                        q[3] = (data_[0][2] + data_[2][0]) / s;
                    } else if ( ad2 >= ad1 && ad2 >= ad3 ) {
                        s = sqrt( 1.0 + data_[1][1] - data_[0][0] - data_[2][2] ) * 2.0;
                        q[0] = (data_[0][2] + data_[2][0]) / s;
                        q[1] = (data_[0][1] + data_[1][0]) / s;
                        q[2] = 0.5 / s;
                        q[3] = (data_[1][2] + data_[2][1]) / s;
                    } else {

                        s = sqrt( 1.0 + data_[2][2] - data_[0][0] - data_[1][1] ) * 2.0;
                        q[0] = (data_[0][1] + data_[1][0]) / s;
                        q[1] = (data_[0][2] + data_[2][0]) / s;
                        q[2] = (data_[1][2] + data_[2][1]) / s;
                        q[3] = 0.5 / s;
                    }
                }             

                return q;
                
            }

            /**
             * Returns the euler angles from this rotation matrix
             * @return the euler angles in a vector 
             * @exception invalid rotation matrix
             * We use so-called "x-convention", which is the most common definition. 
             * In this convention, the rotation given by Euler angles (phi, theta, psi), where the first 
             * rotation is by an angle phi about the z-axis, the second is by an angle  
             * theta (0 <= theta <= 180)about the x-axis, and thethird is by an angle psi about the
             * z-axis (again). 
            */            
            Vector3<Real> toEulerAngles() {
                Vector3<Real> myEuler;
                Real phi,theta,psi,eps;
                Real ctheta,stheta;
                
                // set the tolerance for Euler angles and rotation elements

                theta = acos(std::min(1.0, std::max(-1.0,data_[2][2])));
                ctheta = data_[2][2]; 
                stheta = sqrt(1.0 - ctheta * ctheta);

                // when sin(theta) is close to 0, we need to consider singularity
                // In this case, we can assign an arbitary value to phi (or psi), and then determine 
                // the psi (or phi) or vice-versa. We'll assume that phi always gets the rotation, and psi is 0
                // in cases of singularity.  
                // we use atan2 instead of atan, since atan2 will give us -Pi to Pi. 
                // Since 0 <= theta <= 180, sin(theta) will be always non-negative. Therefore, it never
                // change the sign of both of the parameters passed to atan2.

                if (fabs(stheta) <= oopse::epsilon){
                    psi = 0.0;
                    phi = atan2(-data_[1][0], data_[0][0]);  
                }
                // we only have one unique solution
                else{    
                    phi = atan2(data_[2][0], -data_[2][1]);
                    psi = atan2(data_[0][2], data_[1][2]);
                }

                //wrap phi and psi, make sure they are in the range from 0 to 2*Pi
                if (phi < 0)
                  phi += M_PI;

                if (psi < 0)
                  psi += M_PI;

                myEuler[0] = phi;
                myEuler[1] = theta;
                myEuler[2] = psi;

                return myEuler;
            }
            
            /** Returns the determinant of this matrix. */
            Real determinant() const {
                Real x,y,z;

                x = data_[0][0] * (data_[1][1] * data_[2][2] - data_[1][2] * data_[2][1]);
                y = data_[0][1] * (data_[1][2] * data_[2][0] - data_[1][0] * data_[2][2]);
                z = data_[0][2] * (data_[1][0] * data_[2][1] - data_[1][1] * data_[2][0]);

                return(x + y + z);
            }            
            
            /**
             * Sets the value of this matrix to  the inversion of itself. 
             * @note since simple algorithm can be applied to inverse the 3 by 3 matrix, we hide the 
             * implementation of inverse in SquareMatrix class
             */
            SquareMatrix3<Real>  inverse() {
                SquareMatrix3<Real> m;
                double det = determinant();
                if (fabs(det) <= oopse::epsilon) {
                //"The method was called on a matrix with |determinant| <= 1e-6.",
                //"This is a runtime or a programming error in your application.");
                }

                m(0, 0) = data_[1][1]*data_[2][2] - data_[1][2]*data_[2][1];
                m(1, 0) = data_[1][2]*data_[2][0] - data_[1][0]*data_[2][2];
                m(2, 0) = data_[1][0]*data_[2][1] - data_[1][1]*data_[2][0];
                m(0, 1) = data_[2][1]*data_[0][2] - data_[2][2]*data_[0][1];
                m(1, 1) = data_[2][2]*data_[0][0] - data_[2][0]*data_[0][2];
                m(2, 1) = data_[2][0]*data_[0][1] - data_[2][1]*data_[0][0];
                m(0, 2) = data_[0][1]*data_[1][2] - data_[0][2]*data_[1][1];
                m(1, 2) = data_[0][2]*data_[1][0] - data_[0][0]*data_[1][2];
                m(2, 2) = data_[0][0]*data_[1][1] - data_[0][1]*data_[1][0];

                m /= det;
                return m;
            }
            /**
             * Extract the eigenvalues and eigenvectors from a 3x3 matrix.
             * The eigenvectors (the columns of V) will be normalized. 
             * The eigenvectors are aligned optimally with the x, y, and z
             * axes respectively.
             * @param a symmetric matrix whose eigenvectors are to be computed. On return, the matrix is
             *     overwritten             
             * @param w will contain the eigenvalues of the matrix On return of this function
             * @param v the columns of this matrix will contain the eigenvectors. The eigenvectors are 
             *    normalized and mutually orthogonal.              
             * @warning a will be overwritten
             */
            static void diagonalize(SquareMatrix3<Real>& a, Vector3<Real>& w, SquareMatrix3<Real>& v); 
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
    void SquareMatrix3<Real>::diagonalize(SquareMatrix3<Real>& a, Vector3<Real>& w, 
                                                                           SquareMatrix3<Real>& v) {
        int i,j,k,maxI;
        Real tmp, maxVal;
        Vector3<Real> v_maxI, v_k, v_j;

        // diagonalize using Jacobi
        jacobi(a, w, v);
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
    typedef SquareMatrix3<double> Mat3x3d;
    typedef SquareMatrix3<double> RotMat3x3d;

} //namespace oopse
#endif // MATH_SQUAREMATRIX_HPP

