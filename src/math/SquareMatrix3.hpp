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
#ifndef MATH_SQUAREMATRIX_HPP
#define  MATH_SQUAREMATRIX_HPP

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
                *this = q.toRotationMatrix3();
            }

            SquareMatrix3(Real w, Real x, Real y, Real z) {
                Quaternion<Real> q(w, x, y, z);
                *this = q.toRotationMatrix3();
            }
            
            /** copy assignment operator */
            SquareMatrix3<Real>& operator =(const SquareMatrix<Real, 3>& m) {
                if (this == &m)
                    return *this;
                 SquareMatrix<Real, 3>::operator=(m);
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
                *this = quat.toRotationMatrix3();
            }

            /**
             * Sets this matrix to a rotation matrix by quaternion
             * @param w the first element 
             * @param x the second element
             * @param y the third element
             * @parma z the fourth element
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
                Vector<Real> myEuler;
                Real phi,theta,psi,eps;
                Real ctheta,stheta;
                
                // set the tolerance for Euler angles and rotation elements

                theta = acos(min(1.0,max(-1.0,data_[2][2])));
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
            
            /**
             * Sets the value of this matrix to  the inversion of itself. 
             * @note since simple algorithm can be applied to inverse the 3 by 3 matrix, we hide the 
             * implementation of inverse in SquareMatrix class
             */
            void  inverse();

            void diagonalize();

    };

    typedef template SquareMatrix3<double> Mat3x3d
    typedef template SquareMatrix3<double> RotMat3x3d;

} //namespace oopse
#endif // MATH_SQUAREMATRIX_HPP
