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
 * @file Quaternion.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */

#ifndef MATH_QUATERNION_HPP
#define MATH_QUATERNION_HPP

#include "math/Vector.hpp"

namespace oopse{

    /**
     * @class Quaternion Quaternion.hpp "math/Quaternion.hpp"
     * Quaternion is a sort of a higher-level complex number.
     * It is defined as Q = w + x*i + y*j + z*k,
     * where w, x, y, and z are numbers of type T (e.g. double), and
     * i*i = -1; j*j = -1; k*k = -1;
     * i*j = k; j*k = i; k*i = j;
     */
    template<typename Real>
    class Quaternion : public Vector<Real, 4> {
        public:
            Quaternion();

            /** Constructs and initializes a Quaternion from w, x, y, z values */     
            Quaternion(Real w, Real x, Real y, Real z) {
                data_[0] = w;
                data_[1] = x;
                data_[2] = y;
                data_[3] = z;                
            }
            
            /**
             *  
             */
            Quaternion(const Vector<Real,4>& v) 
                : Vector<Real, 4>(v){
            }

            /** */
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
                return data_[0];
            }

            /**
             * Returns the reference of the first element of this quaternion.
             * @return the reference of the first element of this quaternion
             */
            Real& w() {
                return data_[0];    
            }

            /**
             * Returns the value of the first element of this quaternion.
             * @return the value of the first element of this quaternion
             */
            Real x() const {
                return data_[1];
            }

            /**
             * Returns the reference of the second element of this quaternion.
             * @return the reference of the second element of this quaternion
             */
            Real& x() {
                return data_[1];    
            }

            /**
             * Returns the value of the thirf element of this quaternion.
             * @return the value of the third element of this quaternion
             */
            Real y() const {
                return data_[2];
            }

            /**
             * Returns the reference of the third element of this quaternion.
             * @return the reference of the third element of this quaternion
             */           
            Real& y() {
                return data_[2];    
            }

            /**
             * Returns the value of the fourth element of this quaternion.
             * @return the value of the fourth element of this quaternion
             */
            Real z() const {
                return data_[3];
            }
            /**
             * Returns the reference of the fourth element of this quaternion.
             * @return the reference of the fourth element of this quaternion
             */
            Real& z() {
                return data_[3];    
            }

            /**
             * Returns the inverse of this quaternion
             * @return inverse
             * @note since quaternion is a complex number, the inverse of quaternion
             * q = w + xi + yj+ zk is inv_q = (w -xi - yj - zk)/(|q|^2)
             */
            Quaternion<Real> inverse(){
                Quaternion<Real> q;
                Real d = this->lengthSquared();
                
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

                Real a0( (z() - y()) * (q.y() - q.z()) );
                Real a1( (w() + x()) * (q.w() + q.x()) );
                Real a2( (w() - x()) * (q.y() + q.z()) );
                Real a3( (y() + z()) * (q.w() - q.x()) );
                Real b0( -(x() - z()) * (q.x() - q.y()) );
                Real b1( -(x() + z()) * (q.x() + q.y()) );
                Real b2( (w() + y()) * (q.w() - q.z()) );
                Real b3( (w() - y()) * (q.w() + q.z()) );

                data_[0] = a0 + 0.5*(b0 + b1 + b2 + b3),;
                data_[1] = a1 + 0.5*(b0 + b1 - b2 - b3);
                data_[2] = a2 + 0.5*(b0 - b1 + b2 - b3),
                data_[3] = a3 + 0.5*(b0 - b1 - b2 + b3) );
            }


            /** Set the value of this quaternion to the division of itself by another quaternion */
            void div(const Quaternion<Real>& q) {
                mul(q.inverse());
            }
            
            Quaternion<Real>& operator *=(const Quaternion<Real>& q) {
                mul(q);
                return *this;
            }
                        
            Quaternion<Real>& operator /=(const Quaternion<Real>& q) {
                mul(q.inverse());
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
            SquareMatrix<Real, 3, 3> toRotationMatrix3() {
                SquareMatrix<Real, 3, 3> rotMat3;

                Real w2;
                Real x2;
                Real y2;
                Real z2;

                if (!isNormalized())
                    normalize();
                
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
            }

    };//end Quaternion

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
    inline Quaternion<Real> operator /(const Quaternion<Real>& q1, const Quaternion<Real>& q2) {
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
    Quaternion<Real> operator /(const Quaternion<Real>& s, const Quaternion<Real>& q) {

        Quaternion<Real> x = q.inv();
        return x * s;
    }

    typedef Quaternion<double> Quat4d;
}
#endif //MATH_QUATERNION_HPP 
