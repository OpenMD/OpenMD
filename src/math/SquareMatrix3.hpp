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
#ifndef MATH_SQUAREMATRIX#_HPP
#define  MATH_SQUAREMATRIX#_HPP

#include "SquareMatrix.hpp"
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
            void setupRotMat(const Vector3d& euler);

            /**
             * Sets this matrix to a rotation matrix by three euler angles
             * @param phi
             * @param theta
             * @psi theta
             */
            void setupRotMat(double phi, double theta, double psi);


            /**
             * Sets this matrix to a rotation matrix by quaternion
             * @param quat
            */
            void setupRotMat(const Vector4d& quat);

            /**
             * Sets this matrix to a rotation matrix by quaternion
             * @param q0
             * @param q1
             * @param q2
             * @parma q3
            */
            void setupRotMat(double q0, double q1, double q2, double q4);

            /**
             * Returns the quaternion from this rotation matrix
             * @return the quaternion from this rotation matrix
             * @exception invalid rotation matrix
            */            
            Quaternion rotMatToQuat();

            /**
             * Returns the euler angles from this rotation matrix
             * @return the quaternion from this rotation matrix
             * @exception invalid rotation matrix
            */            
            Vector3d rotMatToEuler();
            
            /**
             * Sets the value of this matrix to  the inversion of itself. 
             * @note since simple algorithm can be applied to inverse the 3 by 3 matrix, we hide the 
             * implementation of inverse in SquareMatrix class
             */
            void  inverse();

            void diagonalize();

    }

    };

}
#endif // MATH_SQUAREMATRIX#_HPP
