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
 * @file Matrix3x3d.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */
#ifndef MATH_MATRIX3X3D_HPP 
#define MATH_MATRIX3X3D_HPP 

#include "math/Vector3.hpp"

namespace oopse {

    /**
     * @class Matrix3x3d Matrix3x3d.hpp "math/Matrix3x3d.hpp"
     * @brief A 3x3 matrix class
     */
    class Matrix3x3d{
        public:

        /** default constructor */
        Matrix3x3d();

        /** Constructs and initializes every element of this matrix to a scalar */ 
        Matrix3x3d(double s);
        
        /**
         * Return the reference of a single element of this matrix.
         * @return the reference of a single element of this matrix 
         * @param i row index
         * @param j colum index
         */
        double& operator()(unsigned int i, unsigned int j) {  return data_[i][j];  }

        /**
         * Return the value of a single element of this matrix.
         * @return the value of a single element of this matrix 
         * @param i row index
         * @param j colum index
         */        
        double operator()(unsigned int i, unsigned int j) const  {  return data_[i][j];  }

        /**
         * Returns a row of  this matrix as a vector.
         * @return a row of  this matrix as a vector 
         * @param i the row index
         */                
        Vector3d getRow(unsigned int i);

        /**
         * Sets a row of  this matrix
         * @param i the row index
         * @param v the vector to be set
         */                
         void setRow(unsigned int i, const Vector3d& v);

        /**
         * Returns a column of  this matrix as a vector.
         * @return a column of  this matrix as a vector 
         * @param i the column index
         */                
        Vector3d getColum(unsigned int i);

        /**
         * Sets a column of  this matrix
         * @param i the column index
         * @param v the vector to be set
         */                
         void setColum(unsigned int i, const Vector3d& v);

        /**
         *  Adds a scalar into every element of this matrix.
         * @param s the scalar value
         */
        Matrix3x3d& operator +=(const double s);

        /**
         *  Subtracts a scalar from every element of this matrix.
         * @param s the scalar value
         */
        Matrix3x3d& operator -=(const double s);

        /**
         *  Multiples a scalar into every element of this matrix.
         * @param s the scalar value
         */
        Matrix3x3d& operator *=(const double s);

        Matrix3x3d& operator /=(const double s);

        /**
         * Sets the value of this matrix to the sum of the other matrix and itself (*this += m).
         * @param m the other matrix
         */
        Matrix3x3d& operator += (const Matrix3x3d& m);

        /**
         * Sets the value of this matrix to the differerence of itself and the other matrix (*this -= m) 
         * @param m the other matrix
         */
        Matrix3x3d& operator -= (const Matrix3x3d& m);

        /** Returns the inverse of this matrix. */
        Matrix3x3d inverse();

        /** Returns the transpose of this matrix. */
        Matrix3x3d transpose();

        /** Returns the determinant of this matrix. */
        double determinant() const;

        /** Returns the trace of this matrix. */
        double trace() const;

        /** Tests if this matrix is symmetrix. */            
        bool isSymmetric() const;

        /** Tests if this matrix is orthogona. */            
        bool isOrthogonal() const;

        /** Tests if this matrix is diagonal. */
        bool isDiagonal() const;

        /** Tests if this matrix is the unit matrix. */
        bool isUnitMatrix() const;
        
        /** Negate the value of every element of this matrix. */
        friend inline Matrix3x3d operator -(const Matrix3x3d& m);
        
        /**
        * Return the sum of two matrixes  (m1 + m2). 
        * @return the sum of two matrixes
        * @param m1 the first matrix
        * @param m2 the second matrix
        */ 
        friend inline Matrix3x3d operator + (const Matrix3x3d& m1, const Matrix3x3d& m2);

        /**
        * Return the difference of two matrixes  (m1 - m2). 
        * @return the sum of two matrixes
        * @param m1 the first matrix
        * @param m2 the second matrix
        */
        friend inline Matrix3x3d operator - (const Matrix3x3d& m1, const Matrix3x3d& m2);

        /**
        * Return the multiplication of two matrixes  (m1 * m2). 
        * @return the multiplication of two matrixes
        * @param m1 the first matrix
        * @param m2 the second matrix
        */
        friend inline Matrix3x3d operator * (const Matrix3x3d& m1, const Matrix3x3d& m2);

        /**
        * Return the multiplication of  matrixes and vector (m1 * v1). 
        * @return the multiplication of matrixes and vector
        * @param m1 the matrix
        * @param v1 the vector
        */
        friend inline Matrix3x3d operator * (const Matrix3x3d& m1, const Matrix3x3d& v1);
        
        protected:
            double data_[3][3]; /**< matrix element */            

    }
}
#endif //MATH_MATRIX3X3D_HPP 
