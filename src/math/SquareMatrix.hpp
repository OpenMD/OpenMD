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
 * @file SquareMatrix.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */
#ifndef MATH_SQUAREMATRIX_HPP 
#define MATH_SQUAREMATRIX_HPP 

#include "Vector3d.hpp"

namespace oopse {

    /**
     * @class SquareMatrix SquareMatrix.hpp "math/SquareMatrix.hpp"
     * @brief A square matrix class
     * @template Real the element type
     * @template Dim the dimension of the square matrix
     */
    template<typename Real, int Dim>
    class SquareMatrix{
        public:

        /** default constructor */
        SquareMatrix() {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    data_[i][j] = 0.0;
         }

        /** Constructs and initializes every element of this matrix to a scalar */ 
        SquareMatrix(double s) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    data_[i][j] = s;
        }

        /** copy constructor */
        SquareMatrix(const SquareMatrix<Real, Dim>& m) {
            *this = m;
        }
        
        /** destructor*/
        ~SquareMatrix() {}

        /** copy assignment operator */
        SquareMatrix<Real, Dim>& operator =(const SquareMatrix<Real, Dim>& m) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    data_[i][j] = m.data_[i][j];
        }
        
        /**
         * Return the reference of a single element of this matrix.
         * @return the reference of a single element of this matrix 
         * @param i row index
         * @param j colum index
         */
        double& operator()(unsigned int i, unsigned int j) {
            return data_[i][j];
        }

        /**
         * Return the value of a single element of this matrix.
         * @return the value of a single element of this matrix 
         * @param i row index
         * @param j colum index
         */        
        double operator()(unsigned int i, unsigned int j) const  {
            return data_[i][j];  
        }

        /**
         * Returns a row of  this matrix as a vector.
         * @return a row of  this matrix as a vector 
         * @param row the row index
         */                
        Vector<Real, Dim> getRow(unsigned int row) {
            Vector<Real, Dim> v;

            for (unsigned int i = 0; i < Dim; i++)
                v[i] = data_[row][i];

            return v;
        }

        /**
         * Sets a row of  this matrix
         * @param row the row index
         * @param v the vector to be set
         */                
         void setRow(unsigned int row, const Vector<Real, Dim>& v) {
            Vector<Real, Dim> v;

            for (unsigned int i = 0; i < Dim; i++)
                data_[row][i] = v[i];
         }

        /**
         * Returns a column of  this matrix as a vector.
         * @return a column of  this matrix as a vector 
         * @param col the column index
         */                
        Vector<Real, Dim> getColum(unsigned int col) {
            Vector<Real, Dim> v;

            for (unsigned int i = 0; i < Dim; i++)
                v[i] = data_[i][col];

            return v;
        }

        /**
         * Sets a column of  this matrix
         * @param col the column index
         * @param v the vector to be set
         */                
         void setColum(unsigned int col, const Vector<Real, Dim>& v){
            Vector<Real, Dim> v;

            for (unsigned int i = 0; i < Dim; i++)
                data_[i][col] = v[i];
         }         

        /** Negates the value of this matrix in place. */           
        inline void negate() {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    data_[i][j] = -data_[i][j];
        }
        
        /**
        * Sets the value of this matrix to the negation of matrix m.
        * @param m the source matrix
        */
        inline void negate(const SquareMatrix<Real, Dim>& m) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    data_[i][j] = -m.data_[i][j];        
        }
        
        /**
        * Sets the value of this matrix to the sum of itself and m (*this += m).
        * @param m the other matrix
        */
        inline void add( const SquareMatrix<Real, Dim>& m ) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)        
                data_[i][j] += m.data_[i][j];
            }
        
        /**
        * Sets the value of this matrix to the sum of m1 and m2 (*this = m1 + m2).
        * @param m1 the first matrix
        * @param m2 the second matrix
        */
        inline void add( const SquareMatrix<Real, Dim>& m1, const SquareMatrix<Real, Dim>& m2 ) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)        
                data_[i][j] = m1.data_[i][j] + m2.data_[i][j];
        }
        
        /**
        * Sets the value of this matrix to the difference  of itself and m (*this -= m).
        * @param m the other matrix
        */
        inline void sub( const SquareMatrix<Real, Dim>& m ) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)        
                data_[i][j] -= m.data_[i][j];
        }
        
        /**
        * Sets the value of this matrix to the difference of matrix m1 and m2 (*this = m1 - m2).
        * @param m1 the first matrix
        * @param m2 the second matrix
        */
        inline void sub( const SquareMatrix<Real, Dim>& m1, const Vector  &m2){
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)        
                data_[i][j] = m1.data_[i][j] - m2.data_[i][j];
        }
        
        /**
        * Sets the value of this matrix to the scalar multiplication of itself (*this *= s).
        * @param s the scalar value
        */
        inline void mul( double s ) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)  
                    data_[i][j] *= s;
        }

        /**
        * Sets the value of this matrix to the scalar multiplication of matrix m  (*this = s * m).
        * @param s the scalar value
        * @param m the matrix
        */
        inline void mul( double s, const SquareMatrix<Real, Dim>& m ) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)  
                    data_[i][j] = s * m.data_[i][j];
        }

        /**
        * Sets the value of this matrix to the  multiplication of this matrix and matrix m
        * (*this = *this * m).
        * @param m the matrix
        */
        inline void mul(const SquareMatrix<Real, Dim>& m ) {
            SquareMatrix<Real, Dim> tmp(*this);
            
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++) {  
                    
                    data_[i][j] = 0.0;
                    for (unsigned int k = 0; k < Dim; k++)
                        data_[i][j]  = tmp.data_[i][k] * m.data_[k][j]
                }
        }
        
        /**
        * Sets the value of this matrix to the  left multiplication of matrix m into itself
        * (*this = m *  *this).
        * @param m the matrix
        */
        inline void leftmul(const SquareMatrix<Real, Dim>& m ) {
            SquareMatrix<Real, Dim> tmp(*this);
            
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++) {  
                    
                    data_[i][j] = 0.0;
                    for (unsigned int k = 0; k < Dim; k++)
                        data_[i][j]  = m.data_[i][k] * tmp.data_[k][j]
                }
        }

        /**
        * Sets the value of this matrix to the  multiplication of matrix m1 and matrix m2
        * (*this = m1 * m2).
        * @param m1 the first  matrix
        * @param m2 the second matrix
        */
        inline void mul(const SquareMatrix<Real, Dim>& m1, 
                                  const SquareMatrix<Real, Dim>& m2 ) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++) {  
                    
                    data_[i][j] = 0.0;
                    for (unsigned int k = 0; k < Dim; k++)
                        data_[i][j]  = m1.data_[i][k] * m2.data_[k][j]
                }

        }
        
        /**
        * Sets the value of this matrix to the scalar division of itself  (*this /= s ).
        * @param s the scalar value
        */             
        inline void div( double s) {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)  
                    data_[i][j] /= s;
        }
        
        inline SquareMatrix<Real, Dim>& operator=(const SquareMatrix<Real, Dim>& v) {
            if (this == &v)
                return *this;
            
            for (unsigned int i = 0; i < Dim; i++)            
                data_[i] = v[i];
            
            return *this;
        }
        
        /**
        * Sets the value of this matrix to the scalar division of matrix v1  (*this = v1 / s ).
        * @paran v1 the source matrix
        * @param s the scalar value
        */                         
        inline void div( const SquareMatrix<Real, Dim>& v1, double s ) {
            for (unsigned int i = 0; i < Dim; i++)
                data_[i] = v1.data_[i] / s;
        }

        /**
         *  Multiples a scalar into every element of this matrix.
         * @param s the scalar value
         */
        SquareMatrix<Real, Dim>& operator *=(const double s) {
            this->mul(s);
            return *this;
        }

        /**
         *  Divides every element of this matrix by a scalar.
         * @param s the scalar value
         */
        SquareMatrix<Real, Dim>& operator /=(const double s) {
            this->div(s);
            return *this;
        }

        /**
         * Sets the value of this matrix to the sum of the other matrix and itself (*this += m).
         * @param m the other matrix
         */
        SquareMatrix<Real, Dim>& operator += (const SquareMatrix<Real, Dim>& m) {
            add(m);
            return *this;
         }

        /**
         * Sets the value of this matrix to the differerence of itself and the other matrix (*this -= m) 
         * @param m the other matrix
         */
        SquareMatrix<Real, Dim>& operator -= (const SquareMatrix<Real, Dim>& m){
            sub(m);
            return *this;
        }

        /** set this matrix to an identity matrix*/

       void identity() {
            for (unsigned int i = 0; i < Dim; i++) 
                for (unsigned int i = 0; i < Dim; i++) 
                    if (i == j)
                        data_[i][j] = 1.0;
                    else
                        data_[i][j] = 0.0;
        }

        /** Sets the value of this matrix to  the inversion of itself. */
        void  inverse() {
            inverse(*this);
        }

        /**
         * Sets the value of this matrix to  the inversion of other matrix.
         * @ param m the source matrix
         */        
        void inverse(const SquareMatrix<Real, Dim>& m);
        
        /** Sets the value of this matrix to  the transpose of itself. */
        void transpose() {
            for (unsigned int i = 0; i < Dim - 1; i++)
                for (unsigned int j = i; j < Dim; j++)
                    std::swap(data_[i][j], data_[j][i]);
        }

        /**
         * Sets the value of this matrix to  the transpose of other matrix.
         * @ param m the source matrix
         */        
        void transpose(const SquareMatrix<Real, Dim>& m) {
            
            if (this == &m) {
                transpose();
            } else {
                for (unsigned int i = 0; i < Dim; i++)
                    for (unsigned int j =0; j < Dim; j++)
                        data_[i][j] = m.data_[i][j];
            }
        }

        /** Returns the determinant of this matrix. */
        double determinant() const {

        }

        /** Returns the trace of this matrix. */
        double trace() const {
           double tmp = 0;
           
            for (unsigned int i = 0; i < Dim ; i++)
                tmp += data_[i][i];

            return tmp;
        }

        /** Tests if this matrix is symmetrix. */            
        bool isSymmetric() const {
            for (unsigned int i = 0; i < Dim - 1; i++)
                for (unsigned int j = i; j < Dim; j++)
                    if (fabs(data_[i][j] - data_[j][i]) > epsilon) 
                        return false;
                    
            return true;
        }

        /** Tests if this matrix is orthogona. */            
        bool isOrthogonal() const {
            SquareMatrix<Real, Dim> t(*this);

            t.transpose();

            return isUnitMatrix(*this * t);
        }

        /** Tests if this matrix is diagonal. */
        bool isDiagonal() const {
            for (unsigned int i = 0; i < Dim ; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    if (i !=j && fabs(data_[i][j]) > epsilon) 
                        return false;
                    
            return true;
        }

        /** Tests if this matrix is the unit matrix. */
        bool isUnitMatrix() const {
            if (!isDiagonal())
                return false;
            
            for (unsigned int i = 0; i < Dim ; i++)
                if (fabs(data_[i][i] - 1) > epsilon)
                    return false;
                
            return true;
        }
        
        protected:
            double data_[Dim][Dim]; /**< matrix element */            

    };//end SquareMatrix

    
    /** Negate the value of every element of this matrix. */
    template<typename Real, int Dim>
    inline SquareMatrix<Real, Dim> operator -(const SquareMatrix& m) {
        SquareMatrix<Real, Dim> result(m);

        result.negate();

        return result;
    }
    
    /**
    * Return the sum of two matrixes  (m1 + m2). 
    * @return the sum of two matrixes
    * @param m1 the first matrix
    * @param m2 the second matrix
    */ 
    template<typename Real, int Dim>
    inline SquareMatrix<Real, Dim> operator + (const SquareMatrix<Real, Dim>& m1,
                                                                                         const SquareMatrix<Real, Dim>& m2) {
        SquareMatrix<Real, Dim>result;

        result.add(m1, m2);

        return result;
    }
    
    /**
    * Return the difference of two matrixes  (m1 - m2). 
    * @return the sum of two matrixes
    * @param m1 the first matrix
    * @param m2 the second matrix
    */
    template<typename Real, int Dim>
    inline SquareMatrix<Real, Dim> operator - (const SquareMatrix<Real, Dim>& m1, 
                                                                                        const SquareMatrix<Real, Dim>& m2) {
        SquareMatrix<Real, Dim>result;

        result.sub(m1, m2);

        return result;
    }
    
    /**
    * Return the multiplication of two matrixes  (m1 * m2). 
    * @return the multiplication of two matrixes
    * @param m1 the first matrix
    * @param m2 the second matrix
    */
    template<typename Real, int Dim>
    inline SquareMatrix<Real, Dim> operator *(const SquareMatrix<Real, Dim>& m1,
                                                                                       const SquareMatrix<Real, Dim>& m2) {
        SquareMatrix<Real, Dim> result;

        result.mul(m1, m2);

        return result;
    }
    
    /**
    * Return the multiplication of  matrixes m  and vector v (m * v). 
    * @return the multiplication of matrixes and vector
    * @param m the matrix
    * @param v the vector
    */
    template<typename Real, int Dim>
    inline Vector<Real, Dim> operator *(const SquareMatrix<Real, Dim>& m, 
                                                                 const SquareMatrix<Real, Dim>& v) {
        Vector<Real, Dim> result;

        for (unsigned int i = 0; i < Dim ; i++)
            for (unsigned int j = 0; j < Dim ; j++)            
                result[i] += m(i, j) * v[j];
            
        return result;                                                                 
    }
}
#endif //MATH_SQUAREMATRIX_HPP 
