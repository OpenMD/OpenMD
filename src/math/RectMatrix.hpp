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
 * @file RectMatrix.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */

#ifndef MATH_RECTMATRIX_HPP
#define MATH_RECTMATRIX_HPP

#include <cmath>
#include "Vector.hpp"

namespace oopse {

    /**
     * @class RectMatrix RectMatrix.hpp "math/RectMatrix.hpp"
     * @brief rectangular matrix class
     */
    template<typename Real, unsigned int Row, unsigned int Col> 
    class RectMatrix {
        public:

        /** default constructor */
        RectMatrix() {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)
                    data_[i][j] = 0.0;
         }

        /** Constructs and initializes every element of this matrix to a scalar */ 
        RectMatrix(Real s) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)
                    data_[i][j] = s;
        }

        /** copy constructor */
        RectMatrix(const RectMatrix<Real, Row, Col>& m) {
            *this = m;
        }
        
        /** destructor*/
        ~RectMatrix() {}

        /** copy assignment operator */
        RectMatrix<Real, Row, Col>& operator =(const RectMatrix<Real, Row, Col>& m) {
            if (this == &m)
                return *this;
            
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)
                    data_[i][j] = m.data_[i][j];
            return *this;
        }
        
        /**
         * Return the reference of a single element of this matrix.
         * @return the reference of a single element of this matrix 
         * @param i row index
         * @param j colum index
         */
        double& operator()(unsigned int i, unsigned int j) {
            //assert( i < Row && j < Col);
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
        Vector<Real, Row> getRow(unsigned int row) {
            Vector<Real, Row> v;

            for (unsigned int i = 0; i < Row; i++)
                v[i] = data_[row][i];

            return v;
        }

        /**
         * Sets a row of  this matrix
         * @param row the row index
         * @param v the vector to be set
         */                
         void setRow(unsigned int row, const Vector<Real, Row>& v) {

            for (unsigned int i = 0; i < Row; i++)
                data_[row][i] = v[i];
         }

        /**
         * Returns a column of  this matrix as a vector.
         * @return a column of  this matrix as a vector 
         * @param col the column index
         */                
        Vector<Real, Col> getColum(unsigned int col) {
            Vector<Real, Col> v;

            for (unsigned int j = 0; j < Col; j++)
                v[j] = data_[j][col];

            return v;
        }

        /**
         * Sets a column of  this matrix
         * @param col the column index
         * @param v the vector to be set
         */                
         void setColum(unsigned int col, const Vector<Real, Col>& v){

            for (unsigned int j = 0; j < Col; j++)
                data_[j][col] = v[j];
         }         

        /**
         * swap two rows of this matrix
         * @param i the first row
         * @param j the second row
         */
        void swapRow(unsigned int i, unsigned int j){
                assert(i < Row && j < Row);

                for (unsigned int k = 0; k < Col; k++)
                    std::swap(data_[i][k], data_[j][k]);
        }

       /**
         * swap two colums of this matrix
         * @param i the first colum
         * @param j the second colum
         */
        void swapColum(unsigned int i, unsigned int j){
                assert(i < Col && j < Col);
                
                for (unsigned int k = 0; k < Row; k++)
                    std::swap(data_[k][i], data_[k][j]);
        }

        /**
         * Tests if this matrix is identical to matrix m
         * @return true if this matrix is equal to the matrix m, return false otherwise
         * @m matrix to be compared
         *
         * @todo replace operator == by template function equal
         */
        bool operator ==(const RectMatrix<Real, Row, Col>& m) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)
                    if (!equal(data_[i][j], m.data_[i][j]))
                        return false;

            return true;
        }

        /**
         * Tests if this matrix is not equal to matrix m
         * @return true if this matrix is not equal to the matrix m, return false otherwise
         * @m matrix to be compared
         */
        bool operator !=(const RectMatrix<Real, Row, Col>& m) {
            return !(*this == m);
        }

        /** Negates the value of this matrix in place. */           
        inline void negate() {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)
                    data_[i][j] = -data_[i][j];
        }
        
        /**
        * Sets the value of this matrix to the negation of matrix m.
        * @param m the source matrix
        */
        inline void negate(const RectMatrix<Real, Row, Col>& m) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)
                    data_[i][j] = -m.data_[i][j];        
        }
        
        /**
        * Sets the value of this matrix to the sum of itself and m (*this += m).
        * @param m the other matrix
        */
        inline void add( const RectMatrix<Real, Row, Col>& m ) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)        
                data_[i][j] += m.data_[i][j];
        }
        
        /**
        * Sets the value of this matrix to the sum of m1 and m2 (*this = m1 + m2).
        * @param m1 the first matrix
        * @param m2 the second matrix
        */
        inline void add( const RectMatrix<Real, Row, Col>& m1, const RectMatrix<Real, Row, Col>& m2 ) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)        
                data_[i][j] = m1.data_[i][j] + m2.data_[i][j];
        }
        
        /**
        * Sets the value of this matrix to the difference  of itself and m (*this -= m).
        * @param m the other matrix
        */
        inline void sub( const RectMatrix<Real, Row, Col>& m ) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)        
                    data_[i][j] -= m.data_[i][j];
        }
        
        /**
        * Sets the value of this matrix to the difference of matrix m1 and m2 (*this = m1 - m2).
        * @param m1 the first matrix
        * @param m2 the second matrix
        */
        inline void sub( const RectMatrix<Real, Row, Col>& m1, const RectMatrix<Real, Row, Col>& m2){
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)        
                    data_[i][j] = m1.data_[i][j] - m2.data_[i][j];
        }

        /**
        * Sets the value of this matrix to the scalar multiplication of itself (*this *= s).
        * @param s the scalar value
        */
        inline void mul( double s ) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)  
                    data_[i][j] *= s;
        }

        /**
        * Sets the value of this matrix to the scalar multiplication of matrix m  (*this = s * m).
        * @param s the scalar value
        * @param m the matrix
        */
        inline void mul( double s, const RectMatrix<Real, Row, Col>& m ) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)  
                    data_[i][j] = s * m.data_[i][j];
        }

        /**
        * Sets the value of this matrix to the scalar division of itself  (*this /= s ).
        * @param s the scalar value
        */             
        inline void div( double s) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)  
                    data_[i][j] /= s;
        }

        /**
        * Sets the value of this matrix to the scalar division of matrix m  (*this = m /s).
        * @param s the scalar value
        * @param m the matrix
        */
        inline void div( double s, const RectMatrix<Real, Row, Col>& m ) {
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)  
                    data_[i][j] = m.data_[i][j] / s;
        }

        /**
         *  Multiples a scalar into every element of this matrix.
         * @param s the scalar value
         */
        RectMatrix<Real, Row, Col>& operator *=(const double s) {
            this->mul(s);
            return *this;
        }

        /**
         *  Divides every element of this matrix by a scalar.
         * @param s the scalar value
         */
        RectMatrix<Real, Row, Col>& operator /=(const double s) {
            this->div(s);
            return *this;
        }

        /**
         * Sets the value of this matrix to the sum of the other matrix and itself (*this += m).
         * @param m the other matrix
         */
        RectMatrix<Real, Row, Col>& operator += (const RectMatrix<Real, Row, Col>& m) {
            add(m);
            return *this;
         }

        /**
         * Sets the value of this matrix to the differerence of itself and the other matrix (*this -= m) 
         * @param m the other matrix
         */
        RectMatrix<Real, Row, Col>& operator -= (const RectMatrix<Real, Row, Col>& m){
            sub(m);
            return *this;
        }

        /** Return the transpose of this matrix */
        RectMatrix<Real,  Col, Row> transpose(){
            RectMatrix<Real,  Col, Row> result;
            
            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)              
                    result(j, i) = data_[i][j];

            return result;
        }
        
        protected:
            Real data_[Row][Col];
    };

    /** Negate the value of every element of this matrix. */
    template<typename Real, unsigned int Row, unsigned int Col> 
    inline RectMatrix<Real, Row, Col> operator -(const RectMatrix<Real, Row, Col>& m) {
        RectMatrix<Real, Row, Col> result(m);

        result.negate();

        return result;
    }
    
    /**
    * Return the sum of two matrixes  (m1 + m2). 
    * @return the sum of two matrixes
    * @param m1 the first matrix
    * @param m2 the second matrix
    */ 
    template<typename Real, unsigned int Row, unsigned int Col> 
    inline RectMatrix<Real, Row, Col> operator + (const RectMatrix<Real, Row, Col>& m1,const RectMatrix<Real, Row, Col>& m2) {
        RectMatrix<Real, Row, Col> result;

        result.add(m1, m2);

        return result;
    }
    
    /**
    * Return the difference of two matrixes  (m1 - m2). 
    * @return the sum of two matrixes
    * @param m1 the first matrix
    * @param m2 the second matrix
    */
    template<typename Real, unsigned int Row, unsigned int Col> 
    inline RectMatrix<Real, Row, Col> operator - (const RectMatrix<Real, Row, Col>& m1, const RectMatrix<Real, Row, Col>& m2) {
        RectMatrix<Real, Row, Col> result;

        result.sub(m1, m2);

        return result;
    }

    /**
    * Return the multiplication of scalra and  matrix  (m * s). 
    * @return the multiplication of a scalra and  a matrix 
    * @param m the matrix
    * @param s the scalar
    */
    template<typename Real, unsigned int Row, unsigned int Col> 
    inline RectMatrix<Real, Row, Col> operator *(const RectMatrix<Real, Row, Col>& m, Real s) {
        RectMatrix<Real, Row, Col> result;

        result.mul(s, m);

        return result;
    }

    /**
    * Return the multiplication of a scalra and  a matrix  (s * m). 
    * @return the multiplication of a scalra and  a matrix 
    * @param s the scalar
    * @param m the matrix
    */
    template<typename Real, unsigned int Row, unsigned int Col> 
    inline RectMatrix<Real, Row, Col> operator *(Real s, const RectMatrix<Real, Row, Col>& m) {
        RectMatrix<Real, Row, Col> result;

        result.mul(s, m);

        return result;
    }
    
    /**
    * Return the multiplication of two matrixes  (m1 * m2). 
    * @return the multiplication of two matrixes
    * @param m1 the first matrix
    * @param m2 the second matrix
    */
    template<typename Real, unsigned int Row, unsigned int Col, unsigned int SameDim> 
    inline RectMatrix<Real, Row, Col> operator *(const RectMatrix<Real, Row, SameDim>& m1, const RectMatrix<Real, SameDim, Col>& m2) {
        RectMatrix<Real, Row, Col> result;

            for (unsigned int i = 0; i < Row; i++)
                for (unsigned int j = 0; j < Col; j++)
                    for (unsigned int k = 0; k < SameDim; k++)
                        result(i, j)  += m1(i, k) * m2(k, j);                

        return result;
    }
    
    /**
    * Return the multiplication of  a matrix and a vector  (m * v). 
    * @return the multiplication of a matrix and a vector
    * @param m the matrix
    * @param v the vector
    */
    template<typename Real, unsigned int Row, unsigned int Col>
    inline Vector<Real, Row> operator *(const RectMatrix<Real, Row, Col>& m, const Vector<Real, Col>& v) {
        Vector<Real, Row> result;

        for (unsigned int i = 0; i < Row ; i++)
            for (unsigned int j = 0; j < Col ; j++)            
                result[i] += m(i, j) * v[j];
            
        return result;                                                                 
    }

    /**
    * Return the scalar division of matrix   (m / s). 
    * @return the scalar division of matrix  
    * @param m the matrix
    * @param s the scalar
    */
    template<typename Real, unsigned int Row, unsigned int Col> 
    inline RectMatrix<Real, Row, Col> operator /(const RectMatrix<Real, Row, Col>& m, Real s) {
        RectMatrix<Real, Row, Col> result;

        result.div(s, m);

        return result;
    }    

    /**
     * Write to an output stream
     */
    template<typename Real,  unsigned int Row, unsigned int Col>
    std::ostream &operator<< ( std::ostream& o, const RectMatrix<Real, Row, Col>& m) {
        for (unsigned int i = 0; i < Row ; i++) {
            o << "(";
            for (unsigned int j = 0; j < Col ; j++) {
                o << m(i, j) << "\t";
            }
            o << ")" << std::endl;
        }
        return o;        
    }    
}
#endif //MATH_RECTMATRIX_HPP
