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

#include "math/RectMatrix.hpp"

namespace oopse {

    /**
     * @class SquareMatrix SquareMatrix.hpp "math/SquareMatrix.hpp"
     * @brief A square matrix class
     * @template Real the element type
     * @template Dim the dimension of the square matrix
     */
    template<typename Real, int Dim>
    class SquareMatrix : public RectMatrix<Real, Dim, Dim> {
        public:

        /** default constructor */
        SquareMatrix() {
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    data_[i][j] = 0.0;
         }

        /** copy constructor */
        SquareMatrix(const RectMatrix<Real, Dim, Dim>& m)  : RectMatrix<Real, Dim, Dim>(m) {
        }
        
        /** copy assignment operator */
        SquareMatrix<Real, Dim>& operator =(const RectMatrix<Real, Dim, Dim>& m) {
            RectMatrix<Real, Dim, Dim>::operator=(m);
            return *this;
        }
                               
        /** Retunrs  an identity matrix*/

       static SquareMatrix<Real, Dim> identity() {
            SquareMatrix<Real, Dim> m;
            
            for (unsigned int i = 0; i < Dim; i++) 
                for (unsigned int j = 0; j < Dim; j++) 
                    if (i == j)
                        m(i, j) = 1.0;
                    else
                        m(i, j) = 0.0;

            return m;
        }

        /** Retunrs  the inversion of this matrix. */
         SquareMatrix<Real, Dim>  inverse() {
             SquareMatrix<Real, Dim> result;

             return result;
        }

        

        /** Returns the determinant of this matrix. */
        double determinant() const {
            double det;
            return det;
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
                    if (fabs(data_[i][j] - data_[j][i]) > oopse::epsilon) 
                        return false;
                    
            return true;
        }

        /** Tests if this matrix is orthogona. */            
        bool isOrthogonal() {
            SquareMatrix<Real, Dim> tmp;

            tmp = *this * transpose();

            return tmp.isUnitMatrix();
        }

        /** Tests if this matrix is diagonal. */
        bool isDiagonal() const {
            for (unsigned int i = 0; i < Dim ; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    if (i !=j && fabs(data_[i][j]) > oopse::epsilon) 
                        return false;
                    
            return true;
        }

        /** Tests if this matrix is the unit matrix. */
        bool isUnitMatrix() const {
            if (!isDiagonal())
                return false;
            
            for (unsigned int i = 0; i < Dim ; i++)
                if (fabs(data_[i][i] - 1) > oopse::epsilon)
                    return false;
                
            return true;
        }         

    };//end SquareMatrix

}
#endif //MATH_SQUAREMATRIX_HPP 
