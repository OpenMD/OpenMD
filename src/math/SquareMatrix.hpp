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

        /** 
         * Retunrs  the inversion of this matrix. 
         * @todo
         */
         SquareMatrix<Real, Dim>  inverse() {
             SquareMatrix<Real, Dim> result;

             return result;
        }        

        /**
         * Returns the determinant of this matrix.
         * @todo
         */
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

        /** Tests if this matrix is orthogonal. */            
        bool isOrthogonal() {
            SquareMatrix<Real, Dim> tmp;

            tmp = *this * transpose();

            return tmp.isDiagonal();
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

        /** @todo need implement */
        void diagonalize() {
            //jacobi(m, eigenValues, ortMat);
        }

        /**
         * Finds the eigenvalues and eigenvectors of a symmetric matrix
         * @param eigenvals a reference to a vector3 where the
         * eigenvalues will be stored. The eigenvalues are ordered so
         * that eigenvals[0] <= eigenvals[1] <= eigenvals[2].
         * @return an orthogonal matrix whose ith column is an
         * eigenvector for the eigenvalue eigenvals[i]
         */
        SquareMatrix<Real, Dim>  findEigenvectors(Vector<Real, Dim>& eigenValues) {
            SquareMatrix<Real, Dim> ortMat;
            
            if ( !isSymmetric()){
                throw();
            }
            
            SquareMatrix<Real, Dim> m(*this);
            jacobi(m, eigenValues, ortMat);

            return ortMat;
        }
        /**
         * Jacobi iteration routines for computing eigenvalues/eigenvectors of 
         * real symmetric matrix
         *
         * @return true if success, otherwise return false
         * @param a source matrix
         * @param w output eigenvalues 
         * @param v output eigenvectors 
         */
        bool jacobi(const SquareMatrix<Real, Dim>& a, Vector<Real, Dim>& w, 
                              SquareMatrix<Real, Dim>& v);
    };//end SquareMatrix


#define ROT(a,i,j,k,l) g=a(i, j);h=a(k, l);a(i, j)=g-s*(h+g*tau);a(k, l)=h+s*(g-h*tau)
#define MAX_ROTATIONS 60

template<typename Real, int Dim>
bool SquareMatrix<Real, Dim>::jacobi(const SquareMatrix<Real, Dim>& a, Vector<Real, Dim>& w, 
                              SquareMatrix<Real, Dim>& v) {
    const int N = Dim;                                                                       
    int i, j, k, iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c;
    double tmp;
    Vector<Real, Dim> b, z;

    // initialize
    for (ip=0; ip<N; ip++) {
        for (iq=0; iq<N; iq++)
            v(ip, iq) = 0.0;
        v(ip, ip) = 1.0;
    }
    
    for (ip=0; ip<N; ip++) {
        b(ip) = w(ip) = a(ip, ip);
        z(ip) = 0.0;
    }

    // begin rotation sequence
    for (i=0; i<MAX_ROTATIONS; i++) {
        sm = 0.0;
        for (ip=0; ip<2; ip++) {
            for (iq=ip+1; iq<N; iq++)
                sm += fabs(a(ip, iq));
        }
        
        if (sm == 0.0)
            break;

        if (i < 4)
            tresh = 0.2*sm/(9);
        else
            tresh = 0.0;

        for (ip=0; ip<2; ip++) {
            for (iq=ip+1; iq<N; iq++) {
                g = 100.0*fabs(a(ip, iq));
                if (i > 4 && (fabs(w(ip))+g) == fabs(w(ip))
                    && (fabs(w(iq))+g) == fabs(w(iq))) {
                    a(ip, iq) = 0.0;
                } else if (fabs(a(ip, iq)) > tresh) {
                    h = w(iq) - w(ip);
                    if ( (fabs(h)+g) == fabs(h)) {
                        t = (a(ip, iq)) / h;
                    } else {
                        theta = 0.5*h / (a(ip, iq));
                        t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));

                        if (theta < 0.0)
                            t = -t;
                    }

                    c = 1.0 / sqrt(1+t*t);
                    s = t*c;
                    tau = s/(1.0+c);
                    h = t*a(ip, iq);
                    z(ip) -= h;
                    z(iq) += h;
                    w(ip) -= h;
                    w(iq) += h;
                    a(ip, iq)=0.0;
                    
                    for (j=0;j<ip-1;j++) 
                        ROT(a,j,ip,j,iq);

                    for (j=ip+1;j<iq-1;j++) 
                        ROT(a,ip,j,j,iq);

                    for (j=iq+1; j<N; j++) 
                        ROT(a,ip,j,iq,j);
                    
                    for (j=0; j<N; j++) 
                        ROT(v,j,ip,j,iq);
                }
            }
        }//for (ip=0; ip<2; ip++) 

        for (ip=0; ip<N; ip++) {
            b(ip) += z(ip);
            w(ip) = b(ip);
            z(ip) = 0.0;
        }
        
    } // end for (i=0; i<MAX_ROTATIONS; i++) 

    if ( i >= MAX_ROTATIONS )
        return false;

    // sort eigenfunctions
    for (j=0; j<N; j++) {
        k = j;
        tmp = w(k);
        for (i=j; i<N; i++) {
            if (w(i) >= tmp) {
            k = i;
            tmp = w(k);
            }
        }
    
        if (k != j) {
            w(k) = w(j);
            w(j) = tmp;
            for (i=0; i<N; i++)  {
                tmp = v(i, j);
                v(i, j) = v(i, k);
                v(i, k) = tmp;
            }
        }
    }

    //    insure eigenvector consistency (i.e., Jacobi can compute
    //    vectors that are negative of one another (.707,.707,0) and
    //    (-.707,-.707,0). This can reek havoc in
    //    hyperstreamline/other stuff. We will select the most
    //    positive eigenvector.
    int numPos;
    for (j=0; j<N; j++) {
        for (numPos=0, i=0; i<N; i++) if ( v(i, j) >= 0.0 ) numPos++;
        if ( numPos < 2 ) for(i=0; i<N; i++) v(i, j) *= -1.0;
    }

    return true;
}

#undef ROT
#undef MAX_ROTATIONS

}

#endif //MATH_SQUAREMATRIX_HPP 
