/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
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
#include "utils/NumericConstant.hpp"

namespace OpenMD {

  /**
   * @class SquareMatrix SquareMatrix.hpp "math/SquareMatrix.hpp"
   * @brief A square matrix class
   * @template Real the element type
   * @template Dim the dimension of the square matrix
   */
  template<typename Real, int Dim>
  class SquareMatrix : public RectMatrix<Real, Dim, Dim> {
  public:
    typedef Real ElemType;
    typedef Real* ElemPoinerType;

    /** default constructor */
    SquareMatrix() {
      for (unsigned int i = 0; i < Dim; i++)
	for (unsigned int j = 0; j < Dim; j++)
	  this->data_[i][j] = 0.0;
    }

    /** Constructs and initializes every element of this matrix to a scalar */ 
    SquareMatrix(Real s) : RectMatrix<Real, Dim, Dim>(s){
    }

    /** Constructs and initializes from an array */ 
    SquareMatrix(Real* array) : RectMatrix<Real, Dim, Dim>(array){
    }


    /** copy constructor */
    SquareMatrix(const RectMatrix<Real, Dim, Dim>& m) : RectMatrix<Real, Dim, Dim>(m) {
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
     * @todo need implementation
     */
    SquareMatrix<Real, Dim>  inverse() {
      SquareMatrix<Real, Dim> result;

      return result;
    }        

    /**
     * Returns the determinant of this matrix.
     * @todo need implementation
     */
    Real determinant() const {
      Real det;
      return det;
    }

    /** Returns the trace of this matrix. */
    Real trace() const {
      Real tmp = 0;
               
      for (unsigned int i = 0; i < Dim ; i++)
	tmp += this->data_[i][i];

      return tmp;
    }

    /** Tests if this matrix is symmetrix. */            
    bool isSymmetric() const {
      for (unsigned int i = 0; i < Dim - 1; i++)
	for (unsigned int j = i; j < Dim; j++)
	  if (fabs(this->data_[i][j] - this->data_[j][i]) > epsilon) 
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
	  if (i !=j && fabs(this->data_[i][j]) > epsilon) 
	    return false;
                        
      return true;
    }

    /** Tests if this matrix is the unit matrix. */
    bool isUnitMatrix() const {
      if (!isDiagonal())
	return false;
                
      for (unsigned int i = 0; i < Dim ; i++)
	if (fabs(this->data_[i][i] - 1) > epsilon)
	  return false;
                    
      return true;
    }         

    /** Return the transpose of this matrix */
    SquareMatrix<Real,  Dim> transpose() const{
      SquareMatrix<Real,  Dim> result;
                
      for (unsigned int i = 0; i < Dim; i++)
	for (unsigned int j = 0; j < Dim; j++)              
	  result(j, i) = this->data_[i][j];

      return result;
    }
            
    /** @todo need implementation */
    void diagonalize() {
      //jacobi(m, eigenValues, ortMat);
    }

    /**
     * Jacobi iteration routines for computing eigenvalues/eigenvectors of 
     * real symmetric matrix
     *
     * @return true if success, otherwise return false
     * @param a symmetric matrix whose eigenvectors are to be computed. On return, the matrix is
     *     overwritten
     * @param w will contain the eigenvalues of the matrix On return of this function
     * @param v the columns of this matrix will contain the eigenvectors. The eigenvectors are 
     *    normalized and mutually orthogonal. 
     */
           
    static int jacobi(SquareMatrix<Real, Dim>& a, Vector<Real, Dim>& d, 
		      SquareMatrix<Real, Dim>& v);
  };//end SquareMatrix


  /*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: SquareMatrix.hpp,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

#define VTK_ROTATE(a,i,j,k,l) g=a(i, j);h=a(k, l);a(i, j)=g-s*(h+g*tau); \
    a(k, l)=h+s*(g-h*tau)

#define VTK_MAX_ROTATIONS 20

  // Jacobi iteration for the solution of eigenvectors/eigenvalues of a nxn
  // real symmetric matrix. Square nxn matrix a; size of matrix in n;
  // output eigenvalues in w; and output eigenvectors in v. Resulting
  // eigenvalues/vectors are sorted in decreasing order; eigenvectors are
  // normalized.
  template<typename Real, int Dim>
  int SquareMatrix<Real, Dim>::jacobi(SquareMatrix<Real, Dim>& a, Vector<Real, Dim>& w, 
				      SquareMatrix<Real, Dim>& v) {
    const int n = Dim;  
    int i, j, k, iq, ip, numPos;
    Real tresh, theta, tau, t, sm, s, h, g, c, tmp;
    Real bspace[4], zspace[4];
    Real *b = bspace;
    Real *z = zspace;

    // only allocate memory if the matrix is large
    if (n > 4) {
      b = new Real[n];
      z = new Real[n]; 
    }

    // initialize
    for (ip=0; ip<n; ip++) {
      for (iq=0; iq<n; iq++) {
	v(ip, iq) = 0.0;
      }
      v(ip, ip) = 1.0;
    }
    for (ip=0; ip<n; ip++) {
      b[ip] = w[ip] = a(ip, ip);
      z[ip] = 0.0;
    }

    // begin rotation sequence
    for (i=0; i<VTK_MAX_ROTATIONS; i++) {
      sm = 0.0;
      for (ip=0; ip<n-1; ip++) {
	for (iq=ip+1; iq<n; iq++) {
	  sm += fabs(a(ip, iq));
	}
      }
      if (sm == 0.0) {
	break;
      }

      if (i < 3) {                                // first 3 sweeps
	tresh = 0.2*sm/(n*n);
      } else {
	tresh = 0.0;
      }

      for (ip=0; ip<n-1; ip++) {
	for (iq=ip+1; iq<n; iq++) {
	  g = 100.0*fabs(a(ip, iq));

	  // after 4 sweeps
	  if (i > 3 && (fabs(w[ip])+g) == fabs(w[ip])
	      && (fabs(w[iq])+g) == fabs(w[iq])) {
	    a(ip, iq) = 0.0;
	  } else if (fabs(a(ip, iq)) > tresh) {
	    h = w[iq] - w[ip];
	    if ( (fabs(h)+g) == fabs(h)) {
	      t = (a(ip, iq)) / h;
	    } else {
	      theta = 0.5*h / (a(ip, iq));
	      t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
	      if (theta < 0.0) {
		t = -t;
	      }
	    }
	    c = 1.0 / sqrt(1+t*t);
	    s = t*c;
	    tau = s/(1.0+c);
	    h = t*a(ip, iq);
	    z[ip] -= h;
	    z[iq] += h;
	    w[ip] -= h;
	    w[iq] += h;
	    a(ip, iq)=0.0;

	    // ip already shifted left by 1 unit
	    for (j = 0;j <= ip-1;j++) {
	      VTK_ROTATE(a,j,ip,j,iq);
	    }
	    // ip and iq already shifted left by 1 unit
	    for (j = ip+1;j <= iq-1;j++) {
	      VTK_ROTATE(a,ip,j,j,iq);
	    }
	    // iq already shifted left by 1 unit
	    for (j=iq+1; j<n; j++) {
	      VTK_ROTATE(a,ip,j,iq,j);
	    }
	    for (j=0; j<n; j++) {
	      VTK_ROTATE(v,j,ip,j,iq);
	    }
	  }
	}
      }

      for (ip=0; ip<n; ip++) {
	b[ip] += z[ip];
	w[ip] = b[ip];
	z[ip] = 0.0;
      }
    }

    //// this is NEVER called
    if ( i >= VTK_MAX_ROTATIONS ) {
      std::cout << "vtkMath::Jacobi: Error extracting eigenfunctions" << std::endl;
      return 0;
    }

    // sort eigenfunctions                 these changes do not affect accuracy 
    for (j=0; j<n-1; j++) {                  // boundary incorrect
      k = j;
      tmp = w[k];
      for (i=j+1; i<n; i++) {                // boundary incorrect, shifted already
	if (w[i] >= tmp) {                   // why exchage if same?
	  k = i;
	  tmp = w[k];
	}
      }
      if (k != j) {
	w[k] = w[j];
	w[j] = tmp;
	for (i=0; i<n; i++) {
	  tmp = v(i, j);
	  v(i, j) = v(i, k);
	  v(i, k) = tmp;
	}
      }
    }
    // insure eigenvector consistency (i.e., Jacobi can compute vectors that
    // are negative of one another (.707,.707,0) and (-.707,-.707,0). This can
    // reek havoc in hyperstreamline/other stuff. We will select the most
    // positive eigenvector.
    int ceil_half_n = (n >> 1) + (n & 1);
    for (j=0; j<n; j++) {
      for (numPos=0, i=0; i<n; i++) {
	if ( v(i, j) >= 0.0 ) {
	  numPos++;
	}
      }
      //    if ( numPos < ceil(RealType(n)/RealType(2.0)) )
      if ( numPos < ceil_half_n) {
	for (i=0; i<n; i++) {
	  v(i, j) *= -1.0;
	}
      }
    }

    if (n > 4) {
      delete [] b;
      delete [] z;
    }
    return 1;
  }


  typedef SquareMatrix<RealType, 6> Mat6x6d;
}
#endif //MATH_SQUAREMATRIX_HPP 

