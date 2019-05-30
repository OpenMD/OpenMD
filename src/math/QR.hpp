#ifndef JAMA_QR_H
#define JAMA_QR_H

#include "math/DynamicRectMatrix.hpp"
#include <cmath>

using namespace OpenMD;
using namespace std;

namespace JAMA
{

  /** 
      <p>
      Classical QR Decompisition:
      for an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
      orthogonal matrix Q and an n-by-n upper triangular matrix R so that
      A = Q*R.
      <P>
      The QR decompostion always exists, even if the matrix does not have
      full rank, so the constructor will never fail.  The primary use of the
      QR decomposition is in the least squares solution of nonsquare systems
      of simultaneous linear equations.  This will fail if isFullRank()
      returns 0 (false).
      
      <p>
      The Q and R factors can be retrived via the getQ() and getR()
      methods. Furthermore, a solve() method is provided to find the
      least squares solution of Ax=b using the QR factors.  
      
      <p>
      (Adapted from JAMA, a Java Matrix Library, developed by jointly 
      by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
  */
  
  template <class Real>
  class QR {

    
    /** Array for internal storage of decomposition.
        @serial internal array storage.
    */
   
    DynamicRectMatrix<Real> QR_;
    
    /** Row and column dimensions.
        @serial column dimension.
        @serial row dimension.
    */
    int m, n;
    
    /** Array for internal storage of diagonal of R.
        @serial diagonal of R.
    */
    DynamicVector<Real> Rdiag;
    
    
  public:
    
    /**
       Create a QR factorization object for A.
       
       @param A rectangular (m>=n) matrix.
    */
    QR(const DynamicRectMatrix<Real> &A)		/* constructor */
    {
      QR_ = DynamicRectMatrix<Real>(A);
      m = A.getNRow();
      n = A.getNCol();
      Rdiag = DynamicVector<Real>(n);
      int i=0, j=0, k=0;

      // Main loop.
      for (k = 0; k < n; k++) {
        // Compute 2-norm of k-th column without under/overflow.
        Real nrm = 0;
        for (i = k; i < m; i++) {
          nrm = hypot(nrm,QR_(i,k));
        }

        if (nrm != 0.0) {
          // Form k-th Householder vector.
          if (QR_(k,k) < 0) {
            nrm = -nrm;
          }
          for (i = k; i < m; i++) {
            QR_(i,k) /= nrm;
          }
          QR_(k,k) += 1.0;

          // Apply transformation to remaining columns.
          for (j = k+1; j < n; j++) {
            Real s = 0.0; 
            for (i = k; i < m; i++) {
              s += QR_(i,k)*QR_(i,j);
            }
            s = -s/QR_(k,k);
            for (i = k; i < m; i++) {
              QR_(i,j) += s*QR_(i,k);
            }
          }
        }
        Rdiag(k) = -nrm;
      }
    }


    /**
       Flag to denote the matrix is of full rank.

       @return 1 if matrix is full rank, 0 otherwise.
    */
    int isFullRank() const		
    {
      for (int j = 0; j < n; j++) 
        {
          if (Rdiag(j) == 0)
            return 0;
        }
      return 1;
    }
	
	


    /** 
   
        Retreive the Householder vectors from QR factorization
        @returns lower trapezoidal matrix whose columns define the reflections
    */
    DynamicRectMatrix<Real> getHouseholder (void)  const
    {
      DynamicRectMatrix<Real> H(m,n);

      /* note: H is completely filled in by algorithm, so
         initializaiton of H is not necessary.
      */
      for (int i = 0; i < m; i++) 
        {
          for (int j = 0; j < n; j++) 
            {
              if (i >= j) {
                H(i,j) = QR_(i,j);
              } else {
                H(i,j) = 0.0;
              }
            }
        }
      return H;
    }



    /** Return the upper triangular factor, R, of the QR factorization
        @return     R
    */
    DynamicRectMatrix<Real> getR() const
    {
      DynamicRectMatrix<Real> R(n,n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (i < j) {
            R(i,j) = QR_(i,j);
          } else if (i == j) {
            R(i,j) = Rdiag(i);
          } else {
            R(i,j) = 0.0;
          }
        }
      }
      return R;
    }
	
	



    /** 
   	Generate and return the (economy-sized) orthogonal factor
        @return     Q the (economy-sized) orthogonal factor (Q*R=A).
    */
    DynamicRectMatrix<Real> getQ() const
    {
      int i=0, j=0, k=0;

      DynamicRectMatrix<Real> Q(m,n);
      for (k = n-1; k >= 0; k--) {
        for (i = 0; i < m; i++) {
          Q(i,k) = 0.0;
        }
        Q(k,k) = 1.0;
        for (j = k; j < n; j++) {
          if (QR_(k,k) != 0) {
            Real s = 0.0;
            for (i = k; i < m; i++) {
              s += QR_(i,k)*Q(i,j);
            }
            s = -s/QR_(k,k);
            for (i = k; i < m; i++) {
              Q(i,j) += s*QR_(i,k);
            }
          }
        }
      }
      return Q;
    }


    /** Least squares solution of A*x = b
        @param b     m-length array (vector).
        @return x    n-length array (vector) that minimizes the two norm of Q*R*X-B.
        If B is non-conformant, or if QR.isFullRank() is false,
        the routine returns a null (0-length) vector.
    */
    DynamicVector<Real> solve(const DynamicVector<Real> &b) const
    {
      if (b.size() != (unsigned int)m )		/* arrays must be conformant */
        return DynamicVector<Real>();

      if ( !isFullRank() )		/* matrix is rank deficient */
        {
          return DynamicVector<Real>();
        }

      DynamicVector<Real> x(b);

      // Compute Y = transpose(Q)*b
      for (int k = 0; k < n; k++) 
        {
          Real s = 0.0; 
          for (int i = k; i < m; i++) 
            {
              s += QR_(i,k)*x(i);
            }
          s = -s/QR_(k,k);
          for (int i = k; i < m; i++) 
            {
              x(i) += s*QR_(i,k);
            }
        }
      // Solve R*X = Y;
      for (int k = n-1; k >= 0; k--) 
        {
          x(k) /= Rdiag(k);
          for (int i = 0; i < k; i++) {
            x(i) -= x(k)*QR_(i,k);
          }
        }


      /* return n x nx portion of X */
      DynamicVector<Real> x_(n);
      for (int i=0; i<n; i++)
        x_(i) = x(i);

      return x_;
    }

    /** Least squares solution of A*X = B
        @param B     m x k Array (must conform).
        @return X     n x k Array that minimizes the two norm of Q*R*X-B. If
        B is non-conformant, or if QR.isFullRank() is false,
        the routine returns a null (0x0) array.
    */
    DynamicRectMatrix<Real> solve(const DynamicRectMatrix<Real> &B) const
    {
      if (B.getNRow() != m)		/* arrays must be conformant */
        return DynamicRectMatrix<Real>(0,0);

      if ( !isFullRank() )		/* matrix is rank deficient */
        {
          return DynamicRectMatrix<Real>(0,0);
        }

      int nx = B.getNCol(); 
      DynamicRectMatrix<Real> X(B);
      int i=0, j=0, k=0;

      // Compute Y = transpose(Q)*B
      for (k = 0; k < n; k++) {
        for (j = 0; j < nx; j++) {
          Real s = 0.0; 
          for (i = k; i < m; i++) {
            s += QR_(i,k)*X(i,j);
          }
          s = -s/QR_(k,k);
          for (i = k; i < m; i++) {
            X(i,j) += s*QR_(i,k);
          }
        }
      }
      // Solve R*X = Y;
      for (k = n-1; k >= 0; k--) {
        for (j = 0; j < nx; j++) {
          X(k,j) /= Rdiag(k);
        }
        for (i = 0; i < k; i++) {
          for (j = 0; j < nx; j++) {
            X(i,j) -= X(k,j)*QR_(i,k);
          }
        }
      }


      /* return n x nx portion of X */
      DynamicRectMatrix<Real> X_(n,nx);
      for (i=0; i<n; i++)
        for (j=0; j<nx; j++)
          X_(i,j) = X(i,j);

      return X_;
    }


  };


}
// namespace JAMA

#endif
// JAMA_QR__H

