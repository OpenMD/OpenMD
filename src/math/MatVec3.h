 /*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#ifndef MATH_MATVEC3_HPP
#define MATH_MATVEC3_HPP

#include <math.h>
#include "utils/simError.h"
#ifdef __cplusplus
extern "C" {
#endif

  void   identityMat3(double A[3][3]);
  void   swapVectors3(double v1[3], double v2[3]);
  static double norm3(const double x[3]){ return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); };
  double normalize3(double x[3]);
  void   matMul3(double a[3][3], double b[3][3], double out[3][3]);
  void   matVecMul3(double m[3][3], double inVec[3], double outVec[3]);
  double matDet3(double m[3][3]);
  void   invertMat3(double in[3][3], double out[3][3]);
  void   transposeMat3(double in[3][3], double out[3][3]);
  void   printMat3(double A[3][3]);
  void   printMat9(double A[9]);
  double matTrace3(double m[3][3]);
  void   crossProduct3(double a[3],double b[3], double out[3]);
  double dotProduct3(double a[3], double b[3]);
  void   diagonalize3x3(const double A[3][3],double w[3],double V[3][3]);
  int    JacobiN(double **a, int n, double *w, double **v);

#ifdef __cplusplus
}
#endif

#endif
