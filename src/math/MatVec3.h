#ifndef __MATVEC3_H__
#define __MATVEC3_H__

#include <math.h>
#include "simError.h"
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
