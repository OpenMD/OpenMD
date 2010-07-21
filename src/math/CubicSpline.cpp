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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
#include "math/CubicSpline.hpp"
#include "utils/simError.h"
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace OpenMD;
using namespace std;

CubicSpline::CubicSpline() : generated(false), isUniform(true) {}

void CubicSpline::addPoint(RealType xp, RealType yp) {
  data.push_back(make_pair(xp, yp));
}

void CubicSpline::addPoints(const vector<RealType>& xps, 
                            const vector<RealType>& yps) {
  
  if (xps.size() != yps.size()) {
    printf( painCave.errMsg,
            "CubicSpline::addPoints was passed vectors of different length!\n");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal = 1;
    simError();    
  }

  for (int i = 0; i < xps.size(); i++) 
    data.push_back(make_pair(xps[i], yps[i]));
}

void CubicSpline::generate() { 
  // Calculate coefficients defining a smooth cubic interpolatory spline.
  //
  // class values constructed:
  //   n   = number of data points.
  //   x   = vector of independent variable values 
  //   y   = vector of dependent variable values
  //   b   = vector of S'(x[i]) values.
  //   c   = vector of S"(x[i])/2 values.
  //   d   = vector of S'''(x[i]+)/6 values (i < n).
  // Local variables:   
  
  RealType fp1, fpn, h, p;
  
  // make sure the sizes match
  
  n = data.size();  
  x.resize(n);
  y.resize(n);
  b.resize(n);
  c.resize(n);
  d.resize(n);
  
  // make sure we are monotonically increasing in x:
  
  bool sorted = true;
  
  for (int i = 1; i < n; i++) {
    if ( (data[i].first - data[i-1].first ) <= 0.0 ) sorted = false;
  }
  
  // sort if necessary
  
  if (!sorted) sort(data.begin(), data.end());  
  
  // Copy spline data out to separate arrays:
  
  for (int i = 0; i < n; i++) {
    x[i] = data[i].first;
    y[i] = data[i].second;
  }
  
  // Calculate coefficients for the tridiagonal system: store
  // sub-diagonal in B, diagonal in D, difference quotient in C.  
  
  b[0] = data[1].first - data[0].first;
  c[0] = (data[1].second - data[0].second) / b[0];
  
  if (n == 2) {

    // Assume the derivatives at both endpoints are zero. Another
    // assumption could be made to have a linear interpolant between
    // the two points.  In that case, the b coefficients below would be
    // (data[1].second - data[0].second) / (data[1].first - data[0].first)
    // and the c and d coefficients would both be zero.
    b[0] = 0.0;
    c[0] = -3.0 * pow((data[1].second - data[0].second) /
                      (data[1].first-data[0].first), 2);
    d[0] = -2.0 * pow((data[1].second - data[0].second) / 
                      (data[1].first-data[0].first), 3);
    b[1] = b[0];
    c[1] = 0.0;
    d[1] = 0.0;
    dx = 1.0 / (data[1].first - data[0].first);
    isUniform = true;
    generated = true;
    return;
  }
  
  d[0] = 2.0 * b[0];
  
  for (int i = 1; i < n-1; i++) {
    b[i] = data[i+1].first - data[i].first;
    if ( fabs( b[i] - b[0] ) / b[0] > 1.0e-5) isUniform = false;
    c[i] = (data[i+1].second - data[i].second) / b[i];
    d[i] = 2.0 * (b[i] + b[i-1]);
  }
  
  d[n-1] = 2.0 * b[n-2];
  
  // Calculate estimates for the end slopes using polynomials
  // that interpolate the data nearest the end.
  
  fp1 = c[0] - b[0]*(c[1] - c[0])/(b[0] + b[1]);
  if (n > 3) fp1 = fp1 + b[0]*((b[0] + b[1]) * (c[2] - c[1]) / 
                               (b[1] + b[2]) - 
                               c[1] + c[0]) / (data[3].first - data[0].first);
  
  fpn = c[n-2] + b[n-2]*(c[n-2] - c[n-3])/(b[n-3] + b[n-2]);

  if (n > 3)  fpn = fpn + b[n-2] * 
    (c[n-2] - c[n-3] - (b[n-3] + b[n-2]) * 
     (c[n-3] - c[n-4])/(b[n-3] + b[n-4]))/(data[n-1].first - data[n-4].first);
  
  
  // Calculate the right hand side and store it in C.
  
  c[n-1] = 3.0 * (fpn - c[n-2]);
  for (int i = n-2; i > 0; i--) 
    c[i] = 3.0 * (c[i] - c[i-1]);  
  c[0] = 3.0 * (c[0] - fp1);
  
  // Solve the tridiagonal system.
  
  for (int k = 1; k < n; k++) {
    p = b[k-1] / d[k-1];
    d[k] = d[k] - p*b[k-1];
    c[k] = c[k] - p*c[k-1];
  }
  
  c[n-1] = c[n-1] / d[n-1];
  
  for (int k = n-2; k >= 0; k--) 
    c[k] = (c[k] - b[k] * c[k+1]) / d[k];
  
  // Calculate the coefficients defining the spline.
  
  for (int i = 0; i < n-1; i++) {
    h = data[i+1].first - data[i].first;
    d[i] = (c[i+1] - c[i]) / (3.0 * h);
    b[i] = (data[i+1].second - data[i].second)/h - h * (c[i] + h * d[i]);
  }
  
  b[n-1] = b[n-2] + h * (2.0 * c[n-2] + h * 3.0 * d[n-2]);
  
  if (isUniform) dx = 1.0 / (data[1].first - data[0].first); 
  
  generated = true;
  return;
}

RealType CubicSpline::getValueAt(RealType t) {
  // Evaluate the spline at t using coefficients 
  //
  // Input parameters
  //   t = point where spline is to be evaluated.
  // Output:
  //   value of spline at t.
  
  if (!generated) generate();
  RealType dt;
  
  if ( t < data[0].first || t > data[n-1].first ) {    
    sprintf( painCave.errMsg,
             "CubicSpline::getValueAt was passed a value outside the range of the spline!\n");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal = 1;
    simError();    
  }

  //  Find the interval ( x[j], x[j+1] ) that contains or is nearest
  //  to t.

  int j;

  if (isUniform) {    
    
    j = max(0, min(n-1, int((t - data[0].first) * dx)));   

  } else { 

    j = n-1;
    
    for (int i = 0; i < n; i++) {
      if ( t < data[i].first ) {
        j = i-1;
        break;
      }      
    }
  }
  
  //  Evaluate the cubic polynomial.
  
  dt = t - data[j].first;
  return data[j].second + dt*(b[j] + dt*(c[j] + dt*d[j]));
  
}


pair<RealType, RealType> CubicSpline::getValueAndDerivativeAt(RealType t) {
  // Evaluate the spline and first derivative at t using coefficients 
  //
  // Input parameters
  //   t = point where spline is to be evaluated.
  // Output:
  //   pair containing value of spline at t and first derivative at t

  if (!generated) generate();
  RealType dt;
  
  if ( t < data.front().first || t > data.back().first ) {    
    sprintf( painCave.errMsg,
             "CubicSpline::getValueAndDerivativeAt was passed a value outside the range of the spline!\n");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal = 1;
    simError();    
  }

  //  Find the interval ( x[j], x[j+1] ) that contains or is nearest
  //  to t.

  int j;

  if (isUniform) {    
    
    j = max(0, min(n-1, int((t - data[0].first) * dx)));   

  } else { 

    j = n-1;
    
    for (int i = 0; i < n; i++) {
      if ( t < data[i].first ) {
        j = i-1;
        break;
      }      
    }
  }
  
  //  Evaluate the cubic polynomial.
  
  dt = t - data[j].first;

  RealType yval = data[j].second + dt*(b[j] + dt*(c[j] + dt*d[j]));
  RealType dydx = b[j] + dt*(2.0 * c[j] + 3.0 * dt * d[j]);
  
  return make_pair(yval, dydx);
}
