/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */
 
#include "math/CubicSpline.hpp"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <numeric>

using namespace OpenMD;
using namespace std;

CubicSpline::CubicSpline() : isUniform(true), generated(false) {
  x_.clear();
  y_.clear();
}

void CubicSpline::addPoint(const RealType xp, const RealType yp) {
  x_.push_back(xp);
  y_.push_back(yp);
}

void CubicSpline::addPoints(const vector<RealType>& xps, 
                            const vector<RealType>& yps) {
  
  assert(xps.size() == yps.size());
  
  for (unsigned int i = 0; i < xps.size(); i++){
    x_.push_back(xps[i]);
    y_.push_back(yps[i]);
  }
}

void CubicSpline::generate() { 
  // Calculate coefficients defining a smooth cubic interpolatory spline.
  //
  // class values constructed:
  //   n   = number of data_ points.
  //   x_  = vector of independent variable values 
  //   y_  = vector of dependent variable values
  //   b   = vector of S'(x_[i]) values.
  //   c   = vector of S"(x_[i])/2 values.
  //   d   = vector of S'''(x_[i]+)/6 values (i < n).
  // Local variables:   
 
  RealType fp1, fpn, p;
  RealType h(0.0);
  
  // make sure the sizes match
  
  n = x_.size();  
  b.resize(n);
  c.resize(n);
  d.resize(n);
  
  // make sure we are monotonically increasing in x:
  
  bool sorted = true;
  
  for (int i = 1; i < n; i++) {
    if ( (x_[i] - x_[i-1] ) <= 0.0 ) sorted = false;
  }
  
  // sort if necessary
  
  if (!sorted) {
    vector<int> p = sort_permutation(x_);
    x_ = apply_permutation(x_, p);
    y_ = apply_permutation(y_, p);
  }
  
  
  // Calculate coefficients for the tridiagonal system: store
  // sub-diagonal in B, diagonal in D, difference quotient in C.  
  
  b[0] = x_[1] - x_[0];
  c[0] = (y_[1] - y_[0]) / b[0];
  
  if (n == 2) {

    // Assume the derivatives at both endpoints are zero. Another
    // assumption could be made to have a linear interpolant between
    // the two points.  In that case, the b coefficients below would be
    // (y_[1] - y_[0]) / (x_[1] - x_[0])
    // and the c and d coefficients would both be zero.
    b[0] = 0.0;
    c[0] = -3.0 * pow((y_[1] - y_[0]) / (x_[1] - x_[0]), 2);
    d[0] = -2.0 * pow((y_[1] - y_[0]) / (x_[1] - x_[0]), 3);
    b[1] = b[0];
    c[1] = 0.0;
    d[1] = 0.0;
    dx = 1.0 / (x_[1] - x_[0]);
    isUniform = true;
    generated = true;
    return;
  }
  
  d[0] = 2.0 * b[0];
  
  for (int i = 1; i < n-1; i++) {
    b[i] = x_[i+1] - x_[i];
    if ( fabs( b[i] - b[0] ) / b[0] > 1.0e-5) isUniform = false;
    c[i] = (y_[i+1] - y_[i]) / b[i];
    d[i] = 2.0 * (b[i] + b[i-1]);
  }
  
  d[n-1] = 2.0 * b[n-2];
  
  // Calculate estimates for the end slopes using polynomials
  // that interpolate the data_ nearest the end.
  
  fp1 = c[0] - b[0]*(c[1] - c[0])/(b[0] + b[1]);
  if (n > 3) fp1 = fp1 + b[0]*((b[0] + b[1]) * (c[2] - c[1]) / 
                               (b[1] + b[2]) - 
                               c[1] + c[0]) / (x_[3] - x_[0]);
  
  fpn = c[n-2] + b[n-2]*(c[n-2] - c[n-3])/(b[n-3] + b[n-2]);
  
  if (n > 3)  fpn = fpn + b[n-2] * 
                (c[n-2] - c[n-3] - (b[n-3] + b[n-2]) * 
                 (c[n-3] - c[n-4])/(b[n-3] + b[n-4])) /
                (x_[n-1] - x_[n-4]);
  
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
    h = x_[i+1] - x_[i];
    d[i] = (c[i+1] - c[i]) / (3.0 * h);
    b[i] = (y_[i+1] - y_[i])/h - h * (c[i] + h * d[i]);
  }
  
  b[n-1] = b[n-2] + h * (2.0 * c[n-2] + h * 3.0 * d[n-2]);
  
  if (isUniform) dx = 1.0 / (x_[1] - x_[0]); 
  
  generated = true;
  return;
}

RealType CubicSpline::getValueAt(const RealType& t) {
  // Evaluate the spline at t using coefficients 
  //
  // Input parameters
  //   t = point where spline is to be evaluated.
  // Output:
  //   value of spline at t.
  
  if (!generated) generate();
  
  assert(t >= x_.front());
  assert(t <= x_.back());

  //  Find the interval ( x[j], x[j+1] ) that contains or is nearest
  //  to t.

  if (isUniform) {    
    
    j = max(0, min(n-1, int((t - x_[0]) * dx)));   

  } else { 

    j = n-1;
    
    for (int i = 0; i < n; i++) {
      if ( t < x_[i] ) {
        j = i-1;
        break;
      }      
    }
  }
  
  //  Evaluate the cubic polynomial.
  
  dt = t - x_[j];
  return y_[j] + dt*(b[j] + dt*(c[j] + dt*d[j]));  
}


void CubicSpline::getValueAt(const RealType& t, RealType& v) {
  // Evaluate the spline at t using coefficients 
  //
  // Input parameters
  //   t = point where spline is to be evaluated.
  // Output:
  //   value of spline at t.
  
  if (!generated) generate();
  
  assert(t >= x_.front());
  assert(t <= x_.back());

  //  Find the interval ( x[j], x[j+1] ) that contains or is nearest
  //  to t.

  if (isUniform) {    
    
    j = max(0, min(n-1, int((t - x_[0]) * dx)));   

  } else { 

    j = n-1;
    
    for (int i = 0; i < n; i++) {
      if ( t < x_[i] ) {
        j = i-1;
        break;
      }      
    }
  }
  
  //  Evaluate the cubic polynomial.
  
  dt = t - x_[j];
  v = y_[j] + dt*(b[j] + dt*(c[j] + dt*d[j]));  
}

pair<RealType, RealType> CubicSpline::getLimits(){
  if (!generated) generate();
  return make_pair( x_.front(), x_.back() );
}

RealType CubicSpline::getSpacing(){
  if (!generated) generate();
  assert(isUniform);
  if (isUniform) return 1.0/dx;
  else return 0.0;
}

void CubicSpline::getValueAndDerivativeAt(const RealType& t, RealType& v, 
                                          RealType &dv) {
  // Evaluate the spline and first derivative at t using coefficients 
  //
  // Input parameters
  //   t = point where spline is to be evaluated.

  if (!generated) generate();
  
  assert(t >= x_.front());
  assert(t <= x_.back());

  //  Find the interval ( x[j], x[j+1] ) that contains or is nearest
  //  to t.

  if (isUniform) {    
    
    j = max(0, min(n-1, int((t - x_[0]) * dx)));   

  } else { 

    j = n-1;
    
    for (int i = 0; i < n; i++) {
      if ( t < x_[i] ) {
        j = i-1;
        break;
      }      
    }
  }
  
  //  Evaluate the cubic polynomial.
  
  dt = t - x_[j];

  v = y_[j] + dt*(b[j] + dt*(c[j] + dt*d[j]));
  dv = b[j] + dt*(2.0 * c[j] + 3.0 * dt * d[j]); 
}

std::vector<int> CubicSpline::sort_permutation(std::vector<RealType>& v) {
  std::vector<int> p(v.size());

  // 6 lines to replace std::iota(p.begin(), p.end(), 0);
  int value = 0;
  std::vector<int>::iterator i;
  for (i = p.begin(); i != p.end(); ++i) {
    (*i) = value;
    ++value;
  }
  
  std::sort(p.begin(), p.end(), OpenMD::Comparator(v) );
  return p;
}

std::vector<RealType> CubicSpline::apply_permutation(std::vector<RealType> const& v,
                                                     std::vector<int> const& p) { 
  int n = p.size();
  std::vector<double> sorted_vec(n);
  for (int i = 0; i < n; ++i) {
    sorted_vec[i] = v[p[i]];
  }
  return sorted_vec;
}
