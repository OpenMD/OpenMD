/*
 * Borrowed from OpenMM.
 */

#include <config.h>
#ifndef MATH_ERFC_H
#define MATH_ERFC_H

/*
 * Up to version 11 (VC++ 2012), Microsoft does not support the
 * standard C99 erf() and erfc() functions so we have to fake them
 * here.  These were added in version 12 (VC++ 2013), which sets
 * _MSC_VER=1800 (VC11 has _MSC_VER=1700).
 */

#ifdef _MSC_VER
#if _MSC_VER <= 1700  // 1700 is VC11, 1800 is VC12

/***************************
 *   erf.cpp
 *   author:  Steve Strand
 *   written: 29-Jan-04
 ***************************/

#include <cmath>

static const RealType rel_error = 1E-12;  // calculate 12 significant figures
// you can adjust rel_error to trade off between accuracy and speed
// but don't ask for > 15 figures (assuming usual 52 bit mantissa in a double)

static RealType erfc(RealType x);

static RealType erf(RealType x)
// erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
//       = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
//       = 1-erfc(x)
{
  static const RealType two_sqrtpi = 1.128379167095512574;  // 2/sqrt(pi)
  if (fabs(x) > 2.2) {
    return 1.0 - erfc(x);  // use continued fraction when fabs(x) > 2.2
  }
  RealType sum = x, term = x, xsqr = x * x;
  int j = 1;
  do {
    term *= xsqr / j;
    sum -= term / (2 * j + 1);
    ++j;
    term *= xsqr / j;
    sum += term / (2 * j + 1);
    ++j;
  } while (fabs(term) / sum > rel_error);
  return two_sqrtpi * sum;
}

static RealType erfc(RealType x)
// erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
//        = exp(-x^2)/sqrt(pi) * [1/x+ (1/2)/x+ (2/2)/x+ (3/2)/x+ (4/2)/x+ ...]
//        = 1-erf(x)
// expression inside [] is a continued fraction so '+' means add to denominator
// only
{
  static const RealType one_sqrtpi = 0.564189583547756287;  // 1/sqrt(pi)
  if (fabs(x) < 2.2) {
    return 1.0 - erf(x);  // use series when fabs(x) < 2.2
  }
  // Don't look for x==0 here!
  if (x < 0) {  // continued fraction only valid for x>0
    return 2.0 - erfc(-x);
  }
  RealType a = 1, b = x;            // last two convergent numerators
  RealType c = x, d = x * x + 0.5;  // last two convergent denominators
  RealType q1, q2   = b / d;        // last two convergents (a/c and b/d)
  RealType n = 1.0, t;
  do {
    t = a * n + b * x;
    a = b;
    b = t;
    t = c * n + d * x;
    c = d;
    d = t;
    n += 0.5;
    q1 = q2;
    q2 = b / d;
  } while (fabs(q1 - q2) / q2 > rel_error);
  return one_sqrtpi * exp(-x * x) * q2;
}

#endif  // _MSC_VER <= 1700
#endif  // _MSC_VER
#endif
