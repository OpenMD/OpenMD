/***************************************************************************
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, see <http://www.gnu.org/licenses/>.  *
 ***************************************************************************/

/***
 * This file was imported from qtpie, found at http://code.google.com/p/qtpie
 * and was modified minimally for use in OpenMD.
 *
 * No author attribution was found in the code, but it presumably is
 * the work of J. Chen and Todd J. Martinez.
 *
 * QTPIE (charge transfer with polarization current equilibration) is
 * a new charge model, similar to other charge models like QEq,
 * fluc-q, EEM or ABEEM. Unlike other existing charge models, however,
 * it is capable of describing both charge transfer and polarization
 * phenomena. It is also unique for its ability to describe
 * intermolecular charge transfer at reasonable computational cost.
 *
 * Good references to cite when using this code are:
 *
 * J. Chen and T. J. Martinez, "QTPIE: Charge transfer with
 * polarization current equalization. A fluctuating charge model with
 * correct asymptotics", Chemical Physics Letters 438 (2007),
 * 315-320. DOI: 10.1016/j.cplett.2007.02.065. Erratum: ibid, 463
 * (2008), 288. DOI: 10.1016/j.cplett.2008.08.060
 *
 * J. Chen, D. Hundertmark and T. J. Martinez, "A unified theoretical
 * framework for fluctuating-charge models in atom-space and in
 * bond-space", Journal of Chemical Physics 129 (2008), 214113. DOI:
 * 10.1063/1.3021400.
 *
 * J. Chen and T. J. Martinez, "Charge conservation in
 * electronegativity equalization and its implications for the
 * electrostatic properties of fluctuating-charge models", Journal of
 * Chemical Physics 131 (2009), 044114. DOI: 10.1063/1.3183167
 *
 * J. Chen and T. J. Martinez, "The dissociation catastrophe in
 * fluctuating-charge models and its implications for the concept of
 * atomic electronegativity", Progress in Theoretical Chemistry and
 * Physics, to appear. arXiv:0812.1543
 *
 * J. Chen, "Theory and applications of fluctuating-charge models",
 * PhD (chemical physics) thesis, University of Illinois at
 * Urbana-Champaign, Department of Chemistry, 2009.
 *
 * J. Chen and T. J. Martinez, "Size-extensive polarizabilities with
 * intermolecular charge transfer in a fluctuating-charge model", in
 * preparation. arXiv:0812.1544
 */

#include "config.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "math/Factorials.hpp"
#include "utils/NumericConstant.hpp"

#ifndef NONBONDED_SLATERINTEGRALS_HPP
#define NONBONDED_SLATERINTEGRALS_HPP

template <typename T> inline T sqr(T t) { return t*t; }
template <typename T> inline T mod(T x, T m)
{ return x<0 ? m - 1 - ((-x) - 1)%m : x%m; }

// #include "parameters.h"

/**
 * @brief Computes Rosen's Guillimer-Zener function A
 *  Computes Rosen's A integral, an auxiliary quantity needed to
 *  compute integrals involving Slater-type orbitals of s symmetry.
 * \f[
 *   A_n(\alpha) = \int_1^\infty x^n e^{-\alpha x}dx
 *   = \frac{n! e^{-\alpha}}{\alpha^{n+1}}\sum_{\nu=0}^n
 *   \frac{\alpha^\nu}{\nu!}
 * \f]
 * @param n - principal quantum number
 * @param a - Slater exponent 
 * @return the value of Rosen's A integral
 * @note N. Rosen, Phys. Rev., 38 (1931), 255
 */
inline RealType RosenA(int n, RealType a)
{
  RealType RosenA_ = 0.;
  if (a != 0.)
    {
      RealType Term = 1.;
      RosenA_ = Term;
      for (int nu=1; nu<=n; nu++)
        {
          Term *= a / nu;
          RosenA_ += Term;
        }
      RosenA_ = (RosenA_/Term) * (exp(-a)/a);
    }
  return RosenA_;
}

/**
 * @brief Computes Rosen's Guillimer-Zener function B
 *  Computes Rosen's B integral, an auxiliary quantity needed to
 *  compute integrals involving Slater-type orbitals of s symmetry.
 * \f[
 *  B_n(\alpha) = \int_{-1}^1 x^n e^{-\alpha x} dx
 *              = \frac{n!}{\alpha^{n+1}} 
 *                \sum_{\nu=0}^n \frac{e^\alpha(-\alpha)^\nu
 *                  - e^{-\alpha} \alpha^\nu}{\nu!}
 * \f]
 * @param n - principal quantum number
 * @param alpha - Slater exponent 
 * @return the value of Rosen's B integral
 * @note N. Rosen, Phys. Rev., 38 (1931), 255
 */
inline RealType RosenB(int n, RealType alpha)
{
  RealType TheSum, Term;
  RealType RosenB_, PSinhRosenA, PCoshRosenA, PHyperRosenA;

  if (alpha != 0.)
    {
      Term = 1.;
      bool IsPositive = true;
		
      // These two expressions are (up to constant factors) equivalent
      // to computing the hyperbolic sine and cosine of a respectively
      // The series consists of adding up these terms in an alternating fashion
      PSinhRosenA =  exp(alpha) - exp(-alpha);
      PCoshRosenA = -exp(alpha) - exp(-alpha);
      TheSum = PSinhRosenA;
      for (unsigned nu=1; nu<=n; nu++)
        {
          if (IsPositive)
            {
              PHyperRosenA = PCoshRosenA;
              IsPositive = false;
            }
          else // term to add should be negative
            {
              PHyperRosenA = PSinhRosenA;
              IsPositive = true;
            }
          Term *= alpha / nu;
          TheSum += Term * PHyperRosenA;
        }
      RosenB_ = TheSum / (alpha*Term);
    }
  else // pathological case of a=0
    {
      printf("WARNING, a = 0 in RosenB\n");
      RosenB_ = (1. - pow(-1., n)) / (n + 1.);
    }
  return RosenB_;
}

/** @brief Computes Rosen's D combinatorial factor
 *  Computes Rosen's D factor, an auxiliary quantity needed to
 *  compute integrals involving Slater-type orbitals of s symmetry.
 *  \f[
 *  RosenD^{mn}_p = \sum_k (-1)^k \frac{m! n!}
 *                  {(p-k)!(m-p+k)!(n-k)!k!}
 *  \f]
 * @return the value of Rosen's D factor
 * @note N. Rosen, Phys. Rev., 38 (1931), 255
 */
inline RealType RosenD(int m, int n, int p)
{
  if (m+n+p > maxFact)
    {
      printf("Error, arguments exceed maximum factorial computed %d > %d\n", m+n+p, maxFact);
      ::exit(0);
    }
	
  RealType RosenD_ = 0;
  for (int k=max(p-m,0); k<=min(n,p); k++)
    {
      if (mod(k,2) == 0)
        RosenD_ += (fact[m] / (fact[p-k] * fact[m-p+k])) * (fact[n] / (fact[n-k] * fact[k]));
      else
        RosenD_ -= (fact[m] / ( fact[p-k] * fact[m-p+k])) * (fact[n] / (fact[n-k] * fact[k]));
    }
  return RosenD_;
}

/** @brief Computes Coulomb integral analytically over s-type STOs
 * Computes the two-center Coulomb integral over Slater-type orbitals of s symmetry.
 * @param a : Slater zeta exponent of first atom in inverse Bohr (au)
 * @param b : Slater zeta exponent of second atom in inverse Bohr (au)
 * @param m : principal quantum number of first atom
 * @param n : principal quantum number of second atom
 * @param R : internuclear distance in atomic units (bohr)
 * @return value of the Coulomb potential energy integral
 * @note N. Rosen, Phys. Rev., 38 (1931), 255
 * @note In Rosen's paper, this integral is known as K2.
 */
inline RealType sSTOCoulInt(RealType a, RealType b, int m, int n, RealType R)
{
  RealType x, K2;
  RealType Factor1, Factor2, Term, OneElectronTerm;
  RealType eps, epsi;
	
  // To speed up calculation, we terminate loop once contributions
  // to integral fall below the bound, epsilon
  RealType epsilon = 0.;

  // x is the argument of the auxiliary RosenA and RosenB functions
  x = 2. * a * R;
	
  // First compute the two-electron component
  RealType sSTOCoulInt_ = 0.;
  if (std::fabs(x) < OpenMD::NumericConstant::epsilon) // Pathological case
    {

      // This solution for the one-center coulomb integrals comes from 
      // Yoshiyuki Hase, Computers & Chemistry 9(4), pp. 285-287 (1985).
      
      RealType Term1 = fact[2*m - 1] / pow(2*a, 2*m);
      RealType Term2 = 0.;
      for (int nu = 1; nu <= 2*n; nu++) {
        Term2 += nu * pow(2*b, 2*n - nu) * fact[2*(m+n)-nu-1] / (fact[2*n-nu]*2*n * pow(2*(a+b), 2*(m+n)-nu));
      }
      sSTOCoulInt_ = pow(2*a, 2*m+1) * (Term1 - Term2) / fact[2*m];

      // Original QTPIE code for the one-center case is below.  Doesn't
      // appear to generate the correct one-center results.
      //
      // if ((a==b) && (m==n))
      //   {
      //     for (int nu=0; nu<=2*n-1; nu++)
      //       {
      //         K2 = 0.;
      //         for (unsigned p=0; p<=2*n+m; p++) K2 += 1. / fact[p];
      //         sSTOCoulInt_ += K2 * fact[2*n+m] / fact[m];
      //       }
      //     sSTOCoulInt_ = 2 * a / (n * fact[2*n]) * sSTOCoulInt_;
      //   }
      // else
      //   {
      //     // Not implemented
      //     printf("ERROR, sSTOCoulInt cannot compute from arguments\n");
      //     printf("a = %lf b = %lf m = %d n = %d R = %lf\n",a, b, m, n, R);
      //     exit(0);
      //   }

    }
  else
    {
      OneElectronTerm = 1./R + pow(x, 2*m)/(fact[2*m]*R)*
        ((x-2*m)*RosenA(2*m-1,x)-exp(-x)) + sSTOCoulInt_;
      eps = epsilon / OneElectronTerm;
      if (a == b)
        {
          // Apply Rosen (48)
          Factor1 = -a*pow(a*R, 2*m)/(n*fact[2*m]);
          for (int nu=0; nu<=2*n-1; nu++)
            {
              Factor2 = (2.*n-nu)/fact[nu]*pow(a*R,nu);
              epsi = eps / fabs(Factor1 * Factor2);
              K2 = 0.;
              for (int p=0; p<=m+(nu-1)/2; p++)
                {
                  Term = RosenD(2*m-1, nu, 2*p)/(2.*p+1.) *RosenA(2*m+nu-1-2*p,x);
                  K2 += Term;
                  if ((Term > 0) && (Term < epsi)) goto label1;
                }
              sSTOCoulInt_ +=  K2 * Factor2;
            }
        label1:
          sSTOCoulInt_ *= Factor1;
        }
      else
        {
          Factor1 = -a*pow(a*R,2*m)/(2.*n*fact[2*m]);
          epsi = eps/fabs(Factor1);
          if (b == 0.)
            printf("WARNING: b = 0 in sSTOCoulInt\n");
          else
            {
              // Apply Rosen (54)
              for (int nu=0; nu<=2*n-1; nu++)
                {
                  K2 = 0;
                  for (int p=0; p<=2*m+nu-1; p++)
                    K2=K2+RosenD(2*m-1,nu,p)*RosenB(p,R*(a-b))
                      *RosenA(2*m+nu-1-p,R*(a+b));
                  Term = K2*(2*n-nu)/fact[nu]*pow(b*R, nu);
                  sSTOCoulInt_ += Term;
                  if (fabs(Term) < epsi) goto label2;
                }
            label2:
              sSTOCoulInt_ *= Factor1;
            }
        }
      // Now add the one-electron term from Rosen (47) = Rosen (53)
      sSTOCoulInt_ += OneElectronTerm;
    }
  return sSTOCoulInt_;
}

/**
 * @brief Computes overlap integral analytically over s-type STOs
 *   Computes the overlap integral over two
 *   Slater-type orbitals of s symmetry.
 * @param a : Slater zeta exponent of first atom in inverse Bohr (au)
 * @param b : Slater zeta exponent of second atom in inverse Bohr (au)
 * @param m : principal quantum number of first atom
 * @param n : principal quantum number of second atom
 * @param R : internuclear distance in atomic units (bohr)
 * @return the value of the sSTOOvInt integral
 * @note N. Rosen, Phys. Rev., 38 (1931), 255
 * @note In the Rosen paper, this integral is known as I.
 */
inline RealType sSTOOvInt(RealType a, RealType b, int m, int n, RealType R)
{
  RealType Factor, Term, eps;
	
  // To speed up calculation, we terminate loop once contributions
  // to integral fall below the bound, epsilon
  RealType epsilon = 0.;
  RealType sSTOOvInt_ = 0.;
	
  if (a == b)
    {
      Factor = pow(a*R, m+n+1)/sqrt(fact[2*m]*fact[2*n]);
      eps = epsilon / fabs(Factor);
      for (int q=0; q<=(m+n)/2; q++)
        {
          Term = RosenD(m,n,2*q)/(2.*q+1.)*RosenA(m+n-2*q,a*R);
          sSTOOvInt_ += Term;
          if (fabs(Term) < eps) exit(0);
        }
      sSTOOvInt_ *= Factor;
    }
  else
    {
      Factor = 0.5*pow(a*R, m+0.5)*pow(b*R,n+0.5)
        /sqrt(fact[2*m]*fact[2*n]);
      eps = epsilon / fabs(Factor);
      for (int q=0; q<=m+n; q++)
        {
          Term = RosenD(m,n,q)*RosenB(q, R/2.*(a-b))
            * RosenA(m+n-q,R/2.*(a+b));
          sSTOOvInt_ += Term;
          if (fabs(Term) < eps) exit(0);
        }
      sSTOOvInt_ *= Factor;
    }
  return sSTOOvInt_;
}

/**
 * @brief Computes kinetic energy integral analytically over s-type STOs
 *   Computes the overlap integral over two Slater-type orbitals of s symmetry.
 * @param a : Slater zeta exponent of first atom in inverse Bohr (au)
 * @param b : Slater zeta exponent of second atom in inverse Bohr (au)
 * @param m : principal quantum number of first atom
 * @param n : principal quantum number of second atom
 * @param R : internuclear distance in atomic units (bohr)
 * @return the value of the kinetic energy integral
 * @note N. Rosen, Phys. Rev., 38 (1931), 255
 * @note untested
 */
inline RealType KinInt(RealType a, RealType b, int m, int n,RealType R)
{
  RealType KinInt_ = -0.5*b*b*sSTOOvInt(a, b, m, n, R);
  if (n > 0)
    {
      KinInt_ += b*b*pow(2*b/(2*b-1),0.5) * sSTOOvInt(a, b, m, n-1, R);
      if (n > 1) KinInt_ += pow(n*(n-1)/((n-0.5)*(n-1.5)), 0.5)
                   * sSTOOvInt(a, b, m, n-2, R);
    }
  return KinInt_;
}

/**
 * @brief Computes derivative of Coulomb integral with respect to the interatomic distance
 *   Computes the two-center Coulomb integral over Slater-type orbitals of s symmetry.
 * @param a: Slater zeta exponent of first atom in inverse Bohr (au)
 * @param b: Slater zeta exponent of second atom in inverse Bohr (au)
 * @param m: principal quantum number of first atom
 * @param n: principal quantum number of second atom
 * @param R: internuclear distance in atomic units (bohr)
 * @return the derivative of the Coulomb potential energy integral
 * @note Derived in QTPIE research notes, May 15 2007
 */
inline RealType sSTOCoulIntGrad(RealType a, RealType b, int m, int n, RealType R)
{
  RealType x, y, z, K2, TheSum;
  // x is the argument of the auxiliary RosenA and RosenB functions
  x = 2. * a * R;
	
  // First compute the two-electron component
  RealType sSTOCoulIntGrad_ = 0.;
  if (x==0) // Pathological case
    {
      printf("WARNING: argument given to sSTOCoulIntGrad is 0\n");
      printf("a = %lf R= %lf\n", a, R);
    }
  else
    {
      if (a == b)
        {
          TheSum = 0.;
          for (int nu=0; nu<=2*(n-1); nu++)
            {
              K2 = 0.;
              for (int p=0; p<=(m+nu)/2; p++)
                K2 += RosenD(2*m-1, nu+1, 2*p)/(2*p + 1.) * RosenA(2*m+nu-1-2*p, x);
              TheSum += (2*n-nu-1)/fact[nu]*pow(a*R, nu) * K2;
            }
          sSTOCoulIntGrad_ = -pow(a, 2*m+2)*pow(R, 2*m) /(n*fact[2*m])*TheSum;
          TheSum = 0.;
          for (int nu=0; nu<=2*n-1; nu++)
            {
              K2 = 0.;
              for (int p=0; p<=(m+nu-1)/2; p++)
                K2 += RosenD(2*m-1, nu, 2*p)/(2*p + 1.) * RosenA(2*m+nu-2*p, x);
              TheSum += (2*n-nu)/fact[nu]*pow(a*R,nu) * K2;
            }
          sSTOCoulIntGrad_ += 2*pow(a, 2*m+2)*pow(R, 2*m) /(n*fact[2*m])*TheSum;
        }
      else
        {
          // Slater exponents are different
          // First calculate some useful arguments
          y = R*(a+b);
          z = R*(a-b);
          TheSum = 0.;
          for (int nu=0; nu<=2*n-1; nu++)
            {
              K2 = 0.;
              for (int p=0; p<=2*m+nu; p++)
                K2 += RosenD(2*m-1, nu+1, p)
                  * RosenB(p,z)*RosenA(2*m+nu-p, y);
              TheSum += (2*n-nu-1)/fact[nu]*pow(b*R,nu) * K2;
            }
          sSTOCoulIntGrad_ = -b*pow(a,2*m+1)*pow(R,2*m)/
            (2*n*fact[2*m])*TheSum;
          TheSum = 0.;
          for (int nu=0; nu<=2*n; nu++)
            {
              K2 = 0.;
              for (int p=0; p<=2*m-1+nu; p++)
                K2 += RosenD(2*m-1, nu, p)
                  * ((a-b)*RosenB(p+1,z)*RosenA(2*m+nu-p-1, y)
                     +(a+b)*RosenB(p  ,z)*RosenA(2*m+nu-p  , y));
              TheSum += (2*n-nu)/fact[nu]*pow(b*R,nu) * K2;
            }
          sSTOCoulIntGrad_ += pow(a,2*m+1)*pow(R,2*m)/(2*n*fact[2*m])*TheSum;
        }
      // Now add one-electron terms and common term
      sSTOCoulIntGrad_ = sSTOCoulIntGrad_ - (2.*m+1.)/sqr(R)
        + 2.*m/R * sSTOCoulInt(a,b,m,n,R)
        + pow(x,2*m)/(fact[2*m]*sqr(R)) * ((2.*m+1.)*exp(-x)
                                           + 2.*m*(1.+2.*m-x)*RosenA(2*m-1,x));
    }
  return sSTOCoulIntGrad_;
}

/** 
 * @brief Computes gradient of overlap integral with respect to the interatomic diatance
 *   Computes the derivative of the overlap integral over two Slater-type orbitals of s symmetry.
 * @param a: Slater zeta exponent of first atom in inverse Bohr (au)
 * @param b: Slater zeta exponent of second atom in inverse Bohr (au)
 * @param m: principal quantum number of first atom
 * @param n: principal quantum number of second atom
 * @param R: internuclear distance in atomic units (bohr)
 * @return the derivative of the sSTOOvInt integral
 * @note Derived in QTPIE research notes, May 15 2007
 */
inline RealType sSTOOvIntGrad(RealType a, RealType b, int m, int n, RealType R)
{
  RealType w, x, y, z, TheSum;
	
  // Calculate first term
  RealType sSTOOvIntGrad_ = (m+n+1.)/R * sSTOOvInt(a, b, m, n, R);
	
  // Calculate remaining terms; answers depend on exponents 
  TheSum = 0.;
  x = a * R;
  if (a == b)
    {
      for (int q=0; q<=(m+n)/2; q++)
        TheSum += RosenD(m,n,2*q) / (2*q + 1.) * RosenA(m+n-2*q+1, x);
      sSTOOvIntGrad_ -= a*pow(x,m+n+1)/ sqrt(fact[2*m]*fact[2*n])*TheSum;
    }
  else
    {
      w = b*R;
      y = 0.5*R*(a+b);
      z = 0.5*R*(a-b);
      for (int q=0; q<m+n; q++)
        TheSum = TheSum + RosenD(m,n,q) *
          ((a-b)*RosenB(q+1,z)*RosenA(m+n-q  ,y)
           +(a+b)*RosenB(q  ,z)*RosenA(m+n-q+1,y));
      sSTOOvIntGrad_ -= 0.25*sqrt((pow(x, 2*m+1)*pow(w, 2*n+1))/(fact[2*m]*fact[2*n]))*TheSum;
    }
  return sSTOOvIntGrad_;
}

/**
 * @brief Calculates a Slater-type orbital exponent based on the hardness parameters
 * @param hardness: chemical hardness in atomic units
 * @param        n: principal quantum number
 * @note Modified for use with OpenMD by Gezelter and Michalka.
 */
inline RealType getSTOZeta(int n, RealType hardness)
{
  //  Approximate the exact value of the constant of proportionality
  //  by its value at a very small distance epsilon
  //  since the exact R = 0 case has not be programmed
  RealType epsilon = 1.0e-8;
  
  // Assign orbital exponent  
  return pow(sSTOCoulInt(1., 1., n, n, epsilon) / hardness, -1./(3. + 2.*n));
}

#endif
