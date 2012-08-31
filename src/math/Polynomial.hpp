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
 * @file Polynomial.hpp
 * @author    teng lin
 * @date  11/16/2004
 * @version 1.0
 */ 

#ifndef MATH_POLYNOMIAL_HPP
#define MATH_POLYNOMIAL_HPP

#include <iostream>
#include <list>
#include <map>
#include <utility>
#include <complex>
#include "config.h"
#include "math/Eigenvalue.hpp"

namespace OpenMD {
  
  template<typename Real> Real fastpow(Real x, int N) {
    Real result(1); //or 1.0?

    for (int i = 0; i < N; ++i) {
      result *= x;
    }

    return result;
  }

  /**
   * @class Polynomial Polynomial.hpp "math/Polynomial.hpp"
   * A generic Polynomial class
   */
  template<typename Real>
  class Polynomial {

  public:
    typedef Polynomial<Real> PolynomialType;    
    typedef int ExponentType;
    typedef Real CoefficientType;
    typedef std::map<ExponentType, CoefficientType> PolynomialPairMap;
    typedef typename PolynomialPairMap::iterator iterator;
    typedef typename PolynomialPairMap::const_iterator const_iterator;

    Polynomial() {}
    Polynomial(Real v) {setCoefficient(0, v);}
    /** 
     * Calculates the value of this Polynomial evaluated at the given x value.
     * @return The value of this Polynomial evaluates at the given x value  
     * @param x the value of the independent variable for this
     * Polynomial function
     */
    Real evaluate(const Real& x) {
      Real result = Real();
      ExponentType exponent;
      CoefficientType coefficient;
            
      for (iterator i = polyPairMap_.begin(); i != polyPairMap_.end(); ++i) {
	exponent = i->first;
	coefficient = i->second;
	result  += fastpow(x, exponent) * coefficient;
      }

      return result;
    }

    /**
     * Returns the first derivative of this polynomial.
     * @return the first derivative of this polynomial
     * @param x
     */
    Real evaluateDerivative(const Real& x) {
      Real result = Real();
      ExponentType exponent;
      CoefficientType coefficient;
            
      for (iterator i = polyPairMap_.begin(); i != polyPairMap_.end(); ++i) {
	exponent = i->first;
	coefficient = i->second;
	result  += fastpow(x, exponent - 1) * coefficient * exponent;
      }

      return result;
    }


    /**
     * Set the coefficent of the specified exponent, if the
     * coefficient is already there, it will be overwritten.
     * @param exponent exponent of a term in this Polynomial 
     * @param coefficient multiplier of a term in this Polynomial 
     */        
    void setCoefficient(int exponent, const Real& coefficient) {
      polyPairMap_[exponent] = coefficient;
    }
    
    /**
     * Set the coefficent of the specified exponent. If the
     * coefficient is already there, just add the new coefficient to
     * the old one, otherwise, just call setCoefficent
     * @param exponent exponent of a term in this Polynomial 
     * @param coefficient multiplier of a term in this Polynomial 
     */        
    void addCoefficient(int exponent, const Real& coefficient) {
      iterator i = polyPairMap_.find(exponent);

      if (i != end()) {
	i->second += coefficient;
      } else {
	setCoefficient(exponent, coefficient);
      }
    }

    /**
     * Returns the coefficient associated with the given power for
     * this Polynomial.
     * @return the coefficient associated with the given power for
     * this Polynomial
     * @exponent exponent of any term in this Polynomial
     */
    Real getCoefficient(ExponentType exponent) {
      iterator i = polyPairMap_.find(exponent);

      if (i != end()) {
	return i->second;
      } else {
	return Real(0);
      }
    }

    iterator begin() {
      return polyPairMap_.begin();
    }

    const_iterator begin() const{
      return polyPairMap_.begin();
    }
        
    iterator end() {
      return polyPairMap_.end();
    }

    const_iterator end() const{
      return polyPairMap_.end();
    }

    iterator find(ExponentType exponent) {
      return polyPairMap_.find(exponent);
    }

    size_t size() {
      return polyPairMap_.size();
    }

    int degree() {
      int deg = 0;
      for (iterator i = polyPairMap_.begin(); i != polyPairMap_.end(); ++i) {
        if (i->first > deg)
          deg = i->first;
      }
      return deg;
    }

    PolynomialType& operator = (const PolynomialType& p) {

      if (this != &p)  // protect against invalid self-assignment
        {
          typename Polynomial<Real>::const_iterator i;

          polyPairMap_.clear();  // clear out the old map
      
          for (i =  p.begin(); i != p.end(); ++i) {
            this->setCoefficient(i->first, i->second);
          }
        }
      // by convention, always return *this
      return *this; 
    }

    PolynomialType& operator += (const PolynomialType& p) {
      typename Polynomial<Real>::const_iterator i;

      for (i =  p.begin(); i  != p.end(); ++i) {
        this->addCoefficient(i->first, i->second);
      }

      return *this;        
    }

    PolynomialType& operator -= (const PolynomialType& p) {
      typename Polynomial<Real>::const_iterator i;
      for (i =  p.begin(); i  != p.end(); ++i) {
        this->addCoefficient(i->first, -i->second);
      }        
      return *this;
    }
    
    PolynomialType& operator *= (const PolynomialType& p) {
      typename Polynomial<Real>::const_iterator i;
      typename Polynomial<Real>::const_iterator j;
      Polynomial<Real> p2(*this);
      
      polyPairMap_.clear();  // clear out old map
      for (i = p2.begin(); i !=p2.end(); ++i) {
        for (j = p.begin(); j !=p.end(); ++j) {
          this->addCoefficient( i->first + j->first, i->second * j->second);
        }
      }
      return *this;
    }

    //PolynomialType& operator *= (const Real v)
    PolynomialType& operator *= (const Real v) {
      typename Polynomial<Real>::const_iterator i;
      //Polynomial<Real> result;
      
      for (i = this->begin(); i != this->end(); ++i) {
	this->setCoefficient( i->first, i->second*v);
      }
      
      return *this;
    }

    PolynomialType& operator += (const Real v) {    
      this->addCoefficient( 0, v);
      return *this;
    }

    /**
     * Returns the first derivative of this polynomial.
     * @return the first derivative of this polynomial
     */
    PolynomialType & getDerivative() {
      Polynomial<Real> p;
      
      typename Polynomial<Real>::const_iterator i;
      ExponentType exponent;
      CoefficientType coefficient;
      
      for (i =  this->begin(); i  != this->end(); ++i) {
        exponent = i->first;
        coefficient = i->second;
        p.setCoefficient(exponent-1, coefficient * exponent);
      }
    
      return p;
    }

    // Creates the Companion matrix for a given polynomial
    DynamicRectMatrix<Real> CreateCompanion() {
      int rank = degree();
      DynamicRectMatrix<Real> mat(rank, rank);
      Real majorCoeff = getCoefficient(rank);
      for(int i = 0; i < rank; ++i) {
        for(int j = 0; j < rank; ++j) {
          if(i - j == 1) {
            mat(i, j) = 1;
          } else if(j == rank-1) {
            mat(i, j) = -1 * getCoefficient(i) / majorCoeff;
          }
        }
      }
      return mat;
    }
    
    // Find the Roots of a given polynomial
    std::vector<complex<Real> > FindRoots() {
      int rank = degree();
      DynamicRectMatrix<Real> companion = CreateCompanion();
      JAMA::Eigenvalue<Real> eig(companion);
      DynamicVector<Real> reals, imags;
      eig.getRealEigenvalues(reals);
      eig.getImagEigenvalues(imags);
      
      std::vector<complex<Real> > roots;
      for (int i = 0; i < rank; i++) {
        roots.push_back(complex<Real>(reals(i), imags(i)));
      }

      return roots;
    }

    std::vector<Real> FindRealRoots() {
      
      const Real fEpsilon = 1.0e-8;
      std::vector<Real> roots;
      roots.clear();
      
      const int deg = degree();
      
      switch (deg) {
      case 1: {
        Real fC1 = getCoefficient(1);
        Real fC0 = getCoefficient(0);
        roots.push_back( -fC0 / fC1);
        return roots;
      }
      case 2: {
        Real fC2 = getCoefficient(2);
        Real fC1 = getCoefficient(1);
        Real fC0 = getCoefficient(0);
        Real fDiscr = fC1*fC1 - 4.0*fC0*fC2;
        if (abs(fDiscr) <= fEpsilon) {
          fDiscr = (Real)0.0;
        }
      
        if (fDiscr < (Real)0.0) {  // complex roots only
          return roots;
        }
      
        Real fTmp = ((Real)0.5)/fC2;
      
        if (fDiscr > (Real)0.0) { // 2 real roots
          fDiscr = sqrt(fDiscr);
          roots.push_back(fTmp*(-fC1 - fDiscr));
          roots.push_back(fTmp*(-fC1 + fDiscr));
        } else {
          roots.push_back(-fTmp * fC1);  // 1 real root
        }
      }
        return roots;        
      case 3: {
        Real fC3 = getCoefficient(3);
        Real fC2 = getCoefficient(2);
        Real fC1 = getCoefficient(1);
        Real fC0 = getCoefficient(0);
      
        // make polynomial monic, x^3+c2*x^2+c1*x+c0
        Real fInvC3 = ((Real)1.0)/fC3;
        fC0 *= fInvC3;
        fC1 *= fInvC3;
        fC2 *= fInvC3;
      
        // convert to y^3+a*y+b = 0 by x = y-c2/3
        const Real fThird = (Real)1.0/(Real)3.0;
        const Real fTwentySeventh = (Real)1.0/(Real)27.0;
        Real fOffset = fThird*fC2;
        Real fA = fC1 - fC2*fOffset;
        Real fB = fC0+fC2*(((Real)2.0)*fC2*fC2-((Real)9.0)*fC1)*fTwentySeventh;
        Real fHalfB = ((Real)0.5)*fB;
      
        Real fDiscr = fHalfB*fHalfB + fA*fA*fA*fTwentySeventh;
        if (fabs(fDiscr) <= fEpsilon) {
          fDiscr = (Real)0.0;
        }
      
        if (fDiscr > (Real)0.0) {  // 1 real, 2 complex roots
        
          fDiscr = sqrt(fDiscr);
          Real fTemp = -fHalfB + fDiscr;
          Real root;
          if (fTemp >= (Real)0.0) {
            root = pow(fTemp,fThird);
          } else {
            root = -pow(-fTemp,fThird);
          }
          fTemp = -fHalfB - fDiscr;
          if ( fTemp >= (Real)0.0 ) {
            root += pow(fTemp,fThird);          
          } else {
            root -= pow(-fTemp,fThird);
          }
          root -= fOffset;
        
          roots.push_back(root);
        } else if (fDiscr < (Real)0.0) {
          const Real fSqrt3 = sqrt((Real)3.0);
          Real fDist = sqrt(-fThird*fA);
          Real fAngle = fThird*atan2(sqrt(-fDiscr), -fHalfB);
          Real fCos = cos(fAngle);
          Real fSin = sin(fAngle);
          roots.push_back(((Real)2.0)*fDist*fCos-fOffset);
          roots.push_back(-fDist*(fCos+fSqrt3*fSin)-fOffset);
          roots.push_back(-fDist*(fCos-fSqrt3*fSin)-fOffset);
        } else {
          Real fTemp;
          if (fHalfB >= (Real)0.0) {
            fTemp = -pow(fHalfB,fThird);
          } else {
            fTemp = pow(-fHalfB,fThird);
          }
          roots.push_back(((Real)2.0)*fTemp-fOffset);
          roots.push_back(-fTemp-fOffset);
          roots.push_back(-fTemp-fOffset);
        }
      }
        return roots;

      case 4: {
        Real fC4 = getCoefficient(4);
        Real fC3 = getCoefficient(3);
        Real fC2 = getCoefficient(2);
        Real fC1 = getCoefficient(1);
        Real fC0 = getCoefficient(0);
      
        // make polynomial monic, x^4+c3*x^3+c2*x^2+c1*x+c0
        Real fInvC4 = ((Real)1.0)/fC4;
        fC0 *= fInvC4;
        fC1 *= fInvC4;
        fC2 *= fInvC4;
        fC3 *= fInvC4;
  
        // reduction to resolvent cubic polynomial y^3+r2*y^2+r1*y+r0 = 0
        Real fR0 = -fC3*fC3*fC0 + ((Real)4.0)*fC2*fC0 - fC1*fC1;
        Real fR1 = fC3*fC1 - ((Real)4.0)*fC0;
        Real fR2 = -fC2;
        Polynomial<Real> tempCubic;
        tempCubic.setCoefficient(0, fR0);
        tempCubic.setCoefficient(1, fR1);
        tempCubic.setCoefficient(2, fR2);
        tempCubic.setCoefficient(3, 1.0);
        std::vector<Real> cubeRoots = tempCubic.FindRealRoots(); // always
        // produces
        // at
        // least
        // one
        // root
        Real fY = cubeRoots[0];
      
        Real fDiscr = ((Real)0.25)*fC3*fC3 - fC2 + fY;
        if (fabs(fDiscr) <= fEpsilon) {
          fDiscr = (Real)0.0;
        }
   
        if (fDiscr > (Real)0.0) {
          Real fR = sqrt(fDiscr);
          Real fT1 = ((Real)0.75)*fC3*fC3 - fR*fR - ((Real)2.0)*fC2;
          Real fT2 = (((Real)4.0)*fC3*fC2 - ((Real)8.0)*fC1 - fC3*fC3*fC3) /
            (((Real)4.0)*fR);
      
          Real fTplus = fT1+fT2;
          Real fTminus = fT1-fT2;
          if (fabs(fTplus) <= fEpsilon) {
            fTplus = (Real)0.0;
          }
          if (fabs(fTminus) <= fEpsilon) {
            fTminus = (Real)0.0;
          }
      
          if (fTplus >= (Real)0.0) {
            Real fD = sqrt(fTplus);
            roots.push_back(-((Real)0.25)*fC3+((Real)0.5)*(fR+fD));
            roots.push_back(-((Real)0.25)*fC3+((Real)0.5)*(fR-fD));
          }
          if (fTminus >= (Real)0.0) {
            Real fE = sqrt(fTminus);
            roots.push_back(-((Real)0.25)*fC3+((Real)0.5)*(fE-fR));
            roots.push_back(-((Real)0.25)*fC3-((Real)0.5)*(fE+fR));
          }
        } else if (fDiscr < (Real)0.0) {
          //roots.clear();
        } else {        
          Real fT2 = fY*fY-((Real)4.0)*fC0;
          if (fT2 >= -fEpsilon) {
            if (fT2 < (Real)0.0) { // round to zero
              fT2 = (Real)0.0;
            }
            fT2 = ((Real)2.0)*sqrt(fT2);
            Real fT1 = ((Real)0.75)*fC3*fC3 - ((Real)2.0)*fC2;
            if (fT1+fT2 >= fEpsilon) {
              Real fD = sqrt(fT1+fT2);
              roots.push_back( -((Real)0.25)*fC3+((Real)0.5)*fD);
              roots.push_back( -((Real)0.25)*fC3-((Real)0.5)*fD);
            }
            if (fT1-fT2 >= fEpsilon) {
              Real fE = sqrt(fT1-fT2);
              roots.push_back( -((Real)0.25)*fC3+((Real)0.5)*fE);
              roots.push_back( -((Real)0.25)*fC3-((Real)0.5)*fE);
            }
          }
        }
      }
        return roots;

      default: {
        DynamicRectMatrix<Real> companion = CreateCompanion();
        JAMA::Eigenvalue<Real> eig(companion);
        DynamicVector<Real> reals, imags;
        eig.getRealEigenvalues(reals);
        eig.getImagEigenvalues(imags);
      
        for (int i = 0; i < deg; i++) {
          if (fabs(imags(i)) < fEpsilon) 
            roots.push_back(reals(i));        
        }      
      }
        return roots;
        
      }     
    }
    
  private:
        
    PolynomialPairMap polyPairMap_;
  };

  
  /**
   * Generates and returns the product of two given Polynomials.
   * @return A Polynomial containing the product of the two given Polynomial parameters
   */
  template<typename Real>
  Polynomial<Real> operator *(const Polynomial<Real>& p1, const Polynomial<Real>& p2) {
    typename Polynomial<Real>::const_iterator i;
    typename Polynomial<Real>::const_iterator j;
    Polynomial<Real> p;
    
    for (i = p1.begin(); i !=p1.end(); ++i) {
      for (j = p2.begin(); j !=p2.end(); ++j) {
	p.addCoefficient( i->first + j->first, i->second * j->second);
      }
    }

    return p;
  }

  template<typename Real>
  Polynomial<Real> operator *(const Polynomial<Real>& p, const Real v) {
    typename Polynomial<Real>::const_iterator i;
    Polynomial<Real> result;
    
    for (i = p.begin(); i !=p.end(); ++i) {
 	result.setCoefficient( i->first , i->second * v);
    }

    return result;
  }

  template<typename Real>
  Polynomial<Real> operator *( const Real v, const Polynomial<Real>& p) {
    typename Polynomial<Real>::const_iterator i;
    Polynomial<Real> result;
    
    for (i = p.begin(); i !=p.end(); ++i) {
 	result.setCoefficient( i->first , i->second * v);
    }

    return result;
  }
  
  /**
   * Generates and returns the sum of two given Polynomials.
   * @param p1 the first polynomial
   * @param p2 the second polynomial
   */
  template<typename Real>
  Polynomial<Real> operator +(const Polynomial<Real>& p1, const Polynomial<Real>& p2) {
    Polynomial<Real> p(p1);

    typename Polynomial<Real>::const_iterator i;

    for (i =  p2.begin(); i  != p2.end(); ++i) {
      p.addCoefficient(i->first, i->second);
    }

    return p;

  }

  /**
   * Generates and returns the difference of two given Polynomials. 
   * @return
   * @param p1 the first polynomial
   * @param p2 the second polynomial
   */
  template<typename Real>
  Polynomial<Real> operator -(const Polynomial<Real>& p1, const Polynomial<Real>& p2) {
    Polynomial<Real> p(p1);

    typename Polynomial<Real>::const_iterator i;

    for (i =  p2.begin(); i  != p2.end(); ++i) {
      p.addCoefficient(i->first, -i->second);
    }

    return p;

  }

  /**
   * Returns the first derivative of this polynomial.
   * @return the first derivative of this polynomial
   */
  template<typename Real>
  Polynomial<Real> getDerivative(const Polynomial<Real>& p1) {
    Polynomial<Real> p;
    
    typename Polynomial<Real>::const_iterator i;
    int exponent;
    Real coefficient;
    
    for (i =  p1.begin(); i  != p1.end(); ++i) {
      exponent = i->first;
      coefficient = i->second;
      p.setCoefficient(exponent-1, coefficient * exponent);
    }
    
    return p;
  }

  /**
   * Tests if two polynomial have the same exponents
   * @return true if all of the exponents in these Polynomial are identical
   * @param p1 the first polynomial
   * @param p2 the second polynomial
   * @note this function does not compare the coefficient
   */
  template<typename Real>
  bool equal(const Polynomial<Real>& p1, const Polynomial<Real>& p2) {

    typename Polynomial<Real>::const_iterator i;
    typename Polynomial<Real>::const_iterator j;

    if (p1.size() != p2.size() ) {
      return false;
    }
    
    for (i =  p1.begin(), j = p2.begin(); i  != p1.end() && j != p2.end(); ++i, ++j) {
      if (i->first != j->first) {
	return false;
      }
    }

    return true;
  }



  typedef Polynomial<RealType> DoublePolynomial;

} //end namespace OpenMD
#endif //MATH_POLYNOMIAL_HPP
