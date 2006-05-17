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
#include "config.h"
namespace oopse {

  template<typename ElemType> ElemType pow(ElemType x, int N) {
    ElemType result(1);

    for (int i = 0; i < N; ++i) {
      result *= x;
    }

    return result;
  }

  /**
   * @class Polynomial Polynomial.hpp "math/Polynomial.hpp"
   * A generic Polynomial class
   */
  template<typename ElemType>
  class Polynomial {

  public:
    typedef Polynomial<ElemType> PolynomialType;    
    typedef int ExponentType;
    typedef ElemType CoefficientType;
    typedef std::map<ExponentType, CoefficientType> PolynomialPairMap;
    typedef typename PolynomialPairMap::iterator iterator;
    typedef typename PolynomialPairMap::const_iterator const_iterator;

    Polynomial() {}
    Polynomial(ElemType v) {setCoefficient(0, v);}
    /** 
     * Calculates the value of this Polynomial evaluated at the given x value.
     * @return The value of this Polynomial evaluates at the given x value
     * @param x the value of the independent variable for this Polynomial function
     */
    ElemType evaluate(const ElemType& x) {
      ElemType result = ElemType();
      ExponentType exponent;
      CoefficientType coefficient;
            
      for (iterator i = polyPairMap_.begin(); i != polyPairMap_.end(); ++i) {
	exponent = i->first;
	coefficient = i->second;
	result  += pow(x, exponent) * coefficient;
      }

      return result;
    }

    /**
     * Returns the first derivative of this polynomial.
     * @return the first derivative of this polynomial
     * @param x
     */
    ElemType evaluateDerivative(const ElemType& x) {
      ElemType result = ElemType();
      ExponentType exponent;
      CoefficientType coefficient;
            
      for (iterator i = polyPairMap_.begin(); i != polyPairMap_.end(); ++i) {
	exponent = i->first;
	coefficient = i->second;
	result  += pow(x, exponent - 1) * coefficient * exponent;
      }

      return result;
    }

    /**
     * Set the coefficent of the specified exponent, if the coefficient is already there, it
     * will be overwritten.
     * @param exponent exponent of a term in this Polynomial 
     * @param coefficient multiplier of a term in this Polynomial 
     */
        
    void setCoefficient(int exponent, const ElemType& coefficient) {
      polyPairMap_.insert(typename PolynomialPairMap::value_type(exponent, coefficient));
    }

    /**
     * Set the coefficent of the specified exponent. If the coefficient is already there,  just add the
     * new coefficient to the old one, otherwise,  just call setCoefficent
     * @param exponent exponent of a term in this Polynomial 
     * @param coefficient multiplier of a term in this Polynomial 
     */
        
    void addCoefficient(int exponent, const ElemType& coefficient) {
      iterator i = polyPairMap_.find(exponent);

      if (i != end()) {
	i->second += coefficient;
      } else {
	setCoefficient(exponent, coefficient);
      }
    }


    /**
     * Returns the coefficient associated with the given power for this Polynomial.
     * @return the coefficient associated with the given power for this Polynomial
     * @exponent exponent of any term in this Polynomial
     */
    ElemType getCoefficient(ExponentType exponent) {
      iterator i = polyPairMap_.find(exponent);

      if (i != end()) {
	return i->second;
      } else {
	return ElemType(0);
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

    PolynomialType& operator += (const PolynomialType& p) {
        typename Polynomial<ElemType>::const_iterator i;

        for (i =  p.begin(); i  != p.end(); ++i) {
          this->addCoefficient(i->first, i->second);
        }

        return *this;        
    }

    PolynomialType& operator -= (const PolynomialType& p) {
        typename Polynomial<ElemType>::const_iterator i;
        for (i =  p.begin(); i  != p.end(); ++i) {
          this->addCoefficient(i->first, -i->second);
        }        
        return *this;
    }

    PolynomialType& operator *= (const PolynomialType& p) {
    typename Polynomial<ElemType>::const_iterator i;
    typename Polynomial<ElemType>::const_iterator j;
    
    for (i = this->begin(); i !=this->end(); ++i) {
      for (j = p.begin(); j !=p.end(); ++j) {
	this->addCoefficient( i->first + j->first, i->second * j->second);
      }
    }

    return *this;
    }

   
  private:
        
    PolynomialPairMap polyPairMap_;
  };


  /**
   * Generates and returns the product of two given Polynomials.
   * @return A Polynomial containing the product of the two given Polynomial parameters
   */
  template<typename ElemType>
  Polynomial<ElemType> operator *(const Polynomial<ElemType>& p1, const Polynomial<ElemType>& p2) {
    typename Polynomial<ElemType>::const_iterator i;
    typename Polynomial<ElemType>::const_iterator j;
    Polynomial<ElemType> p;
    
    for (i = p1.begin(); i !=p1.end(); ++i) {
      for (j = p2.begin(); j !=p2.end(); ++j) {
	p.addCoefficient( i->first + j->first, i->second * j->second);
      }
    }

    return p;
  }

  template<typename ElemType>
  Polynomial<ElemType> operator *(const Polynomial<ElemType>& p, const ElemType v) {
    typename Polynomial<ElemType>::const_iterator i;
    Polynomial<ElemType> result;
    
    for (i = p.begin(); i !=p.end(); ++i) {
 	result.addCoefficient( i->first , i->second * v);
    }

    return result;
  }

  template<typename ElemType>
  Polynomial<ElemType> operator *( const ElemType v, const Polynomial<ElemType>& p) {
    typename Polynomial<ElemType>::const_iterator i;
    Polynomial<ElemType> result;
    
    for (i = p.begin(); i !=p.end(); ++i) {
 	result.addCoefficient( i->first , i->second * v);
    }

    return result;
  }
  
  /**
   * Generates and returns the sum of two given Polynomials.
   * @param p1 the first polynomial
   * @param p2 the second polynomial
   */
  template<typename ElemType>
  Polynomial<ElemType> operator +(const Polynomial<ElemType>& p1, const Polynomial<ElemType>& p2) {
    Polynomial<ElemType> p(p1);

    typename Polynomial<ElemType>::const_iterator i;

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
  template<typename ElemType>
  Polynomial<ElemType> operator -(const Polynomial<ElemType>& p1, const Polynomial<ElemType>& p2) {
    Polynomial<ElemType> p(p1);

    typename Polynomial<ElemType>::const_iterator i;

    for (i =  p2.begin(); i  != p2.end(); ++i) {
      p.addCoefficient(i->first, -i->second);
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
  template<typename ElemType>
  bool equal(const Polynomial<ElemType>& p1, const Polynomial<ElemType>& p2) {

    typename Polynomial<ElemType>::const_iterator i;
    typename Polynomial<ElemType>::const_iterator j;

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

} //end namespace oopse
#endif //MATH_POLYNOMIAL_HPP
