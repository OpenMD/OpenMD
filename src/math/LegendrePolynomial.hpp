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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file LegendrePolynomial.hpp
 * @author    teng lin
 * @date  11/16/2004
 * @version 1.0
 */ 

#ifndef MATH_LEGENDREPOLYNOMIALS_HPP
#define MATH_LEGENDREPOLYNOMIALS_HPP

#include <vector>
#include <cassert>

#include "math/Polynomial.hpp"

namespace OpenMD {

  /**
   * @class LegendrePolynomial
   * A collection of Legendre Polynomials.
   * @todo document
   */
  class LegendrePolynomial {
  public:
    LegendrePolynomial(int maxPower);
    virtual ~LegendrePolynomial() {}
    /**
     * Calculates the value of the nth Legendre Polynomial evaluated at the given x value.
     * @return The value of the nth Legendre Polynomial evaluates at the given x value
     * @param n
     * @param x the value of the independent variable for the nth Legendre Polynomial  function
     */
        
    RealType evaluate(int n, RealType x) {
      assert (n <= maxPower_ && n >=0); 
      return polyList_[n].evaluate(x);
    }

    /**
     * Returns the first derivative of the nth Legendre Polynomial.
     * @return the first derivative of the nth Legendre Polynomial
     * @param n
     * @param x the value of the independent variable for the nth Legendre Polynomial  function
     */
    RealType evaluateDerivative(int n, RealType x) {
      assert (n <= maxPower_ && n >=0); 
      return polyList_[n].evaluateDerivative(x);        
    }

    /**
     * Returns the nth Legendre Polynomial 
     * @return the nth Legendre Polynomial
     * @param n
     */
    const DoublePolynomial& getLegendrePolynomial(int n) const {
      assert (n <= maxPower_ && n >=0); 
      return polyList_[n];
    }

  protected:

    std::vector<DoublePolynomial> polyList_;
                
  private:
        
    void GeneratePolynomials(int maxPower);
    virtual void GenerateFirstTwoTerms();
        
    int maxPower_;
  };    


} 
#endif 

