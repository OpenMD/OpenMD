/*
 * Copyright (c) 2008 The University of Notre Dame. All Rights Reserved.
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
 * @file TrappeTorsionType.hpp
 * @author    Dan Gezelter
 * @date  06/11/2008
 * @version 1.0
 */ 

#ifndef TYPES_TRAPPETORSIONTYPE_HPP
#define TYPES_TRAPPETORSIONTYPE_HPP

#include "types/PolynomialTorsionType.hpp"

namespace OpenMD {

  /**
   * @class TrappeTorsionType TrappeTorsionType.hpp "types/TrappeTorsionType.hpp"
   * These torsion types are defined identically with functional form given
   * in the following paper:
   *
   * "Transferable Potentials for Phase Equilibria. 1. United-Atom 
   * Description of n-Alkanes," by
   * Marcus G. Martin and J. Ilja Siepmann,
   * J. Phys. Chem. B; 1998; 102(14) pp 2569 - 2577;
   *
   * Also listed as type A torsions on this page: 
   *
   *    http://siepmann6.chem.umn.edu/trappe/intro.php
   *
   * This torsion potential has the form:
   * 
   *  \f[ 
         V_{tors} = c_0 + c_1*(1+\cos(\phi)) + c_2*(1-\cos(2*\phi)) + 
                    c_3*(1+\cos(3*\phi)) 
       \f]
   *
   * Notes: 
   * 
   * 1) This is very similar to the OplsTorsionType with coefficients
   *    defined without the factor of 1/2, and an extra \f$c_0\f$ term.
   *
   * 2) Coefficients are assumed to be in kcal / mol, although the papers
   *    usually give the parameters in units of K.
   */
  class TrappeTorsionType : public PolynomialTorsionType{
    
  public:
    
    TrappeTorsionType(RealType c0, RealType c1, RealType c2, RealType c3) :  
      PolynomialTorsionType(), c0_(c0), c1_(c1), c2_(c2), c3_(c3) {

      //convert Trappe Torsion Type to Polynomial Torsion type

      RealType b0 = c0 + c1 + 2.0 * c2 + c3;
      RealType b1 = c1 - 3.0 * c3;
      RealType b2 = -2.0 * c2;
      RealType b3 = 4.0 * c3;
      
      setCoefficient(0, b0);
      setCoefficient(1, b1);
      setCoefficient(2, b2);
      setCoefficient(3, b3);

    }
    
    friend std::ostream& operator <<(std::ostream& os, TrappeTorsionType& ttt);
    
  private:
    
    RealType c0_;
    RealType c1_;
    RealType c2_;
    RealType c3_;
    
  };
  
  std::ostream& operator <<(std::ostream& os, TrappeTorsionType& ttt) {
    
    os << "This TrappeTorsionType has below form:" << std::endl;
    os << ttt.c0_ << " + " 
       << ttt.c1_ << "*(1+cos(phi))" << " + "
       << ttt.c2_ << "*(1-cos(2*phi))" << " + "
       << ttt.c3_ << "*(1+cos(3*phi))" << std::endl;
    return os;
  }
  
  
} //end namespace OpenMD
#endif //TYPES_TRAPPETORSIONTYPE_HPP


