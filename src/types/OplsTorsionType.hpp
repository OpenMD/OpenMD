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
 * @file OplsTorsionType.hpp
 * @author    teng lin
 * @date  11/16/2004
 * @version 1.0
 */ 

#ifndef TYPES_OPLSTORSIONTYPE_HPP
#define TYPES_OPLSTORSIONTYPE_HPP

#include "types/PolynomialTorsionType.hpp"

namespace OpenMD {

  /**
   * @class OplsTorsionType OplsTorsionType.hpp "types/OplsTorsionType.hpp"
   * These torsion types are defined identically with functional form given
   * in the following paper:
   *
   * "Development and Testing of the OPLS All-Atom Force Field on 
   * Conformational Energetics and Properties of Organic Liquids,"  by
   * William L. Jorgensen, David S. Maxwell, and Julian Tirado-Rives, 
   * J. Am. Chem. Soc.; 1996; 118(45) pp 11225 - 11236;
   * DOI: 10.1021/ja9621760
   *
   * This torsion potential has the form:
   * 
   *  Vtors = 0.5* (v1*(1+cos(phi)) + v2*(1-cos(2*phi)) + v3*(1+cos(3*phi)))
   *
   * Notes: 
   * 
   * 1) OpenMD converts internally to a Polynomial torsion type because
   *    all of the phase angles are zero in the OPLS paper.
   * 2) Coefficients are assumed to be in kcal / mol, and be careful about
   *    that factor of 1/2 when importing the coefficients!
   */
  class OplsTorsionType : public PolynomialTorsionType{
    
  public:
    
    OplsTorsionType(RealType v1, RealType v2, RealType v3) :  
      PolynomialTorsionType(), v1_(v1), v2_(v2), v3_(v3) {
      
      //convert OPLS Torsion Type to Polynomial Torsion type
      RealType c0 = v2 + 0.5 * (v1 + v3);
      RealType c1 = 0.5 * (v1 - 3.0 * v3);
      RealType c2 = -v2;
      RealType c3 = 2.0 * v3;
      c0 = c0/2;
      c1 = c1/2;
      c2 = c2/2;
      c3 = c3/2;
      // I change the parameter to half to see if this is the problem.
      setCoefficient(0, c0);
      setCoefficient(1, c1);
      setCoefficient(2, c2);
      setCoefficient(3, c3);
    }
    
    friend std::ostream& operator <<(std::ostream& os, OplsTorsionType& ott);
    
  private:
    
    RealType v1_;
    RealType v2_;
    RealType v3_;
    
  };
  
  std::ostream& operator <<(std::ostream& os, OplsTorsionType& ott) {
    
    os << "This OplsTorsionType has below form:" << std::endl;
    os << ott.v1_ << "/2*(1+cos(phi))" << " + "
       << ott.v2_ << "/2*(1-cos(2*phi))" << " + "
       << ott.v3_ << "/2*(1+cos(3*phi))" << std::endl;
    return os;
  }
  
  
} //end namespace OpenMD
#endif //TYPES_OPLSTORSIONTYPE_HPP


