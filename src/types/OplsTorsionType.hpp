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
 * @file OplsTorsionType.hpp
 * @author    teng lin
 * @date  11/16/2004
 * @version 1.0
 */ 

#ifndef TYPES_POLYNOMIALBONDTYPE_HPP
#define TYPES_POLYNOMIALBONDTYPE_HPP

#include "math/Polynomial.hpp"
#include "types/TorsionType.hpp"

namespace oopse {

  /**
   * @class OplsTorsionType OplsTorsionType.hpp "types/OplsTorsionType.hpp"
   * @todo documentation
   */
  class OplsTorsionType : public PolynomialTorsionType{
    
  public:
    
    OplsTorsionType(RealType v0, RealType v1, RealType v2, RealType v3) :  
      PolynomialTorsionType(){
      
      //convert OPLS Torsion Type to Polynomial Torsion type
      RealType c0 = v0 + v2 + 0.5*(v1 + v3);
      RealType c1 = 0.5 *(3*v3- v1);
      RealType c2 = -v2;
      RealType c3 = -2.0* v3;
      
      setCoefficient(0, c0);
      setCoefficient(1, c1);
      setCoefficient(2, c2);
      setCoefficient(3, c3);
    }
    
    friend std::ostream& operator <<(std::ostream& os, OplsTorsionType& ott);
    
  private:
    
    RealType v0_;
    RealType v1_;
    RealType v2_;
    RealType v3_;
    
  };
  
  std::ostream& operator <<(std::ostream& os, OplsTorsionType& ott) {
    
    os << "This OplsTorsionType has below form:" << std::endl;
    os << ott.v0_ << " + " 
       << ott.v1_ << "/2*(1+cos(Omega))" << " + "
       << ott.v2_ << "/2*(1-cos(2*Omega))" << " + "
       << ott.v3_ << "/2*(1+cos(3*Omega))" << std::endl;
    return os;
  }
  
  
} //end namespace oopse
#endif //TYPES_POLYNOMIALBONDTYPE_HPP


