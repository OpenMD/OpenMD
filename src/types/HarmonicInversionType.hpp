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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef TYPES_HARMONICINVERSIONTYPE_HPP
#define TYPES_HARMONICINVERSIONTYPE_HPP

#if defined(_MSC_VER)
#define copysign _copysign
#endif

namespace OpenMD {

  /**
   * @class HarmonicInversionType
   * This inversion potential has the form:
   * 
   *  \f[ 
         V_{inv} = \frac{d_0}{2} \left(\phi - \phi_0\right)^2
       \f]
   *
   */
  class HarmonicInversionType : public InversionType {
    
  public:
        
    HarmonicInversionType(RealType d0, RealType phi0) :  
      InversionType(), d0_(d0), phi0_(phi0) {}


    virtual void calcForce(RealType phi, RealType& V, RealType& dVdPhi) {

      V = 0.5 * d0_ * pow((phi - phi0_), 2);
      dVdPhi = -d0_ * (phi - phi0_);

    }
    
    virtual InversionKey getKey() { return itAngle; }
    friend std::ostream& operator <<(std::ostream& os, HarmonicInversionType& ttt);
    
  private:
    
    RealType d0_;
    RealType phi0_;
    
  };
  
  std::ostream& operator <<(std::ostream& os, HarmonicInversionType& hit) {
    
    os << "This HarmonicInversionType has below form:" << std::endl;
    os << hit.d0_ << "*(phi - " << hit.phi0_ << ")/2" << std::endl; 
    return os;
  }
  
} //end namespace OpenMD
#endif //TYPES_HARMONICINVERSIONTYPE_HPP


