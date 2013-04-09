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
 
#ifndef TYPES_HARMONICTORSIONTYPE_HPP
#define TYPES_HARMONICTORSIONTYPE_HPP

namespace OpenMD {

  /**
   * @class HarmonicTorsionType
   * These torsion types are defined identically with functional form given
   * in equation 5 in the following paper:
   *
   * "Transferable Potentials for Phase Equilibria. 4. United-Atom
   * Description of Linear and Branched Alkenes and Alkylbenzenes" by
   * Collin D. Wick, Marcus G. Martin and J. Ilja Siepmann,
   * J. Phys. Chem. B; 2000; 104(33) pp 8008 - 8016;
   *
   *    http://pubs.acs.org/doi/abs/10.1021/jp001044x
   *
   * This torsion potential has the form:
   * 
   *  \f[ 
         V_{tors} = \frac{d_0}{2} \left(\phi - \phi_0\right)^2
       \f]
   *
   */
  class HarmonicTorsionType : public TorsionType {
    
  public:
    
    HarmonicTorsionType(RealType d0, RealType phi0) :  
      TorsionType(), d0_(d0), phi0_(phi0) {}


    virtual void calcForce(RealType cosPhi, RealType& V, RealType& dVdCosPhi) {
      //check roundoff     
      if (cosPhi > 1.0) {
        cosPhi = 1.0;
      } else if (cosPhi < -1.0) {
        cosPhi = -1.0;
      }

      RealType phi = acos(cosPhi);
      RealType sinPhi = sqrt(1.0 - cosPhi * cosPhi);

      // trick to avoid divergence in angles near 0 and pi:

      if (fabs(sinPhi) < 1.0E-6) {
        sinPhi = copysign(1.0E-6, sinPhi);
      }

      V = 0.5 * d0_ * pow((phi - phi0_), 2);

      dVdCosPhi = -d0_ * (phi - phi0_) / sinPhi;
    }

    
    friend std::ostream& operator <<(std::ostream& os, HarmonicTorsionType& ttt);
    
  private:
    
    RealType d0_;
    RealType phi0_;
   
  };
  
  std::ostream& operator <<(std::ostream& os, HarmonicTorsionType& htt) {
    
    os << "This HarmonicTorsionType has below form:" << std::endl;
    os << htt.d0_ << "*(phi - " << htt.phi0_ << ")/2" << std::endl; 
    return os;
  }
  
} //end namespace OpenMD
#endif //TYPES_HARMONICTORSIONTYPE_HPP


