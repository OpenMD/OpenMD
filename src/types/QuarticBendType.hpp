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
 * @file QuarticBendType.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */ 

#ifndef TYPES_QUARTICBENDTYPE_HPP
#define TYPES_QUARTICBENDTYPE_HPP

#include "types/BendType.hpp"

namespace OpenMD {
  /**
   * @class QuarticBendType 
   * @todo document
   */
  class QuarticBendType : public BendType {

  public:
        
    QuarticBendType(RealType r0, RealType k4, RealType k3, RealType k2, 
                    RealType k1, RealType k0) : BendType(r0), k4_(k4), k3_(k3),
                                                k2_(k2),  k1_(k1), k0_(k0){
    }
    
    void setForceConstant(RealType k4, RealType k3, RealType k2, RealType k1, 
                          RealType k0) {
      k4_ = k4;
      k3_ = k3;
      k2_ = k2;
      k1_ = k1;
      k0_ = k0;
    }

    void getForceConstant(RealType& k4, RealType& k3, RealType& k2, 
                          RealType& k1, RealType& k0) {
      k4 = k4_;
      k3 = k3_;
      k2 = k2_;
      k1 = k1_;
      k0 = k0_;
    }

    virtual void calcForce(RealType theta, RealType& V, RealType& dVdTheta) {
      RealType delta =  theta- theta0_;
      RealType delta2 = delta * delta;
      RealType delta3 = delta2 * delta;
      RealType delta4 = delta3 * delta;
            
      V =k0_ + k1_ * delta + k2_*delta2 + k3_*delta3 + k4_*delta4;
      dVdTheta = k1_ + 2.0*k2_ * delta + 3.0 * k3_*delta2 + 4.0*k4_*delta3;
    }     
        
  private:
    RealType k4_;
    RealType k3_;
    RealType k2_;
    RealType k1_;
    RealType k0_;
  };
}//end namespace OpenMD
#endif //TYPES_QUARTICBENDTYPE_HPP
