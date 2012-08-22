/*
 * Copyright (c) 2009 The University of Notre Dame. All Rights Reserved.
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
 
#ifndef RESTRAINTS_OBJECTRESTRAINT_HPP
#define RESTRAINTS_OBJECTRESTRAINT_HPP

#include "restraints/Restraint.hpp"
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"

namespace OpenMD {
  /**
   * @class ObjectRestraint 
   *
   * ObjectRestraint is the basic harmonic restraint for the 
   * degrees of freedom of a StuntDouble 
   *
   * In the ideal structure:
   * 
   * k_[twist,swing] are the two spring constants of the restraining 
   * potential
   */
  class ObjectRestraint : public Restraint {
    
  public:
    
    ObjectRestraint() : Restraint() { }

    void setReferenceStructure(Vector3d refPos) {
      refPos_ = refPos;
    }
  
    void setReferenceStructure(Vector3d refPos, RotMat3x3d refA) {
      refPos_ = refPos;
      refA_ = refA;
    }

    Vector3d getReferenceStructure(){
      return refPos_;
    }

    
    void calcForce(Vector3d struc);
    void calcForce(Vector3d struc, RotMat3x3d A);    

    Vector3d getRestraintForce() { return force_; }
    Vector3d getRestraintTorque() { return torque_; }
                    
  private:    
    Vector3d refPos_;
    RotMat3x3d refA_;

    Vector3d force_;
    Vector3d torque_;
  };
}
#endif
