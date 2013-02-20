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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file Restraint.hpp
 * @author   cli2
 * @date  06/17/2009
 * @version 1.0
 */ 

#ifndef RESTRAINTS_RESTRAINT_HPP
#define RESTRAINTS_RESTRAINT_HPP

#include "config.h"
#include "utils/GenericData.hpp"
#include "math/Vector3.hpp"
#include <map>

namespace OpenMD {
  
  class Restraint {
  public:

    enum {
      rtDisplacement = 1,
      rtTwist = 2,
      rtSwingX = 4,
      rtSwingY = 8
    };

    typedef std::pair<RealType, RealType> RealPair;

    Restraint() : twist0_(0.0), swingX0_(0.0), swingY0_(0.0), restType_(0) {
    }
    
    virtual ~Restraint() {}

    // these are place-holders.  The subclasses will have different arguments
    // to the two functions.
    void calcForce() {}
    void setReferenceStructure() {}

    RealType getUnscaledPotential() { return pot_; }
    RealType getPotential() { return scaleFactor_ * pot_; }

    void setRestraintName(std::string name) { restName_ = name; }
    std::string getRestraintName() { return restName_; }

    /** Returns the restraint type  */
    int getRestraintType(){ return restType_; }
    /** Sets the restraint type  */
    void setRestraintType(int restType) { restType_ = restType; }

    void setScaleFactor(RealType sf) { scaleFactor_ = sf;}

    void setDisplacementForceConstant(RealType kDisp) { 
      kDisp_ = kDisp; 
      restType_ |= rtDisplacement;
      restInfo_[rtDisplacement] = std::make_pair(0.0, 0.0);
    }
    
    void setTwistForceConstant(RealType kTwist) { 
      kTwist_ = kTwist/4;
      restType_ |= rtTwist; 
      restInfo_[rtTwist] = std::make_pair(0.0, 0.0);
    }
    
    void setSwingXForceConstant(RealType kSwingX) { 
      kSwingX_ = kSwingX; 
      restType_ |= rtSwingX; 
      restInfo_[rtSwingX] = std::make_pair(0.0, 0.0);
    }

    void setSwingYForceConstant(RealType kSwingY) { 
      kSwingY_ = kSwingY; 
      restType_ |= rtSwingY; 
      restInfo_[rtSwingY] = std::make_pair(0.0, 0.0);
    }
    
    /* restraint angles are measured relative to the ideal structure, 
       and are measured in radians.  If you want to restrain to the 
       same structure as the ideal structure, these do not need to be set.
    */    
    void setRestrainedTwistAngle(RealType twist0) { 
      twist0_ = twist0; 
      restType_ |= rtTwist; 
      restInfo_[rtTwist] = std::make_pair(0.0, 0.0);
    }
    
    void setRestrainedSwingXAngle(RealType swingX0) { 
      swingX0_ = swingX0; 
      restType_ |= rtSwingX; 
      restInfo_[rtSwingX] = std::make_pair(0.0, 0.0);
    }

    void setRestrainedSwingYAngle(RealType swingY0) { 
      swingY0_ = swingY0; 
      restType_ |= rtSwingY; 
      restInfo_[rtSwingY] = std::make_pair(0.0, 0.0);
    }
    
    void setPrintRestraint(bool printRest) {
      printRest_ = printRest;
    }

    RealType getDisplacementForceConstant() { return kDisp_; }
    RealType getTwistForceConstant() { return kTwist_; }
    RealType getSwingXForceConstant() { return kSwingX_; }
    RealType getSwingYForceConstant() { return kSwingY_; }
    RealType getRestrainedTwistAngle() { return twist0_; }
    RealType getRestrainedSwingXAngle() { return swingX0_; }
    RealType getRestrainedSwingYAngle() { return swingY0_; }
    std::map<int, RealPair> getRestraintInfo() { return restInfo_; } 
    bool getPrintRestraint() { return printRest_; }

  protected:

    RealType scaleFactor_;
    RealType kDisp_;
    RealType kTwist_;
    RealType kSwingX_;
    RealType kSwingY_;
    RealType pot_;    
    RealType twist0_;
    RealType swingX0_;
    RealType swingY0_;
    bool printRest_;
    
    int restType_;
    std::string restName_;
    std::map<int, RealPair> restInfo_;
  };    
 
  typedef SimpleTypeData<Restraint*> RestraintData;   

  
}
#endif
