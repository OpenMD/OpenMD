/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

/**
 * @file Restraint.hpp
 * @author   cli2
 * @date  06/17/2009
 * @version 1.0
 */

#ifndef RESTRAINTS_RESTRAINT_HPP
#define RESTRAINTS_RESTRAINT_HPP

#include <config.h>

#include <map>

#include "math/Vector3.hpp"
#include "utils/GenericData.hpp"

namespace OpenMD {

  class Restraint {
  public:
    enum {
      rtDisplacement = 1,
      rtAbsoluteZ    = 2,
      rtTwist        = 4,
      rtSwingX       = 8,
      rtSwingY       = 16
    };

    using RealPair = std::pair<RealType, RealType>;

    Restraint() :
        twist0_(0.0), swingX0_(0.0), swingY0_(0.0), posZ0_(0.0),
        printRest_(false), restType_(0) {}

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
    int getRestraintType() { return restType_; }
    /** Sets the restraint type  */
    void setRestraintType(int restType) { restType_ = restType; }

    void setScaleFactor(RealType sf) { scaleFactor_ = sf; }

    void setDisplacementForceConstant(RealType kDisp) {
      kDisp_ = kDisp;
      restType_ |= rtDisplacement;
      if (printRest_) restInfo_[rtDisplacement] = std::make_pair(0.0, 0.0);
    }

    void setAbsoluteForceConstant(RealType kAbs) {
      kAbs_ = kAbs;
      restType_ |= rtAbsoluteZ;
      if (printRest_) restInfo_[rtAbsoluteZ] = std::make_pair(0.0, 0.0);
    }

    void setTwistForceConstant(RealType kTwist) {
      kTwist_ = kTwist / 4;
      restType_ |= rtTwist;
      if (printRest_) restInfo_[rtTwist] = std::make_pair(0.0, 0.0);
    }

    void setSwingXForceConstant(RealType kSwingX) {
      kSwingX_ = kSwingX;
      restType_ |= rtSwingX;
      if (printRest_) restInfo_[rtSwingX] = std::make_pair(0.0, 0.0);
    }

    void setSwingYForceConstant(RealType kSwingY) {
      kSwingY_ = kSwingY;
      restType_ |= rtSwingY;
      if (printRest_) restInfo_[rtSwingY] = std::make_pair(0.0, 0.0);
    }

    void setAbsolutePositionZ(RealType z0) {
      posZ0_ = z0;
      restType_ |= rtAbsoluteZ;
      if (printRest_) restInfo_[rtAbsoluteZ] = std::make_pair(0.0, 0.0);
    }

    /* restraint angles are measured relative to the ideal structure,
       and are measured in radians.  If you want to restrain to the
       same structure as the ideal structure, these do not need to be set.
    */
    void setRestrainedTwistAngle(RealType twist0) {
      twist0_ = twist0;
      restType_ |= rtTwist;
      if (printRest_) restInfo_[rtTwist] = std::make_pair(0.0, 0.0);
    }

    void setRestrainedSwingXAngle(RealType swingX0) {
      swingX0_ = swingX0;
      restType_ |= rtSwingX;
      if (printRest_) restInfo_[rtSwingX] = std::make_pair(0.0, 0.0);
    }

    void setRestrainedSwingYAngle(RealType swingY0) {
      swingY0_ = swingY0;
      restType_ |= rtSwingY;
      if (printRest_) restInfo_[rtSwingY] = std::make_pair(0.0, 0.0);
    }

    void setPrintRestraint(bool printRest) { printRest_ = printRest; }

    RealType getDisplacementForceConstant() { return kDisp_; }
    RealType getAbsoluteForceConstant() { return kAbs_; }
    RealType getAbsolutePositionZ() { return posZ0_; }
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
    RealType kAbs_;
    RealType kTwist_;
    RealType kSwingX_;
    RealType kSwingY_;
    RealType pot_;
    RealType twist0_;
    RealType swingX0_;
    RealType swingY0_;
    RealType posZ0_;
    bool printRest_;

    int restType_;
    std::string restName_;
    std::map<int, RealPair> restInfo_;
  };

  using RestraintData = SimpleTypeData<Restraint*>;
}  // namespace OpenMD

#endif
