/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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
                    RealType k1, RealType k0) :
        BendType(r0),
        k4_(k4), k3_(k3), k2_(k2), k1_(k1), k0_(k0) {}

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
      RealType delta  = theta - theta0_;
      RealType delta2 = delta * delta;
      RealType delta3 = delta2 * delta;
      RealType delta4 = delta3 * delta;

      V = k0_ + k1_ * delta + k2_ * delta2 + k3_ * delta3 + k4_ * delta4;
      dVdTheta =
          k1_ + 2.0 * k2_ * delta + 3.0 * k3_ * delta2 + 4.0 * k4_ * delta3;
    }

  private:
    RealType k4_;
    RealType k3_;
    RealType k2_;
    RealType k1_;
    RealType k0_;
  };
}  // namespace OpenMD

#endif  // TYPES_QUARTICBENDTYPE_HPP
