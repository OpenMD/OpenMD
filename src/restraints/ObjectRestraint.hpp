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

#ifndef RESTRAINTS_OBJECTRESTRAINT_HPP
#define RESTRAINTS_OBJECTRESTRAINT_HPP

#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "restraints/Restraint.hpp"

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
    ObjectRestraint() : Restraint() {}

    void setReferenceStructure(Vector3d refPos) { refPos_ = refPos; }

    void setReferenceStructure(Vector3d refPos, RotMat3x3d refA) {
      refPos_ = refPos;
      refA_   = refA;
    }

    Vector3d getReferenceStructure() { return refPos_; }

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
}  // namespace OpenMD

#endif
