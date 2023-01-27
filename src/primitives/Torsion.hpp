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
 * @file Torsion.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_TORSION_HPP
#define PRIMITIVES_TORSION_HPP

#include <limits>

#include "primitives/Atom.hpp"
#include "primitives/ShortRangeInteraction.hpp"
#include "types/TorsionType.hpp"

namespace OpenMD {
  struct TorsionData {
    RealType angle;
    RealType potential;
  };

  struct TorsionDataSet {
    RealType deltaV;
    TorsionData prev;
    TorsionData curr;
  };

  /**
   * @class Torsion Torsion.hpp "types/Torsion.hpp"
   */
  class Torsion : public ShortRangeInteraction {
  public:
    using ShortRangeInteraction::getPrevValue;
    using ShortRangeInteraction::getValue;

    Torsion(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4,
            TorsionType* tt);
    virtual ~Torsion() {}
    virtual void calcForce(RealType& angle, bool doParticlePot);

    RealType getValue(int snapshotNo) {
      Vector3d pos1 = atoms_[0]->getPos(snapshotNo);
      Vector3d pos2 = atoms_[1]->getPos(snapshotNo);
      Vector3d pos3 = atoms_[2]->getPos(snapshotNo);
      Vector3d pos4 = atoms_[3]->getPos(snapshotNo);

      Vector3d r21 = pos1 - pos2;
      snapshotMan_->getSnapshot(snapshotNo)->wrapVector(r21);
      Vector3d r32 = pos2 - pos3;
      snapshotMan_->getSnapshot(snapshotNo)->wrapVector(r32);
      Vector3d r43 = pos3 - pos4;
      snapshotMan_->getSnapshot(snapshotNo)->wrapVector(r43);

      //  Calculate the cross products and distances
      Vector3d A  = cross(r21, r32);
      RealType rA = A.length();
      Vector3d B  = cross(r32, r43);
      RealType rB = B.length();

      /*
         If either of the two cross product vectors is tiny, that means
         the three atoms involved are colinear, and the torsion angle is
         going to be undefined.  The easiest check for this problem is
         to use the product of the two lengths.
      */
      if (rA * rB < OpenMD::epsilon) return numeric_limits<double>::quiet_NaN();

      A.normalize();
      B.normalize();

      //  Calculate the sin and cos
      RealType cos_phi = dot(A, B);
      if (cos_phi > 1.0) cos_phi = 1.0;
      if (cos_phi < -1.0) cos_phi = -1.0;
      return acos(cos_phi);
    }

    RealType getPotential() { return potential_; }

    Atom* getAtomA() { return atoms_[0]; }

    Atom* getAtomB() { return atoms_[1]; }

    Atom* getAtomC() { return atoms_[2]; }

    Atom* getAtomD() { return atoms_[3]; }

    TorsionType* getTorsionType() { return torsionType_; }

    virtual std::string getName() { return name_; }
    /** Sets the name of this torsion for selections */
    virtual void setName(const std::string& name) { name_ = name; }

    void accept(BaseVisitor* v) { v->visit(this); }

  protected:
    TorsionType* torsionType_;
    std::string name_;

    RealType potential_;
  };
}  // namespace OpenMD

#endif  // PRIMITIVES_TORSION_HPP
