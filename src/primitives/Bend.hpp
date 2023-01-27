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
 * @file Bend.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_BEND_HPP
#define PRIMITIVES_BEND_HPP

#include "primitives/Atom.hpp"
#include "primitives/ShortRangeInteraction.hpp"
#include "types/BendType.hpp"

namespace OpenMD {
  struct BendData {
    RealType angle;
    RealType potential;
  };

  struct BendDataSet {
    RealType deltaV;
    BendData prev;
    BendData curr;
  };

  class Bend : public ShortRangeInteraction {
  public:
    using ShortRangeInteraction::getPrevValue;
    using ShortRangeInteraction::getValue;
    Bend(Atom* atom1, Atom* atom2, Atom* atom3, BendType* bt) :
        ShortRangeInteraction(), bendType_(bt) {
      atoms_.resize(3);
      atoms_[0] = atom1;
      atoms_[1] = atom2;
      atoms_[2] = atom3;
    }

    virtual ~Bend() {}
    virtual void calcForce(RealType& angle, bool doParticlePot);

    RealType getValue(int snapshotNo) {
      Vector3d pos1 = atoms_[0]->getPos(snapshotNo);
      Vector3d pos2 = atoms_[1]->getPos(snapshotNo);
      Vector3d pos3 = atoms_[2]->getPos(snapshotNo);

      Vector3d r21 = pos1 - pos2;
      snapshotMan_->getSnapshot(snapshotNo)->wrapVector(r21);
      RealType d21 = r21.length();

      Vector3d r23 = pos3 - pos2;
      snapshotMan_->getSnapshot(snapshotNo)->wrapVector(r23);
      RealType d23 = r23.length();

      RealType cosTheta = dot(r21, r23) / (d21 * d23);

      // check roundoff
      if (cosTheta > 1.0) {
        cosTheta = 1.0;
      } else if (cosTheta < -1.0) {
        cosTheta = -1.0;
      }

      return acos(cosTheta);
    }

    RealType getPotential() { return potential_; }

    Atom* getAtomA() { return atoms_[0]; }

    Atom* getAtomB() { return atoms_[1]; }

    Atom* getAtomC() { return atoms_[2]; }

    BendType* getBendType() { return bendType_; }

    virtual std::string getName() { return name_; }
    /** Sets the name of this bend for selections */
    virtual void setName(const std::string& name) { name_ = name; }

    void accept(BaseVisitor* v) { v->visit(this); }

  protected:
    RealType potential_ {};
    BendType* bendType_; /**< bend type */
    std::string name_;
  };
}  // namespace OpenMD

#endif  // PRIMITIVES_BEND_HPP
