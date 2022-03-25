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
 * @file Bond.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_BOND_HPP
#define PRIMITIVES_BOND_HPP

#include "primitives/Atom.hpp"
#include "primitives/ShortRangeInteraction.hpp"
#include "types/BondType.hpp"

namespace OpenMD {

  class Bond : public ShortRangeInteraction {
  public:
    using ShortRangeInteraction::getPrevValue;
    using ShortRangeInteraction::getValue;
    Bond(Atom* atom1, Atom* atom2, BondType* bt) :
        ShortRangeInteraction(), bondType_(bt) {
      atoms_.resize(2);
      atoms_[0] = atom1;
      atoms_[1] = atom2;
    }
    virtual ~Bond() {}
    void calcForce(bool doParticlePot) {
      RealType len;
      RealType dvdr;
      Vector3d r12;
      Vector3d force;

      r12 = atoms_[1]->getPos() - atoms_[0]->getPos();
      snapshotMan_->getCurrentSnapshot()->wrapVector(r12);
      len = r12.length();
      bondType_->calcForce(len, potential_, dvdr);

      force = r12 * (-dvdr / len);

      atoms_[0]->addFrc(-force);
      atoms_[1]->addFrc(force);
      if (doParticlePot) {
        atoms_[0]->addParticlePot(potential_);
        atoms_[1]->addParticlePot(potential_);
      }
    }

    RealType getValue(int snap) {
      Vector3d r12 = atoms_[1]->getPos(snap) - atoms_[0]->getPos(snap);
      snapshotMan_->getSnapshot(snap)->wrapVector(r12);
      return r12.length();
    }

    RealType getPotential() { return potential_; }

    Atom* getAtomA() { return atoms_[0]; }

    Atom* getAtomB() { return atoms_[1]; }

    BondType* getBondType() { return bondType_; }

    virtual std::string getName() { return name_; }
    /** Sets the name of this bond for selections */
    virtual void setName(const std::string& name) { name_ = name; }

    void accept(BaseVisitor* v) { v->visit(this); }

  private:
    RealType potential_;
    BondType* bondType_; /**< bond type */
    std::string name_;
  };
}  // namespace OpenMD

#endif  // PRIMITIVES_BOND_HPP
