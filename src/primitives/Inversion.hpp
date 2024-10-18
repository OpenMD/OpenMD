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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

/**
 * @file Inversion.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_INVERSION_HPP
#define PRIMITIVES_INVERSION_HPP

#include "primitives/Atom.hpp"
#include "primitives/ShortRangeInteraction.hpp"
#include "types/InversionType.hpp"

namespace OpenMD {
  struct InversionData {
    RealType angle;
    RealType potential;
  };

  struct InversionDataSet {
    RealType deltaV;
    InversionData prev;
    InversionData curr;
  };

  /**
   * @class Inversion Inversion.hpp "primitives/Inversion.hpp"
   */
  class Inversion : public ShortRangeInteraction {
  public:
    using ShortRangeInteraction::getPrevValue;
    using ShortRangeInteraction::getValue;

    Inversion(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4,
              InversionType* it);
    virtual ~Inversion() {}
    virtual void calcForce(RealType& angle, bool doParticlePot);

    RealType getValue(int snapshotNo) {
      // In OpenMD's version of an inversion, the central atom
      // comes first.  However, to get the planarity in a typical cosine
      // version of this potential (i.e. Amber-style), the central atom
      // is treated as atom *3* in a standard torsion form:

      Vector3d pos1 = atoms_[1]->getPos(snapshotNo);
      Vector3d pos2 = atoms_[2]->getPos(snapshotNo);
      Vector3d pos3 = atoms_[0]->getPos(snapshotNo);
      Vector3d pos4 = atoms_[3]->getPos(snapshotNo);

      Vector3d r31 = pos1 - pos3;
      snapshotMan_->getSnapshot(snapshotNo)->wrapVector(r31);
      Vector3d r23 = pos3 - pos2;
      snapshotMan_->getSnapshot(snapshotNo)->wrapVector(r23);
      Vector3d r43 = pos3 - pos4;
      snapshotMan_->getSnapshot(snapshotNo)->wrapVector(r43);

      //  Calculate the cross products and distances
      Vector3d A = cross(r31, r43);
      Vector3d B = cross(r43, r23);

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

    InversionType* getInversionType() { return inversionType_; }
    virtual std::string getName() { return name_; }
    /** Sets the name of this inversion for selections */
    virtual void setName(const std::string& name) { name_ = name; }

    void accept(BaseVisitor* v) { v->visit(this); }

  protected:
    InversionType* inversionType_;
    InversionKey inversionKey_;
    std::string name_;

    RealType potential_;
  };
}  // namespace OpenMD

#endif  // PRIMITIVES_INVERSION_HPP
