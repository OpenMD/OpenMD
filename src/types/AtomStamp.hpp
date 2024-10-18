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

#ifndef TYPES_ATOMSTAMP_HPP
#define TYPES_ATOMSTAMP_HPP

#include <set>
#include <vector>

#include "types/DataHolder.hpp"

namespace OpenMD {

  class AtomStamp : public DataHolder {
    DeclareParameter(Type, std::string);

  public:
    AtomStamp(int index);

  public:
    bool setPosition(const std::vector<RealType>& pos);
    bool setOrientation(const std::vector<RealType>& ort);
    bool havePosition() { return havePos_; }
    bool haveOrientation() { return haveOrt_; }
    RealType getPosX() { return position_[0]; }
    RealType getPosY() { return position_[1]; }
    RealType getPosZ() { return position_[2]; }
    RealType getEulerPhi() { return orientation_[0]; }
    RealType getEulerTheta() { return orientation_[1]; }
    RealType getEulerPsi() { return orientation_[2]; }
    int getIndex() { return index_; }
    virtual void validate();
    using AtomIter = std::set<int>::iterator;
    using BondIter = std::vector<int>::iterator;
    int getFirstBondedAtom(AtomIter& ai) {
      ai = bondedAtoms_.begin();
      return ai != bondedAtoms_.end() ? *ai : -1;
    }
    int getNextBondedAtom(AtomIter& ai) {
      ++ai;
      return ai != bondedAtoms_.end() ? *ai : -1;
    }
    int getFirstBond(BondIter& bi) {
      bi = bonds_.begin();
      return bi != bonds_.end() ? *bi : -1;
    }
    int getNextBond(BondIter& bi) {
      ++bi;
      return bi != bonds_.end() ? *bi : -1;
    }
    void addBond(int bondIndex) { bonds_.push_back(bondIndex); }
    void addBondedAtom(int atomIndex) { bondedAtoms_.insert(atomIndex); }
    size_t getCoordination() { return bonds_.size(); }
    void overrideCharge(RealType c) {
      orCharge_    = c;
      hasOverride_ = true;
    }

    bool hasOverride() { return hasOverride_; }
    RealType getOverrideCharge() { return orCharge_; }

  private:
    Vector3d position_;
    Vector3d orientation_;
    RealType orCharge_;
    bool havePos_;
    bool haveOrt_;
    bool hasOverride_;
    int index_;
    std::vector<int> bonds_;
    std::set<int> bondedAtoms_;
  };
}  // namespace OpenMD

#endif
