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

#ifndef TYPES_MOLECULESTAMP_HPP
#define TYPES_MOLECULESTAMP_HPP

#include <utility>
#include <vector>

#include "types/AtomStamp.hpp"
#include "types/BendStamp.hpp"
#include "types/BondStamp.hpp"
#include "types/ConstraintStamp.hpp"
#include "types/CutoffGroupStamp.hpp"
#include "types/FragmentStamp.hpp"
#include "types/InversionStamp.hpp"
#include "types/RigidBodyStamp.hpp"
#include "types/TorsionStamp.hpp"

namespace OpenMD {
  class MoleculeStamp : public DataHolder {
    DeclareParameter(Name, std::string);
    DeclareParameter(ConstrainTotalCharge, bool);

  public:
    MoleculeStamp();
    virtual ~MoleculeStamp();

    bool addAtomStamp(AtomStamp* atom);
    bool addBondStamp(BondStamp* bond);
    bool addBendStamp(BendStamp* bend);
    bool addTorsionStamp(TorsionStamp* torsion);
    bool addInversionStamp(InversionStamp* inversion);
    bool addRigidBodyStamp(RigidBodyStamp* rigidbody);
    bool addCutoffGroupStamp(CutoffGroupStamp* cutoffgroup);
    bool addFragmentStamp(FragmentStamp* fragment);
    bool addConstraintStamp(ConstraintStamp* constraint);

    void setIdent(int id) { ident_ = id; }
    int getIdent() { return ident_; }

    void setRegion(int r) { region_ = r; }
    int getRegion() { return region_; }

    std::size_t getNAtoms() { return atomStamps_.size(); }
    std::size_t getNBonds() { return bondStamps_.size(); }
    std::size_t getNBends() { return bendStamps_.size(); }
    std::size_t getNTorsions() { return torsionStamps_.size(); }
    std::size_t getNInversions() { return inversionStamps_.size(); }
    std::size_t getNRigidBodies() { return rigidBodyStamps_.size(); }
    std::size_t getNCutoffGroups() { return cutoffGroupStamps_.size(); }
    std::size_t getNConstraints() { return constraintStamps_.size(); }
    int getNIntegrable() { return nintegrable_; }
    std::size_t getNFreeAtoms() { return freeAtoms_.size(); }
    virtual void validate();

    AtomStamp* getAtomStamp(int index) { return atomStamps_[index]; }
    BondStamp* getBondStamp(int index) { return bondStamps_[index]; }
    BendStamp* getBendStamp(int index) { return bendStamps_[index]; }
    TorsionStamp* getTorsionStamp(int index) { return torsionStamps_[index]; }
    InversionStamp* getInversionStamp(int index) {
      return inversionStamps_[index];
    }
    RigidBodyStamp* getRigidBodyStamp(int index) {
      return rigidBodyStamps_[index];
    }
    CutoffGroupStamp* getCutoffGroupStamp(int index) {
      return cutoffGroupStamps_[index];
    }
    FragmentStamp* getFragmentStamp(int index) {
      return fragmentStamps_[index];
    }
    ConstraintStamp* getConstraintStamp(int index) {
      return constraintStamps_[index];
    }

    bool isBondInSameRigidBody(BondStamp* bond);
    bool isAtomInRigidBody(int atomIndex);
    bool isAtomInRigidBody(int atomIndex, int& whichRigidBody,
                           int& consAtomIndex);
    std::vector<std::pair<int, int>> getJointAtoms(int rb1, int rb2);

  private:
    void fillBondInfo();
    void checkAtoms();
    void checkBonds();
    void checkBends();
    void checkTorsions();
    void checkInversions();
    void checkRigidBodies();
    void checkCutoffGroups();
    void checkFragments();
    void checkConstraints();

    template<class Cont, class T>
    bool addIndexSensitiveStamp(Cont& cont, T* stamp) {
      // typename Cont::iterator i;
      unsigned int index = stamp->getIndex();
      bool ret           = false;
      size_t size        = cont.size();

      if (size >= index + 1) {
        if (cont[index] != NULL) {
          ret = false;
        } else {
          cont[index] = stamp;
          ret         = true;
        }
      } else {
        cont.insert(cont.end(), index - cont.size() + 1, NULL);
        cont[index] = stamp;
        ret         = true;
      }

      return ret;
    }

    std::vector<AtomStamp*> atomStamps_;
    std::vector<int> freeAtoms_;
    std::vector<BondStamp*> bondStamps_;
    std::vector<BendStamp*> bendStamps_;
    std::vector<TorsionStamp*> torsionStamps_;
    std::vector<InversionStamp*> inversionStamps_;
    std::vector<RigidBodyStamp*> rigidBodyStamps_;
    std::vector<CutoffGroupStamp*> cutoffGroupStamps_;
    std::vector<FragmentStamp*> fragmentStamps_;
    std::vector<ConstraintStamp*> constraintStamps_;
    std::vector<int> atom2Rigidbody;
    int nintegrable_;
    int region_;
    int ident_;
  };

  class MoleculeStampCompare {
  public:
    bool operator()(MoleculeStamp* lhs, MoleculeStamp* rhs) const {
      return lhs->getIdent() < rhs->getIdent();
    }
  };

  using MoleculeStampSet = std::set<MoleculeStamp*, MoleculeStampCompare>;
}  // namespace OpenMD

#endif
