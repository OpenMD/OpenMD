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

#ifndef PRIMITIVES_CUTOFFGROUP_HPP

#define PRIMITIVES_CUTOFFGROUP_HPP

#include "math/Vector3.hpp"
#include "primitives/Atom.hpp"

namespace OpenMD {
  class CutoffGroup {
  public:
    CutoffGroup() : globalIndex(-1), localIndex_(-1), snapshotMan_(NULL) {
      storage_      = &Snapshot::cgData;
      haveTotalMass = false;
      totalMass     = 0.0;
    }

    /**
     * Sets the Snapshot Manager of this cutoffGroup
     */
    void setSnapshotManager(SnapshotManager* sman) { snapshotMan_ = sman; }

    void addAtom(Atom* atom) { cutoffAtomList.push_back(atom); }

    Atom* beginAtom(std::vector<Atom*>::iterator& i) {
      i = cutoffAtomList.begin();
      return i != cutoffAtomList.end() ? *i : NULL;
    }

    Atom* nextAtom(std::vector<Atom*>::iterator& i) {
      i++;
      return i != cutoffAtomList.end() ? *i : NULL;
    }

    std::vector<Atom*> getAtoms() { return cutoffAtomList; }
    RealType getMass() {
      if (!haveTotalMass) {
        totalMass = 0.0;

        std::vector<Atom*>::iterator i;
        for (Atom* atom = beginAtom(i); atom != NULL; atom = nextAtom(i)) {
          RealType mass = atom->getMass();
          totalMass += mass;
        }

        haveTotalMass = true;
      }

      return totalMass;
    }

    void updateCOM() {
      DataStorage& data = snapshotMan_->getCurrentSnapshot()->*storage_;
      bool needsVel     = false;
      if (data.getStorageLayout() & DataStorage::dslVelocity) needsVel = true;

      if (cutoffAtomList.size() == 1) {
        data.position[localIndex_] = cutoffAtomList[0]->getPos();
        if (needsVel) data.velocity[localIndex_] = cutoffAtomList[0]->getVel();
      } else {
        std::vector<Atom*>::iterator i;
        Atom* atom;
        RealType totalMass         = getMass();
        data.position[localIndex_] = V3Zero;
        if (needsVel) data.velocity[localIndex_] = V3Zero;

        for (atom = beginAtom(i); atom != NULL; atom = nextAtom(i)) {
          data.position[localIndex_] += atom->getMass() * atom->getPos();
          if (needsVel)
            data.velocity[localIndex_] += atom->getMass() * atom->getVel();
        }
        data.position[localIndex_] /= totalMass;
        if (needsVel) data.velocity[localIndex_] /= totalMass;
      }
    }

    Vector3d getPos() {
      return ((snapshotMan_->getCurrentSnapshot())->*storage_)
          .position[localIndex_];
    }
    Vector3d getVel() {
      return ((snapshotMan_->getCurrentSnapshot())->*storage_)
          .velocity[localIndex_];
    }

    size_t getNumAtom() { return cutoffAtomList.size(); }

    int getGlobalIndex() { return globalIndex; }

    void setGlobalIndex(int id) { this->globalIndex = id; }

    /**
     * Returns the local index of this cutoffGroup
     * @return the local index of this cutoffGroup
     */
    int getLocalIndex() { return localIndex_; }

    /**
     * Sets the local index of this cutoffGroup
     * @param index new index to be set
     */
    void setLocalIndex(int index) { localIndex_ = index; }

  private:
    std::vector<Atom*> cutoffAtomList;
    bool haveTotalMass;
    RealType totalMass;
    int globalIndex;

    int localIndex_;
    DataStoragePointer storage_;
    SnapshotManager* snapshotMan_;
  };
}  // namespace OpenMD

#endif  // PRIMITIVES_CUTOFFGROUP_HPP
