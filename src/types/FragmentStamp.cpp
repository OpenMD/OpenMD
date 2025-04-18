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
#include "types/FragmentStamp.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <tuple>

#include "utils/MemoryUtils.hpp"

using namespace std;

namespace OpenMD {

  template<class ContainerType>
  bool hasDuplicateElement(const ContainerType& cont) {
    ContainerType tmp = cont;
    std::sort(tmp.begin(), tmp.end());
    tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
    return tmp.size() != cont.size();
  }

  FragmentStamp::FragmentStamp() { DefineParameter(Name, "name"); }

  FragmentStamp::~FragmentStamp() {
    Utils::deletePointers(atomStamps_);
    Utils::deletePointers(bondStamps_);
    Utils::deletePointers(bendStamps_);
    Utils::deletePointers(torsionStamps_);
    Utils::deletePointers(inversionStamps_);
    Utils::deletePointers(rigidBodyStamps_);
    Utils::deletePointers(cutoffGroupStamps_);
    Utils::deletePointers(nodesStamps_);
    Utils::deletePointers(constraintStamps_);
  }

  bool FragmentStamp::addAtomStamp(AtomStamp* atom) {
    bool ret = addIndexSensitiveStamp(atomStamps_, atom);
    if (!ret) {
      std::ostringstream oss;
      oss << "Error in Fragment " << getName()
          << ": multiple atoms have the same indices" << atom->getIndex()
          << "\n";
      throw OpenMDException(oss.str());
    }
    return ret;
  }

  bool FragmentStamp::addBondStamp(BondStamp* bond) {
    bondStamps_.push_back(bond);
    return true;
  }

  bool FragmentStamp::addBendStamp(BendStamp* bend) {
    bendStamps_.push_back(bend);
    return true;
  }

  bool FragmentStamp::addTorsionStamp(TorsionStamp* torsion) {
    torsionStamps_.push_back(torsion);
    return true;
  }
  bool FragmentStamp::addInversionStamp(InversionStamp* inversion) {
    inversionStamps_.push_back(inversion);
    return true;
  }

  bool FragmentStamp::addRigidBodyStamp(RigidBodyStamp* rigidbody) {
    bool ret = addIndexSensitiveStamp(rigidBodyStamps_, rigidbody);
    if (!ret) {
      std::ostringstream oss;
      oss << "Error in Fragment " << getName()
          << ": multiple rigidbodies have the same indices: "
          << rigidbody->getIndex() << "\n";
      throw OpenMDException(oss.str());
    }
    return ret;
  }

  bool FragmentStamp::addCutoffGroupStamp(CutoffGroupStamp* cutoffgroup) {
    cutoffGroupStamps_.push_back(cutoffgroup);
    return true;
  }

  bool FragmentStamp::addNodesStamp(NodesStamp* nodes) {
    nodesStamps_.push_back(nodes);
    return true;
  }

  bool FragmentStamp::addConstraintStamp(ConstraintStamp* constraint) {
    constraintStamps_.push_back(constraint);
    return true;
  }

  void FragmentStamp::validate() {
    DataHolder::validate();
    CheckParameter(Name, isNotEmpty());

    atom2Rigidbody.resize(getNAtoms());

    // A negative number means the atom is a free atom, and does not
    // belong to rigidbody. Every element in atom2Rigidbody has unique
    // negative number at the very beginning

    for (unsigned int i = 0; i < atom2Rigidbody.size(); ++i) {
      atom2Rigidbody[i] = -1 - int(i);
    }
    for (std::size_t i = 0; i < getNRigidBodies(); ++i) {
      RigidBodyStamp* rbStamp  = getRigidBodyStamp(i);
      std::vector<int> members = rbStamp->getMembers();
      for (std::vector<int>::iterator j = members.begin(); j != members.end();
           ++j) {
        atom2Rigidbody[*j] = i;
      }
    }
    checkAtoms();
    checkBonds();
    fillBondInfo();
    checkBends();
    checkTorsions();
    checkInversions();
    checkRigidBodies();
    checkCutoffGroups();
    checkConstraints();
    checkNodes();
  }

  void FragmentStamp::checkAtoms() {
    std::vector<AtomStamp*>::iterator ai = std::find(
        atomStamps_.begin(), atomStamps_.end(), static_cast<AtomStamp*>(NULL));
    if (ai != atomStamps_.end()) {
      std::ostringstream oss;
      oss << "Error in Fragment " << getName() << ": atom["
          << ai - atomStamps_.begin() << "] is missing\n";
      throw OpenMDException(oss.str());
    }
  }

  void FragmentStamp::checkBonds() {
    std::ostringstream oss;
    // make sure index is not out of range
    int natoms = getNAtoms();
    for (std::size_t i = 0; i < getNBonds(); ++i) {
      BondStamp* bondStamp = getBondStamp(i);
      if (bondStamp->getA() > natoms - 1 || bondStamp->getA() < 0 ||
          bondStamp->getB() > natoms - 1 || bondStamp->getB() < 0 ||
          bondStamp->getA() == bondStamp->getB()) {
        oss << "Error in Fragment " << getName() << ": bond("
            << bondStamp->getA() << ", " << bondStamp->getB()
            << ") is invalid\n";
        throw OpenMDException(oss.str());
      }
    }

    // make sure bonds are unique
    std::set<std::pair<int, int>> allBonds;
    for (std::size_t i = 0; i < getNBonds(); ++i) {
      BondStamp* bondStamp = getBondStamp(i);
      std::pair<int, int> bondPair(bondStamp->getA(), bondStamp->getB());
      // make sure bondPair.first is always less than or equal to
      // bondPair.third
      if (bondPair.first > bondPair.second) {
        std::swap(bondPair.first, bondPair.second);
      }

      std::set<std::pair<int, int>>::iterator iter = allBonds.find(bondPair);
      if (iter != allBonds.end()) {
        oss << "Error in Fragment " << getName() << ": "
            << "bond(" << iter->first << ", " << iter->second
            << ") appears multiple times\n";
        throw OpenMDException(oss.str());
      } else {
        allBonds.insert(bondPair);
      }
    }

    // make sure atoms belong to same rigidbody do not bond to each other
    for (std::size_t i = 0; i < getNBonds(); ++i) {
      BondStamp* bondStamp = getBondStamp(i);
      if (atom2Rigidbody[bondStamp->getA()] ==
          atom2Rigidbody[bondStamp->getB()]) {
        oss << "Error in Fragment " << getName() << ": "
            << "bond(" << bondStamp->getA() << ", " << bondStamp->getB()
            << ") belong to same rigidbody "
            << atom2Rigidbody[bondStamp->getA()] << "\n";
        throw OpenMDException(oss.str());
      }
    }
  }

  void FragmentStamp::checkBends() {
    std::ostringstream oss;
    for (std::size_t i = 0; i < getNBends(); ++i) {
      BendStamp* bendStamp         = getBendStamp(i);
      std::vector<int> bendAtoms   = bendStamp->getMembers();
      std::vector<int>::iterator j = std::find_if(
          bendAtoms.begin(), bendAtoms.end(),
          std::bind(std::greater<int>(), placeholders::_1, getNAtoms() - 1));
      std::vector<int>::iterator k =
          std::find_if(bendAtoms.begin(), bendAtoms.end(),
                       std::bind(std::less<int>(), placeholders::_1, 0));

      if (j != bendAtoms.end() || k != bendAtoms.end()) {
        oss << "Error in Fragment " << getName() << " : atoms of bend"
            << containerToString(bendAtoms) << " have invalid indices\n";
        throw OpenMDException(oss.str());
      }

      if (hasDuplicateElement(bendAtoms)) {
        oss << "Error in Fragment " << getName() << " : atoms of bend"
            << containerToString(bendAtoms) << " have duplicated indices\n";
        throw OpenMDException(oss.str());
      }

      if (bendAtoms.size() == 2) {
        if (!bendStamp->haveGhostVectorSource()) {
          oss << "Error in Fragment " << getName()
              << ": ghostVectorSouce is missing\n";
          throw OpenMDException(oss.str());
        } else {
          std::size_t ghostIndex = bendStamp->getGhostVectorSource();
          if (ghostIndex < getNAtoms()) {
            if (std::find(bendAtoms.begin(), bendAtoms.end(), ghostIndex) ==
                bendAtoms.end()) {
              oss << "Error in Fragment " << getName() << ": ghostVectorSouce "
                  << ghostIndex << "is invalid\n";
              throw OpenMDException(oss.str());
            }
            if (!getAtomStamp(ghostIndex)->haveOrientation()) {
              oss << "Error in Fragment " << getName()
                  << ": ghost atom must be a directional atom\n";
              throw OpenMDException(oss.str());
            }
          } else {
            oss << "Error in Fragment " << getName() << ": ghostVectorSource "
                << ghostIndex << "  is invalid\n";
            throw OpenMDException(oss.str());
          }
        }
      } else if (bendAtoms.size() == 3 && bendStamp->haveGhostVectorSource()) {
        oss << "Error in Fragment " << getName()
            << ": normal bend should not have ghostVectorSouce\n";
        throw OpenMDException(oss.str());
      }
    }

    for (std::size_t i = 0; i < getNBends(); ++i) {
      BendStamp* bendStamp       = getBendStamp(i);
      std::vector<int> bendAtoms = bendStamp->getMembers();
      std::vector<int> rigidSet(getNRigidBodies(), 0);
      std::vector<int>::iterator j;
      for (j = bendAtoms.begin(); j != bendAtoms.end(); ++j) {
        int rigidbodyIndex = atom2Rigidbody[*j];
        if (rigidbodyIndex >= 0) {
          ++rigidSet[rigidbodyIndex];
          if (rigidSet[rigidbodyIndex] > 2) {
            oss << "Error in Fragment " << getName() << ": bend"
                << containerToString(bendAtoms)
                << "has three atoms on the same rigid body\n";
            throw OpenMDException(oss.str());
          }
        }
      }
    }

    std::set<std::tuple<int, int, int>> allBends;
    std::set<std::tuple<int, int, int>>::iterator iter;
    for (std::size_t i = 0; i < getNBends(); ++i) {
      BendStamp* bendStamp  = getBendStamp(i);
      std::vector<int> bend = bendStamp->getMembers();
      if (bend.size() == 2) {
        // in case we have two ghost bend. For example,
        // bend {
        // members (0, 1);
        //   ghostVectorSource = 0;
        // }
        // and
        // bend {
        //   members (0, 1);
        // ghostVectorSource = 0;
        // }
        // In order to distinguish them. we expand them to Tuple3.
        // the first one is expanded to (0, 0, 1) while the second one
        // is expaned to (0, 1, 1)
        int ghostIndex = bendStamp->getGhostVectorSource();
        std::vector<int>::iterator j =
            std::find(bend.begin(), bend.end(), ghostIndex);
        if (j != bend.end()) { bend.insert(j, ghostIndex); }
      }

      std::tuple<int, int, int> bendTuple {bend[0], bend[1], bend[2]};
      auto& [first, second, third] = bendTuple;

      // make sure bendTuple.first is always less than or equal to
      // bendTuple.third
      if (first > third) { std::swap(first, third); }

      iter = allBends.find(bendTuple);
      if (iter != allBends.end()) {
        oss << "Error in Fragment " << getName() << ": "
            << "Bend" << containerToString(bend) << " appears multiple times\n";
        throw OpenMDException(oss.str());
      } else {
        allBends.insert(bendTuple);
      }
    }
  }

  void FragmentStamp::checkTorsions() {
    std::ostringstream oss;
    for (std::size_t i = 0; i < getNTorsions(); ++i) {
      TorsionStamp* torsionStamp    = getTorsionStamp(i);
      std::vector<int> torsionAtoms = torsionStamp->getMembers();
      std::vector<int>::iterator j  = std::find_if(
          torsionAtoms.begin(), torsionAtoms.end(),
          std::bind(std::greater<int>(), placeholders::_1, getNAtoms() - 1));
      std::vector<int>::iterator k =
          std::find_if(torsionAtoms.begin(), torsionAtoms.end(),
                       std::bind(std::less<int>(), placeholders::_1, 0));

      if (j != torsionAtoms.end() || k != torsionAtoms.end()) {
        oss << "Error in Fragment " << getName() << ": atoms of torsion"
            << containerToString(torsionAtoms) << " have invalid indices\n";
        throw OpenMDException(oss.str());
      }
      if (hasDuplicateElement(torsionAtoms)) {
        oss << "Error in Fragment " << getName() << " : atoms of torsion"
            << containerToString(torsionAtoms) << " have duplicated indices\n";
        throw OpenMDException(oss.str());
      }
    }

    for (std::size_t i = 0; i < getNTorsions(); ++i) {
      TorsionStamp* torsionStamp    = getTorsionStamp(i);
      std::vector<int> torsionAtoms = torsionStamp->getMembers();
      std::vector<int> rigidSet(getNRigidBodies(), 0);
      std::vector<int>::iterator j;
      for (j = torsionAtoms.begin(); j != torsionAtoms.end(); ++j) {
        int rigidbodyIndex = atom2Rigidbody[*j];
        if (rigidbodyIndex >= 0) {
          ++rigidSet[rigidbodyIndex];
          if (rigidSet[rigidbodyIndex] > 3) {
            oss << "Error in Fragment " << getName() << ": torsion"
                << containerToString(torsionAtoms)
                << "has four atoms on the same rigid body\n";
            throw OpenMDException(oss.str());
          }
        }
      }
    }

    std::set<std::tuple<int, int, int, int>> allTorsions;
    std::set<std::tuple<int, int, int, int>>::iterator iter;
    for (std::size_t i = 0; i < getNTorsions(); ++i) {
      TorsionStamp* torsionStamp = getTorsionStamp(i);
      std::vector<int> torsion   = torsionStamp->getMembers();
      if (torsion.size() == 3) {
        int ghostIndex = torsionStamp->getGhostVectorSource();
        std::vector<int>::iterator j =
            std::find(torsion.begin(), torsion.end(), ghostIndex);
        if (j != torsion.end()) { torsion.insert(j, ghostIndex); }
      }

      std::tuple<int, int, int, int> torsionTuple(torsion[0], torsion[1],
                                                  torsion[2], torsion[3]);
      auto& [first, second, third, fourth] = torsionTuple;

      if (first > fourth) {
        std::swap(first, fourth);
        std::swap(second, third);
      }

      iter = allTorsions.find(torsionTuple);
      if (iter == allTorsions.end()) {
        allTorsions.insert(torsionTuple);
      } else {
        oss << "Error in Fragment " << getName() << ": "
            << "Torsion" << containerToString(torsion)
            << " appears multiple times\n";
        throw OpenMDException(oss.str());
      }
    }
  }

  void FragmentStamp::checkInversions() {
    std::ostringstream oss;

    // first we automatically find the other three atoms that
    // are satellites of an inversion center:

    for (std::size_t i = 0; i < getNInversions(); ++i) {
      InversionStamp* inversionStamp = getInversionStamp(i);
      int center                     = inversionStamp->getCenter();
      std::vector<int> satellites;

      // Some inversions come pre-programmed with the satellites.  If
      // so, don't add the satellites as they are already there.

      if (inversionStamp->getNSatellites() != 3) {
        for (std::size_t j = 0; j < getNBonds(); ++j) {
          BondStamp* bondStamp = getBondStamp(j);
          int a                = bondStamp->getA();
          int b                = bondStamp->getB();

          if (a == center) { satellites.push_back(b); }
          if (b == center) { satellites.push_back(a); }
        }

        if (satellites.size() == 3) {
          std::sort(satellites.begin(), satellites.end());
          inversionStamp->setSatellites(satellites);
        } else {
          oss << "Error in Fragment " << getName() << ": found wrong number"
              << " of bonds for inversion center " << center;
          throw OpenMDException(oss.str());
        }
      }
    }

    // then we do some sanity checking on the inversions:

    for (std::size_t i = 0; i < getNInversions(); ++i) {
      InversionStamp* inversionStamp = getInversionStamp(i);

      std::vector<int> inversionAtoms = inversionStamp->getSatellites();
      // add the central atom to the beginning of the list:
      inversionAtoms.insert(inversionAtoms.begin(),
                            inversionStamp->getCenter());

      std::vector<int>::iterator j = std::find_if(
          inversionAtoms.begin(), inversionAtoms.end(),
          std::bind(std::greater<int>(), placeholders::_1, getNAtoms() - 1));
      std::vector<int>::iterator k =
          std::find_if(inversionAtoms.begin(), inversionAtoms.end(),
                       std::bind(std::less<int>(), placeholders::_1, 0));

      if (j != inversionAtoms.end() || k != inversionAtoms.end()) {
        oss << "Error in Fragment " << getName() << ": atoms of inversion"
            << containerToString(inversionAtoms) << " have invalid indices\n";
        throw OpenMDException(oss.str());
      }

      if (hasDuplicateElement(inversionAtoms)) {
        oss << "Error in Fragment " << getName() << " : atoms of inversion"
            << containerToString(inversionAtoms)
            << " have duplicated indices\n";
        throw OpenMDException(oss.str());
      }
    }

    for (std::size_t i = 0; i < getNInversions(); ++i) {
      InversionStamp* inversionStamp  = getInversionStamp(i);
      std::vector<int> inversionAtoms = inversionStamp->getSatellites();
      inversionAtoms.push_back(inversionStamp->getCenter());
      std::vector<int> rigidSet(getNRigidBodies(), 0);
      std::vector<int>::iterator j;
      for (j = inversionAtoms.begin(); j != inversionAtoms.end(); ++j) {
        int rigidbodyIndex = atom2Rigidbody[*j];
        if (rigidbodyIndex >= 0) {
          ++rigidSet[rigidbodyIndex];
          if (rigidSet[rigidbodyIndex] > 3) {
            oss << "Error in Fragment " << getName()
                << ": inversion centered on atom "
                << inversionStamp->getCenter()
                << " has four atoms that belong to same rigidbody "
                << rigidbodyIndex << "\n";
            throw OpenMDException(oss.str());
          }
        }
      }
    }

    std::set<std::tuple<int, int, int, int>> allInversions;
    std::set<std::tuple<int, int, int, int>>::iterator iter;
    for (std::size_t i = 0; i < getNInversions(); ++i) {
      InversionStamp* inversionStamp = getInversionStamp(i);
      int cent                       = inversionStamp->getCenter();
      std::vector<int> inversion     = inversionStamp->getSatellites();

      std::tuple<int, int, int, int> inversionTuple(cent, inversion[0],
                                                    inversion[1], inversion[2]);
      auto& [first, second, third, fourth] = inversionTuple;

      // In OpenMD, the Central atom in an inversion comes first, and
      // has a special position.  The other three atoms can come in
      // random order, and should be sorted in increasing numerical
      // order to check for duplicates.  This requires three pairwise
      // swaps:
      if (third > fourth) std::swap(third, fourth);
      if (second > third) std::swap(second, third);
      if (third > fourth) std::swap(third, fourth);

      iter = allInversions.find(inversionTuple);
      if (iter == allInversions.end()) {
        allInversions.insert(inversionTuple);
      } else {
        oss << "Error in Fragment " << getName() << ": "
            << "Inversion" << containerToString(inversion)
            << " appears multiple times\n";
        throw OpenMDException(oss.str());
      }
    }
  }

  void FragmentStamp::checkRigidBodies() {
    std::ostringstream oss;
    std::vector<RigidBodyStamp*>::iterator ri =
        std::find(rigidBodyStamps_.begin(), rigidBodyStamps_.end(),
                  static_cast<RigidBodyStamp*>(NULL));
    if (ri != rigidBodyStamps_.end()) {
      oss << "Error in Fragment " << getName() << ":rigidBody["
          << ri - rigidBodyStamps_.begin() << "] is missing\n";
      throw OpenMDException(oss.str());
    }

    for (std::size_t i = 0; i < getNRigidBodies(); ++i) {
      RigidBodyStamp* rbStamp      = getRigidBodyStamp(i);
      std::vector<int> rigidAtoms  = rbStamp->getMembers();
      std::vector<int>::iterator j = std::find_if(
          rigidAtoms.begin(), rigidAtoms.end(),
          std::bind(std::greater<int>(), placeholders::_1, getNAtoms() - 1));
      if (j != rigidAtoms.end()) {
        oss << "Error in Fragment " << getName();
        throw OpenMDException(oss.str());
      }
    }
  }

  void FragmentStamp::checkCutoffGroups() {
    std::vector<AtomStamp*>::iterator ai;
    std::vector<int>::iterator fai;

    // add all atoms into freeAtoms_ set
    for (ai = atomStamps_.begin(); ai != atomStamps_.end(); ++ai) {
      freeAtoms_.push_back((*ai)->getIndex());
    }

    for (std::size_t i = 0; i < getNCutoffGroups(); ++i) {
      CutoffGroupStamp* cutoffGroupStamp = getCutoffGroupStamp(i);
      std::vector<int> cutoffGroupAtoms  = cutoffGroupStamp->getMembers();
      std::vector<int>::iterator j       = std::find_if(
          cutoffGroupAtoms.begin(), cutoffGroupAtoms.end(),
          std::bind(std::greater<int>(), placeholders::_1, getNAtoms() - 1));
      if (j != cutoffGroupAtoms.end()) {
        std::ostringstream oss;
        oss << "Error in Fragment " << getName() << ": cutoffGroup"
            << " is out of range\n";
        throw OpenMDException(oss.str());
      }

      for (fai = cutoffGroupAtoms.begin(); fai != cutoffGroupAtoms.end();
           ++fai) {
        // erase the atoms belonging to cutoff groups from freeAtoms_ vector
        freeAtoms_.erase(
            std::remove(freeAtoms_.begin(), freeAtoms_.end(), (*fai)),
            freeAtoms_.end());
      }
    }
  }

  void FragmentStamp::checkConstraints() {
    std::ostringstream oss;
    // make sure index is not out of range
    int natoms = getNAtoms();
    for (std::size_t i = 0; i < getNConstraints(); ++i) {
      ConstraintStamp* constraintStamp = getConstraintStamp(i);
      if (constraintStamp->getA() > natoms - 1 || constraintStamp->getA() < 0 ||
          constraintStamp->getB() > natoms - 1 || constraintStamp->getB() < 0 ||
          constraintStamp->getA() == constraintStamp->getB()) {
        oss << "Error in Fragment " << getName() << ": constraint("
            << constraintStamp->getA() << ", " << constraintStamp->getB()
            << ") is invalid\n";
        throw OpenMDException(oss.str());
      }
    }

    // make sure constraints are unique
    std::set<std::pair<int, int>> allConstraints;
    for (std::size_t i = 0; i < getNConstraints(); ++i) {
      ConstraintStamp* constraintStamp = getConstraintStamp(i);
      std::pair<int, int> constraintPair(constraintStamp->getA(),
                                         constraintStamp->getB());
      // make sure constraintPair.first is always less than or equal to
      // constraintPair.third
      if (constraintPair.first > constraintPair.second) {
        std::swap(constraintPair.first, constraintPair.second);
      }

      std::set<std::pair<int, int>>::iterator iter =
          allConstraints.find(constraintPair);
      if (iter != allConstraints.end()) {
        oss << "Error in Fragment " << getName() << ": "
            << "constraint(" << iter->first << ", " << iter->second
            << ") appears multiple times\n";
        throw OpenMDException(oss.str());
      } else {
        allConstraints.insert(constraintPair);
      }
    }

    // make sure atoms belong to same rigidbody are not constrained to
    // each other
    for (std::size_t i = 0; i < getNConstraints(); ++i) {
      ConstraintStamp* constraintStamp = getConstraintStamp(i);
      if (atom2Rigidbody[constraintStamp->getA()] ==
          atom2Rigidbody[constraintStamp->getB()]) {
        oss << "Error in Fragment " << getName() << ": "
            << "constraint(" << constraintStamp->getA() << ", "
            << constraintStamp->getB() << ") belong to same rigidbody "
            << atom2Rigidbody[constraintStamp->getA()] << "\n";
        throw OpenMDException(oss.str());
      }
    }
  }

  void FragmentStamp::checkNodes() {
    std::ostringstream oss;
    std::vector<NodesStamp*>::iterator ni =
        std::find(nodesStamps_.begin(), nodesStamps_.end(),
                  static_cast<NodesStamp*>(NULL));
    if (ni != nodesStamps_.end()) {
      oss << "Error in Molecule " << getName() << ":nodes["
          << ni - nodesStamps_.begin() << "] is missing\n";
      throw OpenMDException(oss.str());
    }

    for (std::size_t i = 0; i < getNNodes(); ++i) {
      NodesStamp* nStamp           = getNodesStamp(i);
      std::vector<int> nodeAtoms   = nStamp->getMembers();
      std::vector<int>::iterator j = std::find_if(
          nodeAtoms.begin(), nodeAtoms.end(),
          std::bind(std::greater<int>(), placeholders::_1, getNAtoms() - 1));
      if (j != nodeAtoms.end()) {
        oss << "Error in Fragment " << getName();
        throw OpenMDException(oss.str());
      }
    }
  }

  void FragmentStamp::fillBondInfo() {
    for (std::size_t i = 0; i < getNBonds(); ++i) {
      BondStamp* bondStamp = getBondStamp(i);
      int a                = bondStamp->getA();
      int b                = bondStamp->getB();
      AtomStamp* atomA     = getAtomStamp(a);
      AtomStamp* atomB     = getAtomStamp(b);
      atomA->addBond(i);
      atomA->addBondedAtom(b);
      atomB->addBond(i);
      atomB->addBondedAtom(a);
    }
  }

  // Function Name: isBondInSameRigidBody
  // Returns true is both atoms of the bond belong to the same rigid
  // body, otherwise return false
  bool FragmentStamp::isBondInSameRigidBody(BondStamp* bond) {
    int rbA;
    int rbB;
    int consAtomA;
    int consAtomB;

    if (!isAtomInRigidBody(bond->getA(), rbA, consAtomA)) return false;

    if (!isAtomInRigidBody(bond->getB(), rbB, consAtomB)) return false;

    if (rbB == rbA)
      return true;
    else
      return false;
  }

  // Function Name: isAtomInRigidBody
  // Returns false if atom does not belong to a rigid body, otherwise
  // returns true
  bool FragmentStamp::isAtomInRigidBody(int atomIndex) {
    return atom2Rigidbody[atomIndex] >= 0;
  }

  // Function Name: isAtomInRigidBody
  // Returns false if atom does not belong to a rigid body otherwise
  // returns true and sets whichRigidBody and consAtomIndex
  // atomIndex : the index of atom in component
  // whichRigidBody: the index of the rigidbody in the component
  // consAtomIndex:  the position the joint atom appears in the rigidbody's
  //                 definition
  bool FragmentStamp::isAtomInRigidBody(int atomIndex, int& whichRigidBody,
                                        int& consAtomIndex) {
    whichRigidBody = -1;
    consAtomIndex  = -1;

    if (atom2Rigidbody[atomIndex] >= 0) {
      whichRigidBody          = atom2Rigidbody[atomIndex];
      RigidBodyStamp* rbStamp = getRigidBodyStamp(whichRigidBody);
      int numAtom             = rbStamp->getNMembers();
      for (int j = 0; j < numAtom; j++) {
        if (rbStamp->getMemberAt(j) == atomIndex) {
          consAtomIndex = j;
          return true;
        }
      }
    }

    return false;
  }

  // Returns the position of joint atoms apearing in a rigidbody's definition
  // For the time being, we will use the most inefficient algorithm,
  // the complexity is O(N^2).  We could improve the
  // complexity to O(NlogN) by sorting the atom index in rigid body
  // first
  std::vector<std::pair<int, int>> FragmentStamp::getJointAtoms(int rb1,
                                                                int rb2) {
    RigidBodyStamp* rbStamp1;
    RigidBodyStamp* rbStamp2;
    int natomInRb1;
    int natomInRb2;
    int atomIndex1;
    int atomIndex2;
    std::vector<std::pair<int, int>> jointAtomIndexPair;

    rbStamp1   = this->getRigidBodyStamp(rb1);
    natomInRb1 = rbStamp1->getNMembers();

    rbStamp2   = this->getRigidBodyStamp(rb2);
    natomInRb2 = rbStamp2->getNMembers();

    for (int i = 0; i < natomInRb1; i++) {
      atomIndex1 = rbStamp1->getMemberAt(i);

      for (int j = 0; j < natomInRb2; j++) {
        atomIndex2 = rbStamp2->getMemberAt(j);

        if (atomIndex1 == atomIndex2) {
          jointAtomIndexPair.push_back(std::make_pair(i, j));
          break;
        }
      }
    }

    return jointAtomIndexPair;
  }

}  // namespace OpenMD
