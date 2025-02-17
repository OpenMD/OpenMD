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
#include "selection/NameFinder.hpp"

#include "primitives/Molecule.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"
#include "utils/wildcards.hpp"

namespace OpenMD {

  NameFinder::NameFinder(SimInfo* info) : info_(info) {
    nObjects_.push_back(info_->getNGlobalAtoms() +
                        info_->getNGlobalRigidBodies());
    nObjects_.push_back(info_->getNGlobalBonds());
    nObjects_.push_back(info_->getNGlobalBends());
    nObjects_.push_back(info_->getNGlobalTorsions());
    nObjects_.push_back(info_->getNGlobalInversions());
    nObjects_.push_back(info_->getNGlobalMolecules());

    loadNames();
  }

  void NameFinder::loadNames() {
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;

    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;

    root_ = std::make_shared<TreeNode>();
    root_->bs.resize(nObjects_);
    root_->bs.setAll();  //

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      std::string molName = mol->getType();
      TreeNodePtr molNode = createNode(root_, molName);
      molNode->bs.bitsets_[MOLECULE].setBitOn(mol->getGlobalIndex());

      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        std::string atomName = atom->getType();
        TreeNodePtr atomNode = createNode(molNode, atomName);

        molNode->bs.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
        atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
      }

      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        std::string rbName = rb->getType();
        TreeNodePtr rbNode = createNode(molNode, rbName);

        molNode->bs.bitsets_[STUNTDOUBLE].setBitOn(rb->getGlobalIndex());
        rbNode->bs.bitsets_[STUNTDOUBLE].setBitOn(rb->getGlobalIndex());

        // COMMENTED OUT because rigid bodies are IntegrableObjects
        // (e.g. they are independently mobile, so selecting their
        // member atoms will give some odd results if we are computing
        // degrees of freedom elsewhere.

        // //create nodes for atoms belong to this rigidbody
        // for(atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai))
        // {
        //   std::string rbAtomName = atom->getType();
        //   TreeNodePtr rbAtomNode = createNode(rbNode, rbAtomName);

        //   rbAtomNode->bs.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
        // }
      }

      for (bond = mol->beginBond(bondIter); bond != NULL;
           bond = mol->nextBond(bondIter)) {
        std::string bondName = bond->getName();
        TreeNodePtr bondNode = createNode(molNode, bondName);

        molNode->bs.bitsets_[BOND].setBitOn(bond->getGlobalIndex());
        bondNode->bs.bitsets_[BOND].setBitOn(bond->getGlobalIndex());

        std::vector<Atom*> atoms = bond->getAtoms();
        std::vector<Atom*>::iterator ai;

        for (ai = atoms.begin(); ai != atoms.end(); ++ai) {
          std::string atomName = (*ai)->getType();
          TreeNodePtr atomNode = createNode(bondNode, atomName);
          atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn((*ai)->getGlobalIndex());
        }
      }
      for (bend = mol->beginBend(bendIter); bend != NULL;
           bend = mol->nextBend(bendIter)) {
        std::string bendName = bend->getName();
        TreeNodePtr bendNode = createNode(molNode, bendName);

        molNode->bs.bitsets_[BEND].setBitOn(bend->getGlobalIndex());
        bendNode->bs.bitsets_[BEND].setBitOn(bend->getGlobalIndex());

        std::vector<Atom*> atoms = bend->getAtoms();
        std::vector<Atom*>::iterator ai;

        for (ai = atoms.begin(); ai != atoms.end(); ++ai) {
          std::string atomName = (*ai)->getType();
          TreeNodePtr atomNode = createNode(bendNode, atomName);
          atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn((*ai)->getGlobalIndex());
        }
      }
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL;
           torsion = mol->nextTorsion(torsionIter)) {
        std::string torsionName = torsion->getName();
        TreeNodePtr torsionNode = createNode(molNode, torsionName);

        molNode->bs.bitsets_[TORSION].setBitOn(torsion->getGlobalIndex());
        torsionNode->bs.bitsets_[TORSION].setBitOn(torsion->getGlobalIndex());

        std::vector<Atom*> atoms = torsion->getAtoms();
        std::vector<Atom*>::iterator ai;

        for (ai = atoms.begin(); ai != atoms.end(); ++ai) {
          std::string atomName = (*ai)->getType();
          TreeNodePtr atomNode = createNode(torsionNode, atomName);
          atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn((*ai)->getGlobalIndex());
        }
      }
      for (inversion = mol->beginInversion(inversionIter); inversion != NULL;
           inversion = mol->nextInversion(inversionIter)) {
        std::string inversionName = inversion->getName();
        TreeNodePtr inversionNode = createNode(molNode, inversionName);

        molNode->bs.bitsets_[INVERSION].setBitOn(inversion->getGlobalIndex());
        inversionNode->bs.bitsets_[INVERSION].setBitOn(
            inversion->getGlobalIndex());
        std::vector<Atom*> atoms = inversion->getAtoms();
        std::vector<Atom*>::iterator ai;

        for (ai = atoms.begin(); ai != atoms.end(); ++ai) {
          std::string atomName = (*ai)->getType();
          TreeNodePtr atomNode = createNode(inversionNode, atomName);
          atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn((*ai)->getGlobalIndex());
        }
      }
    }
  }

  TreeNodePtr NameFinder::createNode(TreeNodePtr parent,
                                     const std::string& name) {
    TreeNodePtr node;
    std::map<std::string, TreeNodePtr>::iterator foundIter;
    foundIter = parent->children.find(name);
    if (foundIter == parent->children.end()) {
      node       = std::make_shared<TreeNode>();
      node->name = name;
      node->bs.resize(nObjects_);
      parent->children.insert(std::make_pair(name, node));
    } else {
      node = foundIter->second;
    }
    return node;
  }

  SelectionSet NameFinder::match(const std::string& name) {
    SelectionSet bs(nObjects_);

    StringTokenizer tokenizer(name, ".");

    std::vector<std::string> names;
    while (tokenizer.hasMoreTokens()) {
      names.push_back(tokenizer.nextToken());
    }

    int size = names.size();

    switch (size) {
    case 1:
      // could be molecule name, atom name and rigidbody name
      matchMolecule(names[0], bs);
      matchStuntDouble("*", names[0], bs);
      matchBond("*", names[0], bs);
      matchBend("*", names[0], bs);
      matchTorsion("*", names[0], bs);
      matchInversion("*", names[0], bs);

      break;
    case 2:
      // could be molecule.*(include atoms and rigidbodies) or rigidbody.*(atoms
      // belong to rigidbody)

      if (!isInteger(names[1])) {
        matchRigidAtoms("*", names[0], names[1], bs);
        matchStuntDouble(names[0], names[1], bs);
      } else {
        int internalIndex = lexi_cast<int>(names[1]);
        if (internalIndex < 0) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "NameFinder : Name %s.%s is an invalid name.\n",
                   names[0].c_str(), names[1].c_str());
          painCave.severity = OPENMD_WARNING;
          painCave.isFatal  = 0;
          simError();
        } else {
          matchInternalIndex(names[0], internalIndex, bs);
        }
      }

      break;
    case 3:
      // must be molecule.rigidbody.*
      matchRigidAtoms(names[0], names[1], names[2], bs);
      break;
    default:
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "NameFinder : Invalid Name %s.\n", name.c_str());
      painCave.severity = OPENMD_WARNING;
      painCave.isFatal  = 0;
      simError();
      break;
    }
    return bs;
  }

  void NameFinder::matchMolecule(const std::string& molName, SelectionSet& bs) {
    std::vector<TreeNodePtr> molNodes = getMatchedChildren(root_, molName);
    std::vector<TreeNodePtr>::iterator i;
    for (i = molNodes.begin(); i != molNodes.end(); ++i) {
      bs |= (*i)->bs;
    }
  }

  void NameFinder::matchStuntDouble(const std::string& molName,
                                    const std::string& sdName,
                                    SelectionSet& bs) {
    std::vector<TreeNodePtr> molNodes = getMatchedChildren(root_, molName);
    std::vector<TreeNodePtr>::iterator i;
    for (i = molNodes.begin(); i != molNodes.end(); ++i) {
      std::vector<TreeNodePtr> sdNodes = getMatchedChildren(*i, sdName);
      std::vector<TreeNodePtr>::iterator j;
      for (j = sdNodes.begin(); j != sdNodes.end(); ++j) {
        bs |= (*j)->bs;
      }
    }
  }

  void NameFinder::matchBond(const std::string& molName,
                             const std::string& bondName, SelectionSet& bs) {
    std::vector<TreeNodePtr> molNodes = getMatchedChildren(root_, molName);
    std::vector<TreeNodePtr>::iterator i;
    for (i = molNodes.begin(); i != molNodes.end(); ++i) {
      std::vector<TreeNodePtr> bondNodes = getMatchedChildren(*i, bondName);
      std::vector<TreeNodePtr>::iterator j;
      for (j = bondNodes.begin(); j != bondNodes.end(); ++j) {
        bs |= (*j)->bs;
        std::vector<TreeNodePtr> bondAtomNodes = getAllChildren(*j);
        std::vector<TreeNodePtr>::iterator k;
        for (k = bondAtomNodes.begin(); k != bondAtomNodes.end(); ++k) {
          bs |= (*k)->bs;
        }
      }
    }
  }

  void NameFinder::matchBend(const std::string& molName,
                             const std::string& bendName, SelectionSet& bs) {
    std::vector<TreeNodePtr> molNodes = getMatchedChildren(root_, molName);
    std::vector<TreeNodePtr>::iterator i;
    for (i = molNodes.begin(); i != molNodes.end(); ++i) {
      std::vector<TreeNodePtr> bendNodes = getMatchedChildren(*i, bendName);
      std::vector<TreeNodePtr>::iterator j;
      for (j = bendNodes.begin(); j != bendNodes.end(); ++j) {
        std::vector<TreeNodePtr> bendAtomNodes = getAllChildren(*j);
        std::vector<TreeNodePtr>::iterator k;
        for (k = bendAtomNodes.begin(); k != bendAtomNodes.end(); ++k) {
          bs |= (*k)->bs;
        }
      }
    }
  }
  void NameFinder::matchTorsion(const std::string& molName,
                                const std::string& torsionName,
                                SelectionSet& bs) {
    std::vector<TreeNodePtr> molNodes = getMatchedChildren(root_, molName);
    std::vector<TreeNodePtr>::iterator i;
    for (i = molNodes.begin(); i != molNodes.end(); ++i) {
      std::vector<TreeNodePtr> torsionNodes =
          getMatchedChildren(*i, torsionName);
      std::vector<TreeNodePtr>::iterator j;
      for (j = torsionNodes.begin(); j != torsionNodes.end(); ++j) {
        std::vector<TreeNodePtr> torsionAtomNodes = getAllChildren(*j);
        std::vector<TreeNodePtr>::iterator k;
        for (k = torsionAtomNodes.begin(); k != torsionAtomNodes.end(); ++k) {
          bs |= (*k)->bs;
        }
      }
    }
  }
  void NameFinder::matchInversion(const std::string& molName,
                                  const std::string& inversionName,
                                  SelectionSet& bs) {
    std::vector<TreeNodePtr> molNodes = getMatchedChildren(root_, molName);
    std::vector<TreeNodePtr>::iterator i;
    for (i = molNodes.begin(); i != molNodes.end(); ++i) {
      std::vector<TreeNodePtr> inversionNodes =
          getMatchedChildren(*i, inversionName);
      std::vector<TreeNodePtr>::iterator j;
      for (j = inversionNodes.begin(); j != inversionNodes.end(); ++j) {
        std::vector<TreeNodePtr> inversionAtomNodes = getAllChildren(*j);
        std::vector<TreeNodePtr>::iterator k;
        for (k = inversionAtomNodes.begin(); k != inversionAtomNodes.end();
             ++k) {
          bs |= (*k)->bs;
        }
      }
    }
  }

  void NameFinder::matchRigidAtoms(const std::string& molName,
                                   const std::string& rbName,
                                   const std::string& rbAtomName,
                                   SelectionSet& bs) {
    std::vector<TreeNodePtr> molNodes = getMatchedChildren(root_, molName);
    std::vector<TreeNodePtr>::iterator i;
    for (i = molNodes.begin(); i != molNodes.end(); ++i) {
      std::vector<TreeNodePtr> rbNodes = getMatchedChildren(*i, rbName);
      std::vector<TreeNodePtr>::iterator j;
      for (j = rbNodes.begin(); j != rbNodes.end(); ++j) {
        std::vector<TreeNodePtr> rbAtomNodes =
            getMatchedChildren(*j, rbAtomName);
        std::vector<TreeNodePtr>::iterator k;
        for (k = rbAtomNodes.begin(); k != rbAtomNodes.end(); ++k) {
          bs |= (*k)->bs;
        }
      }
    }
  }

  std::vector<TreeNodePtr> NameFinder::getAllChildren(TreeNodePtr node) {
    std::vector<TreeNodePtr> childNodes;
    std::map<std::string, TreeNodePtr>::iterator i;
    for (i = node->children.begin(); i != node->children.end(); ++i) {
      childNodes.push_back(i->second);
    }
    return childNodes;
  }

  std::vector<TreeNodePtr> NameFinder::getMatchedChildren(
      TreeNodePtr node, const std::string& name) {
    std::vector<TreeNodePtr> matchedNodes;
    std::map<std::string, TreeNodePtr>::iterator i;
    for (i = node->children.begin(); i != node->children.end(); ++i) {
      if (isMatched(i->first, name)) { matchedNodes.push_back(i->second); }
    }

    return matchedNodes;
  }

  bool NameFinder::isMatched(const std::string& str,
                             const std::string& wildcard) {
    return Wildcard::wildcardfit(wildcard.c_str(), str.c_str()) > 0 ? true :
                                                                      false;
  }

  void NameFinder::matchInternalIndex(const std::string& name,
                                      int internalIndex, SelectionSet& bs) {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      if (isMatched(mol->getType(), name)) {
        int natoms       = mol->getNAtoms();
        int nrigidbodies = mol->getNRigidBodies();
        if (internalIndex >= natoms + nrigidbodies) {
          continue;
        } else if (internalIndex < natoms) {
          bs.bitsets_[STUNTDOUBLE].setBitOn(
              mol->getAtomAt(internalIndex)->getGlobalIndex());
          continue;
        } else if (internalIndex < natoms + nrigidbodies) {
          bs.bitsets_[STUNTDOUBLE].setBitOn(
              mol->getRigidBodyAt(internalIndex - natoms)->getGlobalIndex());
        }
      }
    }
  }

  bool NameFinder::isInteger(const std::string& str) {
    for (unsigned int i = 0; i < str.size(); ++i) {
      if (!std::isdigit(str[i])) { return false; }
    }
    return true;
  }
}  // namespace OpenMD
