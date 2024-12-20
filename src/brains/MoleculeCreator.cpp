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
 * @file MoleculeCreator.cpp
 * @author tlin
 * @date 11/04/2004
 * @version 1.0
 */

#include "brains/MoleculeCreator.hpp"

#include <cassert>
#include <set>
#include <typeinfo>

#include "primitives/GhostBend.hpp"
#include "primitives/GhostTorsion.hpp"
#include "types/AtomType.hpp"
#include "types/BendTypeParser.hpp"
#include "types/BondTypeParser.hpp"
#include "types/FixedBondType.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/InversionTypeParser.hpp"
#include "types/TorsionTypeParser.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {

  void MoleculeCreator::createOverrideAtomTypes(ForceField* ff,
						MoleculeStamp* molStamp) {

    // Some stamps have overrides of default types in the force field.
    // This runs through atomStamps which can have overrides and makes sure
    // that the types are created, even if this processor doesn't end up
    // owning that molecule:

    AtomStamp* stamp;
    size_t nAtom = molStamp->getNAtoms();
    
    for (size_t i = 0; i < nAtom; ++i) {
      stamp = molStamp->getAtomStamp(i);
      AtomType* atomType = ff->getAtomType(stamp->getType());
    
      if (atomType == NULL) {
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "Can not find Matching Atom Type for[%s]",
		 stamp->getType().c_str());
	
	painCave.isFatal = 1;
	simError();
      }

      if (stamp->hasOverride()) {
	std::string baseType = atomType->getName();
	RealType oc          = stamp->getOverrideCharge();
	
	// Create a new atom type name that builds in the override charge:
	std::ostringstream ss;
	ss << oc;
	std::string atomTypeOverrideName = baseType + "_q=" + ss.str();
	
	// Maybe we've seen this before?
	
	AtomType* atB = ff->getAtomType(atomTypeOverrideName);
	
	if (atB == NULL) {
	  // Nope, we've never seen it before, so make a new one:
	  AtomType* atomTypeOverride = new AtomType();
	  // Base points to the atomType we already found
	  atomTypeOverride->useBase(atomType);
	  int ident = ff->getNAtomType();
	  atomTypeOverride->setIdent(ident);
	  atomTypeOverride->setName(atomTypeOverrideName);
	  ff->addAtomType(atomTypeOverrideName, atomTypeOverride);
	  FixedChargeAdapter fca = FixedChargeAdapter(atomTypeOverride);
	  RealType charge =
            ff->getForceFieldOptions().getChargeUnitScaling() * oc;
	  fca.makeFixedCharge(charge);
	} 
      }
    }    
  }
  
  Molecule* MoleculeCreator::createMolecule(ForceField* ff,
                                            MoleculeStamp* molStamp,
                                            int globalIndex,
                                            LocalIndexManager* localIndexMan) {
    Molecule* mol = new Molecule(globalIndex, molStamp);

    // create atoms
    Atom* atom;
    AtomStamp* currentAtomStamp;
    size_t nAtom = molStamp->getNAtoms();
    for (size_t i = 0; i < nAtom; ++i) {
      currentAtomStamp = molStamp->getAtomStamp(i);
      atom             = createAtom(ff, currentAtomStamp, localIndexMan);
      mol->addAtom(atom);
    }

    // create rigidbodies
    RigidBody* rb;
    RigidBodyStamp* currentRigidBodyStamp;
    size_t nRigidbodies = molStamp->getNRigidBodies();

    for (size_t i = 0; i < nRigidbodies; ++i) {
      currentRigidBodyStamp = molStamp->getRigidBodyStamp(i);
      rb = createRigidBody(molStamp, mol, currentRigidBodyStamp, localIndexMan);
      mol->addRigidBody(rb);
    }

    // create bonds
    Bond* bond;
    BondStamp* currentBondStamp;
    size_t nBonds = molStamp->getNBonds();

    for (size_t i = 0; i < nBonds; ++i) {
      currentBondStamp = molStamp->getBondStamp(i);
      bond             = createBond(ff, mol, currentBondStamp, localIndexMan);
      mol->addBond(bond);
    }

    // create bends
    Bend* bend;
    BendStamp* currentBendStamp;
    size_t nBends = molStamp->getNBends();
    for (size_t i = 0; i < nBends; ++i) {
      currentBendStamp = molStamp->getBendStamp(i);
      bend             = createBend(ff, mol, currentBendStamp, localIndexMan);
      mol->addBend(bend);
    }

    // create torsions
    Torsion* torsion;
    TorsionStamp* currentTorsionStamp;
    size_t nTorsions = molStamp->getNTorsions();
    for (size_t i = 0; i < nTorsions; ++i) {
      currentTorsionStamp = molStamp->getTorsionStamp(i);
      torsion = createTorsion(ff, mol, currentTorsionStamp, localIndexMan);
      mol->addTorsion(torsion);
    }

    // create inversions
    Inversion* inversion;
    InversionStamp* currentInversionStamp;
    size_t nInversions = molStamp->getNInversions();
    for (size_t i = 0; i < nInversions; ++i) {
      currentInversionStamp = molStamp->getInversionStamp(i);
      inversion =
          createInversion(ff, mol, currentInversionStamp, localIndexMan);
      if (inversion != NULL) { mol->addInversion(inversion); }
    }

    // create cutoffGroups
    CutoffGroup* cutoffGroup;
    CutoffGroupStamp* currentCutoffGroupStamp;
    size_t nCutoffGroups = molStamp->getNCutoffGroups();
    for (size_t i = 0; i < nCutoffGroups; ++i) {
      currentCutoffGroupStamp = molStamp->getCutoffGroupStamp(i);
      cutoffGroup =
          createCutoffGroup(mol, currentCutoffGroupStamp, localIndexMan);
      mol->addCutoffGroup(cutoffGroup);
    }

    // every free atom is a cutoff group
    std::vector<Atom*> freeAtoms;
    std::vector<Atom*>::iterator ai;
    std::vector<Atom*>::iterator fai;

    // add all atoms into allAtoms set
    for (atom = mol->beginAtom(fai); atom != NULL; atom = mol->nextAtom(fai)) {
      freeAtoms.push_back(atom);
    }

    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;

    for (cg = mol->beginCutoffGroup(ci); cg != NULL;
         cg = mol->nextCutoffGroup(ci)) {
      for (atom = cg->beginAtom(ai); atom != NULL; atom = cg->nextAtom(ai)) {
        // erase the atoms belong to cutoff groups from freeAtoms vector
        freeAtoms.erase(std::remove(freeAtoms.begin(), freeAtoms.end(), atom),
                        freeAtoms.end());
      }
    }

    // loop over the free atoms and then create one cutoff group for
    // every single free atom

    for (fai = freeAtoms.begin(); fai != freeAtoms.end(); ++fai) {
      cutoffGroup = createCutoffGroup(mol, *fai, localIndexMan);
      mol->addCutoffGroup(cutoffGroup);
    }

    // create bonded constraintPairs:
    createConstraintPair(mol);

    // create non-bonded constraintPairs
    for (std::size_t i = 0; i < molStamp->getNConstraints(); ++i) {
      ConstraintStamp* cStamp = molStamp->getConstraintStamp(i);
      Atom* atomA;
      Atom* atomB;

      atomA = mol->getAtomAt(cStamp->getA());
      atomB = mol->getAtomAt(cStamp->getB());
      assert(atomA && atomB);

      bool printConstraintForce = false;

      if (cStamp->havePrintConstraintForce()) {
        printConstraintForce = cStamp->getPrintConstraintForce();
      }

      if (!cStamp->haveConstrainedDistance()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Constraint Error: A non-bond constraint was specified\n"
                 "\twithout providing a value for the constrainedDistance.\n");
        painCave.isFatal = 1;
        simError();
      } else {
        RealType distance         = cStamp->getConstrainedDistance();
        ConstraintElem* consElemA = new ConstraintElem(atomA);
        ConstraintElem* consElemB = new ConstraintElem(atomB);
        ConstraintPair* cPair     = new ConstraintPair(
            consElemA, consElemB, distance, printConstraintForce);
        mol->addConstraintPair(cPair);
      }
    }

    // now create the constraint elements:

    createConstraintElem(mol);

    // Does this molecule stamp define a total constrained charge value?
    // If so, let the created molecule know about it.
    if (molStamp->haveConstrainTotalCharge()) {
      mol->setConstrainTotalCharge(molStamp->getConstrainTotalCharge());
    }

    // The construction of this molecule is finished:
    mol->complete();

    return mol;
  }

  Atom* MoleculeCreator::createAtom(ForceField* ff, AtomStamp* stamp,
                                    LocalIndexManager* localIndexMan) {
    AtomType* atomType;
    Atom* atom;

    atomType = ff->getAtomType(stamp->getType());
    if (atomType == NULL) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Can not find Matching Atom Type for[%s]",
               stamp->getType().c_str());

      painCave.isFatal = 1;
      simError();
    }

    if (stamp->hasOverride()) {
      std::string baseType = atomType->getName();
      RealType oc          = stamp->getOverrideCharge();

      // Create a new atom type name that builds in the override charge:
      std::ostringstream ss;
      ss << oc;
      std::string atomTypeOverrideName = baseType + "_q=" + ss.str();

      // Maybe we've seen this before?

      AtomType* atB = ff->getAtomType(atomTypeOverrideName);

      if (atB == NULL) {
        // Nope, we've never seen it before, so make a new one:
        AtomType* atomTypeOverride = new AtomType();
        // Base points to the atomType we already found
        atomTypeOverride->useBase(atomType);
        int ident = ff->getNAtomType();
        atomTypeOverride->setIdent(ident);
        atomTypeOverride->setName(atomTypeOverrideName);
        ff->addAtomType(atomTypeOverrideName, atomTypeOverride);
        FixedChargeAdapter fca = FixedChargeAdapter(atomTypeOverride);
        RealType charge =
            ff->getForceFieldOptions().getChargeUnitScaling() * oc;
        fca.makeFixedCharge(charge);
        // officially use override type for this atom
	atomType = atomTypeOverride;
      } else {
        // we've previously created the override type for this atom, so use that
        // one:
        atomType = atB;
      }
    }

    // below code still have some kind of hard-coding smell
    if (atomType->isDirectional()) {
      DirectionalAtom* dAtom;
      dAtom = new DirectionalAtom(atomType);
      atom  = dAtom;
    } else {
      atom = new Atom(atomType);
    }

    atom->setLocalIndex(localIndexMan->getNextAtomIndex());

    return atom;
  }

  RigidBody* MoleculeCreator::createRigidBody(
      MoleculeStamp* molStamp, Molecule* mol, RigidBodyStamp* rbStamp,
      LocalIndexManager* localIndexMan) {
    Atom* atom;
    size_t nAtoms;
    Vector3d refCoor;
    AtomStamp* atomStamp;

    RigidBody* rb = new RigidBody();
    nAtoms        = rbStamp->getNMembers();
    for (std::size_t i = 0; i < nAtoms; ++i) {
      // rbStamp->getMember(i) return the local index of current atom
      // inside the molecule.  It is not the same as local index of
      // atom which is the index of atom at DataStorage class
      atom      = mol->getAtomAt(rbStamp->getMemberAt(i));
      atomStamp = molStamp->getAtomStamp(rbStamp->getMemberAt(i));
      rb->addAtom(atom, atomStamp);
    }

    // after all of the atoms are added, we need to calculate the
    // reference coordinates
    rb->calcRefCoords();

    // set the local index of this rigid body, global index will be set later
    rb->setLocalIndex(localIndexMan->getNextRigidBodyIndex());

    // The rule for naming a rigidbody is: MoleculeName_RB_Integer
    // The first part is the name of the molecule
    // The second part is always fixed as "RB"
    // The third part is the index of the rigidbody defined in meta-data file
    // For example, Butane_RB_0 is a valid rigid body name of butane molecule

    std::string s = OpenMD_itoa(mol->getNRigidBodies(), 10);
    rb->setType(mol->getType() + "_RB_" + s.c_str());
    return rb;
  }


  Bond* MoleculeCreator::createBond(ForceField* ff, Molecule* mol,
                                    BondStamp* stamp,
                                    LocalIndexManager* localIndexMan) {
    BondTypeParser btParser;
    BondType* bondType = NULL;
    Atom* atomA;
    Atom* atomB;

    atomA = mol->getAtomAt(stamp->getA());
    atomB = mol->getAtomAt(stamp->getB());

    assert(atomA && atomB);

    if (stamp->hasOverride()) {
      try {
        bondType = btParser.parseTypeAndPars(stamp->getOverrideType(),
                                             stamp->getOverridePars());
      } catch (OpenMDException& e) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "MoleculeCreator Error: %s "
                 "for molecule %s\n",
                 e.what(), mol->getType().c_str());
        painCave.isFatal = 1;
        simError();
      }

    } else {
      bondType = ff->getBondType(atomA->getType(), atomB->getType());

      if (bondType == NULL) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Can not find Matching Bond Type for[%s, %s]",
                 atomA->getType().c_str(), atomB->getType().c_str());

        painCave.isFatal = 1;
        simError();
      }
    }

    Bond* bond = new Bond(atomA, atomB, bondType);

    // set the local index of this bond, the global index will be set later
    bond->setLocalIndex(localIndexMan->getNextBondIndex());

    // The rule for naming a bond is: MoleculeName_Bond_Integer
    // The first part is the name of the molecule
    // The second part is always fixed as "Bond"
    // The third part is the index of the bond defined in meta-data file
    // For example, Butane_bond_0 is a valid Bond name in a butane molecule

    std::string s = OpenMD_itoa(mol->getNBonds(), 10);
    bond->setName(mol->getType() + "_Bond_" + s.c_str());
    return bond;
  }

  Bend* MoleculeCreator::createBend(ForceField* ff, Molecule* mol,
                                    BendStamp* stamp,
                                    LocalIndexManager* localIndexMan) {
    BendTypeParser btParser;
    BendType* bendType = NULL;
    Bend* bend         = NULL;

    std::vector<int> bendAtoms = stamp->getMembers();
    if (bendAtoms.size() == 3) {
      Atom* atomA = mol->getAtomAt(bendAtoms[0]);
      Atom* atomB = mol->getAtomAt(bendAtoms[1]);
      Atom* atomC = mol->getAtomAt(bendAtoms[2]);

      assert(atomA && atomB && atomC);

      if (stamp->hasOverride()) {
        try {
          bendType = btParser.parseTypeAndPars(stamp->getOverrideType(),
                                               stamp->getOverridePars());
        } catch (OpenMDException& e) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "MoleculeCreator Error: %s "
                   "for molecule %s\n",
                   e.what(), mol->getType().c_str());
          painCave.isFatal = 1;
          simError();
        }
      } else {
        bendType =
            ff->getBendType(atomA->getType().c_str(), atomB->getType().c_str(),
                            atomC->getType().c_str());

        if (bendType == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not find Matching Bend Type for[%s, %s, %s]",
                   atomA->getType().c_str(), atomB->getType().c_str(),
                   atomC->getType().c_str());

          painCave.isFatal = 1;
          simError();
        }
      }

      bend = new Bend(atomA, atomB, atomC, bendType);

    } else if (bendAtoms.size() == 2 && stamp->haveGhostVectorSource()) {
      int ghostIndex = stamp->getGhostVectorSource();
      int normalIndex =
          ghostIndex != bendAtoms[0] ? bendAtoms[0] : bendAtoms[1];
      Atom* normalAtom = mol->getAtomAt(normalIndex);
      DirectionalAtom* ghostAtom =
          dynamic_cast<DirectionalAtom*>(mol->getAtomAt(ghostIndex));
      if (ghostAtom == NULL) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Can not cast Atom to DirectionalAtom");
        painCave.isFatal = 1;
        simError();
      }

      if (stamp->hasOverride()) {
        try {
          bendType = btParser.parseTypeAndPars(stamp->getOverrideType(),
                                               stamp->getOverridePars());
        } catch (OpenMDException& e) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "MoleculeCreator Error: %s "
                   "for molecule %s\n",
                   e.what(), mol->getType().c_str());
          painCave.isFatal = 1;
          simError();
        }
      } else {
        bendType = ff->getBendType(normalAtom->getType(), ghostAtom->getType(),
                                   "GHOST");

        if (bendType == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not find Matching Bend Type for[%s, %s, %s]",
                   normalAtom->getType().c_str(), ghostAtom->getType().c_str(),
                   "GHOST");

          painCave.isFatal = 1;
          simError();
        }
      }

      bend = new GhostBend(normalAtom, ghostAtom, bendType);
    }

    // set the local index of this bend, the global index will be set later
    bend->setLocalIndex(localIndexMan->getNextBendIndex());

    // The rule for naming a bend is: MoleculeName_Bend_Integer
    // The first part is the name of the molecule
    // The second part is always fixed as "Bend"
    // The third part is the index of the bend defined in meta-data file
    // For example, Butane_Bend_0 is a valid Bend name in a butane molecule

    std::string s = OpenMD_itoa(mol->getNBends(), 10);
    bend->setName(mol->getType() + "_Bend_" + s.c_str());
    return bend;
  }

  Torsion* MoleculeCreator::createTorsion(ForceField* ff, Molecule* mol,
                                          TorsionStamp* stamp,
                                          LocalIndexManager* localIndexMan) {
    TorsionTypeParser ttParser;
    TorsionType* torsionType = NULL;
    Torsion* torsion         = NULL;

    std::vector<int> torsionAtoms = stamp->getMembers();
    if (torsionAtoms.size() < 3) { return torsion; }

    Atom* atomA = mol->getAtomAt(torsionAtoms[0]);
    Atom* atomB = mol->getAtomAt(torsionAtoms[1]);
    Atom* atomC = mol->getAtomAt(torsionAtoms[2]);

    if (torsionAtoms.size() == 4) {
      Atom* atomD = mol->getAtomAt(torsionAtoms[3]);

      assert(atomA && atomB && atomC && atomD);

      if (stamp->hasOverride()) {
        try {
          torsionType = ttParser.parseTypeAndPars(stamp->getOverrideType(),
                                                  stamp->getOverridePars());
        } catch (OpenMDException& e) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "MoleculeCreator Error: %s "
                   "for molecule %s\n",
                   e.what(), mol->getType().c_str());
          painCave.isFatal = 1;
          simError();
        }
      } else {
        torsionType = ff->getTorsionType(atomA->getType(), atomB->getType(),
                                         atomC->getType(), atomD->getType());
        if (torsionType == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not find Matching Torsion Type for[%s, %s, %s, %s]",
                   atomA->getType().c_str(), atomB->getType().c_str(),
                   atomC->getType().c_str(), atomD->getType().c_str());

          painCave.isFatal = 1;
          simError();
        }
      }

      torsion = new Torsion(atomA, atomB, atomC, atomD, torsionType);
    } else {
      DirectionalAtom* dAtom = dynamic_cast<DirectionalAtom*>(
          mol->getAtomAt(stamp->getGhostVectorSource()));
      if (dAtom == NULL) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Can not cast Atom to DirectionalAtom");
        painCave.isFatal = 1;
        simError();
      }

      if (stamp->hasOverride()) {
        try {
          torsionType = ttParser.parseTypeAndPars(stamp->getOverrideType(),
                                                  stamp->getOverridePars());
        } catch (OpenMDException& e) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "MoleculeCreator Error: %s "
                   "for molecule %s\n",
                   e.what(), mol->getType().c_str());
          painCave.isFatal = 1;
          simError();
        }
      } else {
        torsionType = ff->getTorsionType(atomA->getType(), atomB->getType(),
                                         atomC->getType(), "GHOST");

        if (torsionType == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not find Matching Torsion Type for[%s, %s, %s, %s]",
                   atomA->getType().c_str(), atomB->getType().c_str(),
                   atomC->getType().c_str(), "GHOST");

          painCave.isFatal = 1;
          simError();
        }
      }

      torsion = new GhostTorsion(atomA, atomB, dAtom, torsionType);
    }

    // set the local index of this torsion, the global index will be set later
    torsion->setLocalIndex(localIndexMan->getNextTorsionIndex());

    // The rule for naming a torsion is: MoleculeName_Torsion_Integer
    // The first part is the name of the molecule
    // The second part is always fixed as "Torsion"
    // The third part is the index of the torsion defined in meta-data file
    // For example, Butane_Torsion_0 is a valid Torsion name in a
    // butane molecule

    std::string s = OpenMD_itoa(mol->getNTorsions(), 10);
    torsion->setName(mol->getType() + "_Torsion_" + s.c_str());
    return torsion;
  }

  Inversion* MoleculeCreator::createInversion(
      ForceField* ff, Molecule* mol, InversionStamp* stamp,
      LocalIndexManager* localIndexMan) {
    InversionTypeParser itParser;
    InversionType* inversionType = NULL;
    Inversion* inversion         = NULL;

    int center                  = stamp->getCenter();
    std::vector<int> satellites = stamp->getSatellites();
    if (satellites.size() != 3) { return inversion; }

    Atom* atomA = mol->getAtomAt(center);
    Atom* atomB = mol->getAtomAt(satellites[0]);
    Atom* atomC = mol->getAtomAt(satellites[1]);
    Atom* atomD = mol->getAtomAt(satellites[2]);

    assert(atomA && atomB && atomC && atomD);

    if (stamp->hasOverride()) {
      try {
        inversionType = itParser.parseTypeAndPars(stamp->getOverrideType(),
                                                  stamp->getOverridePars());
      } catch (OpenMDException& e) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "MoleculeCreator Error: %s "
                 "for molecule %s\n",
                 e.what(), mol->getType().c_str());
        painCave.isFatal = 1;
        simError();
      }
    } else {
      inversionType = ff->getInversionType(atomA->getType(), atomB->getType(),
                                           atomC->getType(), atomD->getType());

      if (inversionType == NULL) {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "No Matching Inversion Type for[%s, %s, %s, %s]\n"
            "\t(May not be a problem: not all inversions are parametrized)\n",
            atomA->getType().c_str(), atomB->getType().c_str(),
            atomC->getType().c_str(), atomD->getType().c_str());

        painCave.isFatal  = 0;
        painCave.severity = OPENMD_INFO;
        simError();
      }
    }
    if (inversionType != NULL) {
      inversion = new Inversion(atomA, atomB, atomC, atomD, inversionType);

      // set the local index of this inversion, the global index will
      // be set later
      inversion->setLocalIndex(localIndexMan->getNextInversionIndex());

      // The rule for naming an inversion is: MoleculeName_Inversion_Integer
      // The first part is the name of the molecule
      // The second part is always fixed as "Inversion"
      // The third part is the index of the inversion defined in meta-data file
      // For example, Benzene_Inversion_0 is a valid Inversion name in a
      // Benzene molecule

      std::string s = OpenMD_itoa(mol->getNInversions(), 10);
      inversion->setName(mol->getType() + "_Inversion_" + s.c_str());
      return inversion;
    } else {
      return NULL;
    }
  }

  CutoffGroup* MoleculeCreator::createCutoffGroup(
      Molecule* mol, CutoffGroupStamp* stamp,
      LocalIndexManager* localIndexMan) {
    size_t nAtoms;
    CutoffGroup* cg;
    Atom* atom;
    cg = new CutoffGroup();

    nAtoms = stamp->getNMembers();
    for (size_t i = 0; i < nAtoms; ++i) {
      atom = mol->getAtomAt(stamp->getMemberAt(i));
      assert(atom);
      cg->addAtom(atom);
    }

    // set the local index of this cutoffGroup, global index will be set later
    cg->setLocalIndex(localIndexMan->getNextCutoffGroupIndex());

    return cg;
  }

  CutoffGroup* MoleculeCreator::createCutoffGroup(
      Molecule*, Atom* atom, LocalIndexManager* localIndexMan) {
    CutoffGroup* cg;
    cg = new CutoffGroup();
    cg->addAtom(atom);

    // set the local index of this cutoffGroup, global index will be set later
    cg->setLocalIndex(localIndexMan->getNextCutoffGroupIndex());

    return cg;
  }

  void MoleculeCreator::createConstraintPair(Molecule* mol) {
    // add bond constraints
    Molecule::BondIterator bi;
    Bond* bond;
    ConstraintPair* cPair;

    for (bond = mol->beginBond(bi); bond != NULL; bond = mol->nextBond(bi)) {
      BondType* bt = bond->getBondType();

      if (typeid(FixedBondType) == typeid(*bt)) {
        FixedBondType* fbt = dynamic_cast<FixedBondType*>(bt);

        ConstraintElem* consElemA = new ConstraintElem(bond->getAtomA());
        ConstraintElem* consElemB = new ConstraintElem(bond->getAtomB());
        cPair                     = new ConstraintPair(consElemA, consElemB,
                                                       fbt->getEquilibriumBondLength(), false);
        mol->addConstraintPair(cPair);
      }
    }

    // rigidbody -- rigidbody constraint is not support yet
  }

  void MoleculeCreator::createConstraintElem(Molecule* mol) {
    ConstraintPair* consPair;
    Molecule::ConstraintPairIterator cpi;
    std::set<StuntDouble*> sdSet;
    for (consPair = mol->beginConstraintPair(cpi); consPair != NULL;
         consPair = mol->nextConstraintPair(cpi)) {
      StuntDouble* sdA = consPair->getConsElem1()->getStuntDouble();
      if (sdSet.find(sdA) == sdSet.end()) {
        sdSet.insert(sdA);
        mol->addConstraintElem(new ConstraintElem(sdA));
      }

      StuntDouble* sdB = consPair->getConsElem2()->getStuntDouble();
      if (sdSet.find(sdB) == sdSet.end()) {
        sdSet.insert(sdB);
        mol->addConstraintElem(new ConstraintElem(sdB));
      }
    }
  }
}  // namespace OpenMD
