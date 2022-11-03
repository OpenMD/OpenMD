/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
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

#include "hydrodynamics/BeadModel.hpp"
#include "types/LennardJonesAdapter.hpp"

namespace OpenMD {
  bool BeadModel::createBeads(std::vector<BeadParam>& beads) {
    if (sd_->isAtom()) {
      if (!createSingleBead(static_cast<Atom*>(sd_), beads)) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                "BeadModel::createBeads Error: GayBerne and other "
                "non-spherical atoms should use the RoughShell model\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
        return false;
      }
    } else if (sd_->isRigidBody()) {
      RigidBody* rb = static_cast<RigidBody*>(sd_);
      std::vector<Atom*>::iterator ai;
      Atom* atom;
      for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
        if (!createSingleBead(atom, beads)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                  "BeadModel::createBeads Error: GayBerne and other "
                  "non-spherical atoms should use the RoughShell model\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
          return false;
        }
      }
    }
    return true;
  }

  bool BeadModel::createSingleBead(Atom* atom, std::vector<BeadParam>& beads) {
    AtomType* atomType      = atom->getAtomType();
    LennardJonesAdapter lja = LennardJonesAdapter(atomType);
    if (atomType->isGayBerne()) {
      return false;
    } else if (lja.isLennardJones()) {
      BeadParam currBead;
      currBead.atomName = atom->getType();
      currBead.pos      = atom->getPos();
      currBead.mass     = atom->getMass();  // to compute center of mass in
                                            // ApproximationModel.cpp
      currBead.radius = lja.getSigma() / 2.0;
      std::cout << "using rLJ = " << currBead.radius << " for atom "
                << currBead.atomName << "\n";
      beads.push_back(currBead);
    } else {
      std::cout << "For atom " << atom->getType() << ", trying ";
      int obanum(0);
      std::vector<AtomType*> atChain = atomType->allYourBase();
      std::vector<AtomType*>::iterator i;
      for (i = atChain.begin(); i != atChain.end(); ++i) {
        std::cout << (*i)->getName() << " -> ";
        obanum = etab.GetAtomicNum((*i)->getName().c_str());
        if (obanum != 0) {
          BeadParam currBead;
          currBead.atomName = atom->getType();
          currBead.pos      = atom->getPos();
          currBead.mass     = atom->getMass();  // to compute center of mass in
                                                // ApproximationModel.cpp
          currBead.radius = etab.GetVdwRad(obanum);
          std::cout << "using rVdW = " << currBead.radius << "\n";
          beads.push_back(currBead);
          break;
        }
      }
      if (obanum == 0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                "Could not find atom type in default element.txt\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    }
    return true;
  }
}  // namespace OpenMD
