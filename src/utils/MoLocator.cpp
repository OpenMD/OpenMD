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

#include "utils/MoLocator.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "types/AtomType.hpp"
#include "utils/simError.h"

namespace OpenMD {
  MoLocator::MoLocator(MoleculeStamp* theStamp, ForceField* theFF) {
    myStamp            = theStamp;
    myFF               = theFF;
    nIntegrableObjects = myStamp->getNIntegrable();
    calcRef();
  }

  void MoLocator::placeMol(const Vector3d& offset, const Vector3d& ort,
                           Molecule* mol) {
    Vector3d newCoor;
    Vector3d curRefCoor;
    RotMat3x3d rotMat = latVec2RotMat(ort);

    if (mol->getNIntegrableObjects() != nIntegrableObjects) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "MoLocator::placeMol error.\n"
               "\tThe number of integrable objects of MoleculeStamp is not\n"
               "\tthe same as that of Molecule\n");
      painCave.isFatal = 1;
      simError();
    }

    Molecule::IntegrableObjectIterator ii;
    StuntDouble* sd;
    int i;
    for (sd = mol->beginIntegrableObject(ii), i = 0; sd != NULL;
         sd = mol->nextIntegrableObject(ii), ++i) {
      newCoor = rotMat * refCoords[i];
      newCoor += offset;

      sd->setPos(newCoor);
      sd->setVel(V3Zero);

      if (sd->isDirectional()) {
        sd->setA(rotMat * sd->getA());
        sd->setJ(V3Zero);
      }
    }
  }

  void MoLocator::calcRef(void) {
    AtomStamp* currAtomStamp;
    RigidBodyStamp* rbStamp;
    std::vector<RealType> mass;
    Vector3d coor;
    Vector3d refMolCom;
    RealType totMassInRb;
    RealType currAtomMass;
    RealType molMass;

    std::size_t nAtoms       = myStamp->getNAtoms();
    std::size_t nRigidBodies = myStamp->getNRigidBodies();

    for (std::size_t i = 0; i < nAtoms; i++) {
      currAtomStamp = myStamp->getAtomStamp(i);

      if (!currAtomStamp->havePosition()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "MoLocator::calcRef error.\n"
                 "\tComponent %s, atom %s does not have a position specified.\n"
                 "\tThis means MoLocator cannot initalize it's position.\n",
                 myStamp->getName().c_str(), currAtomStamp->getType().c_str());

        painCave.isFatal = 1;
        simError();
      }

      // if atom belongs to rigidbody, just skip it
      if (myStamp->isAtomInRigidBody(i)) continue;
      // get mass and the reference coordinate
      else {
        currAtomMass = getAtomMass(currAtomStamp->getType(), myFF);
        mass.push_back(currAtomMass);
        coor.x() = currAtomStamp->getPosX();
        coor.y() = currAtomStamp->getPosY();
        coor.z() = currAtomStamp->getPosZ();
        refCoords.push_back(coor);
      }
    }

    for (std::size_t i = 0; i < nRigidBodies; i++) {
      rbStamp                = myStamp->getRigidBodyStamp(i);
      std::size_t nAtomsInRb = rbStamp->getNMembers();

      coor.x()    = 0.0;
      coor.y()    = 0.0;
      coor.z()    = 0.0;
      totMassInRb = 0.0;

      for (std::size_t j = 0; j < nAtomsInRb; j++) {
        currAtomStamp = myStamp->getAtomStamp(rbStamp->getMemberAt(j));
        currAtomMass  = getAtomMass(currAtomStamp->getType(), myFF);
        totMassInRb += currAtomMass;

        coor.x() += currAtomStamp->getPosX() * currAtomMass;
        coor.y() += currAtomStamp->getPosY() * currAtomMass;
        coor.z() += currAtomStamp->getPosZ() * currAtomMass;
      }

      mass.push_back(totMassInRb);
      coor /= totMassInRb;
      refCoords.push_back(coor);
    }

    // calculate the reference center of mass
    molMass       = 0;
    refMolCom.x() = 0;
    refMolCom.y() = 0;
    refMolCom.z() = 0;

    for (std::size_t i = 0; i < nIntegrableObjects; i++) {
      refMolCom += refCoords[i] * mass[i];
      molMass += mass[i];
    }

    refMolCom /= molMass;

    // move the reference center of mass to (0,0,0) and adjust the
    // reference coordinate of the integrabel objects
    for (std::size_t i = 0; i < nIntegrableObjects; i++)
      refCoords[i] -= refMolCom;
  }

  RealType MoLocator::getAtomMass(const std::string& at, ForceField* myFF) {
    RealType mass;
    AtomType* atomType = myFF->getAtomType(at);
    if (atomType != NULL) {
      mass = atomType->getMass();
    } else {
      mass = 0.0;
      std::cerr << "Can not find AtomType: " << at << std::endl;
    }
    return mass;
  }

  RealType MoLocator::getMolMass(MoleculeStamp* molStamp, ForceField* myFF) {
    unsigned int nAtoms;
    RealType totMass = 0;
    nAtoms           = molStamp->getNAtoms();

    for (std::size_t i = 0; i < nAtoms; i++) {
      AtomStamp* currAtomStamp = molStamp->getAtomStamp(i);
      totMass += getAtomMass(currAtomStamp->getType(), myFF);
    }
    return totMass;
  }

  RotMat3x3d MoLocator::latVec2RotMat(const Vector3d& lv) {
    RealType theta = acos(lv[2]);
    RealType phi   = atan2(lv[1], lv[0]);
    RealType psi   = 0;

    return RotMat3x3d(phi, theta, psi);
  }
}  // namespace OpenMD
