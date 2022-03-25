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
 * @file LangevinDynamics.cpp
 * @author tlin
 * @date 11/08/2004
 * @version 1.0
 */

#include "integrators/LangevinDynamics.hpp"

#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"

namespace OpenMD {

  LangevinDynamics::LangevinDynamics(SimInfo* info) :
      VelocityVerletIntegrator(info) {}

  void LangevinDynamics::moveA() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d vel;
    Vector3d pos;
    Vector3d frc;
    Vector3d Tb;
    Vector3d ji;
    RealType mass;

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        vel  = sd->getVel();
        pos  = sd->getPos();
        frc  = sd->getFrc();
        mass = sd->getMass();

        // velocity half step
        vel += (dt2 / mass * Constants::energyConvert) * frc;

        // position whole step
        pos += dt * vel;

        sd->setVel(vel);
        sd->setPos(pos);

        if (sd->isDirectional()) {
          // get and convert the torque to body frame

          Tb = sd->lab2Body(sd->getTrq());

          // get the angular momentum, and propagate a half step

          ji = sd->getJ();

          ji += (dt2 * Constants::energyConvert) * Tb;

          rotAlgo_->rotate(sd, ji, dt);

          sd->setJ(ji);
        }
      }
    }  // end for(mol = info_->beginMolecule(i))

    flucQ_->moveA();
    rattle_->constraintA();
  }

  void LangevinDynamics::moveB() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d vel;
    Vector3d frc;
    Vector3d Tb;
    Vector3d ji;
    RealType mass;

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        vel  = sd->getVel();
        frc  = sd->getFrc();
        mass = sd->getMass();

        // velocity half step
        vel += (dt2 / mass * Constants::energyConvert) * frc;

        sd->setVel(vel);

        if (sd->isDirectional()) {
          // get and convert the torque to body frame

          Tb = sd->lab2Body(sd->getTrq());

          // get the angular momentum, and propagate a half step

          ji = sd->getJ();

          ji += (dt2 * Constants::energyConvert) * Tb;

          sd->setJ(ji);
        }
      }
    }  // end for(mol = info_->beginMolecule(i))

    flucQ_->moveB();
    rattle_->constraintB();
  }

  RealType LangevinDynamics::calcConservedQuantity() { return 0.0; }
}  // namespace OpenMD
