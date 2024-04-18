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
#include "integrators/LDForceModifier.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>

#include "brains/ForceModifier.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "hydrodynamics/HydroIO.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "math/CholeskyDecomposition.hpp"
#include "math/SquareMatrix3.hpp"
#include "types/GayBerneAdapter.hpp"
#include "types/LennardJonesAdapter.hpp"
#include "utils/Constants.hpp"
#include "utils/ElementsTable.hpp"

using namespace std;
namespace OpenMD {

  LDForceModifier::LDForceModifier(SimInfo* info) :
      ForceModifier {info}, maxIterNum_ {6}, forceTolerance_ {1e-6},
      randNumGen_ {info->getRandomNumberGenerator()},
      simParams_ {info->getSimParams()} {
    RealType dt = simParams_->getDt();
    dt2_        = dt * 0.5;

    veloMunge_ = std::make_unique<Velocitizer>(info_);

    sphericalBoundaryConditions_ = false;
    if (simParams_->getUseSphericalBoundaryConditions()) {
      sphericalBoundaryConditions_ = true;
      if (simParams_->haveLangevinBufferRadius()) {
        langevinBufferRadius_ = simParams_->getLangevinBufferRadius();
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "langevinBufferRadius must be specified "
                 "when useSphericalBoundaryConditions is turned on.\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }

      if (simParams_->haveFrozenBufferRadius()) {
        frozenBufferRadius_ = simParams_->getFrozenBufferRadius();
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "frozenBufferRadius must be specified "
                 "when useSphericalBoundaryConditions is turned on.\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }

      if (frozenBufferRadius_ < langevinBufferRadius_) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "frozenBufferRadius has been set smaller than the "
                 "langevinBufferRadius.  This is probably an error.\n");
        painCave.severity = OPENMD_WARNING;
        painCave.isFatal  = 0;
        simError();
      }
    }

    // Build the hydroProp_ map:
    Molecule* mol;
    StuntDouble* sd;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    bool needHydroPropFile = false;

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        if (sd->isRigidBody()) {
          RigidBody* rb = static_cast<RigidBody*>(sd);
          if (rb->getNumAtoms() > 1) needHydroPropFile = true;
        }
      }
    }

    if (needHydroPropFile) {
      if (simParams_->haveHydroPropFile()) {
        HydroIO* hio  = new HydroIO();
        hydroPropMap_ = hio->parseHydroFile(simParams_->getHydroPropFile());
      } else {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "HydroPropFile must be set to a file name if Langevin Dynamics\n"
            "\tis specified for rigidBodies which contain more than one atom\n"
            "\tTo create a HydroPropFile, run the \"Hydro\" program.\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          map<string, HydroProp*>::iterator iter =
              hydroPropMap_.find(sd->getType());
          if (iter != hydroPropMap_.end()) {
            hydroProps_.push_back(iter->second);
            moments_.push_back(getMomentData(sd));

          } else {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "Can not find resistance tensor for atom [%s]\n",
                     sd->getType().c_str());
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          }
        }
      }

    } else {
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          Shape* currShape = NULL;

          if (sd->isAtom()) {
            Atom* atom          = static_cast<Atom*>(sd);
            AtomType* atomType  = atom->getAtomType();
            GayBerneAdapter gba = GayBerneAdapter(atomType);
            if (gba.isGayBerne()) {
              currShape = new Ellipsoid(V3Zero, gba.getL() / 2.0,
                                        gba.getD() / 2.0, Mat3x3d::identity());
            } else {
              LennardJonesAdapter lja = LennardJonesAdapter(atomType);
              if (lja.isLennardJones()) {
                currShape = new Sphere(V3Zero, lja.getSigma() / 2.0);
              } else {
                int aNum(0);
                vector<AtomType*> atChain = atomType->allYourBase();
                vector<AtomType*>::iterator i;
                for (i = atChain.begin(); i != atChain.end(); ++i) {
                  aNum = etab.GetAtomicNum((*i)->getName().c_str());
                  if (aNum != 0) {
                    currShape = new Sphere(V3Zero, etab.GetVdwRad(aNum));
                    break;
                  }
                }
                if (aNum == 0) {
                  snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                           "Could not find atom type in default element.txt\n");
                  painCave.severity = OPENMD_ERROR;
                  painCave.isFatal  = 1;
                  simError();
                }
              }
            }
          }

          if (!simParams_->haveTargetTemp()) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "You can't use LangevinDynamics without a targetTemp!\n");
            painCave.isFatal  = 1;
            painCave.severity = OPENMD_ERROR;
            simError();
          }

          if (!simParams_->haveViscosity()) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "You can't use LangevinDynamics without a viscosity!\n");
            painCave.isFatal  = 1;
            painCave.severity = OPENMD_ERROR;
            simError();
          }

          HydroProp* currHydroProp =
              currShape->getHydroProp(simParams_->getViscosity());
          map<string, HydroProp*>::iterator iter =
              hydroPropMap_.find(sd->getType());
          if (iter != hydroPropMap_.end()) {
            hydroProps_.push_back(iter->second);
            moments_.push_back(getMomentData(sd));

          } else {
            currHydroProp->complete();
            hydroPropMap_.insert(map<string, HydroProp*>::value_type(
                sd->getType(), currHydroProp));
            hydroProps_.push_back(currHydroProp);
            moments_.push_back(getMomentData(sd));
          }
          delete currShape;
        }
      }
    }

    RealType stdDev =
        std::sqrt(2.0 * Constants::kb * simParams_->getTargetTemp() / dt);

    forceDistribution_ = std::normal_distribution<RealType>(0.0, stdDev);
  }

  MomentData* LDForceModifier::getMomentData(StuntDouble* sd) {
    map<string, MomentData*>::iterator j = momentsMap_.find(sd->getType());
    if (j != momentsMap_.end()) {
      return j->second;
    } else {
      MomentData* moment                  = new MomentData;
      map<string, HydroProp*>::iterator i = hydroPropMap_.find(sd->getType());

      if (i != hydroPropMap_.end()) {
        // Center of Resistance:
        moment->rcr = i->second->getCenterOfResistance();
      } else {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "LDForceManager createMomentData: Couldn't find HydroProp for\n"
            "object type %s!\n",
            sd->getType().c_str());
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      }

      if (sd->isRigidBody())
        moment->rcr -= dynamic_cast<RigidBody*>(sd)->getRefCOM();

      // Parallel axis formula to get the moment of inertia around
      // the center of resistance:
      RealType mass = sd->getMass();
      moment->Icr   = sd->getI();
      moment->Icr +=
          mass * (dot(moment->rcr, moment->rcr) * Mat3x3d::identity() +
                  outProduct(moment->rcr, moment->rcr));
      moment->IcrInv = moment->Icr.inverse();

      momentsMap_.insert(
          map<string, MomentData*>::value_type(sd->getType(), moment));
      return moment;
    }
  }

  void LDForceModifier::modifyForces() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;
    RealType mass;
    Vector3d pos;
    Vector3d frc;
    Mat3x3d A;
    Mat3x3d Atrans;
    Vector3d Tb;
    Vector3d ji;
    unsigned int index = 0;
    bool doLangevinForces;
    bool freezeMolecule;
    int fdf;

    fdf = 0;

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      doLangevinForces = true;
      freezeMolecule   = false;

      if (sphericalBoundaryConditions_) {
        Vector3d molPos = mol->getCom();
        RealType molRad = molPos.length();

        doLangevinForces = false;

        if (molRad > langevinBufferRadius_) {
          doLangevinForces = true;
          freezeMolecule   = false;
        }
        if (molRad > frozenBufferRadius_) {
          doLangevinForces = false;
          freezeMolecule   = true;
        }
      }

      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        if (freezeMolecule) fdf += sd->freeze();

        if (doLangevinForces) {
          mass = sd->getMass();
          if (sd->isDirectional()) {
            // preliminaries for directional objects:

            A      = sd->getA();
            Atrans = A.transpose();

            Vector3d rcrLab = Atrans * moments_[index]->rcr;

            // apply random force and torque at center of resistance

            Vector3d randomForceBody;
            Vector3d randomTorqueBody;
            genRandomForceAndTorque(randomForceBody, randomTorqueBody, index);
            Vector3d randomForceLab  = Atrans * randomForceBody;
            Vector3d randomTorqueLab = Atrans * randomTorqueBody;

            sd->addFrc(randomForceLab);
            sd->addTrq(randomTorqueLab + cross(rcrLab, randomForceLab));

            Vector3d omegaBody;

            // What remains contains velocity explicitly, but the
            // velocity required is at the full step: v(t + h), while
            // we have initially the velocity at the half step: v(t + h/2).
            // We need to iterate to converge the friction
            // force and friction torque vectors.

            // this is the velocity at the half-step:

            Vector3d vel    = sd->getVel();
            Vector3d angMom = sd->getJ();

            // estimate velocity at full-step using everything but
            // friction forces:

            frc = sd->getFrc();

            Vector3d velStep =
                vel + (dt2_ / mass * Constants::energyConvert) * frc;

            Tb = sd->lab2Body(sd->getTrq());
            Vector3d angMomStep =
                angMom + (dt2_ * Constants::energyConvert) * Tb;

            Vector3d omegaLab;
            Vector3d vcdLab;
            Vector3d vcdBody;
            Vector3d frictionForceBody;
            Vector3d frictionForceLab(0.0);
            Vector3d oldFFL;  // used to test for convergence
            Vector3d frictionTorqueBody(0.0);
            Vector3d oldFTB;  // used to test for convergence
            Vector3d frictionTorqueLab;
            RealType fdot;
            RealType tdot;

            // iteration starts here:

            for (int k = 0; k < maxIterNum_; k++) {
              if (sd->isLinear()) {
                int linearAxis = sd->linearAxis();
                int l          = (linearAxis + 1) % 3;
                int m          = (linearAxis + 2) % 3;
                omegaBody[l]   = angMomStep[l] / moments_[index]->Icr(l, l);
                omegaBody[m]   = angMomStep[m] / moments_[index]->Icr(m, m);

              } else {
                omegaBody = moments_[index]->IcrInv * angMomStep;
                // omegaBody[0] = angMomStep[0] /I(0, 0);
                // omegaBody[1] = angMomStep[1] /I(1, 1);
                // omegaBody[2] = angMomStep[2] /I(2, 2);
              }

              omegaLab = Atrans * omegaBody;

              // apply friction force and torque at center of resistance

              vcdLab             = velStep + cross(omegaLab, rcrLab);
              vcdBody            = A * vcdLab;
              frictionForceBody  = -(hydroProps_[index]->getXitt() * vcdBody +
                                    hydroProps_[index]->getXirt() * omegaBody);
              oldFFL             = frictionForceLab;
              frictionForceLab   = Atrans * frictionForceBody;
              oldFTB             = frictionTorqueBody;
              frictionTorqueBody = -(hydroProps_[index]->getXitr() * vcdBody +
                                     hydroProps_[index]->getXirr() * omegaBody);
              frictionTorqueLab  = Atrans * frictionTorqueBody;

              // re-estimate velocities at full-step using friction forces:

              velStep = vel + (dt2_ / mass * Constants::energyConvert) *
                                  (frc + frictionForceLab);
              angMomStep = angMom + (dt2_ * Constants::energyConvert) *
                                        (Tb + frictionTorqueBody);

              // check for convergence (if the vectors have converged, fdot and
              // tdot will both be 1.0):

              fdot = dot(frictionForceLab, oldFFL) /
                     frictionForceLab.lengthSquare();
              tdot = dot(frictionTorqueBody, oldFTB) /
                     frictionTorqueBody.lengthSquare();

              if (fabs(1.0 - fdot) <= forceTolerance_ &&
                  fabs(1.0 - tdot) <= forceTolerance_)
                break;  // iteration ends here
            }

            sd->addFrc(frictionForceLab);
            sd->addTrq(frictionTorqueLab + cross(rcrLab, frictionForceLab));

          } else {
            // spherical atom

            Vector3d systemForce;
            Vector3d randomForce;
            Vector3d randomTorque;
            genRandomForceAndTorque(randomForce, randomTorque, index);
            systemForce = sd->getFrc();
            sd->addFrc(randomForce);

            // What remains contains velocity explicitly, but the
            // velocity required is at the full step: v(t + h), while
            // we have initially the velocity at the half step: v(t + h/2).
            // We need to iterate to converge the friction
            // force vector.

            // this is the velocity at the half-step:

            Vector3d vel = sd->getVel();

            // estimate velocity at full-step using everything but
            // friction forces:

            frc = sd->getFrc();
            Vector3d velStep =
                vel + (dt2_ / mass * Constants::energyConvert) * frc;

            Vector3d frictionForce(0.0);
            Vector3d oldFF;  // used to test for convergence
            RealType fdot;

            // iteration starts here:

            for (int k = 0; k < maxIterNum_; k++) {
              oldFF         = frictionForce;
              frictionForce = -hydroProps_[index]->getXitt() * velStep;

              // re-estimate velocities at full-step using friction forces:

              velStep = vel + (dt2_ / mass * Constants::energyConvert) *
                                  (frc + frictionForce);

              // check for convergence (if the vector has converged,
              // fdot will be 1.0):

              fdot = dot(frictionForce, oldFF) / frictionForce.lengthSquare();

              if (fabs(1.0 - fdot) <= forceTolerance_)
                break;  // iteration ends here
            }

            sd->addFrc(frictionForce);
          }
        }

        ++index;
      }
    }

    info_->setFdf(fdf);
    if(simParams_->getConserveLinearMomentum()) veloMunge_->removeComDrift();
    // Remove angular drift if we are not using periodic boundary conditions.
    if (!simParams_->getUsePeriodicBoundaryConditions() &&
        simParams_->getConserveAngularMomentum())
      veloMunge_->removeAngularDrift();
  }

  void LDForceModifier::genRandomForceAndTorque(Vector3d& force,
                                                Vector3d& torque,
                                                unsigned int index) {
    Vector<RealType, 6> Z;
    Vector<RealType, 6> generalForce;

    Z[0] = forceDistribution_(*randNumGen_);
    Z[1] = forceDistribution_(*randNumGen_);
    Z[2] = forceDistribution_(*randNumGen_);
    Z[3] = forceDistribution_(*randNumGen_);
    Z[4] = forceDistribution_(*randNumGen_);
    Z[5] = forceDistribution_(*randNumGen_);

    generalForce = hydroProps_[index]->getS() * Z;

    force[0]  = generalForce[0];
    force[1]  = generalForce[1];
    force[2]  = generalForce[2];
    torque[0] = generalForce[3];
    torque[1] = generalForce[4];
    torque[2] = generalForce[5];
  }
}  // namespace OpenMD
