/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
#include <fstream>
#include <iostream>
#include "math/SquareMatrix3.hpp"
#include "integrators/LDForceManager.hpp"
#include "math/CholeskyDecomposition.hpp"
#include "utils/Constants.hpp"
#include "utils/MemoryUtils.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "utils/ElementsTable.hpp"
#include "types/LennardJonesAdapter.hpp"
#include "types/GayBerneAdapter.hpp"

using namespace std;
namespace OpenMD {

  LDForceManager::LDForceManager(SimInfo* info) : ForceManager(info),
						  maxIterNum_(4),
						  forceTolerance_(1e-6) {
    simParams = info->getSimParams();

    // Remove in favor of std::MemoryUtils::make_unique<> when we switch to C++14 and above
    veloMunge = MemoryUtils::make_unique<Velocitizer>(info);

    sphericalBoundaryConditions_ = false;
    if (simParams->getUseSphericalBoundaryConditions()) {
      sphericalBoundaryConditions_ = true;
      if (simParams->haveLangevinBufferRadius()) {
        langevinBufferRadius_ = simParams->getLangevinBufferRadius();
      } else {
        sprintf( painCave.errMsg,
                 "langevinBufferRadius must be specified "
                 "when useSphericalBoundaryConditions is turned on.\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      }

      if (simParams->haveFrozenBufferRadius()) {
        frozenBufferRadius_ = simParams->getFrozenBufferRadius();
      } else {
        sprintf( painCave.errMsg,
                 "frozenBufferRadius must be specified "
                 "when useSphericalBoundaryConditions is turned on.\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      }

      if (frozenBufferRadius_ < langevinBufferRadius_) {
        sprintf( painCave.errMsg,
                 "frozenBufferRadius has been set smaller than the "
                 "langevinBufferRadius.  This is probably an error.\n");
        painCave.severity = OPENMD_WARNING;
        painCave.isFatal = 0;
        simError();
      }
    }

    // Build the hydroProp_ map:

    Molecule* mol;
    StuntDouble* sd;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    bool needHydroPropFile = false;

    for (mol = info->beginMolecule(i); mol != NULL;
         mol = info->nextMolecule(i)) {

      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {

        if (sd->isRigidBody()) {
          RigidBody* rb = static_cast<RigidBody*>(sd);
          if (rb->getNumAtoms() > 1) needHydroPropFile = true;
        }

      }
    }

    if (needHydroPropFile) {
      if (simParams->haveHydroPropFile()) {
        hydroPropMap_ = parseFrictionFile(simParams->getHydroPropFile());
      } else {
        sprintf( painCave.errMsg,
                 "HydroPropFile must be set to a file name if Langevin Dynamics\n"
                 "\tis specified for rigidBodies which contain more than one atom\n"
                 "\tTo create a HydroPropFile, run the \"Hydro\" program.\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      }
       
      for (mol = info->beginMolecule(i); mol != NULL;
           mol = info->nextMolecule(i)) {

        for (sd = mol->beginIntegrableObject(j);  sd != NULL;
             sd = mol->nextIntegrableObject(j)) {

          map<string, HydroProp*>::iterator iter = hydroPropMap_.find(sd->getType());
          if (iter != hydroPropMap_.end()) {

            hydroProps_.push_back(iter->second);
            moments_.push_back(getMomentData(sd));

          } else {
            sprintf( painCave.errMsg,
                     "Can not find resistance tensor for atom [%s]\n",
                     sd->getType().c_str());
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();
          }
        }
      }
      
    } else {

      for (mol = info->beginMolecule(i); mol != NULL;
           mol = info->nextMolecule(i)) {

        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {

          Shape* currShape = NULL;

          if (sd->isAtom()){
            Atom* atom = static_cast<Atom*>(sd);
            AtomType* atomType = atom->getAtomType();
            GayBerneAdapter gba = GayBerneAdapter(atomType);
            if (gba.isGayBerne()) {
              currShape = new Ellipsoid(V3Zero, gba.getL() / 2.0,
                                        gba.getD() / 2.0,
                                        Mat3x3d::identity());
            } else {
              LennardJonesAdapter lja = LennardJonesAdapter(atomType);
              if (lja.isLennardJones()){
                currShape = new Sphere(V3Zero, lja.getSigma()/2.0);
              } else {

                int aNum(0);
                vector<AtomType*> atChain = atomType->allYourBase();
                vector<AtomType*>::iterator i;
                for (i = atChain.begin(); i != atChain.end(); ++i) {
                  aNum = etab.GetAtomicNum((*i)->getName().c_str());
                  if (aNum != 0) {
                    currShape = new Sphere(V3Zero,
                                           etab.GetVdwRad(aNum));
                    break;
                  }
                }
                if (aNum == 0) {
                  sprintf( painCave.errMsg,
                           "Could not find atom type in default element.txt\n");
                  painCave.severity = OPENMD_ERROR;
                  painCave.isFatal = 1;
                  simError();
                }
              }
            }
          }

	  if (!simParams->haveTargetTemp()) {
	    sprintf(painCave.errMsg,
                    "You can't use LangevinDynamics without a targetTemp!\n");
	    painCave.isFatal = 1;
	    painCave.severity = OPENMD_ERROR;
	    simError();
	  }

	  if (!simParams->haveViscosity()) {
	    sprintf(painCave.errMsg,
                    "You can't use LangevinDynamics without a viscosity!\n");
	    painCave.isFatal = 1;
	    painCave.severity = OPENMD_ERROR;
	    simError();
	  }


          HydroProp* currHydroProp = currShape->getHydroProp(simParams->getViscosity(),simParams->getTargetTemp());
          map<string, HydroProp*>::iterator iter = hydroPropMap_.find(sd->getType());
          if (iter != hydroPropMap_.end()) {
            
            hydroProps_.push_back(iter->second);
            moments_.push_back(getMomentData(sd));
            
          } else {
            currHydroProp->complete();
            hydroPropMap_.insert(map<string, HydroProp*>::value_type(sd->getType(), currHydroProp));
            hydroProps_.push_back(currHydroProp);
            moments_.push_back(getMomentData(sd));
          }
          delete currShape;
        }
      }
    }
    variance_ = 2.0 * Constants::kb*simParams->getTargetTemp()/simParams->getDt();
  }

  map<string, HydroProp*> LDForceManager::parseFrictionFile(const string& filename) {
    map<string, HydroProp*> props;
    ifstream ifs(filename.c_str());
    if (ifs.is_open()) {

    }

    const unsigned int BufferSize = 65535;
    char buffer[BufferSize];
    while (ifs.getline(buffer, BufferSize)) {
      HydroProp* currProp = new HydroProp(buffer);
      props.insert(map<string, HydroProp*>::value_type(currProp->getName(),
                                                       currProp));
    }

    return props;
  }

  MomentData* LDForceManager::getMomentData(StuntDouble* sd) {

    map<string, MomentData*>::iterator j = momentsMap_.find(sd->getType());
    if (j != momentsMap_.end()) {
      return j->second;
    } else {
    
      MomentData* moment = new MomentData;
      map<string, HydroProp*>::iterator i = hydroPropMap_.find(sd->getType());
      
      if (i != hydroPropMap_.end()) {            
        // Center of Resistance:
        moment->rcr = i->second->getCOR();
      } else {
        sprintf(painCave.errMsg,
                "LDForceManager createMomentData: Couldn't find HydroProp for\n"
                "object type %s!\n", sd->getType().c_str());
        painCave.isFatal = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      }

      if (sd->isRigidBody())
        moment->rcr -= dynamic_cast<RigidBody*>(sd)->getRefCOM();
      
      // Parallel axis formula to get the moment of inertia around
      // the center of resistance:
      RealType mass = sd->getMass();
      moment->Icr = sd->getI();
      moment->Icr += mass*(dot(moment->rcr, moment->rcr) * Mat3x3d::identity() +
                           outProduct(moment->rcr, moment->rcr));
      moment->IcrInv = moment->Icr.inverse();
      
      momentsMap_.insert(map<string, MomentData*>::value_type(sd->getType(),
                                                              moment));
      return moment;
    }
  }
  
  void LDForceManager::postCalculation(){
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
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
      freezeMolecule = false;

      if (sphericalBoundaryConditions_) {

        Vector3d molPos = mol->getCom();
        RealType molRad = molPos.length();

        doLangevinForces = false;

        if (molRad > langevinBufferRadius_) {
          doLangevinForces = true;
          freezeMolecule = false;
        }
        if (molRad > frozenBufferRadius_) {
          doLangevinForces = false;
          freezeMolecule = true;
        }
      }

      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {

        if (freezeMolecule)
          fdf += sd->freeze();

        if (doLangevinForces) {
          mass = sd->getMass();
          if (sd->isDirectional()){

            // preliminaries for directional objects:

            A = sd->getA();
            Atrans = A.transpose();
            
            Vector3d rcrLab = Atrans * moments_[index]->rcr;

            //apply random force and torque at center of resistance

            Vector3d randomForceBody;
            Vector3d randomTorqueBody;
            genRandomForceAndTorque(randomForceBody, randomTorqueBody,
                                    index, variance_);
            Vector3d randomForceLab = Atrans * randomForceBody;
            Vector3d randomTorqueLab = Atrans * randomTorqueBody;
            sd->addFrc(randomForceLab);
            sd->addTrq(randomTorqueLab + cross(rcrLab, randomForceLab ));
           
            Vector3d omegaBody;

            // What remains contains velocity explicitly, but the
            // velocity required is at the full step: v(t + h), while
            // we have initially the velocity at the half step: v(t + h/2).
            // We need to iterate to converge the friction
            // force and friction torque vectors.

            // this is the velocity at the half-step:

            Vector3d vel =sd->getVel();
            Vector3d angMom = sd->getJ();

            // estimate velocity at full-step using everything but
            // friction forces:

            frc = sd->getFrc();
            Vector3d velStep = vel + (dt2_ /mass * Constants::energyConvert) * frc;

            Tb = sd->lab2Body(sd->getTrq());
            Vector3d angMomStep = angMom + (dt2_ * Constants::energyConvert) * Tb;

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

            //iteration starts here:

            for (int k = 0; k < maxIterNum_; k++) {

              if (sd->isLinear()) {
                int linearAxis = sd->linearAxis();
                int l = (linearAxis +1 )%3;
                int m = (linearAxis +2 )%3;
                omegaBody[l] = angMomStep[l] /moments_[index]->Icr(l, l);
                omegaBody[m] = angMomStep[m] /moments_[index]->Icr(m, m);

              } else {
                omegaBody = moments_[index]->IcrInv * angMomStep;
                //omegaBody[0] = angMomStep[0] /I(0, 0);
                //omegaBody[1] = angMomStep[1] /I(1, 1);
                //omegaBody[2] = angMomStep[2] /I(2, 2);
              }

              omegaLab = Atrans * omegaBody;

              // apply friction force and torque at center of resistance

              vcdLab = velStep + cross(omegaLab, rcrLab);
              vcdBody = A * vcdLab;
              frictionForceBody = -(hydroProps_[index]->getXitt() * vcdBody +
                                    hydroProps_[index]->getXirt() * omegaBody);
              oldFFL = frictionForceLab;
              frictionForceLab = Atrans * frictionForceBody;
              oldFTB = frictionTorqueBody;
              frictionTorqueBody = -(hydroProps_[index]->getXitr() * vcdBody +
                                     hydroProps_[index]->getXirr() * omegaBody);
              frictionTorqueLab = Atrans * frictionTorqueBody;

              // re-estimate velocities at full-step using friction forces:

              velStep = vel + (dt2_ / mass * Constants::energyConvert) * (frc + frictionForceLab);
              angMomStep = angMom + (dt2_ * Constants::energyConvert) * (Tb + frictionTorqueBody);

              // check for convergence (if the vectors have converged, fdot and tdot will both be 1.0):

              fdot = dot(frictionForceLab, oldFFL) / frictionForceLab.lengthSquare();
              tdot = dot(frictionTorqueBody, oldFTB) / frictionTorqueBody.lengthSquare();

              if (fabs(1.0 - fdot) <= forceTolerance_ &&
                  fabs(1.0 - tdot) <= forceTolerance_)
                break; // iteration ends here
            }

            sd->addFrc(frictionForceLab);
            sd->addTrq(frictionTorqueLab + cross(rcrLab, frictionForceLab));


          } else {
            //spherical atom

	    Vector3d systemForce;
            Vector3d randomForce;
            Vector3d randomTorque;
            genRandomForceAndTorque(randomForce, randomTorque,
                                    index, variance_);
	    systemForce = sd->getFrc();
            sd->addFrc(randomForce);

            // What remains contains velocity explicitly, but the
            // velocity required is at the full step: v(t + h), while
            // we have initially the velocity at the half step: v(t + h/2).
            // We need to iterate to converge the friction
            // force vector.

            // this is the velocity at the half-step:

            Vector3d vel =sd->getVel();

            // estimate velocity at full-step using everything but
            // friction forces:

            frc = sd->getFrc();
            Vector3d velStep = vel + (dt2_ / mass * Constants::energyConvert) * frc;

            Vector3d frictionForce(0.0);
            Vector3d oldFF;  // used to test for convergence
            RealType fdot;

            //iteration starts here:

            for (int k = 0; k < maxIterNum_; k++) {

              oldFF = frictionForce;
              frictionForce = -hydroProps_[index]->getXitt() * velStep;

              // re-estimate velocities at full-step using friction forces:

              velStep = vel + (dt2_ / mass * Constants::energyConvert) * (frc + frictionForce);

              // check for convergence (if the vector has converged,
              // fdot will be 1.0):

              fdot = dot(frictionForce, oldFF) / frictionForce.lengthSquare();

              if (fabs(1.0 - fdot) <= forceTolerance_)
                break; // iteration ends here
            }

            sd->addFrc(frictionForce);

          }
        }

        ++index;

      }
    }
    info_->setFdf(fdf);
    veloMunge->removeComDrift();
    // Remove angular drift if we are not using periodic boundary conditions.
    if(!simParams->getUsePeriodicBoundaryConditions())
      veloMunge->removeAngularDrift();
    ForceManager::postCalculation();
  }

  void LDForceManager::genRandomForceAndTorque(Vector3d& force,
                                               Vector3d& torque,
                                               unsigned int index,
                                               RealType variance) {
    Vector<RealType, 6> Z;
    Vector<RealType, 6> generalForce;

    Z[0] = randNumGen_.randNorm(0, variance);
    Z[1] = randNumGen_.randNorm(0, variance);
    Z[2] = randNumGen_.randNorm(0, variance);
    Z[3] = randNumGen_.randNorm(0, variance);
    Z[4] = randNumGen_.randNorm(0, variance);
    Z[5] = randNumGen_.randNorm(0, variance);

    generalForce = hydroProps_[index]->getS()*Z;

    force[0] = generalForce[0];
    force[1] = generalForce[1];
    force[2] = generalForce[2];
    torque[0] = generalForce[3];
    torque[1] = generalForce[4];
    torque[2] = generalForce[5];
  }
}
