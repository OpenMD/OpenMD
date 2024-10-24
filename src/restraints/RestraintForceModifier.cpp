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

#include "restraints/RestraintForceModifier.hpp"

#include <config.h>

#include <cmath>
#include <memory>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/ForceModifier.hpp"
#include "io/RestReader.hpp"
#include "restraints/MolecularRestraint.hpp"
#include "restraints/ObjectRestraint.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/Constants.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {

  RestraintForceModifier::RestraintForceModifier(SimInfo* info) :
      ForceModifier {info} {
    // order of affairs:
    //
    // 1) create restraints from the restraintStamps found in the MD
    // file.
    //
    // 2) Create RestraintReader to parse the input files for the ideal
    // structures.  This reader will set reference structures, and will
    // calculate molecular centers of mass, etc.
    //
    // 3) sit around and wait for calcForces to be called.  When it comes,
    // call the normal force manager calcForces, then loop through the
    // restrained objects and do their restraint forces.

    Globals* simParam = info_->getSimParams();
    currSnapshot_     = info_->getSnapshotManager()->getCurrentSnapshot();

    currRestTime_ = currSnapshot_->getTime();

    if (simParam->haveStatusTime()) {
      restTime_ = simParam->getStatusTime();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Restraint warning: If you use restraints without setting\n"
               "\tstatusTime, no restraint data will be written to the rest\n"
               "\tfile.\n");
      painCave.isFatal = 0;
      simError();
      restTime_ = simParam->getRunTime();
    }

    int nRestraintStamps               = simParam->getNRestraintStamps();
    std::vector<RestraintStamp*> stamp = simParam->getRestraintStamps();

    std::vector<int> stuntDoubleIndex;

    for (int i = 0; i < nRestraintStamps; i++) {
      std::string myType = toUpperCopy(stamp[i]->getType());

      if (myType.compare("MOLECULAR") == 0) {
        int molIndex(-1);
        Vector3d refCom;

        if (!stamp[i]->haveMolIndex()) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Restraint Error: A molecular restraint was specified\n"
                   "\twithout providing a value for molIndex.\n");
          painCave.isFatal = 1;
          simError();
        } else {
          molIndex = stamp[i]->getMolIndex();
        }

        if (molIndex < 0) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Restraint Error: A molecular restraint was specified\n"
                   "\twith a molIndex that was less than 0\n");
          painCave.isFatal = 1;
          simError();
        }
        if (molIndex >= info_->getNGlobalMolecules()) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Restraint Error: A molecular restraint was specified with\n"
                   "\ta molIndex that was greater than the total number of "
                   "molecules\n");
          painCave.isFatal = 1;
          simError();
        }

        Molecule* mol = info_->getMoleculeByGlobalIndex(molIndex);

        if (mol == NULL) {
#ifdef IS_MPI
          // getMoleculeByGlobalIndex returns a NULL in parallel if
          // this proc doesn't have the molecule.  Do a quick check to
          // make sure another processor is supposed to have it.

          int myrank;
          MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

          if (info_->getMolToProc(molIndex) == myrank) {
            // If we were supposed to have it but got a null, then freak out.
#endif

            snprintf(
                painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                "Restraint Error: A molecular restraint was specified, but\n"
                "\tno molecule was found with global index %d.\n",
                molIndex);
            painCave.isFatal = 1;
            simError();

#ifdef IS_MPI
          }
#endif
        }

#ifdef IS_MPI
        // only handle this molecular restraint if this processor owns the
        // molecule
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if (info_->getMolToProc(molIndex) == myrank) {
#endif

          MolecularRestraint* rest = new MolecularRestraint();

          std::string restPre("mol_");
          std::stringstream restName;
          restName << restPre << molIndex;
          rest->setRestraintName(restName.str());

          if (stamp[i]->haveDisplacementSpringConstant()) {
            rest->setDisplacementForceConstant(
                stamp[i]->getDisplacementSpringConstant());
          }
          if (stamp[i]->haveAbsoluteSpringConstant()) {
            rest->setAbsoluteForceConstant(
                stamp[i]->getAbsoluteSpringConstant());
          }
          if (stamp[i]->haveTwistSpringConstant()) {
            rest->setTwistForceConstant(stamp[i]->getTwistSpringConstant());
          }
          if (stamp[i]->haveSwingXSpringConstant()) {
            rest->setSwingXForceConstant(stamp[i]->getSwingXSpringConstant());
          }
          if (stamp[i]->haveSwingYSpringConstant()) {
            rest->setSwingYForceConstant(stamp[i]->getSwingYSpringConstant());
          }
          if (stamp[i]->haveAbsolutePositionZ()) {
            rest->setAbsolutePositionZ(stamp[i]->getAbsolutePositionZ());
          }
          if (stamp[i]->haveRestrainedTwistAngle()) {
            rest->setRestrainedTwistAngle(stamp[i]->getRestrainedTwistAngle() *
                                          Constants::PI / 180.0);
          }
          if (stamp[i]->haveRestrainedSwingYAngle()) {
            rest->setRestrainedSwingYAngle(
                stamp[i]->getRestrainedSwingYAngle() * Constants::PI / 180.0);
          }
          if (stamp[i]->haveRestrainedSwingXAngle()) {
            rest->setRestrainedSwingXAngle(
                stamp[i]->getRestrainedSwingXAngle() * Constants::PI / 180.0);
          }
          if (stamp[i]->havePrint()) {
            rest->setPrintRestraint(stamp[i]->getPrint());
          }

          restraints_.push_back(rest);
          mol->addProperty(std::make_shared<RestraintData>("Restraint", rest));
          restrainedMols_.push_back(mol);
#ifdef IS_MPI
        }
#endif
      } else if (myType.compare("OBJECT") == 0) {
        std::string objectSelection;

        if (!stamp[i]->haveObjectSelection()) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Restraint Error: An object restraint was specified\n"
                   "\twithout providing a selection script in the\n"
                   "\tobjectSelection variable.\n");
          painCave.isFatal = 1;
          simError();
        } else {
          objectSelection = stamp[i]->getObjectSelection();
        }

        SelectionEvaluator evaluator(info);
        SelectionManager seleMan(info);

        evaluator.loadScriptString(objectSelection);
        seleMan.setSelectionSet(evaluator.evaluate());
        int selectionCount = seleMan.getSelectionCount();

#ifdef IS_MPI
        MPI_Allreduce(MPI_IN_PLACE, &selectionCount, 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);
#endif

        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Restraint Info: The specified restraint objectSelection,\n"
                 "\t\t%s\n"
                 "\twill result in %d integrable objects being\n"
                 "\trestrained.\n",
                 objectSelection.c_str(), selectionCount);
        painCave.severity = OPENMD_INFO;
        painCave.isFatal  = 0;
        simError();

        int selei;
        StuntDouble* sd;

        for (sd = seleMan.beginSelected(selei); sd != NULL;
             sd = seleMan.nextSelected(selei)) {
          stuntDoubleIndex.push_back(sd->getGlobalIntegrableObjectIndex());

          ObjectRestraint* rest = new ObjectRestraint();

          if (stamp[i]->haveDisplacementSpringConstant()) {
            rest->setDisplacementForceConstant(
                stamp[i]->getDisplacementSpringConstant());
          }
          if (stamp[i]->haveAbsoluteSpringConstant()) {
            rest->setAbsoluteForceConstant(
                stamp[i]->getAbsoluteSpringConstant());
          }
          if (stamp[i]->haveTwistSpringConstant()) {
            rest->setTwistForceConstant(stamp[i]->getTwistSpringConstant());
          }
          if (stamp[i]->haveSwingXSpringConstant()) {
            rest->setSwingXForceConstant(stamp[i]->getSwingXSpringConstant());
          }
          if (stamp[i]->haveSwingYSpringConstant()) {
            rest->setSwingYForceConstant(stamp[i]->getSwingYSpringConstant());
          }
          if (stamp[i]->haveAbsolutePositionZ()) {
            rest->setAbsolutePositionZ(stamp[i]->getAbsolutePositionZ());
          }
          if (stamp[i]->haveRestrainedTwistAngle()) {
            rest->setRestrainedTwistAngle(stamp[i]->getRestrainedTwistAngle());
          }
          if (stamp[i]->haveRestrainedSwingXAngle()) {
            rest->setRestrainedSwingXAngle(
                stamp[i]->getRestrainedSwingXAngle());
          }
          if (stamp[i]->haveRestrainedSwingYAngle()) {
            rest->setRestrainedSwingYAngle(
                stamp[i]->getRestrainedSwingYAngle());
          }
          if (stamp[i]->havePrint()) {
            rest->setPrintRestraint(stamp[i]->getPrint());
          }

          restraints_.push_back(rest);
          sd->addProperty(std::make_shared<RestraintData>("Restraint", rest));
          restrainedObjs_.push_back(sd);
        }
      }
    }

    // ThermodynamicIntegration subclasses RestraintForceManager, and there
    // are times when it won't use restraints at all, so only open the
    // restraint file if we are actually using restraints:
    if (simParam->getUseRestraints()) {
      std::string refFile = simParam->getRestraint_file();
      RestReader* rr      = new RestReader(info, refFile, stuntDoubleIndex);
      rr->readReferenceStructure();
      delete rr;
    }

    restOutput_ = getPrefix(info_->getFinalConfigFileName()) + ".rest";
    restOut     = new RestWriter(info_, restOutput_.c_str(), restraints_);
    if (!restOut) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Restraint error: Failed to create RestWriter\n");
      painCave.isFatal = 1;
      simError();
    }

    // todo: figure out the scale factor.  Right now, just scale them all to 1
    std::vector<Restraint*>::const_iterator resti;
    for (resti = restraints_.begin(); resti != restraints_.end(); ++resti) {
      (*resti)->setScaleFactor(1.0);
    }
  }

  RestraintForceModifier::~RestraintForceModifier() {
    Utils::deletePointers(restraints_);

    delete restOut;
  }

  void RestraintForceModifier::modifyForces() {
    RealType restPot(0.0);

    restPot = doRestraints(1.0);

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &restPot, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    currSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    RealType rp   = currSnapshot_->getRestraintPotential();
    currSnapshot_->setRestraintPotential(rp + restPot);

    RealType pe = currSnapshot_->getPotentialEnergy();
    currSnapshot_->setRawPotential(pe);
    currSnapshot_->setPotentialEnergy(pe + restPot);

    // write out forces and current positions of restrained molecules
    if (currSnapshot_->getTime() >= currRestTime_) {
      restOut->writeRest(restInfo_);
      currRestTime_ += restTime_;
    }
  }

  RealType RestraintForceModifier::doRestraints(RealType scalingFactor) {
    std::vector<Molecule*>::const_iterator rm;
    std::shared_ptr<GenericData> data;
    Molecule::IntegrableObjectIterator ioi;
    MolecularRestraint* mRest = NULL;
    ObjectRestraint* oRest    = NULL;
    StuntDouble* sd;

    std::vector<StuntDouble*>::const_iterator ro;

    std::map<int, Restraint::RealPair> restInfo;

    unscaledPotential_ = 0.0;

    restInfo_.clear();

    for (rm = restrainedMols_.begin(); rm != restrainedMols_.end(); ++rm) {
      // make sure this molecule (*rm) has a generic data for restraints:
      data = (*rm)->getPropertyByName("Restraint");
      if (data != nullptr) {
        // make sure we can reinterpret the generic data as restraint data:
        std::shared_ptr<RestraintData> restData =
            std::dynamic_pointer_cast<RestraintData>(data);
        if (restData != nullptr) {
          // make sure we can reinterpet the restraint data as a
          // MolecularRestraint
          mRest = dynamic_cast<MolecularRestraint*>(restData->getData());
          if (mRest == NULL) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "Can not cast RestraintData to MolecularRestraint\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          }
        } else {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not cast GenericData to RestraintData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Can not find Restraint for RestrainedObject\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }

      // phew.  At this point, we should have the pointer to the
      // correct MolecularRestraint in the variable mRest.

      Vector3d molCom = (*rm)->getCom();

      std::vector<Vector3d> struc;
      std::vector<Vector3d> forces;

      for (sd = (*rm)->beginIntegrableObject(ioi); sd != NULL;
           sd = (*rm)->nextIntegrableObject(ioi)) {
        struc.push_back(sd->getPos());
      }

      mRest->setScaleFactor(scalingFactor);
      mRest->calcForce(struc, molCom);
      forces    = mRest->getRestraintForces();
      int index = 0;

      for (sd = (*rm)->beginIntegrableObject(ioi); sd != NULL;
           sd = (*rm)->nextIntegrableObject(ioi)) {
        sd->addFrc(forces[index]);
        struc.push_back(sd->getPos());
        index++;
      }

      unscaledPotential_ += mRest->getUnscaledPotential();

      // only collect data on restraints that we're going to print:
      if (mRest->getPrintRestraint()) {
        restInfo = mRest->getRestraintInfo();
        restInfo_.push_back(restInfo);
      }
    }

    for (ro = restrainedObjs_.begin(); ro != restrainedObjs_.end(); ++ro) {
      // make sure this object (*ro) has a generic data for restraints:
      data = (*ro)->getPropertyByName("Restraint");
      if (data != NULL) {
        // make sure we can reinterpret the generic data as restraint data:
        std::shared_ptr<RestraintData> restData =
            std::dynamic_pointer_cast<RestraintData>(data);
        if (restData != nullptr) {
          // make sure we can reinterpet the restraint data as a pointer to
          // an ObjectRestraint:
          oRest = dynamic_cast<ObjectRestraint*>(restData->getData());
          if (oRest == NULL) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "Can not cast RestraintData to ObjectRestraint\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          }
        } else {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not cast GenericData to RestraintData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Can not find Restraint for RestrainedObject\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }

      // phew.  At this point, we should have the pointer to the
      // correct Object restraint in the variable oRest.
      oRest->setScaleFactor(scalingFactor);

      Vector3d pos = (*ro)->getPos();

      if ((*ro)->isDirectional()) {
        // directional objects may have orientational restraints as well
        // as positional, so get the rotation matrix first:

        RotMat3x3d A = (*ro)->getA();
        oRest->calcForce(pos, A);
        (*ro)->addFrc(oRest->getRestraintForce());
        (*ro)->addTrq(oRest->getRestraintTorque());

      } else {
        // plain vanilla positional restraints:

        oRest->calcForce(pos);
        (*ro)->addFrc(oRest->getRestraintForce());
      }

      unscaledPotential_ += oRest->getUnscaledPotential();

      // only collect data on restraints that we're going to print:
      if (oRest->getPrintRestraint()) {
        restInfo = oRest->getRestraintInfo();
        restInfo_.push_back(restInfo);
      }
    }

    return unscaledPotential_ * scalingFactor;
  }
}  // namespace OpenMD
