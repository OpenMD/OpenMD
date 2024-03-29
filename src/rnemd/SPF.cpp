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

#include "rnemd/SPF.hpp"

#include <config.h>

#include <cmath>
#include <random>
#include <vector>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/ForceManager.hpp"
#include "brains/SimInfo.hpp"
#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "rnemd/RNEMD.hpp"
#include "rnemd/RNEMDParameters.hpp"
#include "rnemd/SPFForceManager.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/Constants.hpp"
#include "utils/RandNumGen.hpp"
#include "utils/simError.h"

namespace OpenMD::RNEMD {

  SPFMethod::SPFMethod(SimInfo* info, ForceManager* forceMan) :
      RNEMD {info, forceMan}, selectedMoleculeMan_(info),
      selectedMoleculeEvaluator_(info) {
    rnemdMethodLabel_ = "SPF";

    selectedMoleculeStr_ = "select none";
    selectedMoleculeEvaluator_.loadScriptString(selectedMoleculeStr_);
    selectedMoleculeMan_.setSelectionSet(selectedMoleculeEvaluator_.evaluate());

    if (SPFForceManager* spfForceManager =
            dynamic_cast<SPFForceManager*>(forceMan)) {
      forceManager_ = spfForceManager;
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SPF-RNEMD cannot be used with the default ForceManager.\n");
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    RNEMDParameters* rnemdParams = info->getSimParams()->getRNEMDParameters();

    bool hasParticleFlux = rnemdParams->haveParticleFlux();
    bool hasKineticFlux  = rnemdParams->haveKineticFlux();

    bool methodFluxMismatch = false;
    bool hasCorrectFlux     = false;

    switch (rnemdFluxType_) {
    case rnemdParticle:
      hasCorrectFlux = hasParticleFlux;
      break;
    case rnemdParticleKE:
      hasCorrectFlux = hasParticleFlux && hasKineticFlux;
      break;
    default:
      methodFluxMismatch = true;
      break;
    }

    if (methodFluxMismatch) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RNEMD: The current method, SPF\n"
               "\tcannot be used with the current flux type, %s\n",
               rnemdFluxTypeLabel_.c_str());
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (!hasCorrectFlux) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RNEMD: The current method, SPF, and flux type, %s,\n"
               "\tdid not have the correct flux type specified. Options\n"
               "\tinclude: particleFlux and particleFlux + kineticFlux.\n",
               rnemdFluxTypeLabel_.c_str());
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (hasParticleFlux) {
      setParticleFlux(rnemdParams->getParticleFlux());
    } else {
      setParticleFlux(0.0);
    }

    if (hasKineticFlux) {
      setKineticFlux(rnemdParams->getKineticFlux());
    } else {
      setKineticFlux(0.0);
    }

    uniformKineticScaling_ = rnemdParams->getSPFUniformKineticScaling();
    selectMolecule();
  }

  void SPFMethod::doRNEMDImpl(SelectionManager& smanA,
                              SelectionManager& smanB) {
    if (!doRNEMD_) return;

    // Remove selected molecule from the source selection manager
    if (particleTarget_ > 0.0) {
      smanA -= selectedMoleculeMan_;
    } else {
      smanB -= selectedMoleculeMan_;
    }

    failedLastTrial_ = false;

    int selei {}, selej {};

    StuntDouble* sd;

    Vector3d P_a {V3Zero};
    Vector3d P_b {V3Zero};
    RealType M_a {};
    RealType M_b {};
    RealType K_a {};
    RealType K_b {};

    for (sd = smanA.beginSelected(selei); sd != NULL;
         sd = smanA.nextSelected(selei)) {
      RealType mass = sd->getMass();
      Vector3d vel  = sd->getVel();

      P_a += mass * vel;
      M_a += mass;
      K_a += mass * vel.lengthSquare();

      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I       = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          K_a +=
              angMom[j] * angMom[j] / I(j, j) + angMom[k] * angMom[k] / I(k, k);
        } else {
          K_a += angMom[0] * angMom[0] / I(0, 0) +
                 angMom[1] * angMom[1] / I(1, 1) +
                 angMom[2] * angMom[2] / I(2, 2);
        }
      }
    }

    for (sd = smanB.beginSelected(selej); sd != NULL;
         sd = smanB.nextSelected(selej)) {
      RealType mass = sd->getMass();
      Vector3d vel  = sd->getVel();

      P_b += mass * vel;
      M_b += mass;
      K_b += mass * vel.lengthSquare();

      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I       = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          K_b +=
              angMom[j] * angMom[j] / I(j, j) + angMom[k] * angMom[k] / I(k, k);
        } else {
          K_b += angMom[0] * angMom[0] / I(0, 0) +
                 angMom[1] * angMom[1] / I(1, 1) +
                 angMom[2] * angMom[2] / I(2, 2);
        }
      }
    }

    K_a *= 0.5;
    K_b *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &P_a[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &P_b[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &M_a, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &M_b, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &K_a, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &K_b, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

    bool successfulExchange = false;
    bool doExchange         = false;

    if ((M_a > 0.0) && (M_b > 0.0) &&
        forceManager_->getHasSelectedMolecule()) {  // both slabs are not empty
      Vector3d v_a = P_a / M_a;
      Vector3d v_b = P_b / M_b;

      RealType tempParticleTarget =
          (particleTarget_ != deltaLambda_) ? deltaLambda_ : particleTarget_;

      RealType deltaU = forceManager_->getScaledDeltaU(tempParticleTarget) *
                        Constants::energyConvert;
      RealType a {};
      RealType b {};

      if (uniformKineticScaling_) {
        RealType numerator   = deltaU;
        RealType denominator = K_a + K_b;
        denominator -= 0.5 * M_a * v_a.lengthSquare();
        denominator -= 0.5 * M_b * v_b.lengthSquare();

        RealType a2 = (numerator / denominator) + 1.0;

        if (a2 > 0.0) {
          a = std::sqrt(a2);
          b = a;

          // restrict scaling coefficients
          if ((a > 0.999) && (a < 1.001)) { doExchange = true; }
        }
      } else {
        RealType aNumerator   = deltaU - kineticTarget_;
        RealType aDenominator = 2.0 * K_a;
        aDenominator -= M_a * v_a.lengthSquare();

        RealType bNumerator   = deltaU + kineticTarget_;
        RealType bDenominator = 2.0 * K_b;
        bDenominator -= M_b * v_b.lengthSquare();

        RealType a2 = (aNumerator / aDenominator) + 1.0;
        RealType b2 = (bNumerator / bDenominator) + 1.0;

        if (a2 > 0.0 && b2 > 0.0) {
          a = std::sqrt(a2);
          b = std::sqrt(b2);

          // restrict scaling coefficients
          if ((a > 0.999) && (a < 1.001) && (b > 0.999) && (b < 1.001)) {
            doExchange = true;
          }
        }
      }

      if (doExchange) {
        Vector3d vel;

        int selei2 {}, selej2 {};
        StuntDouble* sd2;

        for (sd2 = smanA.beginSelected(selei2); sd2 != NULL;
             sd2 = smanA.nextSelected(selei2)) {
          vel = (sd2->getVel() - v_a) * a + v_a;

          sd2->setVel(vel);

          if (sd2->isDirectional()) {
            Vector3d angMom = sd2->getJ() * a;
            sd2->setJ(angMom);
          }
        }

        for (sd2 = smanB.beginSelected(selei2); sd2 != NULL;
             sd2 = smanB.nextSelected(selei2)) {
          vel = (sd2->getVel() - v_b) * b + v_b;

          sd2->setVel(vel);

          if (sd2->isDirectional()) {
            Vector3d angMom = sd2->getJ() * b;
            sd2->setJ(angMom);
          }
        }

        currentSnap_->hasTranslationalKineticEnergy = false;
        currentSnap_->hasRotationalKineticEnergy    = false;
        currentSnap_->hasKineticEnergy              = false;
        currentSnap_->hasTotalEnergy                = false;

        RealType deltaLambda = particleTarget_;

        bool updateSelectedMolecule =
            forceManager_->updateLambda(deltaLambda, deltaLambda_);

        if (updateSelectedMolecule) selectMolecule();

        successfulExchange = true;
        particleExchange_ += deltaLambda;
        kineticExchange_ += kineticTarget_;
      }
    }

    if (!forceManager_->getHasSelectedMolecule()) {
      selectMolecule();
      deltaLambda_ = particleTarget_;
      failTrialCount_++;
      failedLastTrial_ = true;
      return;
    }

    if (!successfulExchange) {
      //   snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
      //           "SPF exchange NOT performed - roots that solve\n"
      //           "\tthe constraint equations may not exist or there may
      //           be\n"
      //           "\tno selected objects in one or both slabs.\n");
      //   painCave.isFatal  = 0;
      //   painCave.severity = OPENMD_INFO;
      //   simError();
      failTrialCount_++;
      failedLastTrial_ = true;
    }
  }

  void SPFMethod::selectMolecule() {
    bool hasSelectedMolecule = (currentSnap_->getSPFData()->globalID == -1) ?
                                   setSelectedMolecule() :
                                   getSelectedMolecule();

    if (hasSelectedMolecule) {
      selectedMoleculeStr_ =
          "select " + std::to_string(currentSnap_->getSPFData()->globalID);

      selectedMoleculeEvaluator_.loadScriptString(selectedMoleculeStr_);
      selectedMoleculeMan_.setSelectionSet(
          selectedMoleculeEvaluator_.evaluate());
    }

    forceManager_->setHasSelectedMolecule(hasSelectedMolecule);
  }

  bool SPFMethod::getSelectedMolecule() {
    std::shared_ptr<SPFData> spfData = currentSnap_->getSPFData();

    Molecule* selectedMolecule;

    selectedMolecule = info_->getMoleculeByGlobalIndex(spfData->globalID);

    if (selectedMolecule) {
      forceManager_->setSelectedMolecule(selectedMolecule);
    }

    return true;
  }

  bool SPFMethod::setSelectedMolecule() {
    SelectionManager sourceSman {info_};
    RealType targetSlabCenter {};

    // The sign of our flux determines which slab is the source and which is
    // the sink
    if (particleTarget_ > 0.0) {
      sourceSman       = commonA_;
      targetSlabCenter = slabBCenter_;
    } else {
      sourceSman       = commonB_;
      targetSlabCenter = slabACenter_;
    }

    if (sourceSman.getMoleculeSelectionCount() == 0) {
      forceManager_->setSelectedMolecule(nullptr);
      return false;
    }

    int whichSelectedID {-1};
    Molecule* selectedMolecule;

    std::shared_ptr<SPFData> spfData = currentSnap_->getSPFData();
    Utils::RandNumGenPtr randNumGen  = info_->getRandomNumberGenerator();

#ifdef IS_MPI
    int worldRank {}, size {};

    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status ierr;

    if (worldRank == 0) {
#endif
      std::uniform_int_distribution<> selectedMoleculeDistribution {
          0, sourceSman.getMoleculeSelectionCount() - 1};

      whichSelectedID = selectedMoleculeDistribution(*randNumGen);
#ifdef IS_MPI
    }

    MPI_Bcast(&whichSelectedID, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    selectedMolecule = sourceSman.nthSelectedMolecule(whichSelectedID);

    int globalSelectedID {-1};

    if (selectedMolecule) {
      globalSelectedID = selectedMolecule->getGlobalIndex();

#ifdef IS_MPI
      for (int i {}; i < size; ++i) {
        if (i != worldRank) {
          MPI_Send(&globalSelectedID, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
      }
    } else {
      MPI_Recv(&globalSelectedID, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
               &ierr);
#endif
    }

    int axis0 = (rnemdPrivilegedAxis_ + 1) % 3;
    int axis1 = (rnemdPrivilegedAxis_ + 2) % 3;
    int axis2 = rnemdPrivilegedAxis_;

    spfData->globalID = globalSelectedID;

#ifdef IS_MPI
    if (info_->getMolToProc(globalSelectedID) == worldRank) {
#endif
      std::uniform_real_distribution<RealType> distr0 {0, hmat_(axis0, axis0)};
      std::uniform_real_distribution<RealType> distr1 {0, hmat_(axis1, axis1)};
      std::normal_distribution<RealType> distr2 {targetSlabCenter,
                                                 0.25 * slabWidth_};

      spfData->pos[axis0] = distr0(*randNumGen);
      spfData->pos[axis1] = distr1(*randNumGen);
      spfData->pos[axis2] = distr2(*randNumGen);

      forceManager_->setSelectedMolecule(selectedMolecule);

#ifdef IS_MPI
    }
    MPI_Bcast(&spfData->pos[0], 3, MPI_REALTYPE,
              info_->getMolToProc(globalSelectedID), MPI_COMM_WORLD);
#endif

    return true;
  }
}  // namespace OpenMD::RNEMD
