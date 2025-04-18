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
#include "types/FixedChargeAdapter.hpp"
#include "utils/Constants.hpp"
#include "utils/RandNumGen.hpp"
#include "utils/simError.h"

namespace OpenMD::RNEMD {

  SPFMethod::SPFMethod(SimInfo* info, ForceManager* forceMan) :
      RNEMD {info, forceMan}, smanA_ {info}, smanB_ {info}, anionMan_ {info},
      cationMan_ {info}, selectedMoleculeMan_ {info},
      selectedMoleculeEvaluator_ {info} {
    rnemdMethodLabel_ = "SPF";

    selectedMoleculeStr_ = "select none";
    selectedMoleculeEvaluator_.loadScriptString(selectedMoleculeStr_);
    selectedMoleculeMan_.setSelectionSet(selectedMoleculeEvaluator_.evaluate());

    if (SPFForceManager* spfForceManager =
            dynamic_cast<SPFForceManager*>(forceMan)) {
      forceManager_            = spfForceManager;
      forceManager_->spfRNEMD_ = this;
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SPF-RNEMD cannot be used with the default ForceManager.\n");
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    RNEMDParameters* rnemdParams = info->getSimParams()->getRNEMDParameters();

    // Calculate ion fixed charges for use in the Charged-SPF method
    if (useChargedSPF_) {
      SimInfo::MoleculeIterator i;
      Molecule* mol;
      std::vector<RealType> q_tot(objectTypes_.size());
      std::vector<int> molCount(objectTypes_.size());

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (std::size_t i {}; i < objectTypes_.size(); ++i) {
          if (objectTypes_[i] == mol->getMolStamp()) {
            q_tot[i] += mol->getFixedCharge();
            molCount[i]++;
          }
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &q_tot[0], q_tot.size(), MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &molCount[0], molCount.size(), MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);
#endif

      for (std::size_t i {}; i < objectTypes_.size(); ++i) {
        SelectionEvaluator ionEvaluator {info};
        SelectionManager ionManager {info};

        std::string ionStr = "select " + objectTypes_[i]->getName();
        ionEvaluator.loadScriptString(ionStr);
        ionManager.setSelectionSet(ionEvaluator.evaluate());

        if (molCount[i] > 0) {
          q_tot[i] /= molCount[i];

          if (q_tot[i] > 0.0) {
            cationMan_ |= ionManager;
          } else if (q_tot[i] < 0.0) {
            anionMan_ |= ionManager;
          }
        } else {
          q_tot[i] = 0.0;
        }
      }
    }

    bool hasParticleFlux   = rnemdParams->haveParticleFlux();
    bool hasCurrentDensity = rnemdParams->haveCurrentDensity();
    bool hasKineticFlux    = rnemdParams->haveKineticFlux();

    bool methodFluxMismatch = false;
    bool hasCorrectFlux     = false;

    switch (rnemdFluxType_) {
    case rnemdParticle:
      hasCorrectFlux = hasParticleFlux;
      break;
    case rnemdParticleKE:
      hasCorrectFlux = hasParticleFlux && hasKineticFlux;
      break;
    case rnemdCurrentDensity:
      hasCorrectFlux = hasCurrentDensity;
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
               "\tinclude: particleFlux, particleFlux + kineticFlux,\n"
               "\tand currentDensity.\n",
               rnemdFluxTypeLabel_.c_str());
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (hasParticleFlux) {
      setParticleFlux(rnemdParams->getParticleFlux());
    } else if (hasCurrentDensity) {
      setParticleFlux(rnemdParams->getCurrentDensity());
    } else {
      setParticleFlux(0.0);
    }

    if (hasKineticFlux) {
      setKineticFlux(rnemdParams->getKineticFlux());
    } else {
      setKineticFlux(0.0);
    }

    uniformKineticScaling_ = rnemdParams->getSPFUniformKineticScaling();
  }

  SPFMethod::SlabThermodynamics SPFMethod::calculateSlabTherodynamicQuantities(
      SelectionManager& sman) {
    SlabThermodynamics slab {};

    int selei {}, selej {};

    StuntDouble* sd;

    for (sd = sman.beginSelected(selei); sd != NULL;
         sd = sman.nextSelected(selei)) {
      RealType mass = sd->getMass();
      Vector3d vel  = sd->getVel();

      slab.P += mass * vel;
      slab.M += mass;
      slab.K += mass * vel.lengthSquare();

      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I       = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          slab.K +=
              angMom[j] * angMom[j] / I(j, j) + angMom[k] * angMom[k] / I(k, k);
        } else {
          slab.K += angMom[0] * angMom[0] / I(0, 0) +
                    angMom[1] * angMom[1] / I(1, 1) +
                    angMom[2] * angMom[2] / I(2, 2);
        }
      }
    }

    slab.K *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &(slab.P[0]), 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &(slab.M), 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &(slab.K), 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    return slab;
  }

  void SPFMethod::isValidExchange(Vector3d& v_a, Vector3d& v_b, RealType& a,
                                  RealType& b) {
    const auto& [P_a, M_a, K_a] = calculateSlabTherodynamicQuantities(smanA_);
    const auto& [P_b, M_b, K_b] = calculateSlabTherodynamicQuantities(smanB_);

    RealType deltaU =
        Constants::energyConvert * forceManager_->getScaledDeltaU();

    if ((M_a > 0.0) && (M_b > 0.0)) {  // both slabs are not empty
      v_a = P_a / M_a;
      v_b = P_b / M_b;

      if (uniformKineticScaling_) {
        RealType numerator   = deltaU;
        RealType denominator = K_a + K_b;
        denominator -= 0.5 * M_a * v_a.lengthSquare();
        denominator -= 0.5 * M_b * v_b.lengthSquare();

        RealType a2 = (numerator / denominator) + 1.0;

        if (a2 > 0.0) {
          a = std::sqrt(a2);
          b = a;
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
        }
      }
    }
  }

  void SPFMethod::doRNEMDImpl(SelectionManager& smanA,
                              SelectionManager& smanB) {
    if (!doRNEMD_) return;

    if (!forceManager_->getHasSelectedMolecule()) { selectMolecule(); }

    // Remove selected molecule from the source selection manager
    if (spfTarget_ > 0.0) {
      smanA -= selectedMoleculeMan_;
    } else {
      smanB -= selectedMoleculeMan_;
    }

    smanA_ = smanA;
    smanB_ = smanB;

    if (!failedLastTrial_) {
      Vector3d v_a {};
      Vector3d v_b {};
      RealType a {0.0};
      RealType b {0.0};

      isValidExchange(v_a, v_b, a, b);

      Vector3d vel;

      int selei {}, selej {};
      StuntDouble* sd;

      for (sd = smanA.beginSelected(selei); sd != NULL;
           sd = smanA.nextSelected(selei)) {
        vel = (sd->getVel() - v_a) * a + v_a;
        sd->setVel(vel);

        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ() * a;
          sd->setJ(angMom);
        }
      }

      for (sd = smanB.beginSelected(selej); sd != NULL;
           sd = smanB.nextSelected(selej)) {
        vel = (sd->getVel() - v_b) * b + v_b;
        sd->setVel(vel);

        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ() * b;
          sd->setJ(angMom);
        }
      }

      if (useChargedSPF_) {
        RealType fixedChargeOnIon {};

        if (Molecule* selectedMolecule = forceManager_->getSelectedMolecule();
            selectedMolecule) {
          fixedChargeOnIon = selectedMolecule->getFixedCharge();
        }

#ifdef IS_MPI
        MPI_Bcast(&fixedChargeOnIon, 1, MPI_REALTYPE,
                  info_->getMolToProc(currentSnap_->getSPFData()->globalID),
                  MPI_COMM_WORLD);
#endif

        particleExchange_ += spfTarget_ * fixedChargeOnIon;
      } else {
        particleExchange_ += spfTarget_;
      }

      kineticExchange_ += kineticTarget_;
    } else {
      failTrialCount_++;
    }

    currentSnap_->hasTranslationalKineticEnergy = false;
    currentSnap_->hasRotationalKineticEnergy    = false;
    currentSnap_->hasKineticEnergy              = false;
  }

  void SPFMethod::selectMolecule() {
    std::shared_ptr<SPFData> spfData = currentSnap_->getSPFData();

    // Always reset spfTarget_ before selecting a new particle
    spfTarget_ = particleTarget_;

    bool hasSelectedMolecule = (spfData->globalID == -1) ?
                                   setSelectedMolecule(spfData) :
                                   getSelectedMolecule(spfData);

    if (hasSelectedMolecule) {
      selectedMoleculeStr_ = "select " + std::to_string(spfData->globalID);

      selectedMoleculeEvaluator_.loadScriptString(selectedMoleculeStr_);
      selectedMoleculeMan_.setSelectionSet(
          selectedMoleculeEvaluator_.evaluate());

#ifdef IS_MPI
      MPI_Bcast(&spfTarget_, 1, MPI_REALTYPE,
                info_->getMolToProc(spfData->globalID), MPI_COMM_WORLD);
#endif
    } else {
      failedLastTrial_ = true;
    }

    forceManager_->setHasSelectedMolecule(hasSelectedMolecule);
  }

  bool SPFMethod::getSelectedMolecule(std::shared_ptr<SPFData> spfData) {
    Molecule* selectedMolecule;

    selectedMolecule = info_->getMoleculeByGlobalIndex(spfData->globalID);

    if (selectedMolecule) {
      if (useChargedSPF_) {
        if (selectedMolecule->getFixedCharge() < 0.0) { spfTarget_ *= -1; }

        // Scale the particle flux by the charge yielding a current density
        convertParticlesToElectrons(selectedMolecule);
      }

      forceManager_->setSelectedMolecule(selectedMolecule);
    }

    return true;
  }

  bool SPFMethod::setSelectedMolecule(std::shared_ptr<SPFData> spfData) {
    SelectionManager sourceSman {info_}, oppositeIonSman {info_};
    RealType targetSlabCenter {};

    Utils::RandNumGenPtr randNumGen = info_->getRandomNumberGenerator();

    int ion {-1};

#ifdef IS_MPI
    int worldRank {}, size {};

    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (worldRank == 0) {
#endif
      if (useChargedSPF_) {
        std::uniform_int_distribution<> slabDistribution {0, 1};
        ion = slabDistribution(*randNumGen);
      }
#ifdef IS_MPI
    }

    if (useChargedSPF_) { MPI_Bcast(&ion, 1, MPI_INT, 0, MPI_COMM_WORLD); }
#endif

    SelectedIon selectedIon = static_cast<SelectedIon>(ion);

    if (selectedIon == ANION) {
      spfTarget_ *= -1;
      oppositeIonSman = cationMan_;
    } else if (selectedIon == CATION) {
      oppositeIonSman = anionMan_;
    }

    // The sign of our flux determines which slab is the source and which is
    // the sink
    if (spfTarget_ > 0.0) {
      sourceSman       = commonA_;
      targetSlabCenter = slabBCenter_;
    } else {
      sourceSman       = commonB_;
      targetSlabCenter = slabACenter_;
    }

    if (useChargedSPF_) { sourceSman -= oppositeIonSman; }

    if (sourceSman.getMoleculeSelectionCount() == 0) {
      forceManager_->setSelectedMolecule(nullptr);
      return false;
    }

    // Choose a molecule to move from the designated source slab
    int whichSelectedID {-1};
    Molecule* selectedMolecule;

#ifdef IS_MPI
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

      // Scale the particle flux by the charge yielding a current density
      if (useChargedSPF_) { convertParticlesToElectrons(selectedMolecule); }
    }

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &globalSelectedID, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
#endif

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
