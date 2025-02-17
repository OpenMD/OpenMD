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

#include "integrators/Integrator.hpp"

#include <memory>
#include <utility>

#include "brains/Snapshot.hpp"
#include "flucq/FluctuatingChargeDamped.hpp"
#include "flucq/FluctuatingChargeLangevin.hpp"
#include "flucq/FluctuatingChargeNVE.hpp"
#include "flucq/FluctuatingChargeNVT.hpp"
#include "integrators/DLM.hpp"
#include "rnemd/MethodFactory.hpp"
#include "rnemd/RNEMD.hpp"
#include "rnemd/SPFForceManager.hpp"
#include "utils/CI_String.hpp"
#include "utils/CaseConversion.hpp"
#include "utils/simError.h"

namespace OpenMD {

  Integrator::Integrator(SimInfo* info) :
      info_(info), forceMan_(NULL), rotAlgo_(NULL), flucQ_(NULL), rattle_(NULL),
      velocitizer_(nullptr), needPotential(false), needVirial(false),
      needReset(false), needVelocityScaling(false), useRNEMD(false),
      dumpWriter(NULL), statWriter(NULL), thermo(info_),
      snap(info_->getSnapshotManager()->getCurrentSnapshot()) {
    simParams = info->getSimParams();

    if (simParams->haveDt()) {
      dt  = simParams->getDt();
      dt2 = 0.5 * dt;
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Integrator Error: dt is not set\n");
      painCave.isFatal = 1;
      simError();
    }

    if (simParams->haveRunTime()) {
      runTime = simParams->getRunTime();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Integrator Error: runTime is not set\n");
      painCave.isFatal = 1;
      simError();
    }

    // set the status, sample, and thermal kick times
    if (simParams->haveSampleTime()) {
      sampleTime = simParams->getSampleTime();
      statusTime = sampleTime;
    } else {
      sampleTime = simParams->getRunTime();
      statusTime = sampleTime;
    }

    if (simParams->haveStatusTime()) {
      statusTime = simParams->getStatusTime();
    }

    if (simParams->haveThermalTime()) {
      thermalTime = simParams->getThermalTime();
    } else {
      thermalTime = simParams->getRunTime();
    }

    if (!simParams->getUseInitalTime()) { snap->setTime(0.0); }

    if (simParams->haveResetTime()) {
      needReset = true;
      resetTime = simParams->getResetTime();
    }

    // check for the temperature set flag (velocity scaling)
    needVelocityScaling = false;
    if (simParams->haveTempSet()) {
      needVelocityScaling = simParams->getTempSet();
    }

    if (needVelocityScaling) {
      if (simParams->haveTargetTemp()) {
        targetScalingTemp = simParams->getTargetTemp();
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Integrator Error: Target Temperature must be set to turn on "
                 "tempSet\n");
        painCave.isFatal = 1;
        simError();
      }
    }

    // Create a default a velocitizer: If the subclass wants to use
    // a different velocitizer, use setVelocitizer
    velocitizer_ = std::make_unique<Velocitizer>(info);

    if (simParams->getRNEMDParameters()->haveUseRNEMD()) {
      useRNEMD = simParams->getRNEMDParameters()->getUseRNEMD();

      if (useRNEMD) {
        // ForceManager will only be changed if SPF-RNEMD is enabled
        if (Utils::traits_cast<Utils::ci_char_traits>(
                simParams->getRNEMDParameters()->getMethod()) == "SPF") {
          forceMan_ = new RNEMD::SPFForceManager(info);
        }

        RNEMD::MethodFactory rnemdMethod {
            simParams->getRNEMDParameters()->getMethod()};
        rnemd_ = rnemdMethod.create(info, forceMan_);

        if (simParams->getRNEMDParameters()->haveExchangeTime()) {
          RNEMD_exchangeTime =
              simParams->getRNEMDParameters()->getExchangeTime();

          // check to make sure exchange time is a multiple of dt;
          RealType newET = ceil(RNEMD_exchangeTime / dt) * dt;
          if (std::fabs(newET - RNEMD_exchangeTime) > 1e-6) {
            RNEMD_exchangeTime = newET;
          }
        }
      }
    }

    if (forceMan_ == NULL) forceMan_ = new ForceManager(info);

    rotAlgo_ = new DLM();
    rattle_  = new Rattle(info);

    if (simParams->getFluctuatingChargeParameters()->havePropagator()) {
      std::string prop = toUpperCopy(
          simParams->getFluctuatingChargeParameters()->getPropagator());
      if (prop.compare("NVT") == 0) {
        flucQ_ = new FluctuatingChargeNVT(info);
      } else if (prop.compare("NVE") == 0) {
        flucQ_ = new FluctuatingChargeNVE(info);
      } else if (prop.compare("LANGEVIN") == 0) {
        flucQ_ = new FluctuatingChargeLangevin(info);
      } else if (prop.compare("DAMPED") == 0) {
        flucQ_ = new FluctuatingChargeDamped(info);
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Integrator Error: Unknown Fluctuating Charge propagator (%s) "
                 "requested\n",
                 simParams->getFluctuatingChargeParameters()
                     ->getPropagator()
                     .c_str());
        painCave.isFatal = 1;
      }
    }

    flucQ_->setForceManager(forceMan_);
  }

  Integrator::~Integrator() {
    delete forceMan_;
    delete flucQ_;
    delete rotAlgo_;
    delete rattle_;
    delete dumpWriter;
    delete stats;
    delete statWriter;
  }

  void Integrator::updateSizes() {
    doUpdateSizes();
    flucQ_->updateSizes();
  }

  void Integrator::setVelocitizer(std::unique_ptr<Velocitizer> velocitizer) {
    velocitizer_ = std::move(velocitizer);
  }

  void Integrator::setFluctuatingChargePropagator(
      FluctuatingChargePropagator* prop) {
    if (prop != flucQ_ && flucQ_ != NULL) { delete flucQ_; }
    flucQ_ = prop;
    if (forceMan_ != NULL) { flucQ_->setForceManager(forceMan_); }
  }

  void Integrator::setRotationAlgorithm(RotationAlgorithm* algo) {
    if (algo != rotAlgo_ && rotAlgo_ != NULL) { delete rotAlgo_; }

    rotAlgo_ = algo;
  }

  void Integrator::setRNEMD(std::unique_ptr<RNEMD::RNEMD> rnemd) {
    rnemd_ = std::move(rnemd);
  }

  void Integrator::integrate() {
    initialize();

    while (snap->getTime() <= runTime) {
      preStep();
      step();
      postStep();
    }

    finalize();
  }

  void Integrator::saveConservedQuantity() {
    snap->setConservedQuantity(calcConservedQuantity());
  }

  void Integrator::initialize() {
    forceMan_->initialize();

    // remove center of mass drift velocity (in case we passed in a
    // configuration that was drifting)
    if (simParams->getConserveLinearMomentum()) velocitizer_->removeComDrift();

    // find the initial fluctuating charges.
    flucQ_->initialize();
    // initialize the forces before the first step
    calcForce();
    // execute the constraint algorithm to make sure that the system is
    // constrained at the very beginning
    if (info_->getNGlobalConstraints() > 0) {
      rattle_->constraintA();
      calcForce();
      rattle_->constraintB();
      // copy the current snapshot to previous snapshot
      info_->getSnapshotManager()->advance();
    }

    if (needVelocityScaling) { velocitizer_->randomize(targetScalingTemp); }

    dumpWriter = createDumpWriter();
    statWriter = createStatWriter();
    dumpWriter->writeDumpAndEor();

    progressBar = std::make_unique<ProgressBar>();

    // save statistics, before writeStat,  we must save statistics
    saveConservedQuantity();
    stats->collectStats();

    if (useRNEMD) rnemd_->getStarted();

    statWriter->writeStat();

    currSample  = sampleTime + snap->getTime();
    currStatus  = statusTime + snap->getTime();
    currThermal = thermalTime + snap->getTime();
    if (needReset) { currReset = resetTime + snap->getTime(); }
    if (useRNEMD) { currRNEMD = RNEMD_exchangeTime + snap->getTime(); }
    needPotential = false;
    needVirial    = false;
  }

  void Integrator::preStep() {
    RealType difference = snap->getTime() + dt - currStatus;

    if (difference > 0 || std::fabs(difference) <= dtEps) {
      needPotential = true;
      needVirial    = true;
    }
  }

  void Integrator::calcForce() {
    forceMan_->calcForces();
    flucQ_->applyConstraints();
  }

  void Integrator::postStep() {
    RealType difference;

    if (needVelocityScaling) {
      difference = snap->getTime() - currThermal;

      if (difference > 0 || std::fabs(difference) <= dtEps) {
        velocitizer_->randomize(targetScalingTemp);
        currThermal += thermalTime;
      }
    }

    if (useRNEMD) {
      difference = snap->getTime() - currRNEMD;

      if (difference > 0 || std::fabs(difference) <= dtEps) {
        rnemd_->doRNEMD();
        currRNEMD += RNEMD_exchangeTime;
      }

      rnemd_->collectData();
      rnemd_->writeOutputFile();
    }

    saveConservedQuantity();

    difference = snap->getTime() - currSample;

    if (difference > 0 || std::fabs(difference) <= dtEps) {
      dumpWriter->writeDumpAndEor();
      currSample += sampleTime;
    }

    difference = snap->getTime() - currStatus;

    if (difference > 0 || std::fabs(difference) <= dtEps) {
      stats->collectStats();

      statWriter->writeStat();

      progressBar->setStatus(snap->getTime(), runTime);
      progressBar->update();

      needPotential = false;
      needVirial    = false;
      currStatus += statusTime;
    }

    difference = snap->getTime() - currReset;

    if (needReset && (difference > 0 || std::fabs(difference) <= dtEps)) {
      resetIntegrator();
      currReset += resetTime;
    }
    // save snapshot
    info_->getSnapshotManager()->advance();

    // increase time
    snap->increaseTime(dt);
  }

  void Integrator::finalize() {
    dumpWriter->writeEor();
    if (useRNEMD) { rnemd_->writeOutputFile(); }
    progressBar->setStatus(runTime, runTime);
    progressBar->update();

    statWriter->writeStatReport();
  }

  DumpWriter* Integrator::createDumpWriter() { return new DumpWriter(info_); }

  StatWriter* Integrator::createStatWriter() {
    stats      = new Stats(info_);
    statWriter = new StatWriter(info_->getStatFileName(), stats);
    statWriter->setReportFileName(info_->getReportFileName());

    return statWriter;
  }
}  // namespace OpenMD
