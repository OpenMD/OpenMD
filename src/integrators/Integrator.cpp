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

#include <memory>
#include <utility>
 
#include "brains/Snapshot.hpp"
#include "integrators/Integrator.hpp"
#include "integrators/DLM.hpp"
#include "flucq/FluctuatingChargeDamped.hpp"
#include "flucq/FluctuatingChargeLangevin.hpp"
#include "flucq/FluctuatingChargeNVE.hpp"
#include "flucq/FluctuatingChargeNVT.hpp"
#include "utils/simError.h"
#include "utils/MemoryUtils.hpp"

namespace OpenMD {
  Integrator::Integrator(SimInfo* info) 
    : info_(info), forceMan_(NULL), rotAlgo_(NULL), flucQ_(NULL), 
      rattle_(NULL), velocitizer_(nullptr), rnemd_(NULL), 
      needPotential(false), needVirial(false), 
      needReset(false),  needVelocityScaling(false), 
      useRNEMD(false), dumpWriter(NULL), statWriter(NULL), thermo(info_),
      snap(info_->getSnapshotManager()->getCurrentSnapshot()) {
    
    simParams = info->getSimParams();
    
    if (simParams->haveDt()) {
      dt = simParams->getDt();
      dt2 = 0.5 * dt;
    } else {
      sprintf(painCave.errMsg,
              "Integrator Error: dt is not set\n");
      painCave.isFatal = 1;
      simError();
    }
    
    if (simParams->haveRunTime()) {
      runTime = simParams->getRunTime();
    } else {
      sprintf(painCave.errMsg,
              "Integrator Error: runTime is not set\n");
      painCave.isFatal = 1;
      simError();
    }
    
    // set the status, sample, and thermal kick times
    if (simParams->haveSampleTime()){
      sampleTime = simParams->getSampleTime();
      statusTime = sampleTime;
    } else{
      sampleTime = simParams->getRunTime();
      statusTime = sampleTime;
    }
    
    if (simParams->haveStatusTime()){
      statusTime = simParams->getStatusTime();
    }
    
    if (simParams->haveThermalTime()){
      thermalTime = simParams->getThermalTime();
    } else {
      thermalTime = simParams->getRunTime();
    }
    
    if (!simParams->getUseInitalTime()) {
      snap->setTime(0.0);
    }
    
    if (simParams->haveResetTime()) {
      needReset = true;
      resetTime = simParams->getResetTime();
    }
    
    // Create a default ForceManager: If the subclass wants to use 
    // a different ForceManager, use setForceManager
    forceMan_ = new ForceManager(info);
    
    // check for the temperature set flag (velocity scaling)      
    needVelocityScaling = false;
    if (simParams->haveTempSet()) {
      needVelocityScaling = simParams->getTempSet();
    } 

    if (needVelocityScaling) {
      if (simParams->haveTargetTemp()) {
        targetScalingTemp = simParams->getTargetTemp();
      }
      else {
        sprintf(painCave.errMsg,
                "Integrator Error: Target Temperature must be set to turn on tempSet\n");
        painCave.isFatal = 1;
        simError();

      }
    }
    
    // Create a default a velocitizer: If the subclass wants to use 
    // a different velocitizer, use setVelocitizer
    // Remove in favor of std::MemoryUtils::make_unique<> when we switch to C++14 and above
    velocitizer_ = MemoryUtils::make_unique<Velocitizer>(info);
    
    if (simParams->getRNEMDParameters()->haveUseRNEMD()) {
      if (simParams->getRNEMDParameters()->getUseRNEMD()) {
        // Create a default a RNEMD.
        rnemd_ = new RNEMD(info);
        useRNEMD = simParams->getRNEMDParameters()->getUseRNEMD();
        if (simParams->getRNEMDParameters()->haveExchangeTime()) {
          RNEMD_exchangeTime = simParams->getRNEMDParameters()->getExchangeTime();
	  // check to make sure exchange time is a multiple of dt;
	  RealType newET = ceil(RNEMD_exchangeTime / dt) * dt;
	  if (fabs( newET - RNEMD_exchangeTime ) > 1e-6) {
	    RNEMD_exchangeTime = newET;
	  }
        } 
      }
    }
    
    rotAlgo_ = new DLM();
    rattle_ = new Rattle(info);
    
    if (simParams->getFluctuatingChargeParameters()->havePropagator()) {
      std::string prop = toUpperCopy(simParams->getFluctuatingChargeParameters()->getPropagator());
      if (prop.compare("NVT")==0){
        flucQ_ = new FluctuatingChargeNVT(info);
      } else if (prop.compare("NVE")==0) {
        flucQ_ = new FluctuatingChargeNVE(info);
      } else if (prop.compare("LANGEVIN")==0) {
        flucQ_ = new FluctuatingChargeLangevin(info);
      } else if (prop.compare("DAMPED")==0){
        flucQ_ = new FluctuatingChargeDamped(info);         
      } else {
        sprintf(painCave.errMsg,
                "Integrator Error: Unknown Fluctuating Charge propagator (%s) requested\n",
                simParams->getFluctuatingChargeParameters()->getPropagator().c_str());
        painCave.isFatal = 1;
      }
    }
    flucQ_->setForceManager(forceMan_);
  }
  
  Integrator::~Integrator(){
    delete forceMan_;
    delete rnemd_;
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

  void Integrator::setForceManager(ForceManager* forceMan) {

    if (forceMan_ != forceMan && forceMan_  != NULL) {
      delete forceMan_;
    }
    forceMan_ = forceMan;
    // forward this on:
    if (flucQ_ != NULL) {
      flucQ_->setForceManager(forceMan_);
    }
  }

  void Integrator::setVelocitizer(VelocitizerPtr velocitizer) {
    velocitizer_ = std::move(velocitizer);
  }

  void Integrator::setFluctuatingChargePropagator(FluctuatingChargePropagator* prop) {
    if (prop != flucQ_ && flucQ_ != NULL){            
      delete flucQ_;
    }            
    flucQ_ = prop;
    if (forceMan_ != NULL) {
      flucQ_->setForceManager(forceMan_);
    }
  }

  void Integrator::setRotationAlgorithm(RotationAlgorithm* algo) {
    if (algo != rotAlgo_ && rotAlgo_ != NULL){            
      delete rotAlgo_;
    }
            
    rotAlgo_ = algo;
  }

  void Integrator::setRNEMD(RNEMD* rnemd) {
    if (rnemd_ != rnemd && rnemd_  != NULL) {
      delete rnemd_;
    }
    rnemd_ = rnemd;
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
    snap->setConservedQuantity( calcConservedQuantity() );
  }
  
  void Integrator::initialize(){
    
    forceMan_->initialize();

    // remove center of mass drift velocity (in case we passed in a
    // configuration that was drifting)
    velocitizer_->removeComDrift();

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
      //copy the current snapshot to previous snapshot
      info_->getSnapshotManager()->advance();
    }
    
    if (needVelocityScaling) {      
      velocitizer_->randomize(targetScalingTemp);
    }

    dumpWriter = createDumpWriter();    
    statWriter = createStatWriter(); 
    dumpWriter->writeDumpAndEor();

    // Remove in favor of std::MemoryUtils::make_unique<> when we switch to C++14 and above
    progressBar = MemoryUtils::make_unique<ProgressBar>();

    //save statistics, before writeStat,  we must save statistics
    saveConservedQuantity();
    stats->collectStats();

    if (simParams->getRNEMDParameters()->getUseRNEMD())
      rnemd_->getStarted();

    statWriter->writeStat();
    
    currSample = sampleTime + snap->getTime();
    currStatus =  statusTime + snap->getTime();
    currThermal = thermalTime + snap->getTime();
    if (needReset) {
      currReset = resetTime + snap->getTime();
    }
    if (simParams->getRNEMDParameters()->getUseRNEMD()){
      currRNEMD = RNEMD_exchangeTime + snap->getTime();
    }
    needPotential = false;
    needVirial = false;       
  }

  void Integrator::preStep() {
    
    RealType difference = snap->getTime() + dt - currStatus;
  
    if (difference > 0 || fabs(difference) <= OpenMD::epsilon) {
      needPotential = true;
      needVirial = true;   
    }
  }

  void Integrator::calcForce() { 
    forceMan_->calcForces();
    flucQ_->applyConstraints();
  }
  
  void Integrator::postStep() {

    RealType difference;

    saveConservedQuantity();

    if (needVelocityScaling) {
      difference = snap->getTime() - currThermal;

      if (difference > 0 || fabs(difference) <= OpenMD::epsilon) {
	velocitizer_->randomize(targetScalingTemp);
	currThermal += thermalTime;
      }
    }
    
    if (useRNEMD) {
      difference = snap->getTime() - currRNEMD;

      if (difference > 0 || fabs(difference) <= OpenMD::epsilon) {
	rnemd_->doRNEMD();
	currRNEMD += RNEMD_exchangeTime;
      }
      rnemd_->collectData();
    }

    difference = snap->getTime() - currSample;
  
    if (difference > 0 || fabs(difference) <= OpenMD::epsilon) {
      dumpWriter->writeDumpAndEor();      
      currSample += sampleTime;
    }

    difference = snap->getTime() - currStatus;

    if (difference > 0 || fabs(difference) <= OpenMD::epsilon) {
      stats->collectStats();

      if (simParams->getRNEMDParameters()->getUseRNEMD()) {
	rnemd_->writeOutputFile();
      }

      statWriter->writeStat();

      progressBar->setStatus(snap->getTime(), runTime);
      progressBar->update();

      needPotential = false;
      needVirial = false;
      currStatus += statusTime;
    }

    difference = snap->getTime() - currReset;

    if (needReset &&
        (difference > 0 || fabs(difference) <= OpenMD::epsilon)) {    
      resetIntegrator();
      currReset += resetTime;
    }        
    //save snapshot
    info_->getSnapshotManager()->advance();
  
    //increase time
    snap->increaseTime(dt);        

  }

  void Integrator::finalize() {
    dumpWriter->writeEor();
    if (simParams->getRNEMDParameters()->getUseRNEMD()) {
      rnemd_->writeOutputFile();
    }
    progressBar->setStatus(runTime, runTime);
    progressBar->update();

    statWriter->writeStatReport();
 
    // delete dumpWriter;
    // delete statWriter;
  
    // dumpWriter = NULL;
    // statWriter = NULL;
  }

  DumpWriter* Integrator::createDumpWriter() {
    return new DumpWriter(info_);
  }
  
  StatWriter* Integrator::createStatWriter() {
    
    stats = new Stats(info_);
    statWriter = new StatWriter(info_->getStatFileName(), stats);
    statWriter->setReportFileName(info_->getReportFileName());
    
    return statWriter;
  }
}
