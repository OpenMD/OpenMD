/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "brains/Snapshot.hpp"
#include "integrators/Integrator.hpp"
#include "integrators/DLM.hpp"
#include "flucq/FluctuatingChargeLangevin.hpp"
#include "flucq/FluctuatingChargeNVT.hpp"
#include "utils/simError.h"

namespace OpenMD {
  Integrator::Integrator(SimInfo* info) 
    : info_(info), forceMan_(NULL) , needPotential(false), needStress(false), 
      needReset(false), velocitizer_(NULL), needVelocityScaling(false), 
      rnemd_(NULL), useRNEMD(false), rotAlgo_(NULL), flucQ_(NULL), 
      rattle_(NULL), dumpWriter(NULL), statWriter(NULL), thermo(info),
      snap(info->getSnapshotManager()->getCurrentSnapshot()) {
    
    simParams = info->getSimParams();
    
    if (simParams->haveDt()) {
      dt = simParams->getDt();
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
    velocitizer_ = new Velocitizer(info);
    
    if (simParams->getRNEMDParameters()->haveUseRNEMD()) {
      if (simParams->getRNEMDParameters()->getUseRNEMD()) {
        // Create a default a RNEMD.
        rnemd_ = new RNEMD(info);
        useRNEMD = simParams->getRNEMDParameters()->getUseRNEMD();
        if (simParams->getRNEMDParameters()->haveExchangeTime()) {
          RNEMD_exchangeTime = simParams->getRNEMDParameters()->getExchangeTime();
        } 
      }
    }
    
    rotAlgo_ = new DLM();
    rattle_ = new Rattle(info);
    flucQ_ = new FluctuatingChargeLangevin(info);
    flucQ_->setForceManager(forceMan_);
  }
  
  Integrator::~Integrator(){
    delete forceMan_;
    delete velocitizer_;
    delete rnemd_;
    delete flucQ_;
    delete rotAlgo_;
    delete rattle_;
    
    delete dumpWriter;
    delete statWriter;
  }
}

