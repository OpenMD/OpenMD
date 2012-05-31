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
 
/**
 * @file VelocityVerletIntegrator.cpp
 * @author tlin
 * @date 11/09/2004
 * @time 16:16am
 * @version 1.0
 */

#include "integrators/VelocityVerletIntegrator.hpp"
#include "integrators/DLM.hpp"
#include "utils/StringUtils.hpp"
#include "utils/ProgressBar.hpp"

namespace OpenMD {
  VelocityVerletIntegrator::VelocityVerletIntegrator(SimInfo *info) : Integrator(info) { 
    dt2 = 0.5 * dt;
  }
  
  VelocityVerletIntegrator::~VelocityVerletIntegrator() { 
  }
  
  void VelocityVerletIntegrator::initialize(){
    
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
      velocitizer_->velocitize(targetScalingTemp);
    }
    
    dumpWriter = createDumpWriter();
    
    statWriter = createStatWriter();
 
    dumpWriter->writeDumpAndEor();

    progressBar = new ProgressBar();

    //save statistics, before writeStat,  we must save statistics
    thermo.saveStat();
    saveConservedQuantity();
    if (simParams->getRNEMDParameters()->getUseRNEMD())
      rnemd_->getStarted();

    statWriter->writeStat(currentSnapshot_->statData);
    
    currSample = sampleTime + currentSnapshot_->getTime();
    currStatus =  statusTime + currentSnapshot_->getTime();
    currThermal = thermalTime + currentSnapshot_->getTime();
    if (needReset) {
      currReset = resetTime + currentSnapshot_->getTime();
    }
    if (simParams->getRNEMDParameters()->getUseRNEMD()){
      currRNEMD = RNEMD_exchangeTime + currentSnapshot_->getTime();
    }
    needPotential = false;
    needStress = false;       
    
  }
 
  void VelocityVerletIntegrator::doIntegrate() {
  
    initialize();
  
    while (currentSnapshot_->getTime() < runTime) {    
      preStep();    
      integrateStep();    
      postStep();      
    }  
    finalize();
  }


  void VelocityVerletIntegrator::preStep() {
    RealType difference = currentSnapshot_->getTime() + dt - currStatus;
  
    if (difference > 0 || fabs(difference) < OpenMD::epsilon) {
      needPotential = true;
      needStress = true;   
    }
  }

  void VelocityVerletIntegrator::postStep() {

    //save snapshot
    info_->getSnapshotManager()->advance();
  
    //increase time
    currentSnapshot_->increaseTime(dt);        
   
    if (needVelocityScaling) {
      if (currentSnapshot_->getTime() >= currThermal) {
	velocitizer_->velocitize(targetScalingTemp);
	currThermal += thermalTime;
      }
    }
    if (useRNEMD) {
      if (currentSnapshot_->getTime() >= currRNEMD) {
	rnemd_->doRNEMD();
	currRNEMD += RNEMD_exchangeTime;
      }
      rnemd_->collectData();
    }
    
    if (currentSnapshot_->getTime() >= currSample) {
      dumpWriter->writeDumpAndEor();
      
      currSample += sampleTime;
    }
    
    if (currentSnapshot_->getTime() >= currStatus) {
      //save statistics, before writeStat,  we must save statistics
      thermo.saveStat();
      saveConservedQuantity();

      if (simParams->getRNEMDParameters()->getUseRNEMD()) {
	rnemd_->getStatus();
      }

      statWriter->writeStat(currentSnapshot_->statData);

      progressBar->setStatus(currentSnapshot_->getTime(), runTime);
      progressBar->update();

      needPotential = false;
      needStress = false;
      currStatus += statusTime;
    }
    
    if (needReset && currentSnapshot_->getTime() >= currReset) {    
      resetIntegrator();
      currReset += resetTime;
    }        
  }


  void VelocityVerletIntegrator::finalize() {
    dumpWriter->writeEor();
  
    delete dumpWriter;
    delete statWriter;
  
    dumpWriter = NULL;
    statWriter = NULL;
  }

  void VelocityVerletIntegrator::integrateStep() { 
  
    moveA();
    calcForce();
    moveB();
  }


  void VelocityVerletIntegrator::calcForce() { 
    forceMan_->calcForces();
    flucQ_->applyConstraints();
  }

  DumpWriter* VelocityVerletIntegrator::createDumpWriter() {
    return new DumpWriter(info_);
  }

  StatWriter* VelocityVerletIntegrator::createStatWriter() {

    std::string statFileFormatString = simParams->getStatFileFormat();
    StatsBitSet mask = parseStatFileFormat(statFileFormatString);
   
    // if we're doing a thermodynamic integration, we'll want the raw
    // potential as well as the full potential:

 
    if (simParams->getUseThermodynamicIntegration()) 
      mask.set(Stats::VRAW);

    // if we've got restraints turned on, we'll also want a report of the
    // total harmonic restraints
    if (simParams->getUseRestraints()){
      mask.set(Stats::VHARM);
    }

    if (simParams->havePrintPressureTensor() && 
	simParams->getPrintPressureTensor()){
      mask.set(Stats::PRESSURE_TENSOR_XX);
      mask.set(Stats::PRESSURE_TENSOR_XY);
      mask.set(Stats::PRESSURE_TENSOR_XZ);
      mask.set(Stats::PRESSURE_TENSOR_YX);
      mask.set(Stats::PRESSURE_TENSOR_YY);
      mask.set(Stats::PRESSURE_TENSOR_YZ);
      mask.set(Stats::PRESSURE_TENSOR_ZX);
      mask.set(Stats::PRESSURE_TENSOR_ZY);
      mask.set(Stats::PRESSURE_TENSOR_ZZ);
    }
    
    if (simParams->getAccumulateBoxDipole()) {
      mask.set(Stats::BOX_DIPOLE_X);
      mask.set(Stats::BOX_DIPOLE_Y);
      mask.set(Stats::BOX_DIPOLE_Z);
    }

    if (simParams->getPrintHeatFlux()) {
      mask.set(Stats::HEATFLUX_X);
      mask.set(Stats::HEATFLUX_Y);
      mask.set(Stats::HEATFLUX_Z);
    }
   
    if (simParams->haveTaggedAtomPair() && simParams->havePrintTaggedPairDistance()) {
      if (simParams->getPrintTaggedPairDistance()) {
        mask.set(Stats::TAGGED_PAIR_DISTANCE);
      }
    }

    if (simParams->getRNEMDParameters()->getUseRNEMD()) {
      mask.set(Stats::RNEMD_EXCHANGE_TOTAL);
    }
    

    return new StatWriter(info_->getStatFileName(), mask);
  }


} //end namespace OpenMD
