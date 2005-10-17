/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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

namespace oopse {
  VelocityVerletIntegrator::VelocityVerletIntegrator(SimInfo *info) : Integrator(info), rotAlgo(NULL) { 
    dt2 = 0.5 * dt;
    rotAlgo = new DLM();
    rattle = new Rattle(info);
  }
  
  VelocityVerletIntegrator::~VelocityVerletIntegrator() { 
    delete rotAlgo;
    delete rattle;
  }
  
  void VelocityVerletIntegrator::initialize(){
    
    forceMan_->init();
    
    // remove center of mass drift velocity (in case we passed in a configuration
    // that was drifting
    velocitizer_->removeComDrift();
    
    // initialize the forces before the first step
    calcForce(true, true);
    
    //execute constraint algorithm to make sure at the very beginning the system is constrained  
    if (info_->getNGlobalConstraints() > 0) {
      rattle->constraintA();
      calcForce(true, true);
      rattle->constraintB();        
      info_->getSnapshotManager()->advance();//copy the current snapshot to previous snapshot
    }
    
    if (needVelocityScaling) {
      velocitizer_->velocitize(targetScalingTemp);
    }
    
    dumpWriter = createDumpWriter();
    
    statWriter = createStatWriter();
    
    if (simParams->getUseSolidThermInt()) {
      restWriter = createRestWriter();
      restWriter->writeZangle();
    }
    
    dumpWriter->writeDumpAndEor();
    
    
    //save statistics, before writeStat,  we must save statistics
    thermo.saveStat();
    saveConservedQuantity();
    statWriter->writeStat(currentSnapshot_->statData);
    
    currSample = sampleTime + currentSnapshot_->getTime();
    currStatus =  statusTime + currentSnapshot_->getTime();;
    currThermal = thermalTime + currentSnapshot_->getTime();
    if (needReset) {
      currReset = resetTime + currentSnapshot_->getTime();
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
    double difference = currentSnapshot_->getTime() + dt - currStatus;
  
    if (difference > 0 || fabs(difference) < oopse::epsilon) {
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
  
    if (currentSnapshot_->getTime() >= currSample) {
      dumpWriter->writeDumpAndEor();
    
      if (simParams->getUseSolidThermInt())
	restWriter->writeZangle();
    
      currSample += sampleTime;
    }
  
    if (currentSnapshot_->getTime() >= currStatus) {
      //save statistics, before writeStat,  we must save statistics
      thermo.saveStat();
      saveConservedQuantity();
      statWriter->writeStat(currentSnapshot_->statData);
    
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
  
    if (simParams->getUseSolidThermInt()) {
      restWriter->writeZangle();
      delete restWriter;
      restWriter = NULL;
    }
  
    delete dumpWriter;
    delete statWriter;
  
    dumpWriter = NULL;
    statWriter = NULL;
  
  }

  void VelocityVerletIntegrator::integrateStep() { 
  
    moveA();
    calcForce(needPotential, needStress);
    moveB();
  }


  void VelocityVerletIntegrator::calcForce(bool needPotential,
					   bool needStress) { 
    forceMan_->calcForces(needPotential, needStress);
  }

  DumpWriter* VelocityVerletIntegrator::createDumpWriter() {
    return new DumpWriter(info_);
  }

  StatWriter* VelocityVerletIntegrator::createStatWriter() {

    std::string statFileFormatString = simParams->getStatFileFormat();
    StatsBitSet mask = parseStatFileFormat(statFileFormatString);
    
    // if solidThermInt is true, add extra information to the statfile
    if (simParams->getUseSolidThermInt()){
      mask.set(Stats::VRAW);
      mask.set(Stats::VHARM);
    }

    if (simParams->havePrintPressureTensor() && simParams->getPrintPressureTensor()){
        mask.set(Stats::PRESSURE_TENSOR_X);
        mask.set(Stats::PRESSURE_TENSOR_Y);
        mask.set(Stats::PRESSURE_TENSOR_Z);
    }
    
     return new StatWriter(info_->getStatFileName(), mask);
  }

  RestWriter* VelocityVerletIntegrator::createRestWriter(){
    return new RestWriter(info_);
  }


} //end namespace oopse
