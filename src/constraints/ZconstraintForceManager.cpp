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
 
#include <cmath>
#include "constraints/ZconstraintForceManager.hpp"
#include "integrators/Integrator.hpp"
#include "utils/simError.h"
#include "utils/PhysicalConstants.hpp"
#include "utils/StringUtils.hpp"
#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {
  ZconstraintForceManager::ZconstraintForceManager(SimInfo* info): ForceManager(info), infiniteTime(1e31) {
    currSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Globals* simParam = info_->getSimParams();

    if (simParam->haveDt()){
      dt_ = simParam->getDt();
    } else {
      sprintf(painCave.errMsg,
	      "Integrator Error: dt is not set\n");
      painCave.isFatal = 1;
      simError();    
    }

    if (simParam->haveZconsTime()){
      zconsTime_ = simParam->getZconsTime();
    }
    else{
      sprintf(painCave.errMsg,
	      "ZConstraint error: If you use a ZConstraint,\n"
	      "\tyou must set zconsTime.\n");
      painCave.isFatal = 1;
      simError();
    }

    if (simParam->haveZconsTol()){
      zconsTol_ = simParam->getZconsTol();
    }
    else{
      zconsTol_ = 0.01;
      sprintf(painCave.errMsg,
	      "ZConstraint Warning: Tolerance for z-constraint method is not specified.\n"
	      "\tOpenMD will use a default value of %f.\n"
	      "\tTo set the tolerance, use the zconsTol variable.\n",
	      zconsTol_);
      painCave.isFatal = 0;
      simError();      
    }

    //set zcons gap
    if (simParam->haveZconsGap()){
      usingZconsGap_ = true;
      zconsGap_ = simParam->getZconsGap();
    }else {
      usingZconsGap_ = false;
      zconsGap_ = 0.0;
    }

    //set zcons fixtime
    if (simParam->haveZconsFixtime()){
      zconsFixingTime_ = simParam->getZconsFixtime();
    } else {
      zconsFixingTime_ = infiniteTime;
    }

    //set zconsUsingSMD
    if (simParam->haveZconsUsingSMD()){
      usingSMD_ = simParam->getZconsUsingSMD();
    }else {
      usingSMD_ =false;
    }
    
    zconsOutput_ = getPrefix(info_->getFinalConfigFileName()) + ".fz";

    //estimate the force constant of harmonical potential
    Mat3x3d hmat = currSnapshot_->getHmat();
    RealType halfOfLargestBox = std::max(hmat(0, 0), std::max(hmat(1, 1), hmat(2, 2))) /2;	
    RealType targetTemp;
    if (simParam->haveTargetTemp()) {
      targetTemp = simParam->getTargetTemp();
    } else {
      targetTemp = 298.0;
    }
    RealType zforceConstant = PhysicalConstants::kb * targetTemp / (halfOfLargestBox * halfOfLargestBox);
         
    int nZconstraints = simParam->getNZconsStamps();
    std::vector<ZConsStamp*> stamp = simParam->getZconsStamps();
    //
    for (int i = 0; i < nZconstraints; i++){

      ZconstraintParam param;
      int zmolIndex = stamp[i]->getMolIndex();
      if (stamp[i]->haveZpos()) {
	param.zTargetPos = stamp[i]->getZpos();
      } else {
	param.zTargetPos = getZTargetPos(zmolIndex);
      }

      param.kz = zforceConstant * stamp[i]->getKratio();

      if (stamp[i]->haveCantVel()) {
	param.cantVel = stamp[i]->getCantVel();
      } else {
	param.cantVel = 0.0;
      }

      allZMolIndices_.insert(std::make_pair(zmolIndex, param));
    }

    //create fixedMols_, movingMols_ and unconsMols lists 
    update();
    
    //calculate masss of unconstraint molecules in the whole system (never change during the simulation)
    RealType totMassUnconsMols_local = 0.0;    
    std::vector<Molecule*>::iterator j;
    for ( j = unzconsMols_.begin(); j !=  unzconsMols_.end(); ++j) {
      totMassUnconsMols_local += (*j)->getMass();
    }    
#ifndef IS_MPI
    totMassUnconsMols_ = totMassUnconsMols_local;
#else
    MPI::COMM_WORLD.Allreduce(&totMassUnconsMols_local, &totMassUnconsMols_, 1,
                              MPI::REALTYPE, MPI::SUM);
#endif

    // creat zconsWriter  
    fzOut = new ZConsWriter(info_, zconsOutput_.c_str());   

    if (!fzOut){
      sprintf(painCave.errMsg, "Fail to create ZConsWriter\n");
      painCave.isFatal = 1;
      simError();
    }

  }

  ZconstraintForceManager::~ZconstraintForceManager(){

    if (fzOut){
      delete fzOut;
    }

  }

  void ZconstraintForceManager::update(){
    fixedZMols_.clear();
    movingZMols_.clear();
    unzconsMols_.clear();

    for (std::map<int, ZconstraintParam>::iterator i = allZMolIndices_.begin(); i != allZMolIndices_.end(); ++i) {
#ifdef IS_MPI
      if (info_->getMolToProc(i->first) == worldRank) {
#endif
	ZconstraintMol zmol;
	zmol.mol = info_->getMoleculeByGlobalIndex(i->first);
	assert(zmol.mol);
	zmol.param = i->second;
	zmol.cantPos = zmol.param.zTargetPos; /**@todo fixed me when zmol migrate, it is incorrect*/
	Vector3d com = zmol.mol->getCom();
	RealType diff = fabs(zmol.param.zTargetPos - com[whichDirection]);
	if (diff < zconsTol_) {
	  fixedZMols_.push_back(zmol);
	} else {
	  movingZMols_.push_back(zmol);            
	}

#ifdef IS_MPI
      }
#endif
    }

    calcTotalMassMovingZMols();

    std::set<int> zmolSet;
    for (std::list<ZconstraintMol>::iterator i = movingZMols_.begin(); i !=  movingZMols_.end(); ++i) {
      zmolSet.insert(i->mol->getGlobalIndex());
    }    

    for (std::list<ZconstraintMol>::iterator i = fixedZMols_.begin(); i !=  fixedZMols_.end(); ++i) {
      zmolSet.insert(i->mol->getGlobalIndex());
    }

    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    for(mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
      if (zmolSet.find(mol->getGlobalIndex()) == zmolSet.end()) {
	unzconsMols_.push_back(mol);
      }
    }

  }

  bool ZconstraintForceManager::isZMol(Molecule* mol){
    return allZMolIndices_.find(mol->getGlobalIndex()) == allZMolIndices_.end() ? false : true;
  }

  void ZconstraintForceManager::init(){

    //zero out the velocities of center of mass of unconstrained molecules 
    //and the velocities of center of mass of every single z-constrained molecueles
    zeroVelocity();

    currZconsTime_ = currSnapshot_->getTime();
  }

  void ZconstraintForceManager::calcForces(){
    ForceManager::calcForces();
    
    if (usingZconsGap_){
      updateZPos();
    }

    if (checkZConsState()){
      zeroVelocity();    
      calcTotalMassMovingZMols();
    }  

    //do zconstraint force; 
    if (haveFixedZMols()){
      doZconstraintForce();
    }

    //use external force to move the molecules to the specified positions
    if (haveMovingZMols()){
      doHarmonic();
    }

    //write out forces and current positions of z-constraint molecules    
    if (currSnapshot_->getTime() >= currZconsTime_){
      std::list<ZconstraintMol>::iterator i;
      Vector3d com;
      for(i = fixedZMols_.begin(); i != fixedZMols_.end(); ++i) {
	com = i->mol->getCom();
	i->zpos = com[whichDirection];
      }
        
      fzOut->writeFZ(fixedZMols_);
      currZconsTime_ += zconsTime_;
    }
  }

  void ZconstraintForceManager::zeroVelocity(){

    Vector3d comVel;
    Vector3d vel;
    std::list<ZconstraintMol>::iterator i;
    Molecule* mol;
    StuntDouble* sd;
    Molecule::IntegrableObjectIterator ii;

    //zero out the velocities of center of mass of fixed z-constrained molecules
    for(i = fixedZMols_.begin(); i != fixedZMols_.end(); ++i) {

      mol = i->mol;
      comVel = mol->getComVel();

      for(sd = mol->beginIntegrableObject(ii); sd != NULL; 
	  sd = mol->nextIntegrableObject(ii)) {

	vel = sd->getVel();  
	vel[whichDirection] -= comVel[whichDirection];
	sd->setVel(vel);
      }
    }

    // calculate the vz of center of mass of moving molecules(include unconstrained molecules 
    // and moving z-constrained molecules)  
    RealType pzMovingMols_local = 0.0;
    RealType pzMovingMols;
    
    for ( i = movingZMols_.begin(); i !=  movingZMols_.end(); ++i) {
      mol = i->mol;        
      comVel = mol->getComVel();
      pzMovingMols_local +=  mol->getMass() * comVel[whichDirection];   
    }

    std::vector<Molecule*>::iterator j;
    for ( j = unzconsMols_.begin(); j !=  unzconsMols_.end(); ++j) {
      mol =*j;
      comVel = mol->getComVel();
      pzMovingMols_local += mol->getMass() * comVel[whichDirection];
    }
    
#ifndef IS_MPI
    pzMovingMols = pzMovingMols_local;
#else
    MPI::COMM_WORLD.Allreduce(&pzMovingMols_local, &pzMovingMols, 1, 
                              MPI::REALTYPE, MPI::SUM);
#endif

    RealType vzMovingMols = pzMovingMols / (totMassMovingZMols_ + totMassUnconsMols_);

    //modify the velocities of moving z-constrained molecuels
    for ( i = movingZMols_.begin(); i !=  movingZMols_.end(); ++i) {

      mol = i->mol;

      for(sd = mol->beginIntegrableObject(ii); sd != NULL; 
	  sd = mol->nextIntegrableObject(ii)) {

	vel = sd->getVel();
	vel[whichDirection] -= vzMovingMols;
	sd->setVel(vel); 
      }
    }

    //modify the velocites of unconstrained molecules  
    for ( j = unzconsMols_.begin(); j !=  unzconsMols_.end(); ++j) {

      mol =*j;

      for(sd = mol->beginIntegrableObject(ii); sd != NULL; 
	  sd = mol->nextIntegrableObject(ii)) {

	vel = sd->getVel();
	vel[whichDirection] -= vzMovingMols;
	sd->setVel(vel); 
      }
    }
    
  }


  void ZconstraintForceManager::doZconstraintForce(){
    RealType totalFZ; 
    RealType totalFZ_local;
    Vector3d com;
    Vector3d force(0.0);

    //constrain the molecules which do not reach the specified positions  

    //Zero Out the force of z-contrained molecules    
    totalFZ_local = 0;


    //calculate the total z-contrained force of fixed z-contrained molecules
    std::list<ZconstraintMol>::iterator i;
    Molecule* mol;
    StuntDouble* sd;
    Molecule::IntegrableObjectIterator ii;

    for ( i = fixedZMols_.begin(); i !=  fixedZMols_.end(); ++i) {

      mol = i->mol;
      i->fz = 0.0;

      for( sd = mol->beginIntegrableObject(ii); sd != NULL; 
           sd = mol->nextIntegrableObject(ii)) {

	force = sd->getFrc();    
	i->fz += force[whichDirection]; 
      }
      totalFZ_local += i->fz;
    }

    //calculate total z-constraint force
#ifdef IS_MPI
    MPI::COMM_WORLD.Allreduce(&totalFZ_local, &totalFZ, 1, 
                              MPI::REALTYPE, MPI::SUM);
#else
    totalFZ = totalFZ_local;
#endif


    // apply negative to fixed z-constrained molecues;
    for ( i = fixedZMols_.begin(); i !=  fixedZMols_.end(); ++i) {

      mol = i->mol;

      for(sd = mol->beginIntegrableObject(ii); sd != NULL; 
	  sd = mol->nextIntegrableObject(ii)) {

	force[whichDirection] = -getZFOfFixedZMols(mol, sd, i->fz);
	sd->addFrc(force);
      }
    }

    //modify the forces of moving z-constrained molecules
    for ( i = movingZMols_.begin(); i !=  movingZMols_.end(); ++i) {

      mol = i->mol;

      for(sd = mol->beginIntegrableObject(ii); sd != NULL; 
	  sd = mol->nextIntegrableObject(ii)) {

	force[whichDirection] = -getZFOfMovingMols(mol,totalFZ);
	sd->addFrc(force);
      }
    }

    //modify the forces of unconstrained molecules
    std::vector<Molecule*>::iterator j;
    for ( j = unzconsMols_.begin(); j !=  unzconsMols_.end(); ++j) {

      mol =*j;

      for(sd = mol->beginIntegrableObject(ii); sd != NULL; 
	  sd = mol->nextIntegrableObject(ii)) {

	force[whichDirection] = -getZFOfMovingMols(mol, totalFZ);
	sd->addFrc(force);
      }
    }

  }


  void ZconstraintForceManager::doHarmonic(){
    RealType totalFZ;
    Vector3d force(0.0);
    Vector3d com;
    RealType totalFZ_local = 0;
    RealType lrPot;
    std::list<ZconstraintMol>::iterator i;
    StuntDouble* sd;
    Molecule::IntegrableObjectIterator ii;
    Molecule* mol;
    for ( i = movingZMols_.begin(); i !=  movingZMols_.end(); ++i) {
      mol = i->mol;
      com = mol->getCom();   
      RealType resPos = usingSMD_? i->cantPos : i->param.zTargetPos;
      RealType diff = com[whichDirection] - resPos; 
      RealType harmonicU = 0.5 * i->param.kz * diff * diff;
      lrPot = currSnapshot_->getLongRangePotential();
      lrPot += harmonicU;
      currSnapshot_->setLongRangePotential(lrPot);
      RealType harmonicF = -i->param.kz * diff;
      totalFZ_local += harmonicF;

      //adjust force
      for(sd = mol->beginIntegrableObject(ii); sd != NULL; 
	  sd = mol->nextIntegrableObject(ii)) {

	force[whichDirection] = getHFOfFixedZMols(mol, sd, harmonicF);
	sd->addFrc(force);            
      }
    }

#ifndef IS_MPI
    totalFZ = totalFZ_local;
#else
    MPI::COMM_WORLD.Allreduce(&totalFZ_local, &totalFZ, 1, MPI::REALTYPE, 
                              MPI::SUM);
#endif

    //modify the forces of unconstrained molecules
    std::vector<Molecule*>::iterator j;
    for ( j = unzconsMols_.begin(); j !=  unzconsMols_.end(); ++j) {

      mol = *j;

      for(sd = mol->beginIntegrableObject(ii); sd != NULL; 
	  sd = mol->nextIntegrableObject(ii)) {

	force[whichDirection] = getHFOfUnconsMols(mol, totalFZ);
	sd->addFrc(force);            
      }
    }

  }

  bool ZconstraintForceManager::checkZConsState(){
    Vector3d com;
    RealType diff;
    int changed_local = 0;

    std::list<ZconstraintMol>::iterator i;
    std::list<ZconstraintMol>::iterator j;
    
    std::list<ZconstraintMol> newMovingZMols;
    for ( i = fixedZMols_.begin(); i !=  fixedZMols_.end();) {
      com = i->mol->getCom();
      diff = fabs(com[whichDirection] - i->param.zTargetPos);
      if (diff > zconsTol_) {
	if (usingZconsGap_) {
	  i->endFixingTime = infiniteTime;
	}
	j = i++;
	newMovingZMols.push_back(*j);
	fixedZMols_.erase(j);
	    
	changed_local = 1;            
      }else {
	++i;
      }
    }  

    std::list<ZconstraintMol> newFixedZMols;
    for ( i = movingZMols_.begin(); i !=  movingZMols_.end();) {
      com = i->mol->getCom();
      diff = fabs(com[whichDirection] - i->param.zTargetPos);
      if (diff <= zconsTol_) {
	if (usingZconsGap_) {
	  i->endFixingTime = currSnapshot_->getTime() + zconsFixingTime_;
	}
	//this moving zconstraint molecule is about to fixed
	//moved this molecule to
	j = i++;
	newFixedZMols.push_back(*j);
	movingZMols_.erase(j);
	changed_local = 1;            
      }else {
	++i;
      }
    }     

    //merge the lists
    fixedZMols_.insert(fixedZMols_.end(), newFixedZMols.begin(), newFixedZMols.end());
    movingZMols_.insert(movingZMols_.end(), newMovingZMols.begin(), newMovingZMols.end());

    int changed;
#ifndef IS_MPI
    changed = changed_local; 
#else
    MPI::COMM_WORLD.Allreduce(&changed_local, &changed, 1, MPI::INT, MPI::SUM);
#endif

    return (changed > 0);
  }

  bool ZconstraintForceManager::haveFixedZMols(){
    int havingFixed;
    int havingFixed_local = fixedZMols_.empty() ? 0 : 1;

#ifndef IS_MPI
    havingFixed = havingFixed_local;
#else
    MPI::COMM_WORLD.Allreduce(&havingFixed_local, &havingFixed, 1, 
                              MPI::INT, MPI::SUM);
#endif

    return havingFixed > 0;
  }


  bool ZconstraintForceManager::haveMovingZMols(){
    int havingMoving_local;
    int havingMoving;

    havingMoving_local = movingZMols_.empty()? 0 : 1;

#ifndef IS_MPI
    havingMoving = havingMoving_local;
#else
    MPI::COMM_WORLD.Allreduce(&havingMoving_local, &havingMoving, 1, 
                              MPI::INT, MPI::SUM);
#endif

    return havingMoving > 0;
  }

  void ZconstraintForceManager::calcTotalMassMovingZMols(){

    RealType totMassMovingZMols_local = 0.0;
    std::list<ZconstraintMol>::iterator i;
    for ( i = movingZMols_.begin(); i !=  movingZMols_.end(); ++i) {
      totMassMovingZMols_local += i->mol->getMass();
    }
    
#ifdef IS_MPI
    MPI::COMM_WORLD.Allreduce(&totMassMovingZMols_local, &totMassMovingZMols_, 
                              1, MPI::REALTYPE, MPI::SUM);
#else
    totMassMovingZMols_ = totMassMovingZMols_local;
#endif

  }

  RealType ZconstraintForceManager::getZFOfFixedZMols(Molecule* mol, StuntDouble* sd, RealType totalForce){
    return totalForce * sd->getMass() / mol->getMass();
  }

  RealType ZconstraintForceManager::getZFOfMovingMols(Molecule* mol, RealType totalForce){
    return totalForce * mol->getMass() / (totMassUnconsMols_ + totMassMovingZMols_);
  }

  RealType ZconstraintForceManager::getHFOfFixedZMols(Molecule* mol, StuntDouble*sd, RealType totalForce){
    return totalForce * sd->getMass() / mol->getMass();
  }

  RealType ZconstraintForceManager::getHFOfUnconsMols(Molecule* mol, RealType totalForce){
    return totalForce * mol->getMass() / totMassUnconsMols_;
  }

  void ZconstraintForceManager::updateZPos(){
    RealType curTime = currSnapshot_->getTime();
    std::list<ZconstraintMol>::iterator i;
    for ( i = fixedZMols_.begin(); i !=  fixedZMols_.end(); ++i) {
      i->param.zTargetPos += zconsGap_;     
    }  
  }

  void ZconstraintForceManager::updateCantPos(){
    std::list<ZconstraintMol>::iterator i;
    for ( i = movingZMols_.begin(); i !=  movingZMols_.end(); ++i) {
      i->cantPos += i->param.cantVel * dt_;
    }
  }

  RealType ZconstraintForceManager::getZTargetPos(int index){
    RealType zTargetPos;
#ifndef IS_MPI    
    Molecule* mol = info_->getMoleculeByGlobalIndex(index);
    assert(mol);
    Vector3d com = mol->getCom();
    zTargetPos = com[whichDirection];
#else
    int whicProc = info_->getMolToProc(index);
    MPI::COMM_WORLD.Bcast(&zTargetPos, 1, MPI::REALTYPE, whicProc);
#endif
    return zTargetPos;
  }

}
