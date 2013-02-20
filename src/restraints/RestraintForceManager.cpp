/*
 * Copyright (c) 2009 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "config.h"
#include <cmath>

#include "restraints/RestraintForceManager.hpp"
#include "restraints/MolecularRestraint.hpp"
#include "restraints/ObjectRestraint.hpp"
#include "io/RestReader.hpp"
#include "utils/simError.h"
#include "utils/PhysicalConstants.hpp"
#include "utils/StringUtils.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#ifdef IS_MPI
#include <mpi.h>
#endif


namespace OpenMD {

  RestraintForceManager::RestraintForceManager(SimInfo* info): ForceManager(info) {

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

    currSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Globals* simParam = info_->getSimParams();

    if (simParam->haveStatusTime()){      
      restTime_ = simParam->getStatusTime();
    } else {
      sprintf(painCave.errMsg,
              "Restraint warning: If you use restraints without setting\n"
              "\tstatusTime, no restraint data will be written to the rest\n"
              "\tfile.\n");
      painCave.isFatal = 0;
      simError();
      restTime_ = simParam->getRunTime();
    }

    int nRestraintStamps = simParam->getNRestraintStamps();
    std::vector<RestraintStamp*> stamp = simParam->getRestraintStamps();

    std::vector<int> stuntDoubleIndex;

    for (int i = 0; i < nRestraintStamps; i++){

      std::string myType = toUpperCopy(stamp[i]->getType());

      if (myType.compare("MOLECULAR")==0){

        int molIndex;
        std::vector<Vector3d> ref;
        Vector3d refCom;

        if (!stamp[i]->haveMolIndex()) {
          sprintf(painCave.errMsg,
                  "Restraint Error: A molecular restraint was specified\n"
                  "\twithout providing a value for molIndex.\n");
          painCave.isFatal = 1;
          simError();      
        } else {
          molIndex = stamp[i]->getMolIndex();
        }
        
        if (molIndex < 0) {
          sprintf(painCave.errMsg,
                  "Restraint Error: A molecular restraint was specified\n"
                  "\twith a molIndex that was less than 0\n");
          painCave.isFatal = 1;
          simError();      
        }
        if (molIndex >= info_->getNGlobalMolecules()) {
          sprintf(painCave.errMsg,
                  "Restraint Error: A molecular restraint was specified with\n"
                  "\ta molIndex that was greater than the total number of molecules\n");
          painCave.isFatal = 1;
          simError();      
        }
       
        Molecule* mol = info_->getMoleculeByGlobalIndex(molIndex);
        
        if (mol == NULL) {
#ifdef IS_MPI
          // getMoleculeByGlobalIndex returns a NULL in parallel if
          // this proc doesn't have the molecule.  Do a quick check to
          // make sure another processor is supposed to have it.
         
          int myrank = MPI::COMM_WORLD.Get_rank();
          if (info_->getMolToProc(molIndex) == myrank) {
         
            // If we were supposed to have it but got a null, then freak out.
#endif

            sprintf(painCave.errMsg,
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
        int myrank = MPI::COMM_WORLD.Get_rank();
        if (info_->getMolToProc(molIndex) == myrank) {

#endif

        MolecularRestraint* rest = new MolecularRestraint();

        std::string restPre("mol_");
        std::stringstream restName;
        restName << restPre << molIndex;
        rest->setRestraintName(restName.str());
        
        if (stamp[i]->haveDisplacementSpringConstant()) {
          rest->setDisplacementForceConstant(stamp[i]->getDisplacementSpringConstant());
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
        if (stamp[i]->haveRestrainedTwistAngle()) {
          rest->setRestrainedTwistAngle(stamp[i]->getRestrainedTwistAngle() * M_PI/180.0);
        }
        if (stamp[i]->haveRestrainedSwingYAngle()) {
          rest->setRestrainedSwingYAngle(stamp[i]->getRestrainedSwingYAngle() * M_PI/180.0);
        }
        if (stamp[i]->haveRestrainedSwingXAngle()) {
          rest->setRestrainedSwingXAngle(stamp[i]->getRestrainedSwingXAngle() * M_PI/180.0);
        }
        if (stamp[i]->havePrint()) {
          rest->setPrintRestraint(stamp[i]->getPrint());
        }

        restraints_.push_back(rest);
        mol->addProperty(new RestraintData("Restraint", rest));
        restrainedMols_.push_back(mol);
#ifdef IS_MPI
        }
#endif        
      } else if (myType.compare("OBJECT") == 0) {
        
        std::string objectSelection;

        if (!stamp[i]->haveObjectSelection()) {
          sprintf(painCave.errMsg,
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
        
        sprintf(painCave.errMsg,
                "Restraint Info: The specified restraint objectSelection,\n"
                "\t\t%s\n"
                "\twill result in %d integrable objects being\n"
                "\trestrained.\n", objectSelection.c_str(), selectionCount);
        painCave.isFatal = 0;
        simError();                     

        int selei;
        StuntDouble* sd;

        for (sd = seleMan.beginSelected(selei); sd != NULL; 
             sd = seleMan.nextSelected(selei)) {
          stuntDoubleIndex.push_back(sd->getGlobalIntegrableObjectIndex());

          ObjectRestraint* rest = new ObjectRestraint();
          
          if (stamp[i]->haveDisplacementSpringConstant()) {
            rest->setDisplacementForceConstant(stamp[i]->getDisplacementSpringConstant());
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
          if (stamp[i]->haveRestrainedTwistAngle()) {
            rest->setRestrainedTwistAngle(stamp[i]->getRestrainedTwistAngle());
          }
          if (stamp[i]->haveRestrainedSwingXAngle()) {
            rest->setRestrainedSwingXAngle(stamp[i]->getRestrainedSwingXAngle());
          }
          if (stamp[i]->haveRestrainedSwingYAngle()) {
            rest->setRestrainedSwingYAngle(stamp[i]->getRestrainedSwingYAngle());
          }     
          if (stamp[i]->havePrint()) {
            rest->setPrintRestraint(stamp[i]->getPrint());
          }
    
          restraints_.push_back(rest);
          sd->addProperty(new RestraintData("Restraint", rest));
          restrainedObjs_.push_back(sd);                    
        }

      }
    }

    // ThermodynamicIntegration subclasses RestraintForceManager, and there
    // are times when it won't use restraints at all, so only open the
    // restraint file if we are actually using restraints:

    if (simParam->getUseRestraints()) { 
      std::string refFile = simParam->getRestraint_file();
      RestReader* rr = new RestReader(info, refFile, stuntDoubleIndex);
      rr->readReferenceStructure();
    }

    restOutput_ = getPrefix(info_->getFinalConfigFileName()) + ".rest";
    restOut = new RestWriter(info_, restOutput_.c_str(), restraints_);
    if(!restOut){
      sprintf(painCave.errMsg, "Restraint error: Failed to create RestWriter\n");
      painCave.isFatal = 1;
      simError();
    } 

    // todo: figure out the scale factor.  Right now, just scale them all to 1
    std::vector<Restraint*>::const_iterator resti;
    for(resti=restraints_.begin(); resti != restraints_.end(); ++resti){
      (*resti)->setScaleFactor(1.0);
    }
  }

  RestraintForceManager::~RestraintForceManager(){
    if (restOut)
      delete restOut;
  }

  void RestraintForceManager::init() {
    currRestTime_ = currSnapshot_->getTime();
  }

  void RestraintForceManager::calcForces(){

    ForceManager::calcForces();    
    RealType restPot_local, restPot;

    restPot_local = doRestraints(1.0);

#ifdef IS_MPI    
    MPI::COMM_WORLD.Allreduce(&restPot_local, &restPot, 1, 
                              MPI::REALTYPE, MPI::SUM);
#else
    restPot = restPot_local;
#endif
    currSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    RealType pot = currSnapshot_->getLongRangePotential();
    pot += restPot;
    currSnapshot_->setLongRangePotential(pot);
    currSnapshot_->setRestraintPotential(restPot);

    //write out forces and current positions of restrained molecules    
    if (currSnapshot_->getTime() >= currRestTime_){
      restOut->writeRest(restInfo_);
      currRestTime_ += restTime_;
    }
  }

  RealType RestraintForceManager::doRestraints(RealType scalingFactor){
    std::vector<Molecule*>::const_iterator rm;
    GenericData* data;
    Molecule::IntegrableObjectIterator ioi;
    MolecularRestraint* mRest;
    StuntDouble* sd;

    std::vector<StuntDouble*>::const_iterator ro;
    ObjectRestraint* oRest;

    std::map<int, Restraint::RealPair> restInfo;

    unscaledPotential_ = 0.0;

    restInfo_.clear();

    for(rm=restrainedMols_.begin(); rm != restrainedMols_.end(); ++rm){

      // make sure this molecule (*rm) has a generic data for restraints:
      data = (*rm)->getPropertyByName("Restraint");
      if (data != NULL) {
        // make sure we can reinterpret the generic data as restraint data:
        RestraintData* restData= dynamic_cast<RestraintData*>(data);        
        if (restData != NULL) {
          // make sure we can reinterpet the restraint data as a pointer to
          // an MolecularRestraint:
          mRest = dynamic_cast<MolecularRestraint*>(restData->getData());
          if (mRest == NULL) {
            sprintf( painCave.errMsg,
                     "Can not cast RestraintData to MolecularRestraint\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();                      
          }
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to RestraintData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();       
        }
      } else {
        sprintf( painCave.errMsg, "Can not find Restraint for RestrainedObject\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();          
      }

      // phew.  At this point, we should have the pointer to the
      // correct MolecularRestraint in the variable mRest.

      Vector3d molCom = (*rm)->getCom();

      std::vector<Vector3d> struc;
      std::vector<Vector3d> forces;
      
      for(sd = (*rm)->beginIntegrableObject(ioi); sd != NULL; 
          sd = (*rm)->nextIntegrableObject(ioi)) {
        struc.push_back(sd->getPos());
      }

      mRest->setScaleFactor(scalingFactor);
      mRest->calcForce(struc, molCom);
      forces = mRest->getRestraintForces();
      int index = 0;

      for(sd = (*rm)->beginIntegrableObject(ioi); sd != NULL; 
          sd = (*rm)->nextIntegrableObject(ioi)) {        
        sd->addFrc(forces[index]);
        struc.push_back(sd->getPos());
        index++;
      }
      
      unscaledPotential_ += mRest->getUnscaledPotential();

      restInfo = mRest->getRestraintInfo();

      // only collect data on restraints that we're going to print:
      if (mRest->getPrintRestraint()) 
        restInfo_.push_back(restInfo);
    }

    for(ro=restrainedObjs_.begin(); ro != restrainedObjs_.end(); ++ro){
      // make sure this object (*ro) has a generic data for restraints:
      data = (*ro)->getPropertyByName("Restraint");
      if (data != NULL) {
        // make sure we can reinterpret the generic data as restraint data:
        RestraintData* restData= dynamic_cast<RestraintData*>(data);        
        if (restData != NULL) {
          // make sure we can reinterpet the restraint data as a pointer to
          // an ObjectRestraint:
          oRest = dynamic_cast<ObjectRestraint*>(restData->getData());
          if (oRest == NULL) {
            sprintf( painCave.errMsg,
                     "Can not cast RestraintData to ObjectRestraint\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();                      
          }
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to RestraintData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();       
        }
      } else {
        sprintf( painCave.errMsg, "Can not find Restraint for RestrainedObject\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
      
      // phew.  At this point, we should have the pointer to the
      // correct Object restraint in the variable oRest.

      oRest->setScaleFactor(scalingFactor);
      
      Vector3d pos = (*ro)->getPos();

      if ( (*ro)->isDirectional() ) {

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

      restInfo = oRest->getRestraintInfo();

      // only collect data on restraints that we're going to print:
      if (oRest->getPrintRestraint()) 
        restInfo_.push_back(restInfo);
    }

    return unscaledPotential_ * scalingFactor;
  }
}
