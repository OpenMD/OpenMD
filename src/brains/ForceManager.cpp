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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file ForceManager.cpp
 * @author tlin
 * @date 11/09/2004
 * @version 1.0
 */


#include "brains/ForceManager.hpp"
#include "primitives/Molecule.hpp"
#define __OPENMD_C
#include "utils/simError.h"
#include "primitives/Bond.hpp"
#include "primitives/Bend.hpp"
#include "primitives/Torsion.hpp"
#include "primitives/Inversion.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "perturbations/ElectricField.hpp"
#include "parallel/ForceMatrixDecomposition.hpp"

#include <cstdio>
#include <iostream>
#include <iomanip>

using namespace std;
namespace OpenMD {
  
  ForceManager::ForceManager(SimInfo * info) : info_(info), switcher_(NULL),
                                               initialized_(false) {
    forceField_ = info_->getForceField();
    interactionMan_ = new InteractionManager();
    fDecomp_ = new ForceMatrixDecomposition(info_, interactionMan_);
    thermo = new Thermo(info_);
  }

  ForceManager::~ForceManager() {
    perturbations_.clear();
    
    delete switcher_;
    delete interactionMan_;
    delete fDecomp_;
    delete thermo;
  }
  
  /**
   * setupCutoffs
   *
   * Sets the values of cutoffRadius, switchingRadius, cutoffMethod,
   * and cutoffPolicy
   *
   * cutoffRadius : realType
   *  If the cutoffRadius was explicitly set, use that value.
   *  If the cutoffRadius was not explicitly set:
   *      Are there electrostatic atoms?  Use 12.0 Angstroms.
   *      No electrostatic atoms?  Poll the atom types present in the
   *      simulation for suggested cutoff values (e.g. 2.5 * sigma).
   *      Use the maximum suggested value that was found.
   *
   * cutoffMethod : (one of HARD, SWITCHED, SHIFTED_FORCE, TAYLOR_SHIFTED, 
   *                        SHIFTED_POTENTIAL, or EWALD_FULL)
   *      If cutoffMethod was explicitly set, use that choice.
   *      If cutoffMethod was not explicitly set, use SHIFTED_FORCE
   *
   * cutoffPolicy : (one of MIX, MAX, TRADITIONAL)
   *      If cutoffPolicy was explicitly set, use that choice.
   *      If cutoffPolicy was not explicitly set, use TRADITIONAL
   *
   * switchingRadius : realType
   *  If the cutoffMethod was set to SWITCHED:
   *      If the switchingRadius was explicitly set, use that value
   *          (but do a sanity check first).
   *      If the switchingRadius was not explicitly set: use 0.85 *
   *      cutoffRadius_
   *  If the cutoffMethod was not set to SWITCHED:
   *      Set switchingRadius equal to cutoffRadius for safety.
   */
  void ForceManager::setupCutoffs() {
    
    Globals* simParams_ = info_->getSimParams();
    ForceFieldOptions& forceFieldOptions_ = forceField_->getForceFieldOptions();
    int mdFileVersion;
    rCut_ = 0.0; //Needs a value for a later max() call;   
    
    if (simParams_->haveMDfileVersion()) 
      mdFileVersion = simParams_->getMDfileVersion();
    else
      mdFileVersion = 0;
   
    // We need the list of simulated atom types to figure out cutoffs
    // as well as long range corrections.

    set<AtomType*>::iterator i;
    set<AtomType*> atomTypes_;
    atomTypes_ = info_->getSimulatedAtomTypes();

    if (simParams_->haveCutoffRadius()) {
      rCut_ = simParams_->getCutoffRadius();
    } else {      
      if (info_->usesElectrostaticAtoms()) {
        sprintf(painCave.errMsg,
                "ForceManager::setupCutoffs: No value was set for the cutoffRadius.\n"
                "\tOpenMD will use a default value of 12.0 angstroms"
                "\tfor the cutoffRadius.\n");
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();
        rCut_ = 12.0;
      } else {
        RealType thisCut;
        for (i = atomTypes_.begin(); i != atomTypes_.end(); ++i) {
          thisCut = interactionMan_->getSuggestedCutoffRadius((*i));
          rCut_ = max(thisCut, rCut_);
        }
        sprintf(painCave.errMsg,
                "ForceManager::setupCutoffs: No value was set for the cutoffRadius.\n"
                "\tOpenMD will use %lf angstroms.\n",
                rCut_);
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();
      }
    }

    fDecomp_->setUserCutoff(rCut_);
    interactionMan_->setCutoffRadius(rCut_);

    map<string, CutoffMethod> stringToCutoffMethod;
    stringToCutoffMethod["HARD"] = HARD;
    stringToCutoffMethod["SWITCHED"] = SWITCHED;
    stringToCutoffMethod["SHIFTED_POTENTIAL"] = SHIFTED_POTENTIAL;    
    stringToCutoffMethod["SHIFTED_FORCE"] = SHIFTED_FORCE;
    stringToCutoffMethod["TAYLOR_SHIFTED"] = TAYLOR_SHIFTED;
    stringToCutoffMethod["EWALD_FULL"] = EWALD_FULL;
  
    if (simParams_->haveCutoffMethod()) {
      string cutMeth = toUpperCopy(simParams_->getCutoffMethod());
      map<string, CutoffMethod>::iterator i;
      i = stringToCutoffMethod.find(cutMeth);
      if (i == stringToCutoffMethod.end()) {
        sprintf(painCave.errMsg,
                "ForceManager::setupCutoffs: Could not find chosen cutoffMethod %s\n"
                "\tShould be one of: "
                "HARD, SWITCHED, SHIFTED_POTENTIAL, TAYLOR_SHIFTED,\n"
                "\tSHIFTED_FORCE, or EWALD_FULL\n",
                cutMeth.c_str());
        painCave.isFatal = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      } else {
        cutoffMethod_ = i->second;
      }
    } else {
      if (mdFileVersion > 1) {
        sprintf(painCave.errMsg,
                "ForceManager::setupCutoffs: No value was set for the cutoffMethod.\n"
                "\tOpenMD will use SHIFTED_FORCE.\n");
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();
        cutoffMethod_ = SHIFTED_FORCE;        
      } else {
        // handle the case where the old file version was in play
        // (there should be no cutoffMethod, so we have to deduce it
        // from other data).        

        sprintf(painCave.errMsg,
                "ForceManager::setupCutoffs : DEPRECATED FILE FORMAT!\n"
                "\tOpenMD found a file which does not set a cutoffMethod.\n"
                "\tOpenMD will attempt to deduce a cutoffMethod using the\n"
                "\tbehavior of the older (version 1) code.  To remove this\n"
                "\twarning, add an explicit cutoffMethod and change the top\n"
                "\tof the file so that it begins with <OpenMD version=2>\n");
        painCave.isFatal = 0;
        painCave.severity = OPENMD_WARNING;
        simError();            
                
        // The old file version tethered the shifting behavior to the
        // electrostaticSummationMethod keyword.
        
        if (simParams_->haveElectrostaticSummationMethod()) {
          string myMethod = simParams_->getElectrostaticSummationMethod();
          toUpper(myMethod);
        
          if (myMethod == "SHIFTED_POTENTIAL") {
            cutoffMethod_ = SHIFTED_POTENTIAL;
          } else if (myMethod == "SHIFTED_FORCE") {
            cutoffMethod_ = SHIFTED_FORCE;
          } else if (myMethod == "TAYLOR_SHIFTED") {
            cutoffMethod_ = TAYLOR_SHIFTED;
          } else if (myMethod == "EWALD_FULL") {
            cutoffMethod_ = EWALD_FULL;
          }
        
          if (simParams_->haveSwitchingRadius()) 
            rSwitch_ = simParams_->getSwitchingRadius();

          if (myMethod == "SHIFTED_POTENTIAL" || myMethod == "SHIFTED_FORCE" ||
              myMethod == "TAYLOR_SHIFTED" || myMethod == "EWALD_FULL") {
            if (simParams_->haveSwitchingRadius()){
              sprintf(painCave.errMsg,
                      "ForceManager::setupCutoffs : DEPRECATED ERROR MESSAGE\n"
                      "\tA value was set for the switchingRadius\n"
                      "\teven though the electrostaticSummationMethod was\n"
                      "\tset to %s\n", myMethod.c_str());
              painCave.severity = OPENMD_WARNING;
              painCave.isFatal = 1;
              simError();            
            } 
          }
          if (abs(rCut_ - rSwitch_) < 0.0001) {
            if (cutoffMethod_ == SHIFTED_FORCE) {              
              sprintf(painCave.errMsg,
                      "ForceManager::setupCutoffs : DEPRECATED BEHAVIOR\n" 
                      "\tcutoffRadius and switchingRadius are set to the\n"
                      "\tsame value.  OpenMD will use shifted force\n"
                      "\tpotentials instead of switching functions.\n");
              painCave.isFatal = 0;
              painCave.severity = OPENMD_WARNING;
              simError();            
            } else {
              cutoffMethod_ = SHIFTED_POTENTIAL;
              sprintf(painCave.errMsg,
                      "ForceManager::setupCutoffs : DEPRECATED BEHAVIOR\n" 
                      "\tcutoffRadius and switchingRadius are set to the\n"
                      "\tsame value.  OpenMD will use shifted potentials\n"
                      "\tinstead of switching functions.\n");
              painCave.isFatal = 0;
              painCave.severity = OPENMD_WARNING;
              simError();            
            }
          }
        }
      }
    }

    map<string, CutoffPolicy> stringToCutoffPolicy;
    stringToCutoffPolicy["MIX"] = MIX;
    stringToCutoffPolicy["MAX"] = MAX;
    stringToCutoffPolicy["TRADITIONAL"] = TRADITIONAL;    

    string cutPolicy;
    if (forceFieldOptions_.haveCutoffPolicy()){
      cutPolicy = forceFieldOptions_.getCutoffPolicy();
    }else if (simParams_->haveCutoffPolicy()) {
      cutPolicy = simParams_->getCutoffPolicy();
    }

    if (!cutPolicy.empty()){
      toUpper(cutPolicy);
      map<string, CutoffPolicy>::iterator i;
      i = stringToCutoffPolicy.find(cutPolicy);

      if (i == stringToCutoffPolicy.end()) {
        sprintf(painCave.errMsg,
                "ForceManager::setupCutoffs: Could not find chosen cutoffPolicy %s\n"
                "\tShould be one of: "
                "MIX, MAX, or TRADITIONAL\n",
                cutPolicy.c_str());
        painCave.isFatal = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      } else {
        cutoffPolicy_ = i->second;
      }
    } else {
      sprintf(painCave.errMsg,
              "ForceManager::setupCutoffs: No value was set for the cutoffPolicy.\n"
              "\tOpenMD will use TRADITIONAL.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      cutoffPolicy_ = TRADITIONAL;        
    }

    fDecomp_->setCutoffPolicy(cutoffPolicy_);
        
    // create the switching function object:

    switcher_ = new SwitchingFunction();
   
    if (cutoffMethod_ == SWITCHED) {
      if (simParams_->haveSwitchingRadius()) {
        rSwitch_ = simParams_->getSwitchingRadius();
        if (rSwitch_ > rCut_) {        
          sprintf(painCave.errMsg,
                  "ForceManager::setupCutoffs: switchingRadius (%f) is larger "
                  "than the cutoffRadius(%f)\n", rSwitch_, rCut_);
          painCave.isFatal = 1;
          painCave.severity = OPENMD_ERROR;
          simError();
        }
      } else {      
        rSwitch_ = 0.85 * rCut_;
        sprintf(painCave.errMsg,
                "ForceManager::setupCutoffs: No value was set for the switchingRadius.\n"
                "\tOpenMD will use a default value of 85 percent of the cutoffRadius.\n"
                "\tswitchingRadius = %f. for this simulation\n", rSwitch_);
        painCave.isFatal = 0;
        painCave.severity = OPENMD_WARNING;
        simError();
      }
    } else {
      if (mdFileVersion > 1) {
        // throw an error if we define a switching radius and don't need one.
        // older file versions should not do this.
        if (simParams_->haveSwitchingRadius()) {
          map<string, CutoffMethod>::const_iterator it;
          string theMeth;
          for (it = stringToCutoffMethod.begin(); 
               it != stringToCutoffMethod.end(); ++it) {
            if (it->second == cutoffMethod_) {
              theMeth = it->first;
              break;
            }
          }
          sprintf(painCave.errMsg,
                  "ForceManager::setupCutoffs: the cutoffMethod (%s)\n"
                  "\tis not set to SWITCHED, so switchingRadius value\n"
                  "\twill be ignored for this simulation\n", theMeth.c_str());
          painCave.isFatal = 0;
          painCave.severity = OPENMD_WARNING;
          simError();
        }
      }
      rSwitch_ = rCut_;
    }
    
    // Default to cubic switching function.
    sft_ = cubic;
    if (simParams_->haveSwitchingFunctionType()) {
      string funcType = simParams_->getSwitchingFunctionType();
      toUpper(funcType);
      if (funcType == "CUBIC") {
        sft_ = cubic;
      } else {
        if (funcType == "FIFTH_ORDER_POLYNOMIAL") {
          sft_ = fifth_order_poly;
        } else {
          // throw error        
          sprintf( painCave.errMsg,
                   "ForceManager::setupSwitching : Unknown switchingFunctionType. (Input file specified %s .)\n"
                   "\tswitchingFunctionType must be one of: "
                   "\"cubic\" or \"fifth_order_polynomial\".", 
                   funcType.c_str() );
          painCave.isFatal = 1;
          painCave.severity = OPENMD_ERROR;
          simError();
        }           
      }
    }
    switcher_->setSwitchType(sft_);
    switcher_->setSwitch(rSwitch_, rCut_);
  }



  
  void ForceManager::initialize() {

    if (!info_->isTopologyDone()) {

      info_->update();
      interactionMan_->setSimInfo(info_);
      interactionMan_->initialize();

      // We want to delay the cutoffs until after the interaction
      // manager has set up the atom-atom interactions so that we can
      // query them for suggested cutoff values
      setupCutoffs();

      info_->prepareTopology();      

      doParticlePot_ = info_->getSimParams()->getOutputParticlePotential();
      doHeatFlux_ = info_->getSimParams()->getPrintHeatFlux();
      if (doHeatFlux_) doParticlePot_ = true;

      doElectricField_ = info_->getSimParams()->getOutputElectricField();
   
    }

    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    
    // Force fields can set options on how to scale van der Waals and
    // electrostatic interactions for atoms connected via bonds, bends
    // and torsions in this case the topological distance between
    // atoms is:
    // 0 = topologically unconnected
    // 1 = bonded together 
    // 2 = connected via a bend
    // 3 = connected via a torsion
    
    vdwScale_.reserve(4);
    fill(vdwScale_.begin(), vdwScale_.end(), 0.0);

    electrostaticScale_.reserve(4);
    fill(electrostaticScale_.begin(), electrostaticScale_.end(), 0.0);

    vdwScale_[0] = 1.0;
    vdwScale_[1] = fopts.getvdw12scale();
    vdwScale_[2] = fopts.getvdw13scale();
    vdwScale_[3] = fopts.getvdw14scale();
    
    electrostaticScale_[0] = 1.0;
    electrostaticScale_[1] = fopts.getelectrostatic12scale();
    electrostaticScale_[2] = fopts.getelectrostatic13scale();
    electrostaticScale_[3] = fopts.getelectrostatic14scale();    
    
    if (info_->getSimParams()->haveElectricField()) {
      ElectricField* eField = new ElectricField(info_);
      perturbations_.push_back(eField);
    }

    usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();
    
    fDecomp_->distributeInitialData();
    
    initialized_ = true;
    
  }
  
  void ForceManager::calcForces() {
    
    if (!initialized_) initialize();
    
    preCalculation();   
    shortRangeInteractions();
    longRangeInteractions();
    postCalculation();    
  }
  
  void ForceManager::preCalculation() {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;
    
    // forces and potentials are zeroed here, before any are
    // accumulated.
    
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    snap->setBondPotential(0.0);
    snap->setBendPotential(0.0);
    snap->setTorsionPotential(0.0);
    snap->setInversionPotential(0.0);

    potVec zeroPot(0.0);
    snap->setLongRangePotential(zeroPot);
    snap->setExcludedPotentials(zeroPot);

    snap->setRestraintPotential(0.0);
    snap->setRawPotential(0.0);

    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      for(atom = mol->beginAtom(ai); atom != NULL; 
          atom = mol->nextAtom(ai)) {
	atom->zeroForcesAndTorques();
      }
      
      //change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
	rb->zeroForcesAndTorques();
      }        
      
      if(info_->getNGlobalCutoffGroups() != info_->getNGlobalAtoms()){
        for(cg = mol->beginCutoffGroup(ci); cg != NULL; 
            cg = mol->nextCutoffGroup(ci)) {
          //calculate the center of mass of cutoff group
          cg->updateCOM();
        }
      }      
    }
    
    // Zero out the stress tensor
    stressTensor *= 0.0;
    // Zero out the heatFlux
    fDecomp_->setHeatFlux( Vector3d(0.0) );    
  }
  
  void ForceManager::shortRangeInteractions() {
    Molecule* mol;
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;;
    Molecule::BendIterator  bendIter;
    Molecule::TorsionIterator  torsionIter;
    Molecule::InversionIterator  inversionIter;
    RealType bondPotential = 0.0;
    RealType bendPotential = 0.0;
    RealType torsionPotential = 0.0;
    RealType inversionPotential = 0.0;

    //calculate short range interactions    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {

      //change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
        rb->updateAtoms();
      }

      for (bond = mol->beginBond(bondIter); bond != NULL; 
           bond = mol->nextBond(bondIter)) {
        bond->calcForce(doParticlePot_);
        bondPotential += bond->getPotential();
      }

      for (bend = mol->beginBend(bendIter); bend != NULL; 
           bend = mol->nextBend(bendIter)) {
        
        RealType angle;
        bend->calcForce(angle, doParticlePot_);
        RealType currBendPot = bend->getPotential();          
         
        bendPotential += bend->getPotential();
        map<Bend*, BendDataSet>::iterator i = bendDataSets.find(bend);
        if (i == bendDataSets.end()) {
          BendDataSet dataSet;
          dataSet.prev.angle = dataSet.curr.angle = angle;
          dataSet.prev.potential = dataSet.curr.potential = currBendPot;
          dataSet.deltaV = 0.0;
          bendDataSets.insert(map<Bend*, BendDataSet>::value_type(bend, 
                                                                  dataSet));
        }else {
          i->second.prev.angle = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle = angle;
          i->second.curr.potential = currBendPot;
          i->second.deltaV =  fabs(i->second.curr.potential -  
                                   i->second.prev.potential);
        }
      }
      
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL; 
           torsion = mol->nextTorsion(torsionIter)) {
        RealType angle;
        torsion->calcForce(angle, doParticlePot_);
        RealType currTorsionPot = torsion->getPotential();
        torsionPotential += torsion->getPotential();
        map<Torsion*, TorsionDataSet>::iterator i = torsionDataSets.find(torsion);
        if (i == torsionDataSets.end()) {
          TorsionDataSet dataSet;
          dataSet.prev.angle = dataSet.curr.angle = angle;
          dataSet.prev.potential = dataSet.curr.potential = currTorsionPot;
          dataSet.deltaV = 0.0;
          torsionDataSets.insert(map<Torsion*, TorsionDataSet>::value_type(torsion, dataSet));
        }else {
          i->second.prev.angle = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle = angle;
          i->second.curr.potential = currTorsionPot;
          i->second.deltaV =  fabs(i->second.curr.potential -  
                                   i->second.prev.potential);
        }      
      }      
      
      for (inversion = mol->beginInversion(inversionIter); 
	   inversion != NULL; 
           inversion = mol->nextInversion(inversionIter)) {
        RealType angle;
        inversion->calcForce(angle, doParticlePot_);
        RealType currInversionPot = inversion->getPotential();
        inversionPotential += inversion->getPotential();
        map<Inversion*, InversionDataSet>::iterator i = inversionDataSets.find(inversion);
        if (i == inversionDataSets.end()) {
          InversionDataSet dataSet;
          dataSet.prev.angle = dataSet.curr.angle = angle;
          dataSet.prev.potential = dataSet.curr.potential = currInversionPot;
          dataSet.deltaV = 0.0;
          inversionDataSets.insert(map<Inversion*, InversionDataSet>::value_type(inversion, dataSet));
        }else {
          i->second.prev.angle = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle = angle;
          i->second.curr.potential = currInversionPot;
          i->second.deltaV =  fabs(i->second.curr.potential -  
                                   i->second.prev.potential);
        }      
      }      
    }

#ifdef IS_MPI
    // Collect from all nodes.  This should eventually be moved into a
    // SystemDecomposition, but this is a better place than in
    // Thermo to do the collection.
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &bondPotential, 1, MPI::REALTYPE, 
                              MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &bendPotential, 1, MPI::REALTYPE, 
                              MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &torsionPotential, 1, 
                              MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &inversionPotential, 1, 
                              MPI::REALTYPE, MPI::SUM);
#endif

    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();

    curSnapshot->setBondPotential(bondPotential);
    curSnapshot->setBendPotential(bendPotential);
    curSnapshot->setTorsionPotential(torsionPotential);
    curSnapshot->setInversionPotential(inversionPotential);
    
    // RealType shortRangePotential = bondPotential + bendPotential + 
    //   torsionPotential +  inversionPotential;    

    // curSnapshot->setShortRangePotential(shortRangePotential);
  }
  
  void ForceManager::longRangeInteractions() {

    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    DataStorage* config = &(curSnapshot->atomData);
    DataStorage* cgConfig = &(curSnapshot->cgData);


    //calculate the center of mass of cutoff group

    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;

    if(info_->getNCutoffGroups() > 0){      
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        for(cg = mol->beginCutoffGroup(ci); cg != NULL; 
            cg = mol->nextCutoffGroup(ci)) {
          cg->updateCOM();
        }
      }      
    } else {
      // center of mass of the group is the same as position of the atom  
      // if cutoff group does not exist
      cgConfig->position = config->position;
      cgConfig->velocity = config->velocity;
    }

    fDecomp_->zeroWorkArrays();
    fDecomp_->distributeData();
    
    int cg1, cg2, atom1, atom2, topoDist;
    Vector3d d_grp, dag, d, gvel2, vel2;
    RealType rgrpsq, rgrp, r2, r;
    RealType electroMult, vdwMult;
    RealType vij;
    Vector3d fij, fg, f1;
    tuple3<RealType, RealType, RealType> cuts;
    RealType rCut, rCutSq, rListSq;
    bool in_switching_region;
    RealType sw, dswdr, swderiv;
    vector<int> atomListColumn, atomListRow;
    InteractionData idat;
    SelfData sdat;
    RealType mf;
    RealType vpair;
    RealType dVdFQ1(0.0);
    RealType dVdFQ2(0.0);
    potVec longRangePotential(0.0);
    potVec reciprocalPotential(0.0);
    potVec workPot(0.0);
    potVec exPot(0.0);
    Vector3d eField1(0.0);
    Vector3d eField2(0.0);
    vector<int>::iterator ia, jb;

    int loopStart, loopEnd;
    
    idat.rcut = &rCut;
    idat.vdwMult = &vdwMult;
    idat.electroMult = &electroMult;
    idat.pot = &workPot;
    idat.excludedPot = &exPot;
    sdat.pot = fDecomp_->getEmbeddingPotential();
    sdat.excludedPot = fDecomp_->getExcludedSelfPotential();
    idat.vpair = &vpair;
    idat.dVdFQ1 = &dVdFQ1;
    idat.dVdFQ2 = &dVdFQ2;
    idat.eField1 = &eField1;
    idat.eField2 = &eField2;   
    idat.f1 = &f1;
    idat.sw = &sw;
    idat.shiftedPot = (cutoffMethod_ == SHIFTED_POTENTIAL) ? true : false;
    idat.shiftedForce = (cutoffMethod_ == SHIFTED_FORCE || cutoffMethod_ == TAYLOR_SHIFTED) ? true : false;
    idat.doParticlePot = doParticlePot_;
    idat.doElectricField = doElectricField_;
    sdat.doParticlePot = doParticlePot_;
    
    loopEnd = PAIR_LOOP;
    if (info_->requiresPrepair() ) {
      loopStart = PREPAIR_LOOP;
    } else {
      loopStart = PAIR_LOOP;
    }
    for (int iLoop = loopStart; iLoop <= loopEnd; iLoop++) {
    
      if (iLoop == loopStart) {
        bool update_nlist = fDecomp_->checkNeighborList();
        if (update_nlist) {
          if (!usePeriodicBoundaryConditions_)
            Mat3x3d bbox = thermo->getBoundingBox();
          fDecomp_->buildNeighborList(neighborList_);
        }
      }

      for (vector<pair<int, int> >::iterator it = neighborList_.begin(); 
             it != neighborList_.end(); ++it) {
                
        cg1 = (*it).first;
        cg2 = (*it).second;
        
        fDecomp_->getGroupCutoffs(cg1, cg2, rCut, rCutSq, rListSq);

        d_grp  = fDecomp_->getIntergroupVector(cg1, cg2);

        // already wrapped in the getIntergroupVector call:
        // curSnapshot->wrapVector(d_grp);        
        rgrpsq = d_grp.lengthSquare();

        if (rgrpsq < rCutSq) {
          if (iLoop == PAIR_LOOP) {
            vij = 0.0;
            fij.zero();
            eField1.zero();
            eField2.zero();
          }
          
          in_switching_region = switcher_->getSwitch(rgrpsq, sw, dswdr, 
                                                     rgrp); 

          atomListRow = fDecomp_->getAtomsInGroupRow(cg1);
          atomListColumn = fDecomp_->getAtomsInGroupColumn(cg2);

          if (doHeatFlux_)
            gvel2 = fDecomp_->getGroupVelocityColumn(cg2);

          for (ia = atomListRow.begin(); 
               ia != atomListRow.end(); ++ia) {            
            atom1 = (*ia);

            for (jb = atomListColumn.begin(); 
                 jb != atomListColumn.end(); ++jb) {              
              atom2 = (*jb);

              if (!fDecomp_->skipAtomPair(atom1, atom2, cg1, cg2)) {

                vpair = 0.0;
                workPot = 0.0;
                exPot = 0.0;
                f1.zero();
		dVdFQ1 = 0.0;
		dVdFQ2 = 0.0;

                fDecomp_->fillInteractionData(idat, atom1, atom2);

                topoDist = fDecomp_->getTopologicalDistance(atom1, atom2);
                vdwMult = vdwScale_[topoDist];
                electroMult = electrostaticScale_[topoDist];

                if (atomListRow.size() == 1 && atomListColumn.size() == 1) {
                  idat.d = &d_grp;
                  idat.r2 = &rgrpsq;
                  if (doHeatFlux_)
                    vel2 = gvel2;
                } else {
                  d = fDecomp_->getInteratomicVector(atom1, atom2);
                  curSnapshot->wrapVector( d );
                  r2 = d.lengthSquare();
                  idat.d = &d;
                  idat.r2 = &r2;
                  if (doHeatFlux_)
                    vel2 = fDecomp_->getAtomVelocityColumn(atom2);
                }
               
                r = sqrt( *(idat.r2) );
                idat.rij = &r;

                if (iLoop == PREPAIR_LOOP) {
                  interactionMan_->doPrePair(idat);
                } else {
                  interactionMan_->doPair(idat);
                  fDecomp_->unpackInteractionData(idat, atom1, atom2);
                  vij += vpair;
                  fij += f1;
                  stressTensor -= outProduct( *(idat.d), f1);
                  if (doHeatFlux_) 
                    fDecomp_->addToHeatFlux(*(idat.d) * dot(f1, vel2));
                }
              }
            }
          }

          if (iLoop == PAIR_LOOP) {
            if (in_switching_region) {
              swderiv = vij * dswdr / rgrp;
              fg = swderiv * d_grp;
              fij += fg;

              if (atomListRow.size() == 1 && atomListColumn.size() == 1) {
                if (!fDecomp_->skipAtomPair(atomListRow[0], 
                                            atomListColumn[0], 
                                            cg1, cg2)) {
                  stressTensor -= outProduct( *(idat.d), fg);
                  if (doHeatFlux_)
                    fDecomp_->addToHeatFlux(*(idat.d) * dot(fg, vel2));
                }                
              }
          
              for (ia = atomListRow.begin(); 
                   ia != atomListRow.end(); ++ia) {            
                atom1 = (*ia);                
                mf = fDecomp_->getMassFactorRow(atom1);
                // fg is the force on atom ia due to cutoff group's
                // presence in switching region
                fg = swderiv * d_grp * mf;
                fDecomp_->addForceToAtomRow(atom1, fg);
                if (atomListRow.size() > 1) {
                  if (info_->usesAtomicVirial()) {
                    // find the distance between the atom
                    // and the center of the cutoff group:
                    dag = fDecomp_->getAtomToGroupVectorRow(atom1, cg1);
                    stressTensor -= outProduct(dag, fg);
                    if (doHeatFlux_)
                      fDecomp_->addToHeatFlux( dag * dot(fg, vel2));
                  }
                }
              }
              for (jb = atomListColumn.begin(); 
                   jb != atomListColumn.end(); ++jb) {              
                atom2 = (*jb);
                mf = fDecomp_->getMassFactorColumn(atom2);
                // fg is the force on atom jb due to cutoff group's
                // presence in switching region
                fg = -swderiv * d_grp * mf;
                fDecomp_->addForceToAtomColumn(atom2, fg);

                if (atomListColumn.size() > 1) {
                  if (info_->usesAtomicVirial()) {
                    // find the distance between the atom
                    // and the center of the cutoff group:
                    dag = fDecomp_->getAtomToGroupVectorColumn(atom2, cg2);
                    stressTensor -= outProduct(dag, fg);
                    if (doHeatFlux_)
                      fDecomp_->addToHeatFlux( dag * dot(fg, vel2));
                  }
                }
              }
            }
            //if (!info_->usesAtomicVirial()) {
            //  stressTensor -= outProduct(d_grp, fij);
            //  if (doHeatFlux_)
            //     fDecomp_->addToHeatFlux( d_grp * dot(fij, vel2));
            //}
          }
        }
      }

      if (iLoop == PREPAIR_LOOP) {
        if (info_->requiresPrepair()) {

          fDecomp_->collectIntermediateData();

          for (unsigned int atom1 = 0; atom1 < info_->getNAtoms(); atom1++) {
            fDecomp_->fillSelfData(sdat, atom1);
            interactionMan_->doPreForce(sdat);
          }

          fDecomp_->distributeIntermediateData();

        }
      }
    }
    
    // collects pairwise information
    fDecomp_->collectData();
    if (cutoffMethod_ == EWALD_FULL) {
      interactionMan_->doReciprocalSpaceSum(reciprocalPotential);
    }
        
    if (info_->requiresSelfCorrection()) {
      for (unsigned int atom1 = 0; atom1 < info_->getNAtoms(); atom1++) { 
        fDecomp_->fillSelfData(sdat, atom1);
        interactionMan_->doSelfCorrection(sdat);
      }
    }

    // collects single-atom information
    fDecomp_->collectSelfData();

    longRangePotential = *(fDecomp_->getEmbeddingPotential()) + 
      *(fDecomp_->getPairwisePotential()) + reciprocalPotential;

    curSnapshot->setLongRangePotential(longRangePotential);
    
    curSnapshot->setExcludedPotentials(*(fDecomp_->getExcludedSelfPotential()) +
                                         *(fDecomp_->getExcludedPotential()));

  }

  void ForceManager::postCalculation() {

    vector<Perturbation*>::iterator pi;
    for (pi = perturbations_.begin(); pi != perturbations_.end(); ++pi) {
      (*pi)->applyPerturbation();
    }

    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
   
    // collect the atomic forces onto rigid bodies
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) { 
        Mat3x3d rbTau = rb->calcForcesAndTorquesAndVirial();
        stressTensor += rbTau;
      }
    }
    
#ifdef IS_MPI
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, stressTensor.getArrayPointer(), 9, 
                              MPI::REALTYPE, MPI::SUM);
#endif
    curSnapshot->setStressTensor(stressTensor);
    
    if (info_->getSimParams()->getUseLongRangeCorrections()) {
      /*
      RealType vol = curSnapshot->getVolume();
      RealType Elrc(0.0);
      RealType Wlrc(0.0);

      set<AtomType*>::iterator i;
      set<AtomType*>::iterator j;
    
      RealType n_i, n_j;
      RealType rho_i, rho_j;
      pair<RealType, RealType> LRI;
      
      for (i = atomTypes_.begin(); i != atomTypes_.end(); ++i) {
        n_i = RealType(info_->getGlobalCountOfType(*i));
        rho_i = n_i /  vol;
        for (j = atomTypes_.begin(); j != atomTypes_.end(); ++j) {
          n_j = RealType(info_->getGlobalCountOfType(*j));
          rho_j = n_j / vol;
          
          LRI = interactionMan_->getLongRangeIntegrals( (*i), (*j) );

          Elrc += n_i   * rho_j * LRI.first;
          Wlrc -= rho_i * rho_j * LRI.second;
        }
      }
      Elrc *= 2.0 * NumericConstant::PI;
      Wlrc *= 2.0 * NumericConstant::PI;

      RealType lrp = curSnapshot->getLongRangePotential();
      curSnapshot->setLongRangePotential(lrp + Elrc);
      stressTensor += Wlrc * SquareMatrix3<RealType>::identity();
      curSnapshot->setStressTensor(stressTensor);
      */
     
    }
  }
}
