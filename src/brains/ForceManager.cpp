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

/**
 * @file ForceManager.cpp
 * @author tlin
 * @date 11/09/2004
 * @version 1.0
 */

#include "brains/ForceManager.hpp"

#define __OPENMD_C

#include <cstdio>
#include <iomanip>
#include <iostream>

#include "constraints/ZconstraintForceModifier.hpp"
#include "integrators/LDForceModifier.hpp"
#include "integrators/LangevinHullForceModifier.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "parallel/ForceMatrixDecomposition.hpp"
#include "perturbations/Light.hpp"
#include "perturbations/LightParameters.hpp"
#include "perturbations/MagneticField.hpp"
#include "perturbations/UniformField.hpp"
#include "perturbations/UniformGradient.hpp"
#include "primitives/Bend.hpp"
#include "primitives/Bond.hpp"
#include "primitives/Inversion.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/Torsion.hpp"
#include "restraints/RestraintForceModifier.hpp"
#include "restraints/ThermoIntegrationForceModifier.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  ForceManager::ForceManager(SimInfo* info) :
      initialized_(false), info_(info), switcher_(NULL), seleMan_(info),
      evaluator_(info) {
    forceField_     = info_->getForceField();
    interactionMan_ = new InteractionManager();
    fDecomp_        = new ForceMatrixDecomposition(info_, interactionMan_);
    thermo          = new Thermo(info_);
  }

  ForceManager::~ForceManager() {
    Utils::deletePointers(forceModifiers_);

    delete switcher_;
    delete interactionMan_;
    delete fDecomp_;
    delete thermo;
  }

  /**
   * setupCutoffs
   *
   * Sets the values of cutoffRadius, switchingRadius, and cutoffMethod
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
    int mdFileVersion;
    rCut_ = 0.0;  // Needs a value for a later max() call;

    if (simParams_->haveMDfileVersion())
      mdFileVersion = simParams_->getMDfileVersion();
    else
      mdFileVersion = 0;

    // We need the list of simulated atom types to figure out cutoffs
    // as well as long range corrections.

    AtomTypeSet::iterator i;
    AtomTypeSet atomTypes_;
    atomTypes_ = info_->getSimulatedAtomTypes();

    if (simParams_->haveCutoffRadius()) {
      rCut_ = simParams_->getCutoffRadius();
    } else {
      if (info_->usesElectrostaticAtoms()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "ForceManager::setupCutoffs: No value was set for the "
                 "cutoffRadius.\n"
                 "\tOpenMD will use a default value of 12.0 angstroms "
                 "for the cutoffRadius.\n");
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_INFO;
        simError();
        rCut_ = 12.0;
      } else {
        RealType thisCut;
        for (i = atomTypes_.begin(); i != atomTypes_.end(); ++i) {
          thisCut = interactionMan_->getSuggestedCutoffRadius((*i));
          rCut_   = max(thisCut, rCut_);
        }
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "ForceManager::setupCutoffs: No value was set for the "
                 "cutoffRadius.\n"
                 "\tOpenMD will use %lf angstroms.\n",
                 rCut_);
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_INFO;
        simError();
      }
    }

    fDecomp_->setCutoffRadius(rCut_);
    interactionMan_->setCutoffRadius(rCut_);
    rCutSq_ = rCut_ * rCut_;

    map<string, CutoffMethod> stringToCutoffMethod;
    stringToCutoffMethod["HARD"]              = HARD;
    stringToCutoffMethod["SWITCHED"]          = SWITCHED;
    stringToCutoffMethod["SHIFTED_POTENTIAL"] = SHIFTED_POTENTIAL;
    stringToCutoffMethod["SHIFTED_FORCE"]     = SHIFTED_FORCE;
    stringToCutoffMethod["TAYLOR_SHIFTED"]    = TAYLOR_SHIFTED;
    stringToCutoffMethod["EWALD_FULL"]        = EWALD_FULL;

    if (simParams_->haveCutoffMethod()) {
      string cutMeth = toUpperCopy(simParams_->getCutoffMethod());
      map<string, CutoffMethod>::iterator i;
      i = stringToCutoffMethod.find(cutMeth);
      if (i == stringToCutoffMethod.end()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "ForceManager::setupCutoffs: Could not find chosen "
                 "cutoffMethod %s\n"
                 "\tShould be one of: "
                 "HARD, SWITCHED, SHIFTED_POTENTIAL, TAYLOR_SHIFTED,\n"
                 "\tSHIFTED_FORCE, or EWALD_FULL\n",
                 cutMeth.c_str());
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      } else {
        cutoffMethod_ = i->second;
      }
    } else {
      if (mdFileVersion > 1) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "ForceManager::setupCutoffs: No value was set "
                 "for the cutoffMethod.\n"
                 "\tOpenMD will use SHIFTED_FORCE.\n");
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_INFO;
        simError();
        cutoffMethod_ = SHIFTED_FORCE;
      } else {
        // handle the case where the old file version was in play
        // (there should be no cutoffMethod, so we have to deduce it
        // from other data).

        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "ForceManager::setupCutoffs : DEPRECATED FILE FORMAT!\n"
                 "\tOpenMD found a file which does not set a cutoffMethod.\n"
                 "\tOpenMD will attempt to deduce a cutoffMethod using the\n"
                 "\tbehavior of the older (version 1) code.  To remove this\n"
                 "\twarning, add an explicit cutoffMethod and change the top\n"
                 "\tof the file so that it begins with <OpenMD version=2>\n");
        painCave.isFatal  = 0;
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
            cutoffMethod_   = EWALD_FULL;
            useSurfaceTerm_ = true;
          }

          if (simParams_->haveSwitchingRadius())
            rSwitch_ = simParams_->getSwitchingRadius();

          if (myMethod == "SHIFTED_POTENTIAL" || myMethod == "SHIFTED_FORCE" ||
              myMethod == "TAYLOR_SHIFTED" || myMethod == "EWALD_FULL") {
            if (simParams_->haveSwitchingRadius()) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "ForceManager::setupCutoffs : DEPRECATED ERROR MESSAGE\n"
                       "\tA value was set for the switchingRadius\n"
                       "\teven though the electrostaticSummationMethod was\n"
                       "\tset to %s\n",
                       myMethod.c_str());
              painCave.severity = OPENMD_WARNING;
              painCave.isFatal  = 1;
              simError();
            }
          }
          if (abs(rCut_ - rSwitch_) < 0.0001) {
            if (cutoffMethod_ == SHIFTED_FORCE) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "ForceManager::setupCutoffs : DEPRECATED BEHAVIOR\n"
                       "\tcutoffRadius and switchingRadius are set to the\n"
                       "\tsame value.  OpenMD will use shifted force\n"
                       "\tpotentials instead of switching functions.\n");
              painCave.isFatal  = 0;
              painCave.severity = OPENMD_WARNING;
              simError();
            } else {
              cutoffMethod_ = SHIFTED_POTENTIAL;
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "ForceManager::setupCutoffs : DEPRECATED BEHAVIOR\n"
                       "\tcutoffRadius and switchingRadius are set to the\n"
                       "\tsame value.  OpenMD will use shifted potentials\n"
                       "\tinstead of switching functions.\n");
              painCave.isFatal  = 0;
              painCave.severity = OPENMD_WARNING;
              simError();
            }
          }
        }
      }
    }

    // create the switching function object:

    switcher_ = new SwitchingFunction();

    if (cutoffMethod_ == SWITCHED) {
      if (simParams_->haveSwitchingRadius()) {
        rSwitch_ = simParams_->getSwitchingRadius();
        if (rSwitch_ > rCut_) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "ForceManager::setupCutoffs: switchingRadius (%f) is larger "
                   "than the cutoffRadius(%f)\n",
                   rSwitch_, rCut_);
          painCave.isFatal  = 1;
          painCave.severity = OPENMD_ERROR;
          simError();
        }
      } else {
        rSwitch_ = 0.85 * rCut_;
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "ForceManager::setupCutoffs: No value was set for the "
                 "switchingRadius.\n"
                 "\tOpenMD will use a default value of 85 percent of the "
                 "cutoffRadius.\n"
                 "\tswitchingRadius = %f. for this simulation\n",
                 rSwitch_);
        painCave.isFatal  = 0;
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
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "ForceManager::setupCutoffs: the cutoffMethod (%s)\n"
                   "\tis not set to SWITCHED, so switchingRadius value\n"
                   "\twill be ignored for this simulation\n",
                   theMeth.c_str());
          painCave.isFatal  = 0;
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
          snprintf(
              painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "ForceManager::setupSwitching : Unknown switchingFunctionType. "
              "(Input "
              "file specified %s .)\n"
              "\tswitchingFunctionType must be one of: "
              "\"cubic\" or \"fifth_order_polynomial\".",
              funcType.c_str());
          painCave.isFatal  = 1;
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

      useSurfaceTerm_  = false;  // default, can be set by Ewald or directly
      useSlabGeometry_ = false;

      //! We want to delay the cutoffs until after the interaction
      //! manager has set up the atom-atom interactions so that we can
      //! query them for suggested cutoff values
      setupCutoffs();

      info_->prepareTopology();

      doParticlePot_ = info_->getSimParams()->getOutputParticlePotential();
      doHeatFlux_    = info_->getSimParams()->getPrintHeatFlux();
      if (doHeatFlux_) doParticlePot_ = true;

      doElectricField_ = info_->getSimParams()->getOutputElectricField();
      doElectricField_ |=
          info_->getSimParams()->getRNEMDParameters()->requiresElectricField();
      doSitePotential_ = info_->getSimParams()->getOutputSitePotential();
      if (info_->getSimParams()->haveUseSurfaceTerm() &&
          info_->getSimParams()->getUseSurfaceTerm()) {
        if (info_->usesElectrostaticAtoms()) {
          useSurfaceTerm_ = info_->getSimParams()->getUseSurfaceTerm();
          if (info_->getSimParams()->haveUseSlabGeometry()) {
            useSlabGeometry_ = info_->getSimParams()->getUseSlabGeometry();

            string axis =
                toUpperCopy(info_->getSimParams()->getPrivilegedAxis());
            if (axis.compare("X") == 0)
              axis_ = 0;
            else if (axis.compare("Y") == 0)
              axis_ = 1;
            else
              axis_ = 2;
          }
        } else {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "ForceManager::initialize : useSurfaceTerm was set true\n"
                   "\tbut no electrostatic atoms are present. OpenMD will\n"
                   "\tignore this setting.\n");
          painCave.isFatal  = 0;
          painCave.severity = OPENMD_WARNING;
          simError();
        }
      }
    }

    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();

    //! Force fields can set options on how to scale van der Waals and
    //! electrostatic interactions for atoms connected via bonds, bends
    //! and torsions in this case the topological distance between
    //! atoms is:
    //! 0 = topologically unconnected
    //! 1 = bonded together
    //! 2 = connected via a bend
    //! 3 = connected via a torsion

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

    // Initialize the perturbations
    if (info_->getSimParams()->haveUniformField()) {
      UniformField* eField = new UniformField(info_);
      forceModifiers_.push_back(eField);
    }

    if (info_->getSimParams()->haveMagneticField()) {
      MagneticField* mField = new MagneticField(info_);
      forceModifiers_.push_back(mField);
    }

    if (info_->getSimParams()->haveUniformGradientStrength() ||
        info_->getSimParams()->haveUniformGradientDirection1() ||
        info_->getSimParams()->haveUniformGradientDirection2()) {
      UniformGradient* eGrad = new UniformGradient(info_);
      forceModifiers_.push_back(eGrad);
    }

    if (info_->getSimParams()->getLightParameters()->getUseLight()) {
      Perturbations::Light* light = new Perturbations::Light(info_);
      forceModifiers_.push_back(light);
    }

    // Initialize the force modifiers (order matters)
    if (info_->getSimParams()->getUseThermodynamicIntegration()) {
      ThermoIntegrationForceModifier* thermoInt =
          new ThermoIntegrationForceModifier(info_);
      forceModifiers_.push_back(thermoInt);
    } else if (info_->getSimParams()->getUseRestraints()) {
      RestraintForceModifier* restraint = new RestraintForceModifier(info_);
      forceModifiers_.push_back(restraint);
    }

    std::string ensembleParam =
        toUpperCopy(info_->getSimParams()->getEnsemble());

    if (ensembleParam == "LHULL" || ensembleParam == "LANGEVINHULL" ||
        ensembleParam == "SMIPD") {
#if defined(HAVE_QHULL)
      LangevinHullForceModifier* langevinHullFM =
          new LangevinHullForceModifier(info_);
      forceModifiers_.push_back(langevinHullFM);
#else
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "ForceManager: Cannot use the LangevinHull ensembles without qhull.\n"
          "\tPlease rebuild OpenMD with qhull enabled.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
#endif
    } else if (ensembleParam == "LANGEVINDYNAMICS" || ensembleParam == "LD") {
      LDForceModifier* langevinDynamicsFM = new LDForceModifier(info_);
      forceModifiers_.push_back(langevinDynamicsFM);
    }

    if (info_->getSimParams()->getNZconsStamps() > 0) {
      info_->setNZconstraint(info_->getSimParams()->getNZconsStamps());
      ZConstraintForceModifier* zCons = new ZConstraintForceModifier(info_);
      forceModifiers_.push_back(zCons);
    }

    usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    fDecomp_->distributeInitialData();

    doPotentialSelection_ = false;
    if (info_->getSimParams()->havePotentialSelection()) {
      doPotentialSelection_ = true;
      selectionScript_      = info_->getSimParams()->getPotentialSelection();
      evaluator_.loadScriptString(selectionScript_);
      if (!evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }
    }

    initialized_ = true;
  }

  void ForceManager::calcForces() {
    if (!initialized_) initialize();

    preCalculation();
    shortRangeInteractions();
    if (info_->getSimParams()->getSkipPairLoop() == false)
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

    fDecomp_->setSnapshot(snap);

    snap->setBondPotential(0.0);
    snap->setBendPotential(0.0);
    snap->setTorsionPotential(0.0);
    snap->setInversionPotential(0.0);

    potVec zeroPot(0.0);
    snap->setLongRangePotentials(zeroPot);

    snap->setExcludedPotentials(zeroPot);
    if (doPotentialSelection_) snap->setSelectionPotentials(zeroPot);

    snap->setSelfPotentials(0.0);
    snap->setRestraintPotential(0.0);
    snap->setRawPotential(0.0);

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        atom->zeroForcesAndTorques();
      }

      // change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        rb->zeroForcesAndTorques();
      }

      if (info_->getNGlobalCutoffGroups() != info_->getNGlobalAtoms()) {
        for (cg = mol->beginCutoffGroup(ci); cg != NULL;
             cg = mol->nextCutoffGroup(ci)) {
          // calculate the center of mass of cutoff group
          cg->updateCOM();
        }
      }
    }

    // Zero out the virial tensor
    virialTensor *= 0.0;
    // Zero out the heatFlux
    fDecomp_->setHeatFlux(Vector3d(0.0));

    if (doPotentialSelection_) {
      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }
    }
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
    Molecule::BondIterator bondIter;

    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;
    RealType bondPotential      = 0.0;
    RealType bendPotential      = 0.0;
    RealType torsionPotential   = 0.0;
    RealType inversionPotential = 0.0;
    potVec selectionPotential(0.0);

    // calculate short range interactions
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      // change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        rb->updateAtoms();
      }

      for (bond = mol->beginBond(bondIter); bond != NULL;
           bond = mol->nextBond(bondIter)) {
        bond->calcForce(doParticlePot_);
        bondPotential += bond->getPotential();
        if (doPotentialSelection_) {
          if (seleMan_.isSelected(bond->getAtomA()) ||
              seleMan_.isSelected(bond->getAtomB())) {
            selectionPotential[BONDED_FAMILY] += bond->getPotential();
          }
        }
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
          dataSet.deltaV                                  = 0.0;
          bendDataSets.insert(
              map<Bend*, BendDataSet>::value_type(bend, dataSet));
        } else {
          i->second.prev.angle     = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle     = angle;
          i->second.curr.potential = currBendPot;
          i->second.deltaV =
              fabs(i->second.curr.potential - i->second.prev.potential);
        }
        if (doPotentialSelection_) {
          if (seleMan_.isSelected(bend->getAtomA()) ||
              seleMan_.isSelected(bend->getAtomB()) ||
              seleMan_.isSelected(bend->getAtomC())) {
            selectionPotential[BONDED_FAMILY] += bend->getPotential();
          }
        }
      }

      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL;
           torsion = mol->nextTorsion(torsionIter)) {
        RealType angle;
        torsion->calcForce(angle, doParticlePot_);
        RealType currTorsionPot = torsion->getPotential();
        torsionPotential += torsion->getPotential();
        map<Torsion*, TorsionDataSet>::iterator i =
            torsionDataSets.find(torsion);
        if (i == torsionDataSets.end()) {
          TorsionDataSet dataSet;
          dataSet.prev.angle = dataSet.curr.angle = angle;
          dataSet.prev.potential = dataSet.curr.potential = currTorsionPot;
          dataSet.deltaV                                  = 0.0;
          torsionDataSets.insert(
              map<Torsion*, TorsionDataSet>::value_type(torsion, dataSet));
        } else {
          i->second.prev.angle     = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle     = angle;
          i->second.curr.potential = currTorsionPot;
          i->second.deltaV =
              fabs(i->second.curr.potential - i->second.prev.potential);
        }
        if (doPotentialSelection_) {
          if (seleMan_.isSelected(torsion->getAtomA()) ||
              seleMan_.isSelected(torsion->getAtomB()) ||
              seleMan_.isSelected(torsion->getAtomC()) ||
              seleMan_.isSelected(torsion->getAtomD())) {
            selectionPotential[BONDED_FAMILY] += torsion->getPotential();
          }
        }
      }

      for (inversion = mol->beginInversion(inversionIter); inversion != NULL;
           inversion = mol->nextInversion(inversionIter)) {
        RealType angle;
        inversion->calcForce(angle, doParticlePot_);
        RealType currInversionPot = inversion->getPotential();
        inversionPotential += inversion->getPotential();
        map<Inversion*, InversionDataSet>::iterator i =
            inversionDataSets.find(inversion);
        if (i == inversionDataSets.end()) {
          InversionDataSet dataSet;
          dataSet.prev.angle = dataSet.curr.angle = angle;
          dataSet.prev.potential = dataSet.curr.potential = currInversionPot;
          dataSet.deltaV                                  = 0.0;
          inversionDataSets.insert(
              map<Inversion*, InversionDataSet>::value_type(inversion,
                                                            dataSet));
        } else {
          i->second.prev.angle     = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle     = angle;
          i->second.curr.potential = currInversionPot;
          i->second.deltaV =
              fabs(i->second.curr.potential - i->second.prev.potential);
        }
        if (doPotentialSelection_) {
          if (seleMan_.isSelected(inversion->getAtomA()) ||
              seleMan_.isSelected(inversion->getAtomB()) ||
              seleMan_.isSelected(inversion->getAtomC()) ||
              seleMan_.isSelected(inversion->getAtomD())) {
            selectionPotential[BONDED_FAMILY] += inversion->getPotential();
          }
        }
      }
    }

#ifdef IS_MPI
    // Collect from all nodes.  This should eventually be moved into a
    // SystemDecomposition, but this is a better place than in
    // Thermo to do the collection.

    MPI_Allreduce(MPI_IN_PLACE, &bondPotential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &bendPotential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &torsionPotential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &inversionPotential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &selectionPotential[BONDED_FAMILY], 1,
                  MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    fDecomp_->setSnapshot(curSnapshot);

    curSnapshot->setBondPotential(bondPotential);
    curSnapshot->setBendPotential(bendPotential);
    curSnapshot->setTorsionPotential(torsionPotential);
    curSnapshot->setInversionPotential(inversionPotential);
    curSnapshot->setSelectionPotentials(selectionPotential);

    // RealType shortRangePotential = bondPotential + bendPotential +
    //   torsionPotential +  inversionPotential;

    // curSnapshot->setShortRangePotential(shortRangePotential);
  }

  void ForceManager::longRangeInteractions() {
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    fDecomp_->setSnapshot(curSnapshot);

    DataStorage* config   = &(curSnapshot->atomData);
    DataStorage* cgConfig = &(curSnapshot->cgData);

    // calculate the center of mass of cutoff group

    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;

    if (info_->getNCutoffGroups() != info_->getNAtoms()) {
      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (cg = mol->beginCutoffGroup(ci); cg != NULL;
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

    int cg1, cg2, atom1, atom2;
    Vector3d d_grp, dag, gvel2, vel2;
    RealType rgrpsq, rgrp;
    RealType vij(0.0);
    Vector3d fij, fg, f1;
    bool in_switching_region;
    RealType dswdr, swderiv;
    vector<int> atomListColumn, atomListRow;

    potVec longRangePotential(0.0);
    potVec selfPotential(0.0);
    RealType reciprocalPotential(0.0);
    RealType surfacePotential(0.0);
    potVec selectionPotential(0.0);

    RealType mf;
    bool newAtom1;
    int gid1, gid2;

    vector<int>::iterator ia, jb;

    int loopStart, loopEnd;

    idat.rcut       = rCut_;
    idat.shiftedPot = (cutoffMethod_ == SHIFTED_POTENTIAL) ? true : false;
    idat.shiftedForce =
        (cutoffMethod_ == SHIFTED_FORCE || cutoffMethod_ == TAYLOR_SHIFTED) ?
            true :
            false;
    idat.doParticlePot   = doParticlePot_;
    idat.doElectricField = doElectricField_;
    idat.doSitePotential = doSitePotential_;
    sdat.doParticlePot   = doParticlePot_;

    loopEnd = PAIR_LOOP;
    if (info_->requiresPrepair()) {
      loopStart = PREPAIR_LOOP;
    } else {
      loopStart = PAIR_LOOP;
    }
    for (int iLoop = loopStart; iLoop <= loopEnd; iLoop++) {
      if (iLoop == loopStart) {
        bool update_nlist = fDecomp_->checkNeighborList(savedPositions_);

        if (update_nlist) {
          if (!usePeriodicBoundaryConditions_)
            Mat3x3d bbox = thermo->getBoundingBox();
          fDecomp_->buildNeighborList(neighborList_, point_, savedPositions_);
        }
      }

      for (cg1 = 0; cg1 < int(point_.size()) - 1; cg1++) {
        atomListRow = fDecomp_->getAtomsInGroupRow(cg1);
        newAtom1    = true;

        for (int m2 = point_[cg1]; m2 < point_[cg1 + 1]; m2++) {
          cg2 = neighborList_[m2];

          d_grp = fDecomp_->getIntergroupVector(cg1, cg2);

          // already wrapped in the getIntergroupVector call:
          // curSnapshot->wrapVector(d_grp);
          rgrpsq = d_grp.lengthSquare();

          if (rgrpsq < rCutSq_) {
            if (iLoop == PAIR_LOOP) {
              vij = 0.0;
              fij.zero();
              idat.eField1.zero();
              idat.eField2.zero();
              idat.sPot1 = 0.0;
              idat.sPot2 = 0.0;
            }

            in_switching_region =
                switcher_->getSwitch(rgrpsq, idat.sw, dswdr, rgrp);

            atomListColumn = fDecomp_->getAtomsInGroupColumn(cg2);

            if (doHeatFlux_) gvel2 = fDecomp_->getGroupVelocityColumn(cg2);

            for (ia = atomListRow.begin(); ia != atomListRow.end(); ++ia) {
              atom1 = (*ia);
              if (atomListRow.size() > 1) newAtom1 = true;

              if (doPotentialSelection_) {
                gid1            = fDecomp_->getGlobalIDRow(atom1);
                idat.isSelected = seleMan_.isGlobalIDSelected(gid1);
              }

              for (jb = atomListColumn.begin(); jb != atomListColumn.end();
                   ++jb) {
                atom2 = (*jb);

                if (doPotentialSelection_) {
                  gid2 = fDecomp_->getGlobalIDCol(atom2);
                  idat.isSelected |= seleMan_.isGlobalIDSelected(gid2);
                }

                if (!fDecomp_->skipAtomPair(atom1, atom2, cg1, cg2)) {
                  idat.vpair       = 0.0;
                  idat.pot         = 0.0;
                  idat.excludedPot = 0.0;
                  idat.selePot     = 0.0;
                  idat.f1.zero();
                  idat.dVdFQ1 = 0.0;
                  idat.dVdFQ2 = 0.0;

                  fDecomp_->fillInteractionData(idat, atom1, atom2, newAtom1);
                  newAtom1 = false;

                  idat.topoDist =
                      fDecomp_->getTopologicalDistance(atom1, atom2);
                  idat.vdwMult     = vdwScale_[idat.topoDist];
                  idat.electroMult = electrostaticScale_[idat.topoDist];

                  if (atomListRow.size() == 1 && atomListColumn.size() == 1) {
                    idat.d  = d_grp;
                    idat.r2 = rgrpsq;
                    if (doHeatFlux_) vel2 = gvel2;
                  } else {
                    idat.d = fDecomp_->getInteratomicVector(atom1, atom2);
                    curSnapshot->wrapVector(idat.d);
                    idat.r2 = idat.d.lengthSquare();
                    if (doHeatFlux_)
                      vel2 = fDecomp_->getAtomVelocityColumn(atom2);
                  }

                  idat.rij = sqrt(idat.r2);

                  if (iLoop == PREPAIR_LOOP) {
                    interactionMan_->doPrePair(idat);
                    fDecomp_->unpackPrePairData(idat, atom1, atom2);
                  } else {
                    interactionMan_->doPair(idat);
                    fDecomp_->unpackInteractionData(idat, atom1, atom2);
                    vij += idat.vpair;
                    fij += idat.f1;
                    virialTensor -= outProduct(idat.d, idat.f1);
                    if (doHeatFlux_)
                      fDecomp_->addToHeatFlux(idat.d * dot(idat.f1, vel2));
                  }
                }
              }
            }

            if (iLoop == PAIR_LOOP) {
              if (in_switching_region) {
                swderiv = vij * dswdr / rgrp;
                fg      = swderiv * d_grp;
                fij += fg;

                if (atomListRow.size() == 1 && atomListColumn.size() == 1) {
                  if (!fDecomp_->skipAtomPair(atomListRow[0], atomListColumn[0],
                                              cg1, cg2)) {
                    virialTensor -= outProduct(idat.d, fg);
                    if (doHeatFlux_)
                      fDecomp_->addToHeatFlux(idat.d * dot(fg, vel2));
                  }
                }

                for (ia = atomListRow.begin(); ia != atomListRow.end(); ++ia) {
                  atom1 = (*ia);
                  mf    = fDecomp_->getMassFactorRow(atom1);
                  // fg is the force on atom ia due to cutoff group's
                  // presence in switching region
                  fg = swderiv * d_grp * mf;
                  fDecomp_->addForceToAtomRow(atom1, fg);
                  if (atomListRow.size() > 1) {
                    if (info_->usesAtomicVirial()) {
                      // find the distance between the atom
                      // and the center of the cutoff group:
                      dag = fDecomp_->getAtomToGroupVectorRow(atom1, cg1);
                      virialTensor -= outProduct(dag, fg);
                      if (doHeatFlux_)
                        fDecomp_->addToHeatFlux(dag * dot(fg, vel2));
                    }
                  }
                }
                for (jb = atomListColumn.begin(); jb != atomListColumn.end();
                     ++jb) {
                  atom2 = (*jb);
                  mf    = fDecomp_->getMassFactorColumn(atom2);
                  // fg is the force on atom jb due to cutoff group's
                  // presence in switching region
                  fg = -swderiv * d_grp * mf;
                  fDecomp_->addForceToAtomColumn(atom2, fg);

                  if (atomListColumn.size() > 1) {
                    if (info_->usesAtomicVirial()) {
                      // find the distance between the atom
                      // and the center of the cutoff group:
                      dag = fDecomp_->getAtomToGroupVectorColumn(atom2, cg2);
                      virialTensor -= outProduct(dag, fg);
                      if (doHeatFlux_)
                        fDecomp_->addToHeatFlux(dag * dot(fg, vel2));
                    }
                  }
                }
              }
              // if (!info_->usesAtomicVirial()) {
              //  virialTensor -= outProduct(d_grp, fij);
              //  if (doHeatFlux_)
              //     fDecomp_->addToHeatFlux( d_grp * dot(fij, vel2));
              //}
            }
          }
        }
      }

      if (iLoop == PREPAIR_LOOP) {
        if (info_->requiresPrepair()) {
          fDecomp_->collectIntermediateData();

          for (unsigned int atom1 = 0; atom1 < info_->getNAtoms(); atom1++) {
            if (doPotentialSelection_) {
              gid1            = fDecomp_->getGlobalID(atom1);
              sdat.isSelected = seleMan_.isGlobalIDSelected(gid1);
            }
            fDecomp_->fillPreForceData(sdat, atom1);
            interactionMan_->doPreForce(sdat);
            fDecomp_->unpackPreForceData(sdat, atom1);
          }

          fDecomp_->distributeIntermediateData();
        }
      }
    }

    // collects pairwise information
    fDecomp_->collectData();
    if (cutoffMethod_ == EWALD_FULL) {
      interactionMan_->doReciprocalSpaceSum(reciprocalPotential);
      curSnapshot->setReciprocalPotential(reciprocalPotential);
    }

    if (useSurfaceTerm_) {
      interactionMan_->doSurfaceTerm(useSlabGeometry_, axis_, surfacePotential);
      curSnapshot->setSurfacePotential(surfacePotential);
    }

    if (info_->requiresSelfCorrection()) {
      for (unsigned int atom1 = 0; atom1 < info_->getNAtoms(); atom1++) {
        if (doPotentialSelection_) {
          gid1            = fDecomp_->getGlobalID(atom1);
          sdat.isSelected = seleMan_.isGlobalIDSelected(gid1);
        }
        fDecomp_->fillSelfData(sdat, atom1);
        interactionMan_->doSelfCorrection(sdat);
        fDecomp_->unpackSelfData(sdat, atom1);
      }
    }

    // collects single-atom information
    fDecomp_->collectSelfData();

    longRangePotential = fDecomp_->getPairwisePotential();
    curSnapshot->setLongRangePotentials(longRangePotential);

    selfPotential = fDecomp_->getSelfPotential();
    curSnapshot->setSelfPotentials(selfPotential);

    curSnapshot->setExcludedPotentials(fDecomp_->getExcludedSelfPotential() +
                                       fDecomp_->getExcludedPotential());

    if (doPotentialSelection_) {
      selectionPotential = curSnapshot->getSelectionPotentials();
      selectionPotential += fDecomp_->getSelectedSelfPotential();
      selectionPotential += fDecomp_->getSelectedPotential();
      curSnapshot->setSelectionPotentials(selectionPotential);
    }
  }

  void ForceManager::postCalculation() {
    for (auto& forceModifier : forceModifiers_)
      forceModifier->modifyForces();

    // Modify the rigid bodies in response to the applied force modifications
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();

    // Collect the atomic forces onto rigid bodies
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        Mat3x3d rbTau = rb->calcForcesAndTorquesAndVirial();
        virialTensor += rbTau;
      }
    }

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, virialTensor.getArrayPointer(), 9, MPI_REALTYPE,
                  MPI_SUM, MPI_COMM_WORLD);
#endif
    curSnapshot->setVirialTensor(virialTensor);

    /*
      if (info_->getSimParams()->getUseLongRangeCorrections()) {
        RealType vol = curSnapshot->getVolume();
        RealType Elrc(0.0);
        RealType Wlrc(0.0);

        AtomTypeSet::iterator i;
        AtomTypeSet::iterator j;

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
        Elrc *= 2.0 * Constants::PI;
        Wlrc *= 2.0 * Constants::PI;

        RealType lrp = curSnapshot->getLongRangePotential();
        curSnapshot->setLongRangePotential(lrp + Elrc);
        virialTensor += Wlrc * SquareMatrix3<RealType>::identity();
        curSnapshot->setVirialTensor(virialTensor);
      }
    */

    // Sanity check on PE - don't wait for StatWriter or DumpWriter.
    RealType pe = curSnapshot->getPotentialEnergy();
    if (std::isinf(pe) || std::isnan(pe)) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "ForceManager detected a numerical error in the potential energy.\n"
          "\tStopping the simulation.");
      painCave.isFatal = 1;
      simError();
    }
    Vector3d f;
    Vector3d t;

    Molecule::IntegrableObjectIterator ii;
    StuntDouble* sd;

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ii); sd != NULL;
           sd = mol->nextIntegrableObject(ii)) {
        f = sd->getFrc();
        if (std::isinf(f[0]) || std::isnan(f[0]) || std::isinf(f[1]) ||
            std::isnan(f[1]) || std::isinf(f[2]) || std::isnan(f[2])) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "ForceManager detected a numerical error in the force"
                   " for object %d",
                   sd->getGlobalIndex());
          painCave.isFatal = 1;
          simError();
        }
        if (sd->isDirectional()) {
          t = sd->getTrq();
          if (std::isinf(t[0]) || std::isnan(t[0]) || std::isinf(t[1]) ||
              std::isnan(t[1]) || std::isinf(t[2]) || std::isnan(t[2])) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "ForceManager detected a numerical error in the torque"
                     " for object %d",
                     sd->getGlobalIndex());
            painCave.isFatal = 1;
            simError();
          }
        }
      }
    }
    errorCheckPoint();
  }

  void ForceManager::calcSelectedForces(Molecule* mol1, Molecule* mol2) {
    if (!initialized_) initialize();

    selectedPreCalculation(mol1, mol2);
    selectedShortRangeInteractions(mol1, mol2);
    selectedLongRangeInteractions(mol1, mol2);
    selectedPostCalculation(mol1, mol2);
  }

  void ForceManager::selectedPreCalculation(Molecule* mol1, Molecule* mol2) {
    SimInfo::MoleculeIterator mi;
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
    snap->setLongRangePotentials(zeroPot);
    snap->setExcludedPotentials(zeroPot);
    if (doPotentialSelection_) snap->setSelectionPotentials(zeroPot);

    snap->setRestraintPotential(0.0);
    snap->setRawPotential(0.0);

    // First we zero out for mol1
    for (atom = mol1->beginAtom(ai); atom != NULL; atom = mol1->nextAtom(ai)) {
      atom->zeroForcesAndTorques();
    }
    // change the positions of atoms which belong to the rigidbodies
    for (rb = mol1->beginRigidBody(rbIter); rb != NULL;
         rb = mol1->nextRigidBody(rbIter)) {
      rb->zeroForcesAndTorques();
    }
    if (info_->getNGlobalCutoffGroups() != info_->getNGlobalAtoms()) {
      for (cg = mol1->beginCutoffGroup(ci); cg != NULL;
           cg = mol1->nextCutoffGroup(ci)) {
        // calculate the center of mass of cutoff group
        cg->updateCOM();
      }
    }

    // Next we zero out for mol2
    for (atom = mol2->beginAtom(ai); atom != NULL; atom = mol2->nextAtom(ai)) {
      atom->zeroForcesAndTorques();
    }
    // change the positions of atoms which belong to the rigidbodies
    for (rb = mol2->beginRigidBody(rbIter); rb != NULL;
         rb = mol2->nextRigidBody(rbIter)) {
      rb->zeroForcesAndTorques();
    }
    if (info_->getNGlobalCutoffGroups() != info_->getNGlobalAtoms()) {
      for (cg = mol2->beginCutoffGroup(ci); cg != NULL;
           cg = mol2->nextCutoffGroup(ci)) {
        // calculate the center of mass of cutoff group
        cg->updateCOM();
      }
    }

    // Zero out the virial tensor
    virialTensor *= 0.0;
    // Zero out the heatFlux
    fDecomp_->setHeatFlux(Vector3d(0.0));
  }

  void ForceManager::selectedShortRangeInteractions(Molecule* mol1,
                                                    Molecule* mol2) {
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;
    RealType bondPotential      = 0.0;
    RealType bendPotential      = 0.0;
    RealType torsionPotential   = 0.0;
    RealType inversionPotential = 0.0;
    potVec selectionPotential(0.0);

    // First compute for mol1
    for (rb = mol1->beginRigidBody(rbIter); rb != NULL;
         rb = mol1->nextRigidBody(rbIter)) {
      rb->updateAtoms();
    }

    for (bond = mol1->beginBond(bondIter); bond != NULL;
         bond = mol1->nextBond(bondIter)) {
      bond->calcForce(doParticlePot_);
      bondPotential += bond->getPotential();
      if (doPotentialSelection_) {
        if (seleMan_.isSelected(bond->getAtomA()) ||
            seleMan_.isSelected(bond->getAtomB())) {
          selectionPotential[BONDED_FAMILY] += bond->getPotential();
        }
      }
    }

    for (bend = mol1->beginBend(bendIter); bend != NULL;
         bend = mol1->nextBend(bendIter)) {
      RealType angle;
      bend->calcForce(angle, doParticlePot_);
      RealType currBendPot = bend->getPotential();

      bendPotential += bend->getPotential();
      map<Bend*, BendDataSet>::iterator i = bendDataSets.find(bend);
      if (i == bendDataSets.end()) {
        BendDataSet dataSet;
        dataSet.prev.angle = dataSet.curr.angle = angle;
        dataSet.prev.potential = dataSet.curr.potential = currBendPot;
        dataSet.deltaV                                  = 0.0;
        bendDataSets.insert(map<Bend*, BendDataSet>::value_type(bend, dataSet));
      } else {
        i->second.prev.angle     = i->second.curr.angle;
        i->second.prev.potential = i->second.curr.potential;
        i->second.curr.angle     = angle;
        i->second.curr.potential = currBendPot;
        i->second.deltaV =
            fabs(i->second.curr.potential - i->second.prev.potential);
      }
      if (doPotentialSelection_) {
        if (seleMan_.isSelected(bend->getAtomA()) ||
            seleMan_.isSelected(bend->getAtomB()) ||
            seleMan_.isSelected(bend->getAtomC())) {
          selectionPotential[BONDED_FAMILY] += bend->getPotential();
        }
      }
    }

    for (torsion = mol1->beginTorsion(torsionIter); torsion != NULL;
         torsion = mol1->nextTorsion(torsionIter)) {
      RealType angle;
      torsion->calcForce(angle, doParticlePot_);
      RealType currTorsionPot = torsion->getPotential();
      torsionPotential += torsion->getPotential();
      map<Torsion*, TorsionDataSet>::iterator i = torsionDataSets.find(torsion);
      if (i == torsionDataSets.end()) {
        TorsionDataSet dataSet;
        dataSet.prev.angle = dataSet.curr.angle = angle;
        dataSet.prev.potential = dataSet.curr.potential = currTorsionPot;
        dataSet.deltaV                                  = 0.0;
        torsionDataSets.insert(
            map<Torsion*, TorsionDataSet>::value_type(torsion, dataSet));
      } else {
        i->second.prev.angle     = i->second.curr.angle;
        i->second.prev.potential = i->second.curr.potential;
        i->second.curr.angle     = angle;
        i->second.curr.potential = currTorsionPot;
        i->second.deltaV =
            fabs(i->second.curr.potential - i->second.prev.potential);
      }
      if (doPotentialSelection_) {
        if (seleMan_.isSelected(torsion->getAtomA()) ||
            seleMan_.isSelected(torsion->getAtomB()) ||
            seleMan_.isSelected(torsion->getAtomC()) ||
            seleMan_.isSelected(torsion->getAtomD())) {
          selectionPotential[BONDED_FAMILY] += torsion->getPotential();
        }
      }
    }

    for (inversion = mol1->beginInversion(inversionIter); inversion != NULL;
         inversion = mol1->nextInversion(inversionIter)) {
      RealType angle;
      inversion->calcForce(angle, doParticlePot_);
      RealType currInversionPot = inversion->getPotential();
      inversionPotential += inversion->getPotential();
      map<Inversion*, InversionDataSet>::iterator i =
          inversionDataSets.find(inversion);
      if (i == inversionDataSets.end()) {
        InversionDataSet dataSet;
        dataSet.prev.angle = dataSet.curr.angle = angle;
        dataSet.prev.potential = dataSet.curr.potential = currInversionPot;
        dataSet.deltaV                                  = 0.0;
        inversionDataSets.insert(
            map<Inversion*, InversionDataSet>::value_type(inversion, dataSet));
      } else {
        i->second.prev.angle     = i->second.curr.angle;
        i->second.prev.potential = i->second.curr.potential;
        i->second.curr.angle     = angle;
        i->second.curr.potential = currInversionPot;
        i->second.deltaV =
            fabs(i->second.curr.potential - i->second.prev.potential);
      }
      if (doPotentialSelection_) {
        if (seleMan_.isSelected(inversion->getAtomA()) ||
            seleMan_.isSelected(inversion->getAtomB()) ||
            seleMan_.isSelected(inversion->getAtomC()) ||
            seleMan_.isSelected(inversion->getAtomD())) {
          selectionPotential[BONDED_FAMILY] += inversion->getPotential();
        }
      }
    }

    // Next compute for mol2
    for (rb = mol2->beginRigidBody(rbIter); rb != NULL;
         rb = mol2->nextRigidBody(rbIter)) {
      rb->updateAtoms();
    }

    for (bond = mol2->beginBond(bondIter); bond != NULL;
         bond = mol2->nextBond(bondIter)) {
      bond->calcForce(doParticlePot_);
      bondPotential += bond->getPotential();
      if (doPotentialSelection_) {
        if (seleMan_.isSelected(bond->getAtomA()) ||
            seleMan_.isSelected(bond->getAtomB())) {
          selectionPotential[BONDED_FAMILY] += bond->getPotential();
        }
      }
    }

    for (bend = mol2->beginBend(bendIter); bend != NULL;
         bend = mol2->nextBend(bendIter)) {
      RealType angle;
      bend->calcForce(angle, doParticlePot_);
      RealType currBendPot = bend->getPotential();

      bendPotential += bend->getPotential();
      map<Bend*, BendDataSet>::iterator i = bendDataSets.find(bend);
      if (i == bendDataSets.end()) {
        BendDataSet dataSet;
        dataSet.prev.angle = dataSet.curr.angle = angle;
        dataSet.prev.potential = dataSet.curr.potential = currBendPot;
        dataSet.deltaV                                  = 0.0;
        bendDataSets.insert(map<Bend*, BendDataSet>::value_type(bend, dataSet));
      } else {
        i->second.prev.angle     = i->second.curr.angle;
        i->second.prev.potential = i->second.curr.potential;
        i->second.curr.angle     = angle;
        i->second.curr.potential = currBendPot;
        i->second.deltaV =
            fabs(i->second.curr.potential - i->second.prev.potential);
      }
      if (doPotentialSelection_) {
        if (seleMan_.isSelected(bend->getAtomA()) ||
            seleMan_.isSelected(bend->getAtomB()) ||
            seleMan_.isSelected(bend->getAtomC())) {
          selectionPotential[BONDED_FAMILY] += bend->getPotential();
        }
      }
    }

    for (torsion = mol2->beginTorsion(torsionIter); torsion != NULL;
         torsion = mol2->nextTorsion(torsionIter)) {
      RealType angle;
      torsion->calcForce(angle, doParticlePot_);
      RealType currTorsionPot = torsion->getPotential();
      torsionPotential += torsion->getPotential();
      map<Torsion*, TorsionDataSet>::iterator i = torsionDataSets.find(torsion);
      if (i == torsionDataSets.end()) {
        TorsionDataSet dataSet;
        dataSet.prev.angle = dataSet.curr.angle = angle;
        dataSet.prev.potential = dataSet.curr.potential = currTorsionPot;
        dataSet.deltaV                                  = 0.0;
        torsionDataSets.insert(
            map<Torsion*, TorsionDataSet>::value_type(torsion, dataSet));
      } else {
        i->second.prev.angle     = i->second.curr.angle;
        i->second.prev.potential = i->second.curr.potential;
        i->second.curr.angle     = angle;
        i->second.curr.potential = currTorsionPot;
        i->second.deltaV =
            fabs(i->second.curr.potential - i->second.prev.potential);
      }
      if (doPotentialSelection_) {
        if (seleMan_.isSelected(torsion->getAtomA()) ||
            seleMan_.isSelected(torsion->getAtomB()) ||
            seleMan_.isSelected(torsion->getAtomC()) ||
            seleMan_.isSelected(torsion->getAtomD())) {
          selectionPotential[BONDED_FAMILY] += torsion->getPotential();
        }
      }
    }

    for (inversion = mol2->beginInversion(inversionIter); inversion != NULL;
         inversion = mol2->nextInversion(inversionIter)) {
      RealType angle;
      inversion->calcForce(angle, doParticlePot_);
      RealType currInversionPot = inversion->getPotential();
      inversionPotential += inversion->getPotential();
      map<Inversion*, InversionDataSet>::iterator i =
          inversionDataSets.find(inversion);
      if (i == inversionDataSets.end()) {
        InversionDataSet dataSet;
        dataSet.prev.angle = dataSet.curr.angle = angle;
        dataSet.prev.potential = dataSet.curr.potential = currInversionPot;
        dataSet.deltaV                                  = 0.0;
        inversionDataSets.insert(
            map<Inversion*, InversionDataSet>::value_type(inversion, dataSet));
      } else {
        i->second.prev.angle     = i->second.curr.angle;
        i->second.prev.potential = i->second.curr.potential;
        i->second.curr.angle     = angle;
        i->second.curr.potential = currInversionPot;
        i->second.deltaV =
            fabs(i->second.curr.potential - i->second.prev.potential);
      }
      if (doPotentialSelection_) {
        if (seleMan_.isSelected(inversion->getAtomA()) ||
            seleMan_.isSelected(inversion->getAtomB()) ||
            seleMan_.isSelected(inversion->getAtomC()) ||
            seleMan_.isSelected(inversion->getAtomD())) {
          selectionPotential[BONDED_FAMILY] += inversion->getPotential();
        }
      }
    }

#ifdef IS_MPI
    // Collect from all nodes.  This should eventually be moved into a
    // SystemDecomposition, but this is a better place than in
    // Thermo to do the collection.

    MPI_Allreduce(MPI_IN_PLACE, &bondPotential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &bendPotential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &torsionPotential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &inversionPotential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &selectionPotential[BONDED_FAMILY], 1,
                  MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();

    curSnapshot->setBondPotential(bondPotential);
    curSnapshot->setBendPotential(bendPotential);
    curSnapshot->setTorsionPotential(torsionPotential);
    curSnapshot->setInversionPotential(inversionPotential);
    curSnapshot->setSelectionPotentials(selectionPotential);

    // RealType shortRangePotential = bondPotential + bendPotential +
    //   torsionPotential +  inversionPotential;

    // curSnapshot->setShortRangePotential(shortRangePotential);
  }

  void ForceManager::selectedLongRangeInteractions(Molecule* mol1,
                                                   Molecule* mol2) {
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    DataStorage* config   = &(curSnapshot->atomData);
    DataStorage* cgConfig = &(curSnapshot->cgData);

    // calculate the center of mass of cutoff group

    SimInfo::MoleculeIterator mi;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;

    if (info_->getNCutoffGroups() != info_->getNAtoms()) {
      for (cg = mol1->beginCutoffGroup(ci); cg != NULL;
           cg = mol1->nextCutoffGroup(ci)) {
        cg->updateCOM();
      }
      for (cg = mol2->beginCutoffGroup(ci); cg != NULL;
           cg = mol2->nextCutoffGroup(ci)) {
        cg->updateCOM();
      }
    } else {
      // center of mass of the group is the same as position of the atom
      // if cutoff group does not exist
      cgConfig->position = config->position;
      cgConfig->velocity = config->velocity;
    }

    fDecomp_->zeroWorkArrays();
    fDecomp_->distributeData();

    int cg1, cg2, atom1, atom2;
    Vector3d d_grp, dag, gvel2, vel2;
    RealType rgrpsq, rgrp;
    RealType vij(0.0);
    Vector3d fij, fg;
    bool in_switching_region;
    RealType dswdr, swderiv;
    vector<int> atomListColumn, atomListRow;

    potVec longRangePotential(0.0);
    potVec selfPotential(0.0);
    potVec selectionPotential(0.0);
    RealType reciprocalPotential(0.0);
    RealType surfacePotential(0.0);

    RealType mf;
    bool newAtom1;
    int gid1, gid2;

    vector<int>::iterator ia, jb;

    int loopStart, loopEnd;

    idat.rcut       = rCut_;
    idat.shiftedPot = (cutoffMethod_ == SHIFTED_POTENTIAL) ? true : false;
    idat.shiftedForce =
        (cutoffMethod_ == SHIFTED_FORCE || cutoffMethod_ == TAYLOR_SHIFTED) ?
            true :
            false;
    idat.doParticlePot   = doParticlePot_;
    idat.doElectricField = doElectricField_;
    idat.doSitePotential = doSitePotential_;
    sdat.doParticlePot   = doParticlePot_;

    loopEnd = PAIR_LOOP;
    if (info_->requiresPrepair()) {
      loopStart = PREPAIR_LOOP;
    } else {
      loopStart = PAIR_LOOP;
    }
    for (int iLoop = loopStart; iLoop <= loopEnd; iLoop++) {
      if (iLoop == loopStart) {
        bool update_nlist = fDecomp_->checkNeighborList(savedPositions_);

        if (update_nlist) {
          if (!usePeriodicBoundaryConditions_)
            Mat3x3d bbox = thermo->getBoundingBox();
          fDecomp_->buildNeighborList(neighborList_, point_, savedPositions_);
        }
      }

      for (cg1 = 0; cg1 < int(point_.size()) - 1; cg1++) {
        atomListRow = fDecomp_->getAtomsInGroupRow(cg1);
        newAtom1    = true;

        for (int m2 = point_[cg1]; m2 < point_[cg1 + 1]; m2++) {
          cg2 = neighborList_[m2];

          d_grp = fDecomp_->getIntergroupVector(cg1, cg2);

          // already wrapped in the getIntergroupVector call:
          // curSnapshot->wrapVector(d_grp);
          rgrpsq = d_grp.lengthSquare();

          if (rgrpsq < rCutSq_) {
            if (iLoop == PAIR_LOOP) {
              vij = 0.0;
              fij.zero();
              idat.eField1.zero();
              idat.eField2.zero();
              idat.sPot1 = 0.0;
              idat.sPot2 = 0.0;
            }

            in_switching_region =
                switcher_->getSwitch(rgrpsq, idat.sw, dswdr, rgrp);

            atomListColumn = fDecomp_->getAtomsInGroupColumn(cg2);

            if (doHeatFlux_) gvel2 = fDecomp_->getGroupVelocityColumn(cg2);

            for (ia = atomListRow.begin(); ia != atomListRow.end(); ++ia) {
              if (atomListRow.size() > 1) newAtom1 = true;
              atom1 = (*ia);

              if (doPotentialSelection_) {
                gid1            = fDecomp_->getGlobalIDRow(atom1);
                idat.isSelected = seleMan_.isGlobalIDSelected(gid1);
              }

              for (jb = atomListColumn.begin(); jb != atomListColumn.end();
                   ++jb) {
                atom2 = (*jb);

                if (doPotentialSelection_) {
                  gid2 = fDecomp_->getGlobalIDCol(atom2);
                  idat.isSelected |= seleMan_.isGlobalIDSelected(gid2);
                }

                if (!fDecomp_->skipAtomPair(atom1, atom2, cg1, cg2)) {
                  idat.vpair       = 0.0;
                  idat.pot         = 0.0;
                  idat.excludedPot = 0.0;
                  idat.selePot     = 0.0;
                  idat.f1.zero();
                  idat.dVdFQ1 = 0.0;
                  idat.dVdFQ2 = 0.0;

                  fDecomp_->fillInteractionData(idat, atom1, atom2, newAtom1);
                  newAtom1 = false;

                  idat.topoDist =
                      fDecomp_->getTopologicalDistance(atom1, atom2);
                  idat.vdwMult     = vdwScale_[idat.topoDist];
                  idat.electroMult = electrostaticScale_[idat.topoDist];

                  if (atomListRow.size() == 1 && atomListColumn.size() == 1) {
                    idat.d  = d_grp;
                    idat.r2 = rgrpsq;
                    if (doHeatFlux_) vel2 = gvel2;
                  } else {
                    idat.d = fDecomp_->getInteratomicVector(atom1, atom2);
                    curSnapshot->wrapVector(idat.d);
                    idat.r2 = idat.d.lengthSquare();
                    if (doHeatFlux_)
                      vel2 = fDecomp_->getAtomVelocityColumn(atom2);
                  }

                  idat.rij = sqrt(idat.r2);

                  if (iLoop == PREPAIR_LOOP) {
                    interactionMan_->doPrePair(idat);
                    fDecomp_->unpackPrePairData(idat, atom1, atom2);
                  } else {
                    interactionMan_->doPair(idat);
                    fDecomp_->unpackInteractionData(idat, atom1, atom2);
                    vij += idat.vpair;
                    fij += idat.f1;
                    virialTensor -= outProduct(idat.d, idat.f1);
                    if (doHeatFlux_)
                      fDecomp_->addToHeatFlux(idat.d * dot(idat.f1, vel2));
                  }
                }
              }
            }

            if (iLoop == PAIR_LOOP) {
              if (in_switching_region) {
                swderiv = vij * dswdr / rgrp;
                fg      = swderiv * d_grp;
                fij += fg;

                if (atomListRow.size() == 1 && atomListColumn.size() == 1) {
                  if (!fDecomp_->skipAtomPair(atomListRow[0], atomListColumn[0],
                                              cg1, cg2)) {
                    virialTensor -= outProduct(idat.d, fg);
                    if (doHeatFlux_)
                      fDecomp_->addToHeatFlux(idat.d * dot(fg, vel2));
                  }
                }

                for (ia = atomListRow.begin(); ia != atomListRow.end(); ++ia) {
                  atom1 = (*ia);
                  mf    = fDecomp_->getMassFactorRow(atom1);
                  // fg is the force on atom ia due to cutoff group's
                  // presence in switching region
                  fg = swderiv * d_grp * mf;
                  fDecomp_->addForceToAtomRow(atom1, fg);
                  if (atomListRow.size() > 1) {
                    if (info_->usesAtomicVirial()) {
                      // find the distance between the atom
                      // and the center of the cutoff group:
                      dag = fDecomp_->getAtomToGroupVectorRow(atom1, cg1);
                      virialTensor -= outProduct(dag, fg);
                      if (doHeatFlux_)
                        fDecomp_->addToHeatFlux(dag * dot(fg, vel2));
                    }
                  }
                }
                for (jb = atomListColumn.begin(); jb != atomListColumn.end();
                     ++jb) {
                  atom2 = (*jb);
                  mf    = fDecomp_->getMassFactorColumn(atom2);
                  // fg is the force on atom jb due to cutoff group's
                  // presence in switching region
                  fg = -swderiv * d_grp * mf;
                  fDecomp_->addForceToAtomColumn(atom2, fg);

                  if (atomListColumn.size() > 1) {
                    if (info_->usesAtomicVirial()) {
                      // find the distance between the atom
                      // and the center of the cutoff group:
                      dag = fDecomp_->getAtomToGroupVectorColumn(atom2, cg2);
                      virialTensor -= outProduct(dag, fg);
                      if (doHeatFlux_)
                        fDecomp_->addToHeatFlux(dag * dot(fg, vel2));
                    }
                  }
                }
              }
              // if (!info_->usesAtomicVirial()) {
              //  virialTensor -= outProduct(d_grp, fij);
              //  if (doHeatFlux_)
              //     fDecomp_->addToHeatFlux( d_grp * dot(fij, vel2));
              //}
            }
          }
        }
      }

      if (iLoop == PREPAIR_LOOP) {
        if (info_->requiresPrepair()) {
          fDecomp_->collectIntermediateData();

          for (unsigned int atom1 = 0; atom1 < info_->getNAtoms(); atom1++) {
            if (doPotentialSelection_) {
              gid1            = fDecomp_->getGlobalID(atom1);
              sdat.isSelected = seleMan_.isGlobalIDSelected(gid1);
            }

            fDecomp_->fillPreForceData(sdat, atom1);
            interactionMan_->doPreForce(sdat);
            fDecomp_->unpackPreForceData(sdat, atom1);
          }

          fDecomp_->distributeIntermediateData();
        }
      }
    }

    // collects pairwise information
    fDecomp_->collectData();
    if (cutoffMethod_ == EWALD_FULL) {
      interactionMan_->doReciprocalSpaceSum(reciprocalPotential);
      curSnapshot->setReciprocalPotential(reciprocalPotential);
    }

    if (useSurfaceTerm_) {
      interactionMan_->doSurfaceTerm(useSlabGeometry_, axis_, surfacePotential);
      curSnapshot->setSurfacePotential(surfacePotential);
    }

    if (info_->requiresSelfCorrection()) {
      for (unsigned int atom1 = 0; atom1 < info_->getNAtoms(); atom1++) {
        if (doPotentialSelection_) {
          gid1            = fDecomp_->getGlobalID(atom1);
          sdat.isSelected = seleMan_.isGlobalIDSelected(gid1);
        }

        fDecomp_->fillSelfData(sdat, atom1);
        interactionMan_->doSelfCorrection(sdat);
        fDecomp_->unpackSelfData(sdat, atom1);
      }
    }

    // collects single-atom information
    fDecomp_->collectSelfData();

    longRangePotential = fDecomp_->getPairwisePotential();
    curSnapshot->setLongRangePotentials(longRangePotential);

    selfPotential = fDecomp_->getSelfPotential();
    curSnapshot->setSelfPotentials(selfPotential);

    curSnapshot->setExcludedPotentials(fDecomp_->getExcludedSelfPotential() +
                                       fDecomp_->getExcludedPotential());

    if (doPotentialSelection_) {
      selectionPotential = curSnapshot->getSelectionPotentials();
      selectionPotential += fDecomp_->getSelectedSelfPotential();
      selectionPotential += fDecomp_->getSelectedPotential();
      curSnapshot->setSelectionPotentials(selectionPotential);
    }
  }

  void ForceManager::selectedPostCalculation(Molecule* mol1, Molecule* mol2) {
    for (auto& forceModifier : forceModifiers_)
      forceModifier->modifyForces();

    // Modify the rigid bodies in response to the applied force
    // modifications
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();

    // Collect the atomic forces onto rigid bodies
    for (rb = mol1->beginRigidBody(rbIter); rb != NULL;
         rb = mol1->nextRigidBody(rbIter)) {
      Mat3x3d rbTau = rb->calcForcesAndTorquesAndVirial();
      virialTensor += rbTau;
    }

    for (rb = mol2->beginRigidBody(rbIter); rb != NULL;
         rb = mol2->nextRigidBody(rbIter)) {
      Mat3x3d rbTau = rb->calcForcesAndTorquesAndVirial();
      virialTensor += rbTau;
    }

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, virialTensor.getArrayPointer(), 9, MPI_REALTYPE,
                  MPI_SUM, MPI_COMM_WORLD);
#endif
    curSnapshot->setVirialTensor(virialTensor);

    /*
      if (info_->getSimParams()->getUseLongRangeCorrections()) {

        RealType vol = curSnapshot->getVolume();
        RealType Elrc(0.0);
        RealType Wlrc(0.0);

        AtomTypeSet::iterator i;
        AtomTypeSet::iterator j;

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
        Elrc *= 2.0 * Constants::PI;
        Wlrc *= 2.0 * Constants::PI;

        RealType lrp = curSnapshot->getLongRangePotential();
        curSnapshot->setLongRangePotential(lrp + Elrc);
        virialTensor += Wlrc * SquareMatrix3<RealType>::identity();
        curSnapshot->setVirialTensor(virialTensor);
      }
    */
  }
}  // namespace OpenMD
