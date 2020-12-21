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
 
#include "flucq/FluctuatingChargePropagator.hpp"
#include "flucq/FluctuatingChargeObjectiveFunction.hpp"
#include "optimization/Constraint.hpp"
#include "optimization/Problem.hpp"
#include "optimization/EndCriteria.hpp"
#include "optimization/StatusFunction.hpp"
#include "optimization/OptimizationFactory.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

using namespace QuantLib;
namespace OpenMD {

  FluctuatingChargePropagator::FluctuatingChargePropagator(SimInfo* info) : 
    info_(info), forceMan_(NULL), hasFlucQ_(false), initialized_(false) {
    
    Globals* simParams = info_->getSimParams();
    fqParams_ = simParams->getFluctuatingChargeParameters();    
  }

  FluctuatingChargePropagator::~FluctuatingChargePropagator() {
    if (info_->usesFluctuatingCharges() && info_->getNFluctuatingCharges() > 0)
      delete fqConstraints_;
  }

  void FluctuatingChargePropagator::setForceManager(ForceManager* forceMan) {
    forceMan_ = forceMan;
  }

  void FluctuatingChargePropagator::initialize() {
    if (info_->usesFluctuatingCharges()) {
      if (info_->getNFluctuatingCharges() > 0) {
        hasFlucQ_ = true;
	      fqConstraints_ = new FluctuatingChargeConstraints(info_);
	      fqConstraints_->setConstrainRegions(fqParams_->getConstrainRegions());
      }
    }
    
    if (!hasFlucQ_) {
      initialized_ = true;
      return;
    }

    // SimInfo::MoleculeIterator i;
    // Molecule::FluctuatingChargeIterator  j;
    // Molecule* mol;
    // Atom* atom;
    //  
    // For single-minima flucq, this ensures a net neutral system, but
    // for multiple minima, this is no longer the right thing to do:
    //
    // for (mol = info_->beginMolecule(i); mol != NULL; 
    //      mol = info_->nextMolecule(i)) {
    //   for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
    //        atom = mol->nextFluctuatingCharge(j)) {
    //     atom->setFlucQPos(0.0);
    //     atom->setFlucQVel(0.0);
    //   }
    // }


    if (fqParams_->getDoInitialOptimization()) {
    
      FluctuatingChargeObjectiveFunction flucQobjf(info_, forceMan_, 
                                                   fqConstraints_);
      
      DynamicVector<RealType> initCoords = flucQobjf.setInitialCoords();

      NoConstraint noConstraint {};
      NoStatus noStatus {};
      
      Problem problem(flucQobjf, noConstraint, noStatus, initCoords);
      
      
      int maxIter = fqParams_->getMaxIterations();
      RealType tolerance = fqParams_->getTolerance();
      RealType initialStepSize = fqParams_->getInitialStepSize();
      
      EndCriteria endCriteria(maxIter, maxIter, tolerance, tolerance,
                              tolerance);
      std:: string chargeOptMethod = fqParams_->getChargeOptimizationMethod();      
      OptimizationMethod* minim = OptimizationFactory::getInstance().createOptimization(chargeOptMethod, info_);
      
      minim->minimize(problem, endCriteria, initialStepSize);

      delete minim;
    }
    
    initialized_ = true;
  }

  void FluctuatingChargePropagator::applyConstraints() {
    if (!initialized_) initialize();
    if (!hasFlucQ_) return;

    fqConstraints_->applyConstraints();
  }
}
