/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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
    info_(info), hasFlucQ_(false), forceMan_(NULL) {
    
    Globals* simParams = info_->getSimParams();
    fqParams_ = simParams->getFluctuatingChargeParameters();
    fqConstraints_ = new FluctuatingChargeConstraints(info_);
  }

  FluctuatingChargePropagator::~FluctuatingChargePropagator() {
    if (fqConstraints_ != NULL) delete fqConstraints_;
  }

  void FluctuatingChargePropagator::setForceManager(ForceManager* forceMan) {
    forceMan_ = forceMan;
  }

  void FluctuatingChargePropagator::initialize() {

    if (info_->usesFluctuatingCharges()) {
      if (info_->getNFluctuatingCharges() > 0) {
        hasFlucQ_ = true;
      }
    }

    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator  j;
    Molecule* mol;
    Atom* atom;
    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        atom->setFlucQPos(0.0);
        atom->setFlucQVel(0.0);
      }
    }
    
    FluctuatingChargeObjectiveFunction flucQobjf(info_, forceMan_, 
                                                 fqConstraints_);

    DynamicVector<RealType> initCoords = flucQobjf.setInitialCoords();
    Problem problem(flucQobjf, *(new NoConstraint()), *(new NoStatus()), 
                    initCoords);

    EndCriteria endCriteria(1000, 100, 1e-5, 1e-5, 1e-5);       

    OptimizationMethod* minim = OptimizationFactory::getInstance()->createOptimization("SD", info_);

    DumpStatusFunction dsf(info_);  // we want a dump file written
                                    // every iteration

    minim->minimize(problem, endCriteria);
  }

  void FluctuatingChargePropagator::applyConstraints() {
    if (!hasFlucQ_) return;
    fqConstraints_->applyConstraints();
  }
}
