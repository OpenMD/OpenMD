/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "flucq/FluctuatingChargeObjectiveFunction.hpp"
#include "optimization/Constraint.hpp"
#include "optimization/EndCriteria.hpp"
#include "optimization/OptimizationFactory.hpp"
#include "optimization/Problem.hpp"
#include "optimization/StatusFunction.hpp"

using namespace QuantLib;
namespace OpenMD {

  FluctuatingChargePropagator::FluctuatingChargePropagator(SimInfo* info) :
      info_(info), forceMan_(NULL), hasFlucQ_(false), initialized_(false) {
    Globals* simParams = info_->getSimParams();
    fqParams_          = simParams->getFluctuatingChargeParameters();
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
        hasFlucQ_      = true;
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

      int maxIter              = fqParams_->getMaxIterations();
      RealType tolerance       = fqParams_->getTolerance();
      RealType initialStepSize = fqParams_->getInitialStepSize();

      EndCriteria endCriteria(maxIter, maxIter, tolerance, tolerance,
                              tolerance);
      std::string chargeOptMethod = fqParams_->getChargeOptimizationMethod();
      OptimizationMethod* minim =
          OptimizationFactory::getInstance().createOptimization(chargeOptMethod,
                                                                info_);

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
}  // namespace OpenMD
