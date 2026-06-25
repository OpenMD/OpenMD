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
 */

/*! \file integrators/LHDForceModifier.hpp
    \brief Langevin force modifier with intramolecular RPY hydrodynamic
           interactions for flexible bead molecules in a background flow.

    Translation-only (unstructured spherical atom) case. Each molecule is a
    hydrodynamically coupled group of N spherical beads; molecules do not
    interact hydrodynamically with each other. Per step and per molecule:

      - assemble the 3N x 3N RPY mobility and its resistance R = M^{-1}
        (RPYMobility), Cholesky factor S of R;
      - sample the background flow at each bead (VelocityField) and fold in
        the dipolar disturbance to get v_inf_eff;
      - add the FDT-consistent correlated random force sqrt(2 kT/dt) S Z;
      - solve the coupled friction f = R (v_inf_eff - v) self-consistently for
        the full-step velocity (the LDForceModifier iteration, now coupled
        across the molecule's beads).

    Pairs with a velocity-Verlet (LangevinDynamics-style) integrator, which
    performs the actual velocity update; this modifier only injects forces.
*/

#ifndef INTEGRATOR_LHDFORCEMODIFIER_HPP
#define INTEGRATOR_LHDFORCEMODIFIER_HPP

#include <memory>
#include <random>
#include <vector>

#include "brains/ForceModifier.hpp"
#include "brains/Velocitizer.hpp"
#include "hydrodynamics/RPYMobility.hpp"
#include "perturbations/VelocityField.hpp"
#include "primitives/Atom.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/RandNumGen.hpp"

namespace OpenMD {

  class LHDForceModifier : public ForceModifier {
  public:
    LHDForceModifier(SimInfo* info);

    void modifyForces() override;

    int getMaxIterationNumber() const { return maxIterNum_; }
    void setMaxIterationNumber(int n) { maxIterNum_ = n; }
    RealType getForceTolerance() const { return forceTolerance_; }
    void setForceTolerance(RealType t) { forceTolerance_ = t; }

  private:
    //! Hydrodynamic radius of a bead, following the LDForceModifier rules:
    //! Lennard-Jones sigma/2, else the element van der Waals radius.
    RealType beadRadius(Atom* atom) const;

    //! One hydrodynamically coupled group (a molecule).
    struct HydroMolecule {
      std::vector<StuntDouble*> beads;
      std::vector<RealType> masses;
      std::unique_ptr<RPYMobility> mobility;
    };

    std::vector<HydroMolecule> molecules_;

    int maxIterNum_;
    RealType forceTolerance_;
    RealType dt_;
    RealType dt2_;
    RealType viscosity_;
    RealType kT_;  // Constants::kb * targetTemp

    Globals* simParams_;
    Utils::RandNumGenPtr randNumGen_;
    std::normal_distribution<RealType> normal_ {0.0, 1.0};
    std::unique_ptr<VelocityField> velField_;
    std::unique_ptr<Velocitizer> veloMunge_;
  };
}  // namespace OpenMD

#endif  // INTEGRATOR_LHDFORCEMODIFIER_HPP
