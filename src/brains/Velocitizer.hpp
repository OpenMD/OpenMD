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

#ifndef BRAINS_VELOCITIZER_HPP
#define BRAINS_VELOCITIZER_HPP

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "io/Globals.hpp"
#include "utils/RandNumGen.hpp"

namespace OpenMD {

  //! Velocity-modifying routines.
  /*
    Provides functions that re-sample or re-scale the velocities.
    These are particularly useful for setting initial configurations
    to a desired temperature or re-scaling velocities to match a
    desired energy.
  */
  class Velocitizer {
  public:
    Velocitizer(SimInfo* info);
    virtual ~Velocitizer() = default;

    /** @brief Resamples velocities and angular momenta
     * Resamples velocities and angular momenta from a Maxwell-Boltzmann
     * distribution.
     * @param t : Temperature of the new distribution.
     */
    void randomize(RealType ct);

    /** @brief Resamples charge velocities
     * Resamples charge velocities from a Maxwell-Boltzmann
     * distribution.
     * @param ct : Charge temperature of the new distribution.
     */
    void randomizeChargeVelocity(RealType ct);

    /** @brief Scales velocities and angular momenta by a scaling factor
     * Rescales velocity (and angular momenta) by a scaling factor.
     * @param lambda : Scaling factor (must be >= 0)
     */
    void scale(RealType lambda);

    /** @brief Removes Center of Mass Drift Velocity
     * Removes the center of mass drift velocity (required for accurate
     * calculations of diffusion).
     */
    void removeComDrift();

    /** @brief Removes Center of Mass Angular momentum
     * Removes the center of mass angular momentum (particularly useful in
     * non-periodic simulations).
     */
    void removeAngularDrift();

  private:
    SimInfo* info_ {nullptr};
    Thermo thermo_;
    Globals* globals_;
    Utils::RandNumGenPtr randNumGen_ {nullptr};
  };
}  // namespace OpenMD

#endif  // BRAINS_VELOCITIZER_HPP
