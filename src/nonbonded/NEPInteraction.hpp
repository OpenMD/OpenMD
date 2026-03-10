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

#ifndef NONBONDED_NEPINTERACTION_HPP
#define NONBONDED_NEPINTERACTION_HPP

#include <string>
#include <vector>

#include "NEP_CPU/src/nep.h"

namespace OpenMD {

  /**
   * @class NEPInteraction NEPInteraction.hpp "nonbonded/NEPInteraction.hpp"
   *
   * Manages per-atom descriptor storage and the two-pass (EAM-like)
   * evaluation of a Neuroevolution Potential (NEP) using OpenMD's own
   * pair-loop infrastructure.
   *
   * Pass 1 (PREPAIR_LOOP): accumulate radial and angular descriptors
   *   for each atom from the pair loop.
   * Between passes: run the ANN per atom to obtain Fp (dE/dq).
   * Pass 2 (PAIR_LOOP): compute forces and virial contributions.
   *
   * Unit convention: NEP outputs energies in eV and forces in eV/Å.
   * OpenMD uses kcal/mol and kcal/(mol·Å).  The conversion factor
   * eV_to_kcal_per_mol = 23.060541945329334 is applied internally.
   */
  class NEPInteraction {
  public:
    explicit NEPInteraction(const std::string& potentialFile);
    ~NEPInteraction() = default;

    /** Allocate (or reallocate) per-atom working arrays for nAtoms. */
    void allocate(int nAtoms);

    /** Zero descriptor accumulators at the start of PREPAIR_LOOP. */
    void zeroDescriptors();

    /**
     * Accumulate descriptor contributions from pair (atom1, atom2).
     * r12[3]: vector from atom1 to atom2 (Å).
     * d12   : ||r12|| (Å).
     * t1,t2 : NEP type indices (0-based) for atom1 and atom2.
     */
    void accumDescriptor(int atom1, int atom2,
                         const double* r12, double d12,
                         int t1, int t2);

    /**
     * Run the ANN for every atom to produce per-atom energies and
     * descriptor gradients (Fp).  Called once between the two passes.
     * nepTypes[atom] is the 0-based NEP type index for each local atom.
     */
    void runANN(const std::vector<int>& nepTypes);

    /**
     * Compute the combined force on atom1 from pair (atom1, atom2),
     * accounting for both atoms acting as the descriptor center.
     * r12[3]: vector from atom1 to atom2 (Å).
     * d12   : ||r12|| (Å).
     * t1,t2 : NEP type indices for atom1 and atom2.
     * f12[3]: force on atom1 (kcal/(mol·Å)) — incremented in place.
     */
    void calcForce(int atom1, int atom2,
                   const double* r12, double d12,
                   int t1, int t2,
                   double* f12);

    /** Sum of per-atom NEP energies in kcal/mol. */
    double getTotalEnergy() const;

    /** Return max(rc_radial, rc_angular) so ForceManager can set rCut_. */
    double getMaxCutoff() const;

    NEP3 nep3_;

  private:
    int nAtoms_         {0};
    int nMaxRadial_     {0};
    int nMaxAngular_    {0};
    int dim_            {0};
    int dimAngular_     {0};
    int sumFxyzStride_  {0};  // (nMaxAngular_+1)*NUM_OF_ABC

    // Per-atom descriptor storage (atom-major layout):
    //   q_radial_[atom*(nMaxRadial_+1) + n]
    std::vector<double> q_radial_;

    //   sum_fxyz_[atom*sumFxyzStride_ + n*NUM_OF_ABC + abc]
    std::vector<double> sum_fxyz_;

    //   Fp_[atom*dim_ + d]
    std::vector<double> Fp_;

    //   energy_[atom]  (eV)
    std::vector<double> energy_;

    // eV → kcal/mol
    static constexpr double eV_to_kcal_ = 23.060541945329334;
  };

}  // namespace OpenMD

#endif  // NONBONDED_NEPINTERACTION_HPP
