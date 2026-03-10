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

#include "nonbonded/NEPInteraction.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace OpenMD {

  NEPInteraction::NEPInteraction(const std::string& potentialFile)
      : nep3_(potentialFile)
  {
    nMaxRadial_    = nep3_.paramb.n_max_radial;
    nMaxAngular_   = nep3_.paramb.n_max_angular;
    dim_           = nep3_.annmb.dim;
    dimAngular_    = nep3_.paramb.dim_angular;
    sumFxyzStride_ = (nMaxAngular_ + 1) * NEP3::NUM_OF_ABC;
  }

  void NEPInteraction::allocate(int nAtoms)
  {
    if (nAtoms == nAtoms_) return;
    nAtoms_ = nAtoms;
    q_radial_.assign(nAtoms * (nMaxRadial_ + 1), 0.0);
    sum_fxyz_.assign(nAtoms * sumFxyzStride_, 0.0);
    Fp_.assign(nAtoms * dim_, 0.0);
    energy_.assign(nAtoms, 0.0);
  }

  void NEPInteraction::zeroDescriptors()
  {
    std::fill(q_radial_.begin(), q_radial_.end(), 0.0);
    std::fill(sum_fxyz_.begin(), sum_fxyz_.end(), 0.0);
  }

  void NEPInteraction::accumDescriptor(int atom1, int atom2,
                                       const double* r12, double d12,
                                       int t1, int t2)
  {
    // Radial: atom1 as center
    nep3_.accum_radial(t1, t2, d12,
                       &q_radial_[atom1 * (nMaxRadial_ + 1)]);
    // Radial: atom2 as center (same d12, same cutoff radius — symmetric)
    nep3_.accum_radial(t2, t1, d12,
                       &q_radial_[atom2 * (nMaxRadial_ + 1)]);

    // Angular: atom1 as center, r12 = pos[atom2] - pos[atom1]
    nep3_.accum_angular(t1, t2, d12, r12,
                        &sum_fxyz_[atom1 * sumFxyzStride_]);
    // Angular: atom2 as center, r21 = -r12
    const double r21[3] = {-r12[0], -r12[1], -r12[2]};
    nep3_.accum_angular(t2, t1, d12, r21,
                        &sum_fxyz_[atom2 * sumFxyzStride_]);
  }

  void NEPInteraction::runANN(const std::vector<int>& nepTypes)
  {
    for (int i = 0; i < nAtoms_; ++i) {
      int t1 = nepTypes[i];
      if (t1 < 0) {
        energy_[i] = 0.0;
        std::fill(&Fp_[i * dim_], &Fp_[i * dim_] + dim_, 0.0);
        continue;
      }
      nep3_.run_ann_one_atom(t1,
                             &q_radial_[i * (nMaxRadial_ + 1)],
                             &sum_fxyz_[i * sumFxyzStride_],
                             energy_[i],
                             &Fp_[i * dim_]);
    }
  }

  void NEPInteraction::calcForce(int atom1, int atom2,
                                 const double* r12, double d12,
                                 int t1, int t2,
                                 double* f12)
  {
    const double* Fp1     = &Fp_[atom1 * dim_];
    const double* Fp2     = &Fp_[atom2 * dim_];
    const double* sfxyz1  = &sum_fxyz_[atom1 * sumFxyzStride_];
    const double* sfxyz2  = &sum_fxyz_[atom2 * sumFxyzStride_];

    // Radial force (eV/Å); both center contributions combined
    double f12_eV[3] = {0.0, 0.0, 0.0};
    nep3_.calc_radial_force(t1, t2, d12, r12, Fp1, Fp2, f12_eV);

    // Angular force (eV/Å); both center contributions combined
    const double* Fp1_ang = Fp1 + (nMaxRadial_ + 1);
    const double* Fp2_ang = Fp2 + (nMaxRadial_ + 1);
    nep3_.calc_angular_force(t1, t2, d12, r12,
                             Fp1_ang, sfxyz1,
                             Fp2_ang, sfxyz2,
                             f12_eV);

    // Convert eV/Å → kcal/(mol·Å) and accumulate
    for (int d = 0; d < 3; ++d) {
      f12[d] += f12_eV[d] * eV_to_kcal_;
    }
  }

  double NEPInteraction::getTotalEnergy() const
  {
    double total = std::accumulate(energy_.begin(), energy_.end(), 0.0);
    return total * eV_to_kcal_;
  }

  double NEPInteraction::getMaxCutoff() const
  {
    return std::max(nep3_.paramb.rc_radial, nep3_.paramb.rc_angular);
  }

}  // namespace OpenMD
