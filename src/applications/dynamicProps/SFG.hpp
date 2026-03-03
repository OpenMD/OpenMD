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
 */

#ifndef APPLICATIONS_DYNAMICPROPS_SFG_HPP
#define APPLICATIONS_DYNAMICPROPS_SFG_HPP

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "applications/dynamicProps/TimeCorrFunc.hpp"
#include "brains/ForceManager.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  /**
   * Computes the vibrational SFG susceptibility spectrum Im χ²(ω)
   * in the OH-stretch region using the exciton Hamiltonian approach.
   *
   * Reference: Buch et al., J. Chem. Phys. 127, 204710 (2007);
   *            Kananenka group MultiSpec code.
   *
   * Cross-correlates the collective Raman polarizability αpq(t)
   * with the collective transition dipole μr(0) through instantaneous
   * exciton eigenstates built from the Skinner/Auer e-field spectroscopic
   * maps (same maps used in OHFrequencyMap).
   *
   * Supported polarization combinations (interface normal = z):
   *   ssp : ⟨αyy(t) μz(0)⟩
   *   ppp : combination of αzz, αxz, αzx components
   *   sps : ⟨αyz(t) μy(0)⟩
   */
  class SFG : public SystemACF<RealType> {
  public:
    SFG(SimInfo* info, const std::string& filename, const std::string& sele1,
        const std::string& sele2, const std::string& polarization = "ssp");

    virtual ~SFG();

  protected:
    // Called once per frame by SystemACF machinery (currentSnapshot_ is set)
    void computeProperty1(int frame) override;

    // Called for each frame pair (frame1=t0, frame2=t)
    RealType calcCorrVal(int frame1, int frame2) override;

    void postCorrelate() override;
    void writeCorrelate() override;

  private:
    // -----------------------------------------------------------------------
    // Spectroscopic map tables — keyed by atom type name (e.g. "H_SPCE",
    // "H_TIP4P").  Populated in the constructor, matching OHFrequencyMap.
    //
    // w10_      : (a0, a1, a2)  ω₁₀(E) = a0 + a1·E + a2·E²   [cm⁻¹]
    // muprime_  : (c0, c1, c2)  |dμ/dq|(E) = c0 + c1·E + c2·E²
    //                           [first entry is the zero-field value in D/Å·amu^0.5]
    // wintra_   : (d0, d1, d2)  J_intra(E) = d0 + d1·E + d2·E²  [cm⁻¹]
    //             field-dependent intramolecular coupling (Auer/Skinner 2008)
    // alpha_par_, alpha_perp_:  bond polarizability derivatives ∥ and ⊥
    //             stored as pairs (value, field_slope) — zero-field value used
    //             here; field dependence can be added later.
    // -----------------------------------------------------------------------
    std::map<std::string, std::tuple<RealType, RealType, RealType>> w10_;
    std::map<std::string, std::tuple<RealType, RealType, RealType>> muprime_;
    std::map<std::string, std::tuple<RealType, RealType, RealType>> wintra_;
    // Bond polarizability model parameters per atom type
    // stored as (alpha_parallel, alpha_perpendicular) — units Å³/Å·amu^0.5
    std::map<std::string, std::pair<RealType, RealType>> alphaMap_;

    // Applied (external) electric field in kcal mol⁻¹ Å⁻¹ e⁻¹
    Vector3d EF_;

    // Whether the dump file already contains electric fields
    bool dumpHasElectricFields_;

    // ForceManager needed when electric fields are absent from dump
    ForceManager* forceMan_ {nullptr};

    // -----------------------------------------------------------------------
    // Per-frame stored exciton data, indexed by frame number.
    // Populated in computeProperty1(); consumed in calcCorrVal().
    // -----------------------------------------------------------------------
    struct FrameData {
      std::vector<RealType> eigvals;  // exciton frequencies [cm⁻¹], length N
      std::vector<Vector3d> mu_ex;    // collective transition dipoles, length N
      std::vector<Mat3x3d>  alpha_ex; // collective Raman tensors, length N
    };
    std::vector<FrameData> frameData_; // frameData_[frame]

    // Polarization combination: "ssp", "ppp", or "sps"
    std::string polarization_;

    // -----------------------------------------------------------------------
    // Helper methods
    // -----------------------------------------------------------------------

    /**
     * Build and diagonalize the N_OH × N_OH exciton Hamiltonian for one frame.
     *
     * @param ohVecs      Unit OH bond vectors (lab frame), length N
     * @param ohPos       OH bond midpoint positions (Å, lab frame), length N
     * @param localFreqs  Local ω₁₀ frequencies from spectroscopic map [cm⁻¹]
     * @param localMu     Local transition dipole vectors (direction × magnitude)
     * @param molIndex    Molecule index for each OH (used for intra coupling)
     * @param intraJ      Field-dependent intramolecular coupling per OH [cm⁻¹]
     * @param eigvals     Output: eigenvalues [cm⁻¹]
     * @param eigvecs     Output: eigenvectors (row k = mode k coefficients)
     */
    void buildExcitonHamiltonian(
        const std::vector<Vector3d>& ohVecs,
        const std::vector<Vector3d>& ohPos,
        const std::vector<RealType>& localFreqs,
        const std::vector<Vector3d>& localMu,
        const std::vector<int>&      molIndex,
        const std::vector<RealType>& intraJ,
        std::vector<RealType>&       eigvals,
        std::vector<std::vector<RealType>>& eigvecs);

    /**
     * Transition dipole coupling between OH bonds i and j.
     * J_ij = prefactor × (μᵢ·μⱼ − 3(μᵢ·r̂)(μⱼ·r̂)) / r³
     * Units: D, Å → cm⁻¹  (prefactor = 5034.12 cm⁻¹ D⁻² Å³)
     */
    RealType tdcCoupling(const Vector3d& r_ij, const Vector3d& mu_i,
                         const Vector3d& mu_j);

    /**
     * Bond-polarizability Raman tensor for a single OH.
     * α = α_perp·I + (α_par − α_perp)·ê⊗ê
     */
    Mat3x3d bondPolarizability(const Vector3d& ohUnit, RealType alpha_par,
			       RealType alpha_perp);
    
    // Quantum correction Q(ω) = ω·β·ħ/2  (harmonic, classical-limit)
    RealType quantumCorrection(RealType omega) const;
  };

}  // namespace OpenMD
#endif
