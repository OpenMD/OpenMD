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

#include <complex>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "applications/dynamicProps/TimeCorrFunc.hpp"
#include "brains/ForceManager.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/DynamicVector.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  /**
   * Computes the vibrational SFG susceptibility spectrum Im χ²(ω)
   * in the OH-stretch region using the exciton Hamiltonian approach with
   * an ADIABATIC quantum propagator.
   *
   * Algorithm mirrors the MultiSpec code (Kananenka group, exc.cpp):
   *
   *   C(t) = α_pq(t)ᵀ · F(t,0) · μ_r(0)
   *
   * where F(t,0) is the time-ordered adiabatic propagator:
   *
   *   F(t+dt, 0) = exp(i H(t) dt / ℏ) · F(t, 0)
   *
   * exp(i H dt / ℏ) is evaluated by diagonalizing H at each step:
   *   exp(i H dt / ℏ) = V · diag(exp(i ω_k dt / ℏ)) · V†
   * where V contains eigenvectors of H as columns.
   *
   * All quantities are kept in the LOCAL (site) basis throughout.
   * The trajectory is divided into navg windows of length nTimeBins_ frames,
   * separated by nSep_ frames, starting at nStart_; these are all populated
   * by the base-class setWindowingParameters() and preCorrelate(). We override
   * correlation() (which is virtual) to implement the adiabatic loop directly,
   * and computeProperty1() to cache per-frame data during the preCorrelate scan.
   *
   * Spectroscopic maps (Auer & Skinner 2008 for SPC/E;
   *                     Gruenbaum et al. 2013 for TIP4P):
   *   w10_   : ω₁₀(E) = a0 + a1·E + a2·E²        [cm⁻¹; E in a.u.]
   *   muPrime_: μ′(E) = m0 + m1·E + m2·E²        [a.u., dipole deriv.]
   *   x10_   : x₁₀(ω) = x0 + x1·ω₁₀             [a.u., coord. matrix el.]
   *   p10_   : p₁₀(ω) = p0 + p1·ω₁₀              [a.u., momentum matrix el.]
   *            Transition dipole |μ₁₀| = μ′·x₁₀ (converted to Debye).
   *   wintra_: ω^intra = [k0+k1·(Eⱼ+Eₖ)]·xⱼ·xₖ + kp·pⱼ·pₖ  [cm⁻¹]
   *   alphaMap_: (α_par, α_perp)                   [Å² amu⁻⁰·⁵]
   *
   * Supported polarization combinations (interface normal = z):
   *   ssp : Im C_yyz = α_yy(t)ᵀ · F(t,0) · μ_z(0)
   *   ppp : Im [−0.5·C_xxz − 0.5·C_yyz + C_zzz] (near-normal incidence)
   *   sps : Im C_yzy = α_yz(t)ᵀ · F(t,0) · μ_y(0)
   */
  class SFG : public SystemACF<RealType> {
  public:
    SFG(SimInfo* info, const std::string& filename, const std::string& sele1,
        const std::string& sele2, const std::string& polarization = "ssp",
        int privilegedAxis = 2,
        RealType t_apod = 0.0, RealType t_zerofill = 0.0,
        RealType fc = 50.0);

    virtual ~SFG();

  protected:
    // computeProperty1 is called by preCorrelate() for every frame.
    // We use it to build and cache the FrameData for that frame.
    void computeProperty1(int frame) override;

    // calcCorrVal satisfies the pure-virtual requirement but is never
    // called because we override correlation() to bypass correlateFrames().
    RealType calcCorrVal(int frame1, int frame2) override { return 0.0; }

    // Override the correlation loop with the adiabatic windowed propagation.
    // Called by the base doCorrelate() after preCorrelate() has read all
    // frames and populated nTimeBins_, dtMean_, nStart_, nSep_, nStride_,
    // navg_, allFrames_, etc.
    void correlation() override;

    void writeCorrelate() override;

  private:
    // -----------------------------------------------------------------------
    // Spectroscopic map tables, keyed by H atom type name.
    //
    // w10_     : (a0,a1,a2)  ω₁₀(E) = a0 + a1·E + a2·E²       [cm⁻¹]
    // muPrime_ : (m0,m1,m2)  μ′(E)  = m0 + m1·E + m2·E²       [a.u.]
    //            Dipole derivative along OH bond.
    // x10_     : (x0,x1)     x₁₀(ω) = x0 + x1·ω₁₀            [a.u.]
    //            Coordinate matrix element.  |μ₁₀| = μ′·x₁₀.
    // p10_     : (p0,p1)     p₁₀(ω) = p0 + p1·ω₁₀             [a.u.]
    //            Momentum matrix element.
    // wintra_  : (k0,k1,kinetic_coeff)
    //            ω^intra = [k0 + k1·(Eⱼ+Eₖ)]·xⱼ·xₖ
    //                    + kinetic_coeff·pⱼ·pₖ               [cm⁻¹]
    //            (kinetic_coeff < 0; includes −cosφ/m_O factor)
    // alphaMap_: (α_par, α_perp)                                [Å² amu⁻⁰·⁵]
    // tdcLoc_ : distance from O along OH bond for TDC location  [Å]
    //           (Auer 2008: 0.58 Å for SPC/E; Gruenbaum 2013: 0.67 Å for TIP4P)
    //
    // Bend maps (Ni & Skinner JCP 143, 014502 (2015), for H2O):
    //   wb01_   : (b0, b1)  ω_{0→1}(E_b) = b0 + b1·E_b   [cm⁻¹]
    //   wb12_   : (c0, c1)  ω_{1→2}(E_b) = c0 + c1·E_b   [cm⁻¹]
    //   ω₂δ on Hamiltonian diagonal = ω_{0→1} + ω_{1→2}.
    //   E_b is the field projected on the HOH bisector (computed in extractFrame).
    // -----------------------------------------------------------------------
    std::map<std::string, std::tuple<RealType,RealType,RealType>> w10_;
    std::map<std::string, std::tuple<RealType,RealType,RealType>> muPrime_;
    std::map<std::string, std::pair<RealType,RealType>>           x10_;
    std::map<std::string, std::pair<RealType,RealType>>           p10_;
    std::map<std::string, std::tuple<RealType,RealType,RealType>> wintra_;
    std::map<std::string, std::pair<RealType,RealType>>           alphaMap_;
    std::map<std::string, RealType>                              tdcLoc_;

    // Bend maps keyed by water-model name (e.g. "TIP4P", "SPCE", "TIP4P-Ice")
    std::map<std::string, std::pair<RealType,RealType>>           wb01_;
    std::map<std::string, std::pair<RealType,RealType>>           wb12_;

    // Fermi coupling [cm⁻¹] between OH stretch and HOH bend overtone
    // on the same molecule.  Default 50 cm⁻¹ following MultiSpec / Ni 2015.
    // If 0, Fermi resonance is effectively disabled.
    RealType fc_ {50.0};

    // True if a bend map is available for the water model(s) selected.
    // Set during the first extractFrame() call.
    bool useFermi_ {false};

    // Applied external field [kcal mol⁻¹ Å⁻¹ e⁻¹], pre-converted in ctor
    Vector3d EF_;

    // Whether the dump already contains electric fields
    bool dumpHasElectricFields_;
    ForceManager* forceMan_ {nullptr};

    // Polarization combination label ("ssp", "ppp", "sps")
    std::string polarization_;

    // Interface normal axis: 0=x, 1=y, 2=z (default 2).
    // s1_ and s2_ are the surface-parallel axes in cyclic order:
    //   n=2 (z): s1=0 (x), s2=1 (y)   [standard water/vapor interface]
    //   n=0 (x): s1=1 (y), s2=2 (z)
    //   n=1 (y): s1=2 (z), s2=0 (x)
    int privilegedAxis_ {2};
    int s1_             {0};
    int s2_             {1};

    // -----------------------------------------------------------------------
    // Per-frame data in the LOCAL (site) basis.
    //
    // When Fermi resonance is enabled, the Hamiltonian and chromophore arrays
    // include bend-overtone entries appended after the OH stretches:
    //   indices [0 .. N_stretch-1]            : OH-stretch chromophores
    //   indices [N_stretch .. N_stretch+N_bend-1] : HOH bend overtones
    //
    // Total chromophore count is fd.N = N_stretch + N_bend.
    // The bend overtones have zero transition dipole and polarizability
    // (their intensity is borrowed from the stretches via Fermi mixing).
    // -----------------------------------------------------------------------
    struct FrameData {
      int N {0};                        // total chromophores (stretches+bends)
      int nStretch {0};                 // number of OH stretches
      int nBend    {0};                 // number of bend overtones
      DynamicRectMatrix<RealType> H;    // NxN exciton Hamiltonian [cm-1]
      std::vector<Vector3d>  mu;        // transition dipole vectors [D]
      std::vector<Mat3x3d>   alpha;     // Raman polarizability tensors
      std::vector<int>       globalIDs; // global H-atom index per chromophore
                                         // (for bend overtones: O atom global ID)
      std::map<int,int>      idToIndex; // reverse map: globalID -> local index
    };

    // All frames from the dump, populated by computeProperty1() during
    // preCorrelate(). Indexed by frame number [0 .. nFrames_-1].
    std::vector<FrameData> allFrames_;

    // Accumulated TCF (complex; length nTimeBins_; averaged over navg_ windows)
    std::vector<std::complex<double>> tcf_ssp_;
    std::vector<std::complex<double>> tcf_ppp_;
    std::vector<std::complex<double>> tcf_sps_;

    // T1 vibrational relaxation time [ps]; 0 = no relaxation
    RealType T1_ps_ {0.0};

    // Apodization time constant [ps].
    // TCF is multiplied by exp(-t/t_apod_ps) before FFT.
    // 0 = no apodization.  Adds ~33.36/(pi*t_apod_ps) cm-1 of Lorentzian
    // broadening; a value of ~T_corr/3 is a good starting point.
    RealType t_apod_ps_ {0.0};

    // Zero-filling duration [ps].
    // The apodized TCF is padded with zeros out to
    // (T_corr + t_zerofill_ps) before FFT, interpolating the spectrum
    // onto a finer frequency grid without adding new information.
    // 0 = no zero-filling.  A value of 2-4x T_corr is typical.
    RealType t_zerofill_ps_ {0.0};

    // Average OH frequency computed from actual per-frame diagonal entries
    // across all frames during preCorrelate(). Finalized at the start of
    // correlation() and subtracted from every stored H diagonal before
    // propagation. Restored on the frequency axis in writeCorrelate().
    RealType wAvg_      {0.0};
    RealType wAvgSum_   {0.0};  // running sum of all ω₁₀ values seen
    long     wAvgCount_ {0};    // total number of ω₁₀ values accumulated

    // Frame interval [ps], set at start of correlation() and reused
    // in writeCorrelate()
    RealType dt_ps_ {0.001};

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /**
     * Extract spectroscopic-map quantities from currentSnapshot_ and build
     * a FrameData. Also fills ohPos and molIndex needed for buildHamiltonian.
     */
    FrameData extractFrame(std::vector<Vector3d>& ohPos,
                           std::vector<int>&      molIndex,
                           std::vector<RealType>& intraJ);

    /**
     * Fill fd.H with the N×N exciton Hamiltonian.
     * Diagonal: local ω₁₀. Off-diagonal: intraJ (same mol) or TDC.
     */
    void buildHamiltonian(FrameData&                   fd,
                          const std::vector<Vector3d>& ohPos,
                          const std::vector<int>&      molIndex,
                          const std::vector<RealType>& intraJ);

    /**
     * Propagate F by one time step:
     *   F ← exp(i H dt/ℏ) · F    (adiabatic step)
     *
     * @param F      N²-element complex vector (row-major N×N)
     * @param H      N×N real symmetric Hamiltonian [cm⁻¹]
     * @param dt_ps  time step [ps]
     * @param N      system size
     */
    void propagateF(std::vector<std::complex<double>>& F,
                    const DynamicRectMatrix<RealType>&  H,
                    RealType dt_ps, int N);

    /**
     * Accumulate one lag point into tcf_ssp_, tcf_ppp_, tcf_sps_.
     *
     * For each polarization element pqr:
     *   C_pqr(t) += Σᵢ α_pq,i(t) * [F · μ_r(0)]ᵢ  × exp(−t/2T1)
     *
     * @param tt     lag index [0 .. ncor_-1]
     * @param F      current N×N propagator (row-major complex vector)
     * @param fd0    FrameData at t=0  (μ₀ components)
     * @param fdt    FrameData at t    (α components)
     * @param N      system size
     * @param exptc  T1 decay: exp(−tt·dt_ps / (2·T1_ps_)); 1.0 if T1=0
     */
    void accumulateTCF(int tt,
                       const std::vector<std::complex<double>>& F,
                       const FrameData& fd0,
                       const FrameData& fdt,
                       int N, RealType exptc,
                       std::vector<std::complex<double>>& tgt_ssp,
                       std::vector<std::complex<double>>& tgt_ppp,
                       std::vector<std::complex<double>>& tgt_sps);

    /**
     * TDC coupling between sites i and j (MultiSpec convention).
     * J = AU_TO_WN × geom × (μ′·x₁₀)ᵢ × (μ′·x₁₀)ⱼ / r_bohr³
     * where geom = (eᵢ·eⱼ − 3(eᵢ·r̂)(eⱼ·r̂)) using unit OH vectors.
     */
    RealType tdcCoupling(const Vector3d& r_ij,
                         const Vector3d& e_i, const Vector3d& e_j,
                         RealType mx_i, RealType mx_j);

    /**
     * Build a FrameData containing only the atoms in refIDs, in that order,
     * drawn from src. Atoms absent from src get zero μ/α and their diagonal
     * H entry is taken from refFreqs. Returns false if more than half the
     * reference atoms are missing.
     */
    bool remapFrame(const FrameData&             src,
                    const std::vector<int>&      refIDs,
                    const std::vector<RealType>& refFreqs,
                    FrameData&                   out,
                    int                          refNStretch,
                    int                          refNBend);

    /**
     * Bond-polarizability tensor for one OH bond.
     * α = [α_perp·I + (α_par − α_perp)·ê⊗ê] × (x10/x10_gas)
     * The x10/x10_gas factor captures the field-dependent coordinate
     * matrix element following MultiSpec's trPol() convention.
     */
    Mat3x3d bondPolarizability(const Vector3d& ohUnit,
                               RealType alpha_par, RealType alpha_perp,
                               RealType x10, RealType x10_gas);
  };

}  // namespace OpenMD
#endif
