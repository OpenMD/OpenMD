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

#include "applications/dynamicProps/SFG.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>

#include "brains/ForceManager.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "io/DumpReader.hpp"
#include "math/Eigenvalue.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

#if defined(HAVE_FFTW3_H)
#include <fftw3.h>
#endif

// LAPACK dsyevd (divide-and-conquer symmetric eigensolver)
#if defined(HAVE_LAPACK)
extern "C" {
  void dsyevd_(const char* jobz, const char* uplo, const int* n,
               double* a, const int* lda, double* w,
               double* work, const int* lwork,
               int* iwork, const int* liwork, int* info);
}
#endif

// CBLAS for optimised dgemm (optional, used alongside LAPACK)
#if defined(HAVE_CBLAS)
#include <cblas.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace OpenMD {

  // ℏ in cm⁻¹·ps  (= ℏ / (hc × 100) × 10¹²)
  // Value from MultiSpec constants.hpp; derivation:
  //   ℏ = 1.054571817e-34 J·s
  //   1 cm⁻¹ = 1.9864458e-23 J
  //   1 ps = 1e-12 s
  //   → ℏ/(1 cm⁻¹ × 1 ps) = 1.054571817e-34 / (1.9864458e-23 × 1e-12) = 5.3088...
  static constexpr double HBAR_CM_PS = 5.3088373;

  // =========================================================================
  // Constructor
  // =========================================================================
  SFG::SFG(SimInfo* info, const std::string& filename,
           const std::string& sele1, const std::string& sele2,
           const std::string& polarization, int privilegedAxis,
           RealType t_apod, RealType t_zerofill,
           RealType fc) :
    SystemACF<RealType>(info, filename, sele1, sele2),
    polarization_(polarization),
    privilegedAxis_(privilegedAxis),
    s1_((privilegedAxis + 1) % 3),
    s2_((privilegedAxis + 2) % 3),
    t_apod_ps_(t_apod * 1e-3),
    t_zerofill_ps_(t_zerofill * 1e-3),
    fc_(fc) {

    setCorrFuncType("SFG");
    setOutputName(getPrefix(dumpFilename_) + ".sfg");
    setLabelString("Im chi2(omega)");

    // -----------------------------------------------------------------------
    // Spectroscopic maps.
    //
    // w10_     : ω₁₀(E) = a0 + a1·E + a2·E²            [cm⁻¹; E in a.u.]
    // muPrime_ : μ′(E)  = m0 + m1·E + m2·E²            [a.u.]
    //            Dipole derivative along OH bond direction.
    // x10_     : x₁₀(ω) = x0 + x1·ω₁₀                 [a.u.]
    //            Coordinate matrix element.  |μ₁₀| = μ′ × x₁₀  [a.u.]
    //            Converted to Debye·amu⁻⁰·⁵ by × 4.8032 × 7.2906.
    // p10_     : p₁₀(ω) = p0 + p1·ω₁₀                  [a.u.]
    //            Momentum matrix element.
    // wintra_  : ω^intra_jk = [k0 + k1·(Eⱼ+Eₖ)]·xⱼ·xₖ
    //                        + kinetic_coeff · pⱼ·pₖ    [cm⁻¹]
    // alphaMap_: (α_par, α_perp)                         [Å² amu⁻⁰·⁵]
    //
    // Sources: Auer & Skinner, J. Chem. Phys. 128, 224511 (2008) — SPC/E
    //          Gruenbaum et al., JCTC 9, 3109 (2013)     — TIP4P
    //          Takayama et al., J. Chem. Phys. 158, 136101 (2023) — TIP4P-Ice
    // -----------------------------------------------------------------------

    // SPC/E  (Auer & Skinner 2008, Table I)
    //   μ′/μ′_g = 0.7112 + 75.59·E  →  μ′(E) = 0.1874*(0.7112 + 75.59·E)
    //           = 0.1333 + 14.17·E + 0·E²
    w10_["H_SPCE"]      = std::make_tuple(3761.6,  -5060.4,  -86225.0);
    muPrime_["H_SPCE"] = std::make_tuple(0.1333, 14.17, 0.0);
    x10_["H_SPCE"]      = std::make_pair(0.1934,  -1.75e-5);
    p10_["H_SPCE"]      = std::make_pair(1.611,    5.893e-4);
    wintra_["H_SPCE"]   = std::make_tuple(-1789.0,  23852.0,  -1.966);
    alphaMap_["H_SPCE"] = std::make_pair(0.185, 0.033);
    tdcLoc_["H_SPCE"]   = 0.58;   // Å from O along OH bond

    // TIP4P  (Gruenbaum et al. 2013, Table 1)
    w10_["H_TIP4P"]      = std::make_tuple(3760.2,  -3541.7, -152677.0);
    muPrime_["H_TIP4P"]  = std::make_tuple(0.1646,  11.39,    63.41);
    x10_["H_TIP4P"]      = std::make_pair(0.19285, -1.7261e-5);
    p10_["H_TIP4P"]      = std::make_pair(1.6466,   5.7692e-4);
    wintra_["H_TIP4P"]   = std::make_tuple(-1361.0,  27165.0,  -1.887);
    alphaMap_["H_TIP4P"] = std::make_pair(0.185, 0.033);
    tdcLoc_["H_TIP4P"]   = 0.67;  // Å from O along OH bond

    // TIP4P-Ice (transferability from TIP4P)
    w10_["H_TIP4P-Ice"]      = w10_["H_TIP4P"];
    muPrime_["H_TIP4P-Ice"]  = muPrime_["H_TIP4P"];
    x10_["H_TIP4P-Ice"]      = x10_["H_TIP4P"];
    p10_["H_TIP4P-Ice"]      = p10_["H_TIP4P"];
    wintra_["H_TIP4P-Ice"]   = wintra_["H_TIP4P"];
    alphaMap_["H_TIP4P-Ice"] = alphaMap_["H_TIP4P"];
    tdcLoc_["H_TIP4P-Ice"]   = tdcLoc_["H_TIP4P"];

    // -----------------------------------------------------------------------
    // HOH bend overtone maps (Ni & Skinner JCP 143, 014502 (2015)).
    //
    //   ω_HOH_10(E_b) = b0_10 + b1_10 · E_b      [cm⁻¹]
    //   ω_HOH_21(E_b) = b0_21 + b1_21 · E_b      [cm⁻¹]
    //   ω_2δ on Hamiltonian diagonal = ω_HOH_10 + ω_HOH_21
    //
    // Field-free limit: ω_2δ ≈ 1581 + 1551 = 3133 cm⁻¹, right in the
    // H-bonded OH stretch region as expected for Fermi resonance.
    //
    // The map was originally fit for TIP4P; we apply it to all model
    // variants for which we have a stretch map.
    // -----------------------------------------------------------------------
    wb01_["TIP4P"]     = std::make_pair(1581.46, 2938.51);
    wb12_["TIP4P"]     = std::make_pair(1551.32, 3147.80);
    wb01_["TIP4P-Ice"] = wb01_["TIP4P"];
    wb12_["TIP4P-Ice"] = wb12_["TIP4P"];
    // Apply same map to SPC/E as well — Auer/Skinner 2008 didn't publish
    // a separate bend map for SPC/E, so we use the Ni 2015 TIP4P map.
    // The bend frequency depends weakly on the water model (compared to
    // the stretch), so this transfer is a reasonable approximation.
    wb01_["SPCE"]      = wb01_["TIP4P"];
    wb12_["SPCE"]      = wb12_["TIP4P"];

    // -----------------------------------------------------------------------
    // Applied external electric field.
    // EF_ stored in kcal mol⁻¹ Å⁻¹ e⁻¹ (same units as getElectricField()).
    // SimParams gives fields in V/Å; chrgToKcal converts to kcal/mol/Å/e.
    // -----------------------------------------------------------------------
    EF_ = V3Zero;
    const RealType chrgToKcal = 23.060548;
    std::vector<RealType> ef;
    bool efSpec = false;
    if (info_->getSimParams()->haveElectricField()) {
      efSpec = true; ef = info_->getSimParams()->getElectricField();
    }
    if (info_->getSimParams()->haveUniformField()) {
      efSpec = true; ef = info_->getSimParams()->getUniformField();
    }
    if (efSpec) {
      if (ef.size() != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "SFG: uniformField must have 3 components (%zu given).\n",
                 ef.size());
        painCave.isFatal = 1; simError();
      }
      EF_.x() = ef[0] * chrgToKcal;
      EF_.y() = ef[1] * chrgToKcal;
      EF_.z() = ef[2] * chrgToKcal;
    }

    // -----------------------------------------------------------------------
    // Electric field storage setup (mirroring OHFrequencyMap).
    // -----------------------------------------------------------------------
    dumpHasElectricFields_ = info_->getSimParams()->getOutputElectricField();

    if (!dumpHasElectricFields_) {
      int asl  = info_->getAtomStorageLayout()      | DataStorage::dslElectricField;
      int rbsl = info_->getRigidBodyStorageLayout() | DataStorage::dslElectricField;
      int cgsl = info_->getCutoffGroupStorageLayout();
      info_->setAtomStorageLayout(asl);
      info_->setRigidBodyStorageLayout(rbsl);
      info_->setCutoffGroupStorageLayout(cgsl);
      info_->setSnapshotManager(new SimSnapshotManager(info_, asl, rbsl, cgsl));
      forceMan_ = new ForceManager(info_);
      forceMan_->setDoElectricField(true);
      forceMan_->initialize();
    }
    info_->getSimParams()->setOutputElectricField(true);
  }

  SFG::~SFG() { delete forceMan_; }

  // =========================================================================
  // computeProperty1
  //
  // Called by preCorrelate() via computeFrame() for every frame istep.
  // We build and cache the FrameData for that frame into allFrames_[istep].
  // preCorrelate() has already loaded the snapshot into currentSnapshot_.
  // If electric fields are not in the dump, we recalculate them here.
  // =========================================================================
  void SFG::computeProperty1(int frame) {
    // Grow allFrames_ to cover this frame index (preCorrelate calls frames
    // 0..nFrames_-1 in order, so a simple push_back would also work, but
    // resize+assign is safer if frames ever come out of order).
    if (static_cast<int>(allFrames_.size()) <= frame)
      allFrames_.resize(frame + 1);

    if (!dumpHasElectricFields_) {
      forceMan_->setDoElectricField(true);
      forceMan_->calcForces();
    }

    std::vector<Vector3d> ohPos;
    std::vector<int>      molIndex;
    std::vector<RealType> intraJ;
    allFrames_[frame] = extractFrame(ohPos, molIndex, intraJ);
    buildHamiltonian(allFrames_[frame], ohPos, molIndex, intraJ);

    // Accumulate raw diagonal frequencies for wAvg_ computation.
    // buildHamiltonian leaves H diagonal as unshifted ω₁₀ (stretch) and
    // ω_2δ (bend overtone) values at this point.  We only average stretch
    // frequencies — including bend overtones (at ~3130 cm⁻¹) would skew
    // wAvg_ and shift the spectrum's frequency axis incorrectly.
    const FrameData& fd = allFrames_[frame];
    for (int i = 0; i < fd.nStretch; ++i) {
      wAvgSum_   += fd.H(i, i);
      wAvgCount_ += 1;
    }
  }

  // =========================================================================
  // correlation
  //
  // Overrides TimeCorrFunc::correlation().  Called by the base doCorrelate()
  // after preCorrelate() has finished, so by the time we get here:
  //   - allFrames_[0..nFrames_-1]  are fully populated
  //   - nTimeBins_  = correlation window length in frames
  //   - nStart_     = frames to skip at start
  //   - nSep_       = gap frames between windows
  //   - nStride_    = nTimeBins_ + nSep_
  //   - navg_       = number of windows that fit
  //   - dtMean_     = mean frame interval [fs]
  //
  // Implements the MultiSpec adiabatic windowed propagation:
  //   for each window iwin (origin frame = nStart_ + iwin * nStride_):
  //     F = I (identity, N×N complex)
  //     for tt = 0 .. nTimeBins_-1:
  //       if tt > 0: F ← exp(i H(t-1) dt/ℏ) · F
  //       C(tt) += αᵀ(t) · F · μ(0)
  // =========================================================================
  void SFG::correlation() {
    const std::complex<double> czero(0.0, 0.0);

    int    ncor     = static_cast<int>(nTimeBins_);
    RealType dt_ps  = dtMean_ * 1.0e-3;  // fs → ps
    dt_ps_          = dt_ps;

    // Finalize wAvg_ from all per-frame diagonal entries accumulated during
    // preCorrelate(). This is the true mean ω₁₀ across the trajectory and
    // the selected chromophores, regardless of water model.
    wAvg_ = (wAvgCount_ > 0)
      ? wAvgSum_ / static_cast<RealType>(wAvgCount_)
      : 0.0;

    // Apply the shift to every stored Hamiltonian diagonal in one pass.
    // The off-diagonals (TDC couplings) are unaffected.
    for (auto& fd : allFrames_) {
      for (int i = 0; i < fd.N; ++i)
        fd.H(i, i) -= wAvg_;
    }

    tcf_ssp_.assign(ncor, czero);
    tcf_ppp_.assign(ncor, czero);
    tcf_sps_.assign(ncor, czero);

    progressBar_->clear();

#if defined(_OPENMP)
    // -----------------------------------------------------------------
    // OpenMP parallel window loop.
    //
    // Each thread accumulates into its own local TCF vectors, which are
    // merged into the shared tcf_ssp_/ppp_/sps_ at the end.  All per-
    // window state (F, refFreqs, fd0, fdt) is already thread-local.
    // The propagateF scratch (Hbuf, Vbuf, etc.) is thread_local static.
    // allFrames_ is read-only and shared.
    // -----------------------------------------------------------------
    int completedWindows = 0;

#pragma omp parallel
    {
      // Per-thread TCF accumulators
      std::vector<std::complex<double>> loc_ssp(ncor, czero);
      std::vector<std::complex<double>> loc_ppp(ncor, czero);
      std::vector<std::complex<double>> loc_sps(ncor, czero);

#pragma omp for schedule(dynamic)
      for (int iwin = 0; iwin < navg_; ++iwin) {
        int istart = nStart_ + iwin * nStride_;
        if (istart >= nFrames_) continue;

        const FrameData& fd0raw = allFrames_[istart];
        int N = fd0raw.N;
        if (N == 0) continue;

        const std::vector<int>& refIDs = fd0raw.globalIDs;

        std::vector<RealType> refFreqs(N);
        for (int i = 0; i < N; ++i)
          refFreqs[i] = fd0raw.H(i, i);

        // F = identity
        std::vector<std::complex<double>> F(N * N, czero);
        for (int i = 0; i < N; ++i) F[i*N+i] = std::complex<double>(1.0, 0.0);

        FrameData fd0;
        remapFrame(fd0raw, refIDs, refFreqs, fd0,
                   fd0raw.nStretch, fd0raw.nBend);

        for (int tt = 0; tt < ncor; ++tt) {
          int iframe = istart + tt;
          if (iframe >= nFrames_) break;

          FrameData fdt;
          if (!remapFrame(allFrames_[iframe], refIDs, refFreqs, fdt,
                          fd0raw.nStretch, fd0raw.nBend))
            break;

          for (int i = 0; i < N; ++i) {
            auto it = allFrames_[iframe].idToIndex.find(refIDs[i]);
            if (it != allFrames_[iframe].idToIndex.end())
              refFreqs[i] = allFrames_[iframe].H(it->second, it->second);
          }

          if (tt > 0)
            propagateF(F, fdt.H, dt_ps, N);

          RealType exptc = (T1_ps_ > 0.0)
            ? std::exp(-static_cast<RealType>(tt) * dt_ps / (2.0 * T1_ps_))
            : 1.0;

          accumulateTCF(tt, F, fd0, fdt, N, exptc,
                        loc_ssp, loc_ppp, loc_sps);
        }

        // Progress update (atomic to avoid garbled output)
#pragma omp atomic
        completedWindows++;

        if (omp_get_thread_num() == 0) {
          progressBar_->setStatus(completedWindows, navg_);
          progressBar_->update();
        }
      }

      // Merge per-thread accumulators into shared TCF vectors
#pragma omp critical
      {
        for (int tt = 0; tt < ncor; ++tt) {
          tcf_ssp_[tt] += loc_ssp[tt];
          tcf_ppp_[tt] += loc_ppp[tt];
          tcf_sps_[tt] += loc_sps[tt];
        }
      }
    }  // end omp parallel

#else
    // -----------------------------------------------------------------
    // Serial fallback (no OpenMP).
    // -----------------------------------------------------------------
    for (int iwin = 0; iwin < navg_; ++iwin) {
      int istart = nStart_ + iwin * nStride_;
      if (istart >= nFrames_) break;

      progressBar_->setStatus(iwin + 1, navg_);
      progressBar_->update();
      
      const FrameData& fd0raw = allFrames_[istart];
      int N = fd0raw.N;
      if (N == 0) continue;

      const std::vector<int>& refIDs = fd0raw.globalIDs;

      std::vector<RealType> refFreqs(N);
      for (int i = 0; i < N; ++i)
        refFreqs[i] = fd0raw.H(i, i);

      // F = identity
      std::vector<std::complex<double>> F(N * N, czero);
      for (int i = 0; i < N; ++i) F[i*N+i] = std::complex<double>(1.0, 0.0);

      FrameData fd0;
      remapFrame(fd0raw, refIDs, refFreqs, fd0,
                 fd0raw.nStretch, fd0raw.nBend);

      for (int tt = 0; tt < ncor; ++tt) {
        int iframe = istart + tt;
        if (iframe >= nFrames_) break;

        FrameData fdt;
        if (!remapFrame(allFrames_[iframe], refIDs, refFreqs, fdt,
                        fd0raw.nStretch, fd0raw.nBend))
          break;

        for (int i = 0; i < N; ++i) {
          auto it = allFrames_[iframe].idToIndex.find(refIDs[i]);
          if (it != allFrames_[iframe].idToIndex.end())
            refFreqs[i] = allFrames_[iframe].H(it->second, it->second);
        }

        if (tt > 0)
          propagateF(F, fdt.H, dt_ps, N);

        RealType exptc = (T1_ps_ > 0.0)
          ? std::exp(-static_cast<RealType>(tt) * dt_ps / (2.0 * T1_ps_))
          : 1.0;

        accumulateTCF(tt, F, fd0, fdt, N, exptc,
                      tcf_ssp_, tcf_ppp_, tcf_sps_);
      }
    }
#endif

    // Normalize by number of windows
    if (navg_ > 0) {
      RealType inv = 1.0 / static_cast<RealType>(navg_);
      for (int tt = 0; tt < ncor; ++tt) {
        tcf_ssp_[tt] *= inv;
        tcf_ppp_[tt] *= inv;
        tcf_sps_[tt] *= inv;
      }
    }
  }



  // =========================================================================
  // extractFrame
  //
  // Reads the current snapshot and builds a FrameData in the local (site)
  // basis. The Hamiltonian matrix is NOT filled here; call buildHamiltonian
  // afterward. ohPos and molIndex are returned for Hamiltonian construction.
  // =========================================================================
  SFG::FrameData SFG::extractFrame(std::vector<Vector3d>& ohPos,
				   std::vector<int>&      molIndex,
				   std::vector<RealType>& intraJ) {
    const RealType kcalToAU = 0.0008432975573;

    FrameData fd;
    ohPos.clear();
    molIndex.clear();
    intraJ.clear();

    std::vector<RealType> freqs;

    // Per-chromophore electric field, x₁₀, p₁₀ needed for
    // intramolecular coupling: ω^intra = [k0+k1(Ei+Ej)] xi xj + kp pi pj
    std::vector<RealType> eFields;   // projected E in a.u.
    std::vector<RealType> x10vals;   // coordinate matrix element [a.u.]
    std::vector<RealType> p10vals;   // momentum matrix element [a.u.]
    std::vector<RealType> mxvals;    // μ′·x₁₀ in a.u. (for TDC)
    std::vector<RealType> intraJ_k0; // per-site k0 from wintra_ map
    std::vector<RealType> intraJ_kp; // per-site kp from wintra_ map

    // Per-molecule data collected during the OH loop, used afterward to
    // append bend overtone chromophores.  Indexed by molID.
    struct MolData {
      Vector3d Opos;
      Vector3d ohUnit1, ohUnit2;     // unit OH vectors (in-plane basis)
      Vector3d Efield1, Efield2;     // field at H1, H2 in a.u.
      RealType rOH1 = 0.0, rOH2 = 0.0;
      int      ohCount = 0;          // 1 or 2 stretches contributed
      int      molOgID = -1;         // global O atom ID (used as bend ID)
      std::string modelName;         // water-model type name (e.g. "TIP4P")
    };
    std::map<int, MolData> molMap;   // molID -> MolData

    int ii;
    StuntDouble* sd1;

    for (sd1 = seleMan1_.beginSelected(ii); sd1 != nullptr;
	 sd1 = seleMan1_.nextSelected(ii)) {

      if (!sd1->isAtom()) continue;

      Atom* atom      = dynamic_cast<Atom*>(sd1);
      AtomType* at    = atom->getAtomType();
      std::string name = at->getName();

      auto wIt = w10_.find(name);
      if (wIt == w10_.end()) continue;
      auto mpIt = muPrime_.find(name);
      if (mpIt == muPrime_.end()) continue;

      int sdID  = sd1->getGlobalIndex();
      int molID = info_->getGlobalMolMembership(sdID);
      Molecule* mol = info_->getMoleculeByGlobalIndex(molID);

      Vector3d Opos;
      int      OgID = -1;
      bool foundO = false;
      Molecule::AtomIterator ai;
      for (Atom* oatom = mol->beginAtom(ai); oatom != nullptr;
           oatom = mol->nextAtom(ai)) {
        if (oatom->getAtomType()->getName()[0] == 'O') {
          Opos   = oatom->getPos();
          OgID   = oatom->getGlobalIndex();
          foundO = true;
          break;
        }
      }
      if (!foundO) continue;

      Vector3d Hpos    = sd1->getPos();
      Vector3d rOH_mic = Hpos - Opos;
      currentSnapshot_->wrapVector(rOH_mic);   // minimum-image OH vector
      RealType rOH_len = rOH_mic.length();
      Vector3d ohUnit  = rOH_mic;
      ohUnit.normalize();

      // Electric field [kcal mol⁻¹ Å⁻¹ e⁻¹] → AU
      Vector3d sE = sd1->getElectricField() + EF_;
      sE *= kcalToAU;
      RealType E = dot(ohUnit, sE);

      // Local frequency ω₁₀ [cm⁻¹]
      auto [a0, a1, a2] = wIt->second;
      RealType freq = a0 + a1*E + a2*E*E;

      auto [m0, m1, m2] = mpIt->second;
      RealType muPrime = m0 + m1*E + m2*E*E;        // [a.u.]

      auto [x0, x1] = x10_.at(name);
      RealType x10  = x0 + x1*freq;                 // [a.u.]

      // Momentum matrix element for intramolecular coupling
      auto [pv0, pv1] = p10_.at(name);
      RealType p10  = pv0 + pv1*freq;               // [a.u.]

      // Transition dipole magnitude [Debye]
      const RealType eaBohr_to_Debye = 2.5417464;
      RealType muMag = muPrime * x10 * eaBohr_to_Debye;

      // TDC dipole location: a fitted distance along the OH bond from O.
      // Auer & Skinner 2008: 0.58 Å (SPC/E); Gruenbaum et al. 2013: 0.67 Å (TIP4P).
      RealType tdcDist = 0.5 * rOH_len;  // fallback: midpoint
      auto tdcIt = tdcLoc_.find(name);
      if (tdcIt != tdcLoc_.end())
        tdcDist = tdcIt->second;
      Vector3d dipolePos = Opos + tdcDist * ohUnit;

      // Bond polarizability tensor.
      auto [apar, aperp] = alphaMap_.at(name);
      RealType x10_gas = x10_.at(name).first;

      fd.N++;
      fd.globalIDs.push_back(sd1->getGlobalIndex());
      fd.mu.push_back(muMag * ohUnit);
      fd.alpha.push_back(bondPolarizability(ohUnit, apar, aperp,
					    x10, x10_gas));

      // Auxiliary data for Hamiltonian construction
      ohPos.push_back(dipolePos);
      molIndex.push_back(molID);
      freqs.push_back(freq);
      eFields.push_back(E);
      x10vals.push_back(x10);
      p10vals.push_back(p10);
      mxvals.push_back(muPrime * x10);   // a.u., for TDC coupling

      auto [k0, k1, kp] = wintra_.at(name);
      intraJ.push_back(k0 + k1*E);
      intraJ_k0.push_back(k0);
      intraJ_kp.push_back(kp);

      // -------------------------------------------------------------------
      // Stash this stretch's per-molecule info for the bend pass
      // -------------------------------------------------------------------
      MolData& md = molMap[molID];
      md.Opos    = Opos;
      md.molOgID = OgID;
      md.modelName = mol->getType();
      if (md.ohCount == 0) {
        md.ohUnit1 = ohUnit;
        md.Efield1 = sE;
        md.rOH1    = rOH_len;
      } else if (md.ohCount == 1) {
        md.ohUnit2 = ohUnit;
        md.Efield2 = sE;
        md.rOH2    = rOH_len;
      }
      md.ohCount++;
    }

    int nStretch = fd.N;
    fd.nStretch  = nStretch;

    // -----------------------------------------------------------------------
    // Bend overtone chromophores.
    //
    // Add one bend chromophore per molecule that contributed BOTH OH
    // stretches (we need both bond vectors to define the bend coordinate).
    // The bend electric field is projected on the in-plane perpendicular
    // directions to each OH bond:
    //
    //   E_b = ê_⊥1 · E(H1) / r_OH1  +  ê_⊥2 · E(H2) / r_OH2
    //
    // where ê_⊥1 lies in the molecular plane perpendicular to OH1, in the
    // direction that opens the HOH angle (and similarly for ê_⊥2).  This
    // matches MultiSpec's calcEf() construction:
    //   v_p   = (e1 × e2) / |e1 × e2|       (out-of-plane unit vector)
    //   ê_⊥1 = v_p × e1                      (in-plane, perpendicular to e1)
    //   ê_⊥2 = e2 × v_p                      (in-plane, perpendicular to e2)
    //
    // Bend overtones have zero transition dipole and polarizability (their
    // intensity is borrowed from the stretches via Fermi mixing in the
    // exciton eigenvectors).
    // -----------------------------------------------------------------------
    int nBend = 0;
    if (fc_ != 0.0) {
      for (auto& [molID, md] : molMap) {
        if (md.ohCount != 2) continue;   // need both OH for bend coordinate
        auto bIt = wb01_.find(md.modelName);
        if (bIt == wb01_.end()) continue;
        auto cIt = wb12_.find(md.modelName);
        if (cIt == wb12_.end()) continue;

        // In-plane perpendiculars to each OH (opens HOH angle direction)
        Vector3d vp = cross(md.ohUnit1, md.ohUnit2);
        RealType vpLen = vp.length();
        if (vpLen < 1.0e-9) continue;     // degenerate (parallel OH bonds)
        vp /= vpLen;

        Vector3d ePerp1 = cross(vp, md.ohUnit1);  // ⊥ OH1, in-plane
        Vector3d ePerp2 = cross(md.ohUnit2, vp);  // ⊥ OH2, in-plane

        // Bend field projection in a.u.
        RealType Eb = dot(ePerp1, md.Efield1) / md.rOH1
                    + dot(ePerp2, md.Efield2) / md.rOH2;

        // Diagonal bend overtone frequency: ω_2δ = ω_10 + ω_21
        auto [b0_10, b1_10] = bIt->second;
        auto [b0_21, b1_21] = cIt->second;
        RealType w_2d = (b0_10 + b1_10*Eb) + (b0_21 + b1_21*Eb);

        // Append bend chromophore.  Use the molecule's O global ID as a
        // unique identifier (distinct from any H atom's ID).
        fd.N++;
        fd.globalIDs.push_back(md.molOgID);
        fd.mu.push_back(V3Zero);                // zero transition dipole
        fd.alpha.push_back(Mat3x3d(0.0));       // zero polarizability

        // For the buildHamiltonian loop, intraJ packing must continue
        // for the bend block too (6 entries per chromophore).  Bends
        // have no intramolecular OH-OH coupling and no inter TDC, so
        // we fill placeholders.  buildHamiltonian recognizes bend
        // chromophores by index >= fd.nStretch and applies Fermi
        // coupling instead of the stretch formulas.
        ohPos.push_back(md.Opos);     // dummy position; bend has no TDC
        molIndex.push_back(molID);
        freqs.push_back(w_2d);
        eFields.push_back(Eb);
        x10vals.push_back(0.0);
        p10vals.push_back(0.0);
        mxvals.push_back(0.0);        // mx=0 -> bend cannot TDC-couple
        intraJ.push_back(0.0);
        intraJ_k0.push_back(0.0);
        intraJ_kp.push_back(0.0);

        nBend++;
      }
    }

    fd.nBend  = nBend;
    useFermi_ = (nBend > 0);

    // Build H once with the final size; set diagonal to ω₁₀ frequencies
    // for stretches and ω_2δ for bend overtones.
    int N = fd.N;
    fd.H = DynamicRectMatrix<RealType>(N, N, 0.0);
    for (int i = 0; i < N; ++i) {
      fd.H(i, i) = freqs[i];
      fd.idToIndex[fd.globalIDs[i]] = i;
    }

    // Stash all per-site values for buildHamiltonian via packed intraJ.
    // Each chromophore contributes 6 consecutive entries:
    //   [0] Ki = k0 + k1·E      (intra potential prefactor)
    //   [1] x₁₀                  (a.u.)
    //   [2] p₁₀                  (a.u.)
    //   [3] k0                   (per atom type, for intra)
    //   [4] kp                   (per atom type, for intra)
    //   [5] mx = μ′·x₁₀          (a.u., for inter TDC)
    // For bend chromophores all of these are 0.0, which makes the
    // intra-stretch coupling formula evaluate to 0 and the inter TDC
    // formula also evaluate to 0 (mx = 0).  Fermi coupling between
    // each bend and its two parent stretches is added separately in
    // buildHamiltonian using the molIndex array.
    {
      std::vector<RealType> packed;
      packed.reserve(6 * N);
      for (int i = 0; i < N; ++i) {
        packed.push_back(intraJ[i]);
        packed.push_back(x10vals[i]);
        packed.push_back(p10vals[i]);
        packed.push_back(intraJ_k0[i]);
        packed.push_back(intraJ_kp[i]);
        packed.push_back(mxvals[i]);
      }
      intraJ = std::move(packed);
    }

    return fd;
  }

  // =========================================================================
  // buildHamiltonian
  //
  // Fills fd.H off-diagonal elements only.
  // Diagonal (ω₁₀) was already set by extractFrame.
  //   same molecule  → J from Gruenbaum eq 3.17:
  //     ω^intra_jk = [k0 + k1·(Eⱼ+Eₖ)]·xⱼ·xₖ + kp·pⱼ·pₖ
  //     where k0, kp, and the per-chromophore potential prefactors,
  //     x₁₀, p₁₀ are packed in intraJ (per atom type from wintra_).
  //   different mols → J = TDC coupling computed in MultiSpec convention:
  //     J = AU_TO_WN · geom · (μ′·x₁₀)ᵢ · (μ′·x₁₀)ⱼ / r_bohr³
  //
  // intraJ packing (from extractFrame): 6 entries per chromophore:
  //   [6*i+0] = Ki = k0 + k1·Eᵢ  (potential coupling prefactor)
  //   [6*i+1] = x₁₀,ᵢ            (coordinate matrix element, a.u.)
  //   [6*i+2] = p₁₀,ᵢ            (momentum matrix element, a.u.)
  //   [6*i+3] = k0                (from wintra_ map, per atom type)
  //   [6*i+4] = kp                (from wintra_ map, per atom type)
  //   [6*i+5] = μ′·x₁₀            (a.u., for inter TDC)
  // =========================================================================
  void SFG::buildHamiltonian(FrameData&                   fd,
			     const std::vector<Vector3d>& ohPos,
			     const std::vector<int>&      molIndex,
			     const std::vector<RealType>& intraJ) {
    int N = fd.N;
    if (N == 0) return;

    // H diagonal already holds ω₁₀ (stretch) and ω_2δ (bend) values
    // (set during extractFrame).  Fill off-diagonal elements only.
    //
    // Coupling logic:
    //   same molecule + both stretches  -> Gruenbaum intramolecular formula
    //   same molecule + stretch ⟷ bend  -> Fermi coupling (constant fc_)
    //   different molecules + stretches -> TDC coupling
    //   any pair involving bends across molecules -> zero (mx=0 makes
    //                                                       TDC vanish)
    const int nStr = fd.nStretch;

    for (int i = 0; i < N; ++i) {
      bool i_is_bend = (i >= nStr);
      for (int j = i+1; j < N; ++j) {
	bool j_is_bend = (j >= nStr);
	RealType J = 0.0;
	if (molIndex[i] == molIndex[j]) {
	  if (!i_is_bend && !j_is_bend) {
	    // Stretch-stretch on same molecule: Gruenbaum eq 3.17
	    //   ω^intra_jk = [k0 + k1·(Eⱼ+Eₖ)]·xⱼ·xₖ + kp·pⱼ·pₖ
	    RealType Ki     = intraJ[6*i + 0];
	    RealType xi     = intraJ[6*i + 1];
	    RealType pi     = intraJ[6*i + 2];
	    RealType k0_val = intraJ[6*i + 3];
	    RealType kp_val = intraJ[6*i + 4];
	    RealType Kj     = intraJ[6*j + 0];
	    RealType xj     = intraJ[6*j + 1];
	    RealType pj     = intraJ[6*j + 2];

	    RealType potentialTerm = (Ki + Kj - k0_val) * xi * xj;
	    RealType kineticTerm   = kp_val * pi * pj;
	    J = potentialTerm + kineticTerm;
	  } else if (i_is_bend != j_is_bend) {
	    // Stretch ⟷ bend overtone Fermi coupling on same molecule.
	    // Constant fc_ (default 50 cm⁻¹, Ni & Skinner 2015 / Kananenka 2018).
	    J = fc_;
	  }
	  // else: both bends → impossible (one bend per molecule), J=0
	} else {
	  // Different molecules: TDC coupling (vanishes naturally for
	  // bend chromophores since mx=0 was packed for them).
	  Vector3d r_ij = ohPos[j] - ohPos[i];
	  currentSnapshot_->wrapVector(r_ij);

	  Vector3d e_i = fd.mu[i];
	  RealType m_i = e_i.length();
	  if (m_i > 1.0e-12) e_i /= m_i;

	  Vector3d e_j = fd.mu[j];
	  RealType m_j = e_j.length();
	  if (m_j > 1.0e-12) e_j /= m_j;

	  RealType mx_i = intraJ[6*i + 5];
	  RealType mx_j = intraJ[6*j + 5];

	  J = tdcCoupling(r_ij, e_i, e_j, mx_i, mx_j);
	}
	fd.H(i, j) = J;
	fd.H(j, i) = J;
      }
    }
  }

  // =========================================================================
  // propagateF
  //
  // Adiabatic propagation step (MultiSpec moveF()):
  //   F <- exp(i H dt/hbar) . F
  //
  // exp(i H dt/hbar) = V . diag(exp(i w_k dt/hbar)) . V^T
  //   (V is real orthogonal since H is real symmetric)
  //
  // Steps:
  //   1. Diagonalize H -> eigenvalues W, eigenvectors V
  //      Prefers dsyevd (LAPACK divide-and-conquer) when available;
  //      falls back to JAMA otherwise.
  //   2. Phase factors: u_k = exp(i w_k dt / hbar)
  //   3. tmp  = V^T . F          (real x complex, via cblas_dgemm if available)
  //   4. tmp_k* *= u_k           (row-wise phase scaling)
  //   5. F    = V  . tmp         (real x complex, via cblas_dgemm if available)
  //
  // Work arrays are thread_local static and resized only when N grows,
  // so per-call heap allocation is zero in the steady state.
  // =========================================================================
  void SFG::propagateF(std::vector<std::complex<double>>& F,
		       const DynamicRectMatrix<RealType>&  H,
		       double dt_ps, int N) {

    using cdbl = std::complex<double>;

    // -----------------------------------------------------------------------
    // Thread-local scratch: grown as needed, never shrunk.
    // Vbuf: eigenvectors, ROW-MAJOR: V(i,k) = Vbuf[i*N + k].
    //       Column k is the k-th eigenvector.
    //       This matches MultiSpec's convention (LAPACK_ROW_MAJOR).
    // -----------------------------------------------------------------------
    thread_local std::vector<double> Hbuf;
    thread_local std::vector<double> Wbuf;
    thread_local std::vector<double> Vbuf;
    thread_local std::vector<double> work_d;
    thread_local std::vector<int>    iwork_d;
    thread_local std::vector<double> Fr, Fi, Tr, Ti;
    thread_local int lastN = 0;  // tracks when workspace query is needed

    const int N2 = N * N;
    if (static_cast<int>(Hbuf.size()) < N2) {
      Hbuf.resize(N2); Wbuf.resize(N); Vbuf.resize(N2);
      Fr.resize(N2);   Fi.resize(N2);
      Tr.resize(N2);   Ti.resize(N2);
    }

    // -----------------------------------------------------------------------
    // Step 1: Diagonalize H -> Wbuf (eigenvalues),
    //                          Vbuf (eigenvectors row-major)
    // -----------------------------------------------------------------------
#if defined(HAVE_LAPACK)
    // Copy H's contiguous row-major storage into Hbuf.
    // H is symmetric so row-major == column-major for dsyevd_.
    // dsyevd_ will overwrite Hbuf with eigenvectors (column-major).
    {
      const RealType* src = H.getArrayPointer();
      for (int i = 0; i < N2; ++i)
        Hbuf[i] = static_cast<double>(src[i]);
    }

    // Workspace query — only needed when N changes.
    if (N != lastN) {
      int lwork_q = -1, liwork_q = -1, info_q = 0;
      double  wq  = 0.0;
      int     wiq = 0;
      std::vector<double> Htmp(N2, 0.0);
      std::vector<double> Wtmp(N, 0.0);
      dsyevd_("V", "U", &N, Htmp.data(), &N, Wtmp.data(),
              &wq, &lwork_q, &wiq, &liwork_q, &info_q);
      int lw  = static_cast<int>(wq);
      int liw = wiq;
      if (static_cast<int>(work_d.size())  < lw)  work_d.resize(lw);
      if (static_cast<int>(iwork_d.size()) < liw) iwork_d.resize(liw);
      lastN = N;
    }

    {
      int lw  = static_cast<int>(work_d.size());
      int liw = static_cast<int>(iwork_d.size());
      int info = 0;
      dsyevd_("V", "U", &N, Hbuf.data(), &N, Wbuf.data(),
              work_d.data(), &lw, iwork_d.data(), &liw, &info);
      if (info == 0) {
        // dsyevd_ (Fortran) returns eigenvectors in column-major:
        //   Hbuf[i + k*N] = V(i,k).
        // Transpose into row-major Vbuf: Vbuf[i*N + k] = V(i,k).
        for (int i = 0; i < N; ++i)
          for (int k = 0; k < N; ++k)
            Vbuf[i*N + k] = Hbuf[i + k*N];
        goto after_diag;
      }
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SFG::propagateF: dsyevd returned info=%d; using JAMA.\n", info);
      painCave.isFatal  = 0;
      painCave.severity = OPENMD_WARNING;
      simError();
    }
#endif  // HAVE_LAPACK

    // JAMA fallback
    {
      JAMA::Eigenvalue<RealType>  eig(H);
      DynamicVector<RealType>     evals_j(N);
      DynamicRectMatrix<RealType> evecs_j(N, N);
      eig.getRealEigenvalues(evals_j);
      eig.getV(evecs_j);
      for (int k = 0; k < N; ++k) {
        Wbuf[k] = static_cast<double>(evals_j[k]);
        for (int i = 0; i < N; ++i)
          Vbuf[i*N + k] = static_cast<double>(evecs_j(i, k));
      }
    }

#if defined(HAVE_LAPACK)
  after_diag:
#endif

    // -----------------------------------------------------------------------
    // Step 2: Phase factors u_k = exp(i w_k dt / hbar)
    // -----------------------------------------------------------------------
    std::vector<cdbl> phase(N);
    for (int k = 0; k < N; ++k) {
      double arg = Wbuf[k] * dt_ps / HBAR_CM_PS;
      phase[k] = std::exp(cdbl(0.0, arg));
    }

    // -----------------------------------------------------------------------
    // Steps 3-5: F <- V · diag(phase) · V^T · F
    //
    // All arrays are ROW-MAJOR: M(i,j) = buf[i*N + j].
    // Vbuf is row-major: V(i,k) = Vbuf[i*N + k].
    // F is row-major:    F(i,j) = F[i*N + j].
    //
    // Split F into real/imag parts.
    //
    // Step 3: T = V^T · F      →  T(k,j) = Σ_i V(i,k)·F(i,j)
    // Step 4: T(k,j) *= phase[k]
    // Step 5: F = V · T        →  F(i,j) = Σ_k V(i,k)·T(k,j)
    //
    // For CBLAS with CblasRowMajor:
    //   Step 3: T = V^T * F  → dgemm(RowMajor, Trans, NoTrans, ...)
    //   Step 5: F = V  * T   → dgemm(RowMajor, NoTrans, NoTrans, ...)
    // -----------------------------------------------------------------------

    for (int idx = 0; idx < N2; ++idx) {
      Fr[idx] = F[idx].real();
      Fi[idx] = F[idx].imag();
    }

#if defined(HAVE_CBLAS)
#  if defined(__APPLE__)
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wdeprecated-declarations"
#  endif
    // Step 3: T = V^T * F  (row-major throughout)
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                N, N, N, 1.0, Vbuf.data(), N, Fr.data(), N, 0.0, Tr.data(), N);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                N, N, N, 1.0, Vbuf.data(), N, Fi.data(), N, 0.0, Ti.data(), N);
#  if defined(__APPLE__)
#    pragma clang diagnostic pop
#  endif
#else
    // Manual fallback: T(k,j) = Σ_i V(i,k) * F(i,j)
    std::fill(Tr.begin(), Tr.begin()+N2, 0.0);
    std::fill(Ti.begin(), Ti.begin()+N2, 0.0);
    for (int k = 0; k < N; ++k)
      for (int i = 0; i < N; ++i) {
        double vik = Vbuf[i*N + k];   // V(i,k) row-major
        for (int j = 0; j < N; ++j) {
          Tr[k*N+j] += vik * Fr[i*N+j];
          Ti[k*N+j] += vik * Fi[i*N+j];
        }
      }
#endif

    // Step 4: Apply phase row-wise: T(k,:) *= phase[k]
    for (int k = 0; k < N; ++k) {
      double pr = phase[k].real(), pi = phase[k].imag();
      for (int j = 0; j < N; ++j) {
        double tr = Tr[k*N+j], ti = Ti[k*N+j];
        Tr[k*N+j] = pr*tr - pi*ti;
        Ti[k*N+j] = pr*ti + pi*tr;
      }
    }

#if defined(HAVE_CBLAS)
#  if defined(__APPLE__)
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wdeprecated-declarations"
#  endif
    // Step 5: F = V * T  (row-major throughout)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N, 1.0, Vbuf.data(), N, Tr.data(), N, 0.0, Fr.data(), N);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N, 1.0, Vbuf.data(), N, Ti.data(), N, 0.0, Fi.data(), N);
#  if defined(__APPLE__)
#    pragma clang diagnostic pop
#  endif
#else
    // Manual fallback: F(i,j) = Σ_k V(i,k) * T(k,j)
    std::fill(Fr.begin(), Fr.begin()+N2, 0.0);
    std::fill(Fi.begin(), Fi.begin()+N2, 0.0);
    for (int i = 0; i < N; ++i)
      for (int k = 0; k < N; ++k) {
        double vik = Vbuf[i*N + k];   // V(i,k) row-major
        for (int j = 0; j < N; ++j) {
          Fr[i*N+j] += vik * Tr[k*N+j];
          Fi[i*N+j] += vik * Ti[k*N+j];
        }
      }
#endif

    for (int idx = 0; idx < N2; ++idx)
      F[idx] = cdbl(Fr[idx], Fi[idx]);
  }
  // =========================================================================
  // accumulateTCF
  //
  // Computes the SFG TCF contributions at lag tt for all polarizations:
  //
  //   For pqr (e.g. yyz, ssp):
  //     C_pqr(tt) += Σᵢ  α_pq,i(t) × [F · μ_r(0)]ᵢ  ×  exp_tc
  //
  // where  [F · μ_r(0)]ᵢ = Σⱼ F_ij · μ_r,j(0)   (propagated reference dipole)
  //
  // Mirroring MultiSpec calcSFG():
  //   work = F · mu0           (N-vector, complex)
  //   C_yyz(t) += conj(α_yy(t)) · work_z   ← note MultiSpec uses conj(α)
  //
  // Actually MultiSpec:
  //   tpt[ii] = conj(plz_yy[ii])     ← conjugate of α at time t
  //   result  = tpt · (F · mu0)       ← dot product
  // So:  C_yyz = Σᵢ α_yy,i*(t) · Σⱼ F_ij · μ_z,j(0)
  //
  // For real α (as in our bond-polarizability model), conj(α) = α.
  // But we keep the conjugate for generality.
  // =========================================================================
  void SFG::accumulateTCF(int tt,
			  const std::vector<std::complex<double>>& F,
			  const FrameData& fd0,
			  const FrameData& fdt,
			  int N, double exptc,
			  std::vector<std::complex<double>>& tgt_ssp,
			  std::vector<std::complex<double>>& tgt_ppp,
			  std::vector<std::complex<double>>& tgt_sps) {

    using cdbl = std::complex<double>;
    const int n  = privilegedAxis_;   // interface normal axis index
    const int s1 = s1_;               // first surface-parallel axis
    const int s2 = s2_;               // second surface-parallel axis

    // Fmu[r][i] = sum_j F[i*N+j] * mu_r,j(0)  for r in {n, s1, s2}
    auto matvec = [&](int r) -> std::vector<cdbl> {
      std::vector<cdbl> res(N, cdbl(0.0, 0.0));
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
          res[i] += F[i*N+j] * static_cast<double>(fd0.mu[j][r]);
      return res;
    };

    std::vector<cdbl> Fmu_n  = matvec(n);
    std::vector<cdbl> Fmu_s1 = matvec(s1);
    std::vector<cdbl> Fmu_s2 = matvec(s2);

    // conj(alpha_pq(t)) dotted against propagated dipole component.
    auto dot_alpha = [&](int p, int q, const std::vector<cdbl>& Fmu) -> cdbl {
      cdbl s(0.0, 0.0);
      for (int i = 0; i < N; ++i)
        s += std::conj(cdbl(fdt.alpha[i](p, q))) * Fmu[i];
      return s;
    };

    cdbl decay(exptc, 0.0);

    // SSP: average chi2_{s1,s1,n} and chi2_{s2,s2,n}
    tgt_ssp[tt] += 0.5 * (dot_alpha(s1, s1, Fmu_n)
			  + dot_alpha(s2, s2, Fmu_n)) * decay;

    // SPS: average s1 and s2 channels
    tgt_sps[tt] += 0.5 * (dot_alpha(s1, n, Fmu_s1)
			  + dot_alpha(s2, n, Fmu_s2)) * decay;

    // PPP: near-normal approximation (Buch et al. 2007)
    // MultiSpec: pppt += (-0.5*cxxz - 0.5*cyyz + czzz) * exptc
    tgt_ppp[tt] += (-0.5 * dot_alpha(s1, s1, Fmu_n)
		    - 0.5 * dot_alpha(s2, s2, Fmu_n)
		    + dot_alpha(n,  n,  Fmu_n)) * decay;
  }

  // =========================================================================
  // remapFrame
  //
  // Build a FrameData of size refIDs.size() by looking up each atom in src.
  // Atoms present in src are copied directly.
  // Atoms absent from src (left the z-selection) get:
  //   - zero μ and α  (no contribution to TCF numerator)
  //   - H diagonal from refFreqs (keeps F propagation sensible)
  //   - zero off-diagonal couplings to/from that row/column
  //
  // Returns false if more than half the reference atoms are missing,
  // signalling the window has drifted too far to be meaningful.
  // =========================================================================
  bool SFG::remapFrame(const FrameData&             src,
                       const std::vector<int>&      refIDs,
                       const std::vector<RealType>& refFreqs,
                       FrameData&                   out,
                       int                          refNStretch,
                       int                          refNBend) {
    int N = static_cast<int>(refIDs.size());
    out.N = N;
    // Use the reference-frame's stretch/bend partition.  refIDs is
    // ordered with stretches first, then bends.
    out.nStretch = refNStretch;
    out.nBend    = refNBend;
    out.globalIDs = refIDs;

    out.mu.assign(N, Vector3d(0.0, 0.0, 0.0));
    out.alpha.assign(N, Mat3x3d(0.0));
    out.H = DynamicRectMatrix<RealType>(N, N, 0.0);

    int nMissing = 0;
    for (int i = 0; i < N; ++i) {
      auto it = src.idToIndex.find(refIDs[i]);
      if (it == src.idToIndex.end()) {
        out.H(i, i) = refFreqs[i];
        ++nMissing;
      } else {
        int j = it->second;
        out.mu[i]    = src.mu[j];
        out.alpha[i] = src.alpha[j];
        out.H(i, i) = src.H(j, j);
        for (int k = 0; k < N; ++k) {
          if (k == i) continue;
          auto kt = src.idToIndex.find(refIDs[k]);
          if (kt != src.idToIndex.end()) {
            out.H(i, k) = src.H(j, kt->second);
            out.H(k, i) = src.H(kt->second, j);
          }
        }
      }
    }
    return nMissing <= N / 2;
  }

  // =========================================================================
  // tdcCoupling
  //
  // Transition-dipole coupling between two OH chromophores, computed in
  // MultiSpec's convention (water.cpp::waterTDC + intermC).  The signature
  // takes UNIT OH vectors (eᵢ, eⱼ), the dipole-location separation r_ij
  // (in Å), and the per-site matrix elements (μ′·x₁₀)ᵢ and (μ′·x₁₀)ⱼ
  // in atomic units.
  //
  //   J [cm⁻¹] = AU_TO_WN × [eᵢ·eⱼ − 3(eᵢ·r̂)(eⱼ·r̂)] / r_bohr³
  //              × (μ′·x₁₀)ᵢ × (μ′·x₁₀)ⱼ
  //
  // The geometric part uses unit OH vectors (no dipole magnitudes).  The
  // four field-dependent matrix elements then convert the geometric
  // coupling into actual cm⁻¹.  This keeps maps, transition dipoles, and
  // TDC couplings all in a consistent atomic-unit framework — which is
  // what Auer/Skinner/Gruenbaum/Kananenka use throughout.
  //
  // The previous Debye-based formula (5034.12 × μ_D²/r_Å³) is numerically
  // equivalent in principle, but accumulates conversion factors at each
  // site and is harder to verify against published values.
  // =========================================================================
  RealType SFG::tdcCoupling(const Vector3d& r_ij,
			    const Vector3d& e_i,
			    const Vector3d& e_j,
			    RealType mx_i, RealType mx_j) {
    RealType r_A = r_ij.length();
    if (r_A < 1.0e-6) return 0.0;
    Vector3d rhat = r_ij / r_A;

    // Convert separation to bohr
    const RealType A_to_bohr = 1.0 / 0.52917721067;
    RealType r_bohr = r_A * A_to_bohr;
    RealType r3 = r_bohr * r_bohr * r_bohr;

    // 1 hartree in cm⁻¹
    const RealType AU_TO_WN = 219474.6313705;

    RealType geom = dot(e_i, e_j) - 3.0 * dot(e_i, rhat) * dot(e_j, rhat);

    return AU_TO_WN * geom * mx_i * mx_j / r3;
  }

  // =========================================================================
  // bondPolarizability
  //
  // Transition polarizability tensor for one OH bond.
  //   α_ab = [α_perp·δ_ab + (α_par − α_perp)·ê_a·ê_b] × (x10 / x10_gas)
  //
  // The published α_par/α_perp values incorporate the gas-phase x10.
  // Following MultiSpec's trPol(), we scale by the local x10 (which
  // depends on the OH frequency / electric field) divided by x10_gas
  // (the constant term of the x10 map) so the Raman tensor reflects
  // the field-dependent coordinate matrix element.
  // =========================================================================
  Mat3x3d SFG::bondPolarizability(const Vector3d& ohUnit,
				  RealType alpha_par, RealType alpha_perp,
				  RealType x10, RealType x10_gas) {
    RealType scale = (x10_gas > 0.0) ? x10 / x10_gas : 1.0;
    RealType da = alpha_par - alpha_perp;
    Mat3x3d a(0.0);
    for (int p = 0; p < 3; ++p) {
      a(p, p) = (alpha_perp + da * ohUnit[p] * ohUnit[p]) * scale;
      for (int q = p+1; q < 3; ++q) {
        a(p, q) = da * ohUnit[p] * ohUnit[q] * scale;
        a(q, p) = a(p, q);
      }
    }
    return a;
  }

  // =========================================================================
  // writeCorrelate
  //
  // Writes:
  //   <prefix>.sfg.tcf   — real and imaginary parts of all TCFs vs time
  //   <prefix>.sfg       — SFG spectrum Im χ²(ω) with quantum correction
  //                        for the polarization selected on the command line
  //
  // FFT convention: r2c with normalization by N; frequency axis in cm⁻¹.
  // Quantum correction (Bose-Einstein):
  //   Q(ω) = (ω/kT) / (1 − exp(−ω/kT))
  // =========================================================================
  void SFG::writeCorrelate() {
    Revision r;
  
    // -----------------------------------------------------------------------
    // Write time-domain TCFs for diagnostics
    // -----------------------------------------------------------------------
  
    std::string tcfName = outputFilename_ + ".tcf";
    std::ofstream tcf(tcfName.c_str());
    if (!tcf.is_open()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	       "SFG::writeCorrelate: cannot open %s\n", tcfName.c_str());
      painCave.isFatal = 1; simError();
    }
    tcf << "# SFG time-domain TCFs\n";
    tcf << "# OpenMD " << r.getFullRevision() << "\n";
    tcf << "# " << r.getBuildDate() << "\n";
    tcf << "# selection: \"" << selectionScript1_ << "\"\n";
    tcf << "# navg = " << navg_ << "  ncor = " << nTimeBins_ << "\n";
    tcf << "#t(ps)\tRe_ssp\tIm_ssp\tRe_ppp\tIm_ppp\tRe_sps\tIm_sps\n";

    for (int tt = 0; tt < static_cast<int>(nTimeBins_); ++tt) {
      tcf << tt * dt_ps_ << "\t"
	  << tcf_ssp_[tt].real() << "\t" << tcf_ssp_[tt].imag() << "\t"
	  << tcf_ppp_[tt].real() << "\t" << tcf_ppp_[tt].imag() << "\t"
	  << tcf_sps_[tt].real() << "\t" << tcf_sps_[tt].imag() << "\n";
    }
    tcf.close();
  
#if defined(HAVE_FFTW3_H)

    const double ps_to_invcm = 33.3564;

    // Select TCF based on polarization argument
    const std::vector<std::complex<double>>* pTCF = &tcf_ssp_;
    if (polarization_ == "ppp") pTCF = &tcf_ppp_;
    else if (polarization_ == "sps") pTCF = &tcf_sps_;

    // N_corr: number of real TCF points (determines intrinsic resolution).
    // N_fft:  FFT size after zero-filling (determines grid spacing).
    // The Hermitian extension requires N_fft >= 2*N_corr.
    int N_corr = static_cast<int>(nTimeBins_);
    if (N_corr < 2) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	       "SFG::writeCorrelate: too few TCF points (%d) for FFT.\n", N_corr);
      painCave.isFatal = 1; simError();
    }

    // Zero-filling: extend to T_corr + t_zerofill_ps_, rounded up to next
    // power of two for FFTW efficiency.  Minimum is 2*N_corr for the
    // Hermitian extension (C(-t) = C*(t) occupies the upper half).
    int N_zf = (t_zerofill_ps_ > 0.0)
      ? static_cast<int>(std::ceil((static_cast<double>(N_corr)
				    + t_zerofill_ps_ / dt_ps_)))
      : N_corr;
    N_zf = std::max(N_zf, 2 * N_corr);
    // Round up to next power of two
    int N_fft = 1;
    while (N_fft < N_zf) N_fft <<= 1;

    fftw_complex* in  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_fft);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_fft);
    if (!in || !out) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	       "SFG: FFTW malloc failed.\n");
      painCave.isFatal = 1; simError();
    }

    // Pre-apply apodization to the TCF data and store in a work vector.
    std::vector<std::complex<double>> apodTCF(N_corr);
    for (int i = 0; i < N_corr; ++i) {
      double t = static_cast<double>(i) * dt_ps_;
      double w = (t_apod_ps_ > 0.0) ? std::exp(-t / t_apod_ps_) : 1.0;
      apodTCF[i] = (*pTCF)[i] * w;
    }

    // Result vector for the two-pass Hermitian FFT
    std::vector<std::complex<double>> spectrum(N_fft);

    fftw_plan plan = fftw_plan_dft_1d(N_fft, in, out,
				      FFTW_FORWARD, FFTW_ESTIMATE);

    // -----------------------------------------------------------------------
    // Two-pass Hermitian FFT (following MultiSpec dofft.cpp).
    //
    // The one-sided TCF C(t) for t >= 0 is extended to negative times
    // using the physical symmetry C(-t) = C*(t).  This Hermitian
    // extension ensures that the real and imaginary parts of the
    // Fourier transform are cleanly separated into dispersive (Re)
    // and absorptive (Im) lineshapes with no phase mixing.
    //
    // Pass 1:  input = C(t) for t>0, C*(-t) for t<0  →  Re(spectrum)
    // Pass 2:  input = i·C(t) for t>0, conj(i·C(t)) for t<0  →  Im(spectrum)
    // -----------------------------------------------------------------------

    // Pass 1: get Re(chi2)
    for (int i = 0; i < N_fft; ++i) { in[i][0] = 0.0; in[i][1] = 0.0; }
    for (int i = 0; i < N_corr; ++i) {
      in[i][0] = apodTCF[i].real();
      in[i][1] = apodTCF[i].imag();
      if (i > 0) {
	// Hermitian mirror: C(-t) = C*(t) at index N_fft - i
	in[N_fft - i][0] =  apodTCF[i].real();
	in[N_fft - i][1] = -apodTCF[i].imag();
      }
    }
    fftw_execute(plan);
    // The Hermitian input guarantees the output is purely real;
    // take the real part as Re(spectrum).
    for (int i = 0; i < N_fft; ++i)
      spectrum[i].real(out[i][0]);

    // Pass 2: get Im(chi2)
    // Multiply TCF by -i: -i·C = Im(C) - i·Re(C)
    // The sign convention ensures that Im(chi2) > 0 for the H-bonded
    // OH stretch region at the water/vapor interface (SSP), consistent
    // with experiment and the one-sided FT convention.
    for (int i = 0; i < N_fft; ++i) { in[i][0] = 0.0; in[i][1] = 0.0; }
    for (int i = 0; i < N_corr; ++i) {
      double re_ic =  apodTCF[i].imag();
      double im_ic = -apodTCF[i].real();
      in[i][0] = re_ic;
      in[i][1] = im_ic;
      if (i > 0) {
	// Hermitian mirror of -i·C(t)
	in[N_fft - i][0] =  re_ic;
	in[N_fft - i][1] = -im_ic;
      }
    }
    fftw_execute(plan);
    for (int i = 0; i < N_fft; ++i)
      spectrum[i].imag(out[i][0]);

    fftw_destroy_plan(plan);

    // Normalization: we normalize by 2*N_corr — the factor of N_corr keeps
    // amplitudes independent of zero-filling, and the factor of 2 accounts
    // for the Hermitian extension which places C(t) at both positive and
    // negative times, doubling the effective signal relative to the one-sided
    // transform.
    const double norm = 1.0 / (2.0 * static_cast<double>(N_corr));
    for (int i = 0; i < N_fft; ++i)
      spectrum[i] *= norm;

    // Frequency grid
    const double dOmega      = ps_to_invcm / (static_cast<double>(N_fft) * dt_ps_);
    const double dOmega_res  = ps_to_invcm / (static_cast<double>(N_corr) * dt_ps_);

    // Temperature
    RealType T = 300.0;
    if (info_->getSimParams()->haveTargetTemp())
      T = info_->getSimParams()->getTargetTemp();
    const double kT_invcm = 0.6950356 * static_cast<double>(T);

    std::ofstream ofs(outputFilename_.c_str());
    if (!ofs.is_open()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	       "SFG: cannot open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1; simError();
    }

    ofs << "# SFG spectrum (" << polarization_ << ")\n";
    ofs << "# OpenMD " << r.getFullRevision() << "\n";
    ofs << "# " << r.getBuildDate() << "\n";
    ofs << "# selection: \"" << selectionScript1_ << "\"\n";
    ofs << "# T = " << T << " K,  kT = " << kT_invcm << " cm-1\n";
    ofs << "# wAvg = " << wAvg_ << " cm-1\n";
    ofs << "# navg = " << navg_ << "  nTimeBins = " << nTimeBins_ << "\n";
    ofs << "# intrinsic resolution = " << dOmega_res << " cm-1";
    if (t_apod_ps_ > 0.0)
      ofs << "  t_apod = " << t_apod_ps_ << " ps"
          << "  (added broadening ~"
          << ps_to_invcm / (M_PI * t_apod_ps_) << " cm-1)";
    ofs << "\n";
    ofs << "# dOmega = " << dOmega << " cm-1";
    if (t_zerofill_ps_ > 0.0)
      ofs << "  (zero-filled: N_corr=" << N_corr
          << " -> N_fft=" << N_fft << ")";
    ofs << "\n";
    ofs << "#omega(cm-1)\tRe(chi2)\tIm(chi2)\t|chi2|^2\n";

    // Output frequency range: wAvg_ ± (N_corr/2)*dOmega_res from wAvg.
    // Convert intrinsic bandwidth to bin count on the zero-filled grid.
    const double halfBW_freq = (N_corr / 2) * dOmega_res;
    const int halfBW = static_cast<int>(std::round(halfBW_freq / dOmega));

    auto writeRow = [&](double omega, int k) {
      if (omega <= 0.0) return;
      double re = spectrum[k].real();
      double im = spectrum[k].imag();
      double x  = omega / kT_invcm;
      double qCorr = (x > 1.0e-6) ?
	((x < 50.0) ? x / (1.0 - std::exp(-x)) : x) : 1.0;
      re *= qCorr;
      im *= qCorr;
      ofs << omega << "\t" << re << "\t" << im << "\t"
	  << re*re + im*im << "\n";
    };

    // negative half → ω below wAvg_ (ascending)
    for (int k = N_fft - halfBW; k < N_fft; ++k)
      writeRow(wAvg_ + (k - N_fft) * dOmega, k);
    // positive half → ω at/above wAvg_
    for (int k = 0; k <= halfBW; ++k)
      writeRow(wAvg_ + k * dOmega, k);
    
    ofs.close();
    fftw_free(in);
    fftw_free(out);

#else
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	     "SFG: FFTW3 required for SFG spectra.\n");
    painCave.isFatal = 1; simError();
#endif
  }
}  // namespace OpenMD
