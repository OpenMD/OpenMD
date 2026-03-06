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
           const std::string& polarization, int privilegedAxis) :
    SystemACF<RealType>(info, filename, sele1, sele2),
    polarization_(polarization),
    privilegedAxis_(privilegedAxis),
    s1_((privilegedAxis + 1) % 3),
    s2_((privilegedAxis + 2) % 3) {

    setCorrFuncType("SFG");
    setOutputName(getPrefix(dumpFilename_) + ".sfg");
    setLabelString("Im chi2(omega)");

    // -----------------------------------------------------------------------
    // Spectroscopic maps.
    // w10_   : ω₁₀(E) = a0 + a1·E + a2·E²    [cm⁻¹;  E in a.u.]
    // p10_   : |μ₁₀|(E) = c0 + c1·E           [Debye]
    //          Transition dipole MOMENT for TDC coupling + collective μ.
    // wintra_: J_intra(E) = d0 + d1·E + d2·E² [cm⁻¹]
    // alphaMap_: (α_par, α_perp)               [Å² amu⁻⁰·⁵]
    //
    // Sources: Auer & Skinner, J. Chem. Phys. 128, 224511 (2008) — SPC/E
    //          Gruenbaum et al., JCTC 9, 3109 (2013) — TIP4P
    //          Takayama et al., J. Chem. Phys. 158, 136101 (2023) — TIP4P-Ice
    // -----------------------------------------------------------------------

    // SPC/E
    w10_["H_SPCE"]      = std::make_tuple(3761.6,  -5060.4,  -86225.0);
    p10_["H_SPCE"]      = std::make_pair(1.611,    5.893e-4);
    wintra_["H_SPCE"]   = std::make_tuple(-1789.0,  23852.0,    -1.966);
    alphaMap_["H_SPCE"] = std::make_pair(0.185, 0.033);

    // TIP4P
    w10_["H_TIP4P"]      = std::make_tuple(3760.2,  -3541.7, -152677.0);
    p10_["H_TIP4P"]      = std::make_pair(1.6466,   5.7692e-4);
    wintra_["H_TIP4P"]   = std::make_tuple(-1361.0,  27165.0,    -1.887);
    alphaMap_["H_TIP4P"] = std::make_pair(0.185, 0.033);

    // TIP4P-Ice (transferability from TIP4P)
    w10_["H_TIP4P-Ice"]      = w10_["H_TIP4P"];
    p10_["H_TIP4P-Ice"]      = p10_["H_TIP4P"];
    wintra_["H_TIP4P-Ice"]   = wintra_["H_TIP4P"];
    alphaMap_["H_TIP4P-Ice"] = alphaMap_["H_TIP4P"];

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
      int asl  = info_->getAtomStorageLayout()        | DataStorage::dslElectricField;
      int rbsl = info_->getRigidBodyStorageLayout()   | DataStorage::dslElectricField;
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
    // buildHamiltonian leaves H diagonal as unshifted ω₁₀ values at this point.
    const FrameData& fd = allFrames_[frame];
    for (int i = 0; i < fd.N; ++i) {
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

    int    ncor   = static_cast<int>(nTimeBins_);
    double dt_ps  = static_cast<double>(dtMean_) * 1.0e-3;  // fs → ps
    dt_ps_        = dt_ps;

    // Finalize wAvg_ from all per-frame diagonal entries accumulated during
    // preCorrelate(). This is the true mean ω₁₀ across the trajectory and
    // the selected chromophores, regardless of water model.
    wAvg_ = (wAvgCount_ > 0)
      ? wAvgSum_ / static_cast<double>(wAvgCount_)
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

    for (int iwin = 0; iwin < navg_; ++iwin) {
      int istart = nStart_ + iwin * nStride_;
      if (istart >= nFrames_) break;

      const FrameData& fd0raw = allFrames_[istart];
      int N = fd0raw.N;
      if (N == 0) continue;

      // Fix the atom set for this entire window at the origin frame's
      // selection. As molecules diffuse across the z-boundary, remapFrame()
      // maps every subsequent frame onto this fixed set so F stays valid.
      const std::vector<int>& refIDs = fd0raw.globalIDs;

      // Seed reference frequencies from the origin frame's shifted diagonal.
      // Updated each lag step so atoms that temporarily leave the selection
      // get a fresh fallback frequency when they return.
      std::vector<RealType> refFreqs(N);
      for (int i = 0; i < N; ++i)
        refFreqs[i] = fd0raw.H(i, i);

      // F = identity
      std::vector<std::complex<double>> F(N * N, czero);
      for (int i = 0; i < N; ++i) F[i*N+i] = std::complex<double>(1.0, 0.0);

      FrameData fd0;
      remapFrame(fd0raw, refIDs, refFreqs, fd0);

      for (int tt = 0; tt < ncor; ++tt) {
        int iframe = istart + tt;
        if (iframe >= nFrames_) break;

        FrameData fdt;
        if (!remapFrame(allFrames_[iframe], refIDs, refFreqs, fdt))
          break;  // too many atoms missing; abandon this window

        // Update refFreqs for atoms present in this frame.
        for (int i = 0; i < N; ++i) {
          auto it = allFrames_[iframe].idToIndex.find(refIDs[i]);
          if (it != allFrames_[iframe].idToIndex.end())
            refFreqs[i] = allFrames_[iframe].H(it->second, it->second);
        }

        // Propagate F using H at the current frame (matches MultiSpec:
        // readHf advances to current frame, moveF() uses it, then calcSFG
        // accumulates).
        if (tt > 0)
          propagateF(F, fdt.H, dt_ps, N);

        double exptc = (T1_ps_ > 0.0)
          ? std::exp(-static_cast<double>(tt) * dt_ps / (2.0 * T1_ps_))
          : 1.0;

        accumulateTCF(tt, F, fd0, fdt, N, exptc);
      }
    }

    // Normalize by number of windows
    if (navg_ > 0) {
      double inv = 1.0 / static_cast<double>(navg_);
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
  // basis.  The Hamiltonian matrix is NOT filled here; call buildHamiltonian
  // afterward.  ohPos and molIndex are returned for Hamiltonian construction.
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
      auto pIt = p10_.find(name);
      if (pIt == p10_.end()) continue;

      int sdID  = sd1->getGlobalIndex();
      int molID = info_->getGlobalMolMembership(sdID);
      Molecule* mol = info_->getMoleculeByGlobalIndex(molID);

      Vector3d Opos;
      bool foundO = false;
      Molecule::AtomIterator ai;
      for (Atom* oatom = mol->beginAtom(ai); oatom != nullptr;
           oatom = mol->nextAtom(ai)) {
        if (oatom->getAtomType()->getName()[0] == 'O') {
          Opos   = oatom->getPos();
          foundO = true;
          break;
        }
      }
      if (!foundO) continue;

      Vector3d Hpos    = sd1->getPos();
      Vector3d rOH_mic = Hpos - Opos;
      currentSnapshot_->wrapVector(rOH_mic);   // minimum-image OH vector
      Vector3d ohUnit  = rOH_mic;
      ohUnit.normalize();

      // Electric field [kcal mol⁻¹ Å⁻¹ e⁻¹] → AU
      Vector3d sE = sd1->getElectricField() + EF_;
      sE *= kcalToAU;
      RealType E = dot(ohUnit, sE);

      // Local frequency ω₁₀ [cm⁻¹]
      auto [a0, a1, a2] = wIt->second;
      RealType freq = a0 + a1*E + a2*E*E;

      // Transition dipole moment |μ₁₀| [D] from p10_ map
      auto [c0, c1] = pIt->second;
      RealType muMag = c0 + c1*E;

      // Intramolecular coupling [cm⁻¹]
      auto [d0, d1, d2] = wintra_.at(name);
      RealType Jintra = d0 + d1*E + d2*E*E;

      // OH midpoint using minimum-image vector (correct across PBC)
      Vector3d midpoint = Opos + 0.5 * rOH_mic;

      // Bond polarizability tensor
      auto [apar, aperp] = alphaMap_.at(name);

      fd.N++;
      fd.globalIDs.push_back(sd1->getGlobalIndex());
      fd.mu.push_back(muMag * ohUnit);
      fd.alpha.push_back(bondPolarizability(ohUnit, apar, aperp));

      // Auxiliary data for Hamiltonian construction
      ohPos.push_back(midpoint);
      molIndex.push_back(molID);
      intraJ.push_back(Jintra);
      freqs.push_back(freq);
    }

    // Build H once with the final size; set diagonal to ω₁₀ frequencies.
    // Off-diagonals are filled by buildHamiltonian.
    int N = fd.N;
    fd.H = DynamicRectMatrix<RealType>(N, N, 0.0);
    for (int i = 0; i < N; ++i) {
      fd.H(i, i) = freqs[i];
      fd.idToIndex[fd.globalIDs[i]] = i;
    }

    return fd;
  }

  // =========================================================================
  // buildHamiltonian
  //
  // Fills fd.H off-diagonal elements only.
  // Diagonal (ω₁₀) was already set by extractFrame.
  //   same molecule  → J = average of the two intramolecular couplings
  //   different mols → J = TDC coupling between OH midpoints
  // =========================================================================
  void SFG::buildHamiltonian(FrameData&                   fd,
			     const std::vector<Vector3d>& ohPos,
			     const std::vector<int>&      molIndex,
			     const std::vector<RealType>& intraJ) {
    int N = fd.N;
    if (N == 0) return;

    // H diagonal already holds ω₁₀ values (set during extractFrame).
    // Fill off-diagonal elements only.
    for (int i = 0; i < N; ++i) {
      for (int j = i+1; j < N; ++j) {
	RealType J = 0.0;
	if (molIndex[i] == molIndex[j]) {
	  J = 0.5 * (intraJ[i] + intraJ[j]);
	} else {
	  Vector3d r_ij = ohPos[j] - ohPos[i];
	  currentSnapshot_->wrapVector(r_ij);
	  J = tdcCoupling(r_ij, fd.mu[i], fd.mu[j]);
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
  //   3. tmp  = V^T . F                  (real x complex, via cblas_dgemm if available)
  //   4. tmp_k* *= u_k                   (row-wise phase scaling)
  //   5. F    = V  . tmp                 (real x complex, via cblas_dgemm if available)
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
    // Vbuf: eigenvectors, column-major (V(i,k) = Vbuf[i + k*N]).
    // -----------------------------------------------------------------------
    thread_local std::vector<double> Hbuf;
    thread_local std::vector<double> Wbuf;
    thread_local std::vector<double> Vbuf;
    thread_local std::vector<double> work_d;
    thread_local std::vector<int>    iwork_d;
    thread_local std::vector<double> Fr, Fi, Tr, Ti;

    const int N2 = N * N;
    if (static_cast<int>(Hbuf.size()) < N2) {
      Hbuf.resize(N2); Wbuf.resize(N); Vbuf.resize(N2);
      Fr.resize(N2);   Fi.resize(N2);
      Tr.resize(N2);   Ti.resize(N2);
    }

    // -----------------------------------------------------------------------
    // Step 1: Diagonalize H -> Wbuf (eigenvalues), Vbuf (eigenvectors col-major)
    // -----------------------------------------------------------------------
#if defined(HAVE_LAPACK)
    // Copy H into Hbuf column-major (symmetric, so row == col order)
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        Hbuf[i + j*N] = static_cast<double>(H(i, j));

    // Workspace query
    {
      int lwork_q = -1, liwork_q = -1, info = 0;
      double  wq = 0.0; int wiq = 0;
      dsyevd_("V", "U", &N, Hbuf.data(), &N, Wbuf.data(),
              &wq, &lwork_q, &wiq, &liwork_q, &info);
      int lw = static_cast<int>(wq), liw = wiq;
      if (static_cast<int>(work_d.size())  < lw)  work_d.resize(lw);
      if (static_cast<int>(iwork_d.size()) < liw) iwork_d.resize(liw);
    }

    {
      int lw = static_cast<int>(work_d.size());
      int liw = static_cast<int>(iwork_d.size());
      int info = 0;
      dsyevd_("V", "U", &N, Hbuf.data(), &N, Wbuf.data(),
              work_d.data(), &lw, iwork_d.data(), &liw, &info);
      if (info == 0) {
        Vbuf = Hbuf;   // dsyevd_ overwrites Hbuf with eigenvectors
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
      JAMA::Eigenvalue<RealType> eig(H);
      DynamicVector<RealType>     evals_j(N);
      DynamicRectMatrix<RealType> evecs_j(N, N);
      eig.getRealEigenvalues(evals_j);
      eig.getV(evecs_j);
      for (int k = 0; k < N; ++k) {
        Wbuf[k] = static_cast<double>(evals_j[k]);
        for (int i = 0; i < N; ++i)
          Vbuf[i + k*N] = static_cast<double>(evecs_j(i, k));
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
    // Steps 3-5: F <- V . diag(phase) . V^T . F
    //
    // Split F into real/imag parts so that real BLAS dgemm can be used.
    // Vbuf is column-major: V(i,k) = Vbuf[i + k*N].
    //
    // Step 3: [Tr,Ti] = V^T * [Fr,Fi]
    //   In column-major BLAS: T = V^T * F  =>  cblas_dgemm(Trans, NoTrans)
    // Step 4: row k of [Tr,Ti] *= phase[k]
    // Step 5: [Fr,Fi] = V * [Tr,Ti]
    //   In column-major BLAS: F = V * T   =>  cblas_dgemm(NoTrans, NoTrans)
    // -----------------------------------------------------------------------

    for (int idx = 0; idx < N2; ++idx) {
      Fr[idx] = F[idx].real();
      Fi[idx] = F[idx].imag();
    }

#if defined(HAVE_CBLAS)
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                N, N, N, 1.0, Vbuf.data(), N, Fr.data(), N, 0.0, Tr.data(), N);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                N, N, N, 1.0, Vbuf.data(), N, Fi.data(), N, 0.0, Ti.data(), N);
#else
    std::fill(Tr.begin(), Tr.begin()+N2, 0.0);
    std::fill(Ti.begin(), Ti.begin()+N2, 0.0);
    for (int k = 0; k < N; ++k)
      for (int i = 0; i < N; ++i) {
        double vik = Vbuf[i + k*N];
        for (int a = 0; a < N; ++a) {
          Tr[k*N+a] += vik * Fr[i*N+a];
          Ti[k*N+a] += vik * Fi[i*N+a];
        }
      }
#endif

    // Apply phase row-wise
    for (int k = 0; k < N; ++k) {
      double pr = phase[k].real(), pi = phase[k].imag();
      for (int a = 0; a < N; ++a) {
        double tr = Tr[k*N+a], ti = Ti[k*N+a];
        Tr[k*N+a] = pr*tr - pi*ti;
        Ti[k*N+a] = pr*ti + pi*tr;
      }
    }

#if defined(HAVE_CBLAS)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                N, N, N, 1.0, Vbuf.data(), N, Tr.data(), N, 0.0, Fr.data(), N);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                N, N, N, 1.0, Vbuf.data(), N, Ti.data(), N, 0.0, Fi.data(), N);
#else
    std::fill(Fr.begin(), Fr.begin()+N2, 0.0);
    std::fill(Fi.begin(), Fi.begin()+N2, 0.0);
    for (int i = 0; i < N; ++i)
      for (int k = 0; k < N; ++k) {
        double vik = Vbuf[i + k*N];
        for (int a = 0; a < N; ++a) {
          Fr[i*N+a] += vik * Tr[k*N+a];
          Fi[i*N+a] += vik * Ti[k*N+a];
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
			  int N, double exptc) {

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
    tcf_ssp_[tt] += 0.5 * (dot_alpha(s1, s1, Fmu_n)
                         + dot_alpha(s2, s2, Fmu_n)) * decay;

    // SPS: average s1 and s2 channels
    tcf_sps_[tt] += 0.5 * (dot_alpha(s1, n, Fmu_s1)
                         + dot_alpha(s2, n, Fmu_s2)) * decay;

    // PPP: near-normal approximation (Buch et al. 2007)
    tcf_ppp_[tt] += (-dot_alpha(s1, s1, Fmu_n)
                    - dot_alpha(s2, s2, Fmu_n)
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
                       FrameData&                   out) {
    int N = static_cast<int>(refIDs.size());
    out.N = N;
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
  // J_ij = 5034.12 [cm⁻¹ D⁻² Å³] × (μᵢ·μⱼ − 3(μᵢ·r̂)(μⱼ·r̂)) / r³
  // =========================================================================
  RealType SFG::tdcCoupling(const Vector3d& r_ij,
			    const Vector3d& mu_i, const Vector3d& mu_j) {
    RealType r = r_ij.length();
    if (r < 1.0e-6) return 0.0;
    Vector3d rhat = r_ij / r;
    const RealType prefactor = 5034.12;
    return prefactor * (dot(mu_i, mu_j)
			- 3.0 * dot(mu_i, rhat) * dot(mu_j, rhat))
      / (r * r * r);
  }

  // =========================================================================
  // bondPolarizability
  // α_ab = α_perp·δ_ab + (α_par − α_perp)·ê_a·ê_b
  // =========================================================================
  Mat3x3d SFG::bondPolarizability(const Vector3d& ohUnit,
                                   RealType alpha_par, RealType alpha_perp) {
    RealType da = alpha_par - alpha_perp;
    Mat3x3d a(0.0);
    for (int p = 0; p < 3; ++p) {
      a(p, p) = alpha_perp + da * ohUnit[p] * ohUnit[p];
      for (int q = p+1; q < 3; ++q) {
        a(p, q) = da * ohUnit[p] * ohUnit[q];
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

    // Select TCF based on polarization argument
    const std::vector<std::complex<double>>* pTCF = &tcf_ssp_;
    if (polarization_ == "ppp") pTCF = &tcf_ppp_;
    else if (polarization_ == "sps") pTCF = &tcf_sps_;

    int N = static_cast<int>(nTimeBins_);
    if (N < 2) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	       "SFG::writeCorrelate: too few TCF points (%d) for FFT.\n", N);
      painCave.isFatal = 1; simError();
    }

    // FFTW: real input (take Im part of the complex TCF, which is the
    // odd/causal part; alternatively take the full complex TCF with c2c).
    // MultiSpec writes Im of the complex FFT as the spectrum.
    // Use c2c FFT on the full complex TCF.
    fftw_complex* in  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    if (!in || !out) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	       "SFG: FFTW malloc failed.\n");
      painCave.isFatal = 1; simError();
    }

    for (int i = 0; i < N; ++i) {
      in[i][0] = (*pTCF)[i].real();
      in[i][1] = (*pTCF)[i].imag();
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Frequency axis using stored dt_ps_
    const double ps_to_invcm = 33.3564;
    const double dOmega = ps_to_invcm / (static_cast<double>(N) * dt_ps_);

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
    ofs << "# navg = " << navg_ << "  nTimeBins = " << nTimeBins_
	<< "  dOmega = " << dOmega << " cm-1\n";
    ofs << "#omega(cm-1)\tRe(chi2)\tIm(chi2)\t|chi2|^2\n";

    // The FFT of the shifted TCF C'(t) is centered at k=0 = wAvg_.
    // Bin k → absolute frequency: wAvg_ + k*dOmega
    // Bins N/2+1..N-1 are the negative side (ω < wAvg_).
    // Output in ascending frequency order: negative half first, then positive.
    auto writeRow = [&](double omega, int k) {
      if (omega <= 0.0) return;
      double re = out[k][0] / static_cast<double>(N);
      double im = out[k][1] / static_cast<double>(N);
      double x  = omega / kT_invcm;
      double qCorr = (x > 1.0e-6) ?
        ((x < 50.0) ? x / (1.0 - std::exp(-x)) : x) : 1.0;
      double im_corr = im * qCorr;
      ofs << omega << "\t" << re << "\t" << im_corr << "\t"
	  << re*re + im_corr*im_corr << "\n";
    };

    for (int k = N/2 + 1; k < N; ++k)   // negative half → ω below wAvg_
      writeRow(wAvg_ + (k - N) * dOmega, k);
    for (int k = 0; k <= N/2; ++k)       // positive half → ω at/above wAvg_
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
