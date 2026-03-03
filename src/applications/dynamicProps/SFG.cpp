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
#include <fstream>
#include <cmath>

#include "brains/ForceManager.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "math/Eigenvalue.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)

#ifdef HAVE_FFTW_H
#include <fftw.h>
#endif

#ifdef HAVE_DFFTW_H
#include <dfftw.h>
#endif

#ifdef HAVE_FFTW3_H
#include <fftw3.h>
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
#endif
#endif

namespace OpenMD {

  // =========================================================================
  // Constructor
  // =========================================================================
  SFG::SFG(SimInfo* info, const std::string& filename,
           const std::string& sele1, const std::string& sele2,
           const std::string& polarization) :
    SystemACF<RealType>(info, filename, sele1, sele2),
    polarization_(polarization) {

    setCorrFuncType("SFG");
    setOutputName(getPrefix(dumpFilename_) + ".sfg");
    setLabelString("Im chi2(omega)");

    // -----------------------------------------------------------------------
    // Spectroscopic maps — identical parameter sets to OHFrequencyMap.
    // Source: Auer & Skinner, J. Chem. Phys. 128, 224511 (2008) for SPC/E;
    //         Gruenbaum et al., J. Chem. Theory Comput. 9, 3109 (2013) for TIP4P.
    //
    // w10_     : ω₁₀(E) = a0 + a1·E + a2·E²    [cm⁻¹; E in atomic units]
    // muprime_ : |dμ/dq|(E) = c0 + c1·E + c2·E² [D Å⁻¹ amu⁻⁰·⁵]
    // wintra_  : J_intra(E) [cm⁻¹] — field-dependent intramolecular coupling
    // alphaMap_: (α_par, α_perp) bond polarizability derivatives [Å² amu⁻⁰·⁵]
    //            Values from Auer & Skinner 2008 Table I.
    // -----------------------------------------------------------------------

    // SPC/E
    w10_["H_SPCE"]      = std::make_tuple(3761.6,  -5060.4,  -86225.0);
    muprime_["H_SPCE"]  = std::make_tuple(0.71116,  75.591,       0.0);
    wintra_["H_SPCE"]   = std::make_tuple(-1789.0,  23852.0,    -1.966);
    alphaMap_["H_SPCE"] = std::make_pair(0.185, 0.033);   // Å² amu⁻⁰·⁵

    // TIP4P
    w10_["H_TIP4P"]      = std::make_tuple(3760.2,  -3541.7, -152677.0);
    muprime_["H_TIP4P"]  = std::make_tuple(0.1646,     11.39,     63.41);
    wintra_["H_TIP4P"]   = std::make_tuple(-1361.0,  27165.0,    -1.887);
    alphaMap_["H_TIP4P"] = std::make_pair(0.185, 0.033);

    // TIP4P-Ice (transferability: Takayama et al., J. Chem. Phys. 158, 136101 (2023))
    w10_["H_TIP4P-Ice"]      = w10_["H_TIP4P"];
    muprime_["H_TIP4P-Ice"]  = muprime_["H_TIP4P"];
    wintra_["H_TIP4P-Ice"]   = wintra_["H_TIP4P"];
    alphaMap_["H_TIP4P-Ice"] = alphaMap_["H_TIP4P"];

    // -----------------------------------------------------------------------
    // Applied external electric field
    // -----------------------------------------------------------------------
    EF_ = V3Zero;
    const RealType chrgToKcal = 23.060548;
    std::vector<RealType> ef;
    bool efSpec = false;

    if (info_->getSimParams()->haveElectricField()) {
      efSpec = true;
      ef     = info_->getSimParams()->getElectricField();
    }
    if (info_->getSimParams()->haveUniformField()) {
      efSpec = true;
      ef     = info_->getSimParams()->getUniformField();
    }
    if (efSpec) {
      if (ef.size() != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "SFG: uniformField must have 3 components (%zu given).\n",
                 ef.size());
        painCave.isFatal = 1;
        simError();
      }
      EF_.x() = ef[0] * chrgToKcal;
      EF_.y() = ef[1] * chrgToKcal;
      EF_.z() = ef[2] * chrgToKcal;
    }

    // -----------------------------------------------------------------------
    // Electric field storage — mirror OHFrequencyMap logic exactly.
    // If the dump lacks electric fields, expand storage and build a
    // ForceManager so we can compute them on-the-fly.
    // -----------------------------------------------------------------------
    dumpHasElectricFields_ = info_->getSimParams()->getOutputElectricField();

    if (!dumpHasElectricFields_) {
      int atomStorageLayout        = info_->getAtomStorageLayout();
      int rigidBodyStorageLayout   = info_->getRigidBodyStorageLayout();
      int cutoffGroupStorageLayout = info_->getCutoffGroupStorageLayout();

      atomStorageLayout        |= DataStorage::dslElectricField;
      rigidBodyStorageLayout   |= DataStorage::dslElectricField;

      info_->setAtomStorageLayout(atomStorageLayout);
      info_->setRigidBodyStorageLayout(rigidBodyStorageLayout);
      info_->setCutoffGroupStorageLayout(cutoffGroupStorageLayout);

      info_->setSnapshotManager(
				new SimSnapshotManager(info_, atomStorageLayout,
						       rigidBodyStorageLayout,
						       cutoffGroupStorageLayout));

      forceMan_ = new ForceManager(info_);
      forceMan_->setDoElectricField(true);
      forceMan_->initialize();
    }

    info_->getSimParams()->setOutputElectricField(true);

    // -----------------------------------------------------------------------
    // Pre-allocate frameData_ to nFrames_ so that computeProperty1 can
    // index directly.  nFrames_ is set by TimeCorrFunc before preCorrelate()
    // calls computeProperty1, but we resize lazily inside computeProperty1
    // to be safe.
    // -----------------------------------------------------------------------
  }

  SFG::~SFG() {
    delete forceMan_;
  }

  // =========================================================================
  // computeProperty1(int frame)
  //
  // Called by the TimeCorrFunc framework for every stored frame, with
  // currentSnapshot_ already pointing at the correct snapshot.
  //
  // Workflow:
  //   1. Optionally recompute electric fields via ForceManager.
  //   2. Iterate over selected H atoms (sele1 = H-atom selection).
  //   3. For each H: project e-field onto OH axis → local frequency,
  //      transition dipole, and intramolecular coupling via spectroscopic maps.
  //   4. Collect per-OH vectors for Hamiltonian construction.
  //   5. Build and diagonalize the exciton Hamiltonian.
  //   6. Project local μ and α into exciton basis.
  //   7. Store FrameData for this frame.
  // =========================================================================
  void SFG::computeProperty1(int frame) {

    // Grow frameData_ on demand (safe if called out-of-order or if
    // nFrames_ was not yet available at construction time).
    if (frame >= static_cast<int>(frameData_.size()))
      frameData_.resize(frame + 1);

    // -----------------------------------------------------------------------
    // 1. Electric fields: compute via ForceManager if not in dump.
    // -----------------------------------------------------------------------
    if (!dumpHasElectricFields_) {
      forceMan_->setDoElectricField(true);
      forceMan_->calcForces();
    }

    // -----------------------------------------------------------------------
    // 2. Update dynamic selection if needed.
    // -----------------------------------------------------------------------
    if (evaluator1_.isDynamic())
      seleMan1_.setSelectionSet(evaluator1_.evaluate());

    // -----------------------------------------------------------------------
    // 3. Iterate over selected H atoms and collect per-OH data.
    //
    // Convention (matching OHFrequencyMap):
    //   - sele1 selects H atoms of water molecules.
    //   - Oxygen is mol->getAtomAt(0); hydrogens are getAtomAt(1) and (2).
    //   - OH unit vector: rOH = (Hpos − Opos).normalize()
    //   - Projected e-field: E = dot(rOH, sE)  with sE in atomic units.
    // -----------------------------------------------------------------------
    const RealType kcalToAU = 0.0008432975573; // kcal mol⁻¹ Å⁻¹ e⁻¹ → au

    // Per-OH data collected before Hamiltonian build
    std::vector<Vector3d> ohVecs;       // unit OH vectors (lab frame)
    std::vector<Vector3d> ohPos;        // OH midpoint positions (Å)
    std::vector<RealType> localFreqs;   // ω₁₀ from spectroscopic map [cm⁻¹]
    std::vector<Vector3d> localMu;      // transition dipole vectors
    std::vector<int>      molIndex;     // which molecule owns this OH
    std::vector<RealType> intraJ;       // field-dependent intramolecular coupling
    std::vector<std::pair<RealType,RealType>> alphaParams; // (alpha_par, alpha_perp) per OH
    
    int ii;
    StuntDouble* sd1;

    for (sd1 = seleMan1_.beginSelected(ii); sd1 != nullptr;
         sd1 = seleMan1_.nextSelected(ii)) {

      if (!sd1->isAtom()) continue;

      Atom* atom   = dynamic_cast<Atom*>(sd1);
      AtomType* at = atom->getAtomType();
      std::string name = at->getName();

      // Look up map entry for this H atom type
      auto wIt = w10_.find(name);
      if (wIt == w10_.end()) continue;   // not a mapped H type; skip

      // Retrieve parent molecule and oxygen position
      int sdID  = sd1->getGlobalIndex();
      int molID = info_->getGlobalMolMembership(sdID);
      Molecule* mol = info_->getMoleculeByGlobalIndex(molID);

      Vector3d Opos = mol->getAtomAt(0)->getPos();
      Vector3d Hpos = sd1->getPos();

      // Minimum-image OH vector and unit vector
      Vector3d rOH_raw = Hpos - Opos;
      Vector3d rOH = rOH_raw;
      currentSnapshot_->wrapVector(rOH);
      rOH.normalize();

      // Electric field on H atom (kcal mol⁻¹ Å⁻¹ e⁻¹) + applied field
      Vector3d sE = sd1->getElectricField() + EF_;
      sE *= kcalToAU;   // convert to atomic units

      // Projected e-field component along OH bond
      RealType E = dot(rOH, sE);

      // ω₁₀ frequency from quadratic map
      auto [a0, a1, a2] = wIt->second;
      RealType freq = a0 + a1 * E + a2 * E * E;

      // Transition dipole magnitude from muprime_ map
      auto [c0, c1, c2] = muprime_.at(name);
      RealType muMag = c0 + c1 * E + c2 * E * E;

      // Transition dipole vector: magnitude × unit OH direction
      Vector3d muVec = muMag * rOH;

      // Field-dependent intramolecular coupling
      auto [d0, d1, d2] = wintra_.at(name);
      RealType Jintra = d0 + d1 * E + d2 * E * E;

      // OH bond midpoint (used for TDC distance calculation)
      Vector3d midpoint = Opos + 0.5 * rOH_raw;

      // Store
      ohVecs.push_back(rOH);
      ohPos.push_back(midpoint);
      localFreqs.push_back(freq);
      localMu.push_back(muVec);
      molIndex.push_back(molID);
      intraJ.push_back(Jintra);
      alphaParams.push_back(alphaMap_.at(name));
    }

    int N = static_cast<int>(ohVecs.size());
    if (N == 0) {
      // No selected H atoms in this frame — store empty FrameData
      frameData_[frame] = FrameData{};
      return;
    }

    // -----------------------------------------------------------------------
    // 4. Build and diagonalize exciton Hamiltonian → eigvals, eigvecs
    // -----------------------------------------------------------------------
    std::vector<RealType>              eigvals;
    std::vector<std::vector<RealType>> eigvecs;

    buildExcitonHamiltonian(ohVecs, ohPos, localFreqs, localMu,
                            molIndex, intraJ, eigvals, eigvecs);

    // -----------------------------------------------------------------------
    // 5. Project local μ and α into exciton basis.
    //
    // For exciton mode k:
    //   μ_k  = Σᵢ eigvecs[k][i] · localMu[i]
    //   α_k  = Σᵢ eigvecs[k][i] · bondPolarizability(ohVecs[i])
    //
    // These are the collective (delocalized) quantities entering the TCF.
    // -----------------------------------------------------------------------
    FrameData fd;
    fd.eigvals.resize(N);
    fd.mu_ex.assign(N, V3Zero);
    fd.alpha_ex.assign(N, M3Zero);  // M3Zero is the predefined OpenMD
				    // zero matrix
    for (int k = 0; k < N; ++k) {
      fd.eigvals[k] = eigvals[k];

      for (int i = 0; i < N; ++i) {
        RealType c = eigvecs[k][i];

        // Collective transition dipole
        fd.mu_ex[k] += c * localMu[i];
	fd.alpha_ex[k] += c * bondPolarizability(ohVecs[i],
						 alphaParams[i].first,
						 alphaParams[i].second);
      }
    }
    frameData_[frame] = std::move(fd);
  }

  // =========================================================================
  // buildExcitonHamiltonian
  // =========================================================================
  void SFG::buildExcitonHamiltonian(
				    const std::vector<Vector3d>& ohVecs,
				    const std::vector<Vector3d>& ohPos,
				    const std::vector<RealType>& localFreqs,
				    const std::vector<Vector3d>& localMu,
				    const std::vector<int>&      molIndex,
				    const std::vector<RealType>& intraJ,
				    std::vector<RealType>&       eigvals,
				    std::vector<std::vector<RealType>>& eigvecs) {

    int N = static_cast<int>(localFreqs.size());
    DynamicRectMatrix<RealType> H(N, N, 0.0);

    // Diagonal: local ω₁₀ frequencies
    for (int i = 0; i < N; ++i)
      H(i, i) = localFreqs[i];

    // Off-diagonal couplings
    for (int i = 0; i < N; ++i) {
      for (int j = i + 1; j < N; ++j) {
        RealType J = 0.0;

        if (molIndex[i] == molIndex[j]) {
          // Intramolecular coupling — field-dependent average of the two OHs
          J = 0.5 * (intraJ[i] + intraJ[j]);
        } else {
          // Intermolecular TDC coupling
          Vector3d r_ij = ohPos[j] - ohPos[i];
          // Apply minimum image convention
          currentSnapshot_->wrapVector(r_ij);
          J = tdcCoupling(r_ij, localMu[i], localMu[j]);
        }

        H(i, j) = J;
        H(j, i) = J;
      }
    }

    // Diagonalize with JAMA symmetric eigensolver
    JAMA::Eigenvalue<RealType> eigensystem(H);
    DynamicRectMatrix<RealType> evects(N, N);
    DynamicVector<RealType>     evals(N);

    eigensystem.getRealEigenvalues(evals);
    eigensystem.getV(evects);

    eigvals.resize(N);
    eigvecs.assign(N, std::vector<RealType>(N));

    for (int k = 0; k < N; ++k) {
      eigvals[k] = evals[k];
      for (int i = 0; i < N; ++i)
        eigvecs[k][i] = evects(i, k);   // evects column k = mode k
    }
  }

  // =========================================================================
  // tdcCoupling
  // J_ij = prefactor × (μᵢ·μⱼ − 3(μᵢ·r̂)(μⱼ·r̂)) / r³
  // prefactor = 5034.12 cm⁻¹ D⁻² Å³  (standard conversion)
  // =========================================================================
  RealType SFG::tdcCoupling(const Vector3d& r_ij, const Vector3d& mu_i,
			    const Vector3d& mu_j) {
    RealType r = r_ij.length();
    if (r < 1.0e-6) return 0.0;

    Vector3d rhat = r_ij / r;

    RealType dot_ij = dot(mu_i, mu_j);
    RealType dot_ir = dot(mu_i, rhat);
    RealType dot_jr = dot(mu_j, rhat);

    const RealType prefactor = 5034.12;
    return prefactor * (dot_ij - 3.0 * dot_ir * dot_jr) / (r * r * r);
  }

  // =========================================================================
  // bondPolarizability
  // α_tensor = α_perp · I + (α_par − α_perp) · ê⊗ê
  // =========================================================================
  Mat3x3d SFG::bondPolarizability(const Vector3d& ohUnit, RealType alpha_par,
				  RealType alpha_perp) {
    Mat3x3d a_tens = Mat3x3d::identity();
    a_tens *= alpha_perp;
    a_tens += (alpha_par - alpha_perp) * outProduct(ohUnit, ohUnit);
    return a_tens;
  }

  // =========================================================================
  // calcCorrVal
  // Polarization-selected cross-correlation at lag (frame2 − frame1).
  //
  // Non-adiabatic approximation (MultiSpec default):
  //   C(t) = [Σₖ αpq,k(t)] × [Σₗ μr,l(0)]
  // where k,l index exciton modes independently at their respective frames.
  // =========================================================================
  RealType SFG::calcCorrVal(int frame1, int frame2) {
    const FrameData& fd0 = frameData_[frame1];  // t = 0
    const FrameData& fd1 = frameData_[frame2];  // t

    if (fd0.mu_ex.empty() || fd1.alpha_ex.empty()) return 0.0;

    RealType val = 0.0;

    if (polarization_ == "ssp") {
      // ssp: α_yy(t) × μ_z(0)
      RealType alpha_yy = 0.0;
      for (const auto& a : fd1.alpha_ex) alpha_yy += a(1, 1);

      RealType mu_z = 0.0;
      for (const auto& m : fd0.mu_ex) mu_z += m[2];

      val = alpha_yy * mu_z;

    } else if (polarization_ == "sps") {
      // sps: α_yz(t) × μ_y(0)
      RealType alpha_yz = 0.0;
      for (const auto& a : fd1.alpha_ex) alpha_yz += a(1, 2);

      RealType mu_y = 0.0;
      for (const auto& m : fd0.mu_ex) mu_y += m[1];

      val = alpha_yz * mu_y;

    } else if (polarization_ == "ppp") {
      // ppp: combination — see Buch et al. 2007 eq. 10
      // χ²_ppp ∝ −cos β_SF sin β_IR sin β_vis · χ²_xxz
      //          −sin β_SF cos β_IR cos β_vis · χ²_zxx
      //          +sin β_SF sin β_IR sin β_vis · χ²_zzz
      // For the near-normal incidence limit, retain χ²_zzz and χ²_xxz:
      RealType alpha_zz = 0.0, alpha_xz = 0.0;
      for (const auto& a : fd1.alpha_ex) {
        alpha_zz += a(2, 2);
        alpha_xz += a(0, 2);
      }
      RealType mu_z = 0.0, mu_x = 0.0;
      for (const auto& m : fd0.mu_ex) {
        mu_z += m[2];
        mu_x += m[0];
      }
      val = alpha_zz * mu_z - alpha_xz * mu_x;
    }

    return val;
  }

  // =========================================================================
  // postCorrelate
  //
  // The base class divides histogram_[i] by count_[i], producing the
  // normalized time-domain TCF C(t).  Nothing more is needed here.
  // =========================================================================
  void SFG::postCorrelate() {
    SystemACF<RealType>::postCorrelate();
  }

  // =========================================================================
  // writeCorrelate
  //
  // Writes two output files:
  //
  //   .sfg      — the SFG spectrum Im χ²(ω), obtained by FFT of the TCF
  //               followed by a quantum correction
  //   .sfg.tcf  — the raw time-domain TCF C(t), matching the format of
  //               the base-class writeCorrelate for diagnostics
  //
  // FFT details:
  //   - histogram_[i] = C(i * dtMean_) is real-valued (RealType scalar).
  //   - We use FFTW's real-to-complex r2c plan: N real inputs →
  //     N/2+1 complex outputs covering 0 … Nyquist.
  //   - The transform is unnormalized; divide by N for correct scaling.
  //   - We use dtMean_ (the measured mean inter-frame interval, in fs)
  //     rather than deltaTime_ (the nominal sampleTime from simParams),
  //     because the base class uses dtMean_ for all time-bin arithmetic
  //     and it may differ slightly from deltaTime_ when allowTimeFuzz_
  //     is set.
  //
  // Frequency axis:
  //   dt_ps    = dtMean_ * 1e-3              [fs → ps]
  //   dOmega   = 33.3564 / (N * dt_ps)      [cm⁻¹ per bin]
  //   omega_k  = k * dOmega                 [cm⁻¹]
  //   33.3564 cm⁻¹ ps = 1e12 / (100 * c_SI)
  //
  // Quantum correction (full Bose-Einstein form):
  //   x      = omega_k / kT_invcm           [dimensionless]
  //   Q(ω)   = x / (1 − exp(−x))
  //   kT_invcm = 0.6950356 * T [K]          [cm⁻¹]
  //   T is read from SimParams (targetTemp); default 300 K.
  //   Guard: for x > 50 (deep quantum limit), exp(−x) ≈ 0 so Q ≈ x.
  //   At ω → 0, Q → 1 (L'Hôpital).
  //
  // SFG convention:
  //   Im χ²(ω) ∝ Im[ FT[C(t)] ] * Q(ω)
  //   The quantum correction is applied to Im only; Re is written
  //   uncorrected.  |χ²|² uses the corrected Im.
  // =========================================================================
  void SFG::writeCorrelate() {
    Revision r;

    // -----------------------------------------------------------------------
    // Write the time-domain TCF for diagnostics, matching base-class format.
    // -----------------------------------------------------------------------
    {
      std::string tcfName = outputFilename_ + ".tcf";
      std::ofstream tcf(tcfName.c_str());
      if (!tcf.is_open()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "SFG::writeCorrelate: cannot open %s\n", tcfName.c_str());
        painCave.isFatal = 1;
        simError();
      }
      tcf << "# " << getCorrFuncType() << " time-domain TCF\n";
      tcf << "# OpenMD " << r.getFullRevision() << "\n";
      tcf << "# " << r.getBuildDate() << "\n";
      tcf << "# selection script1: \"" << selectionScript1_
          << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      tcf << "# polarization: " << polarization_ << "\n";
      tcf << "#time(fs)\tC(t)\n";
      for (unsigned int i = 0; i < nTimeBins_; ++i) {
        if (count_[i] > 0)
          tcf << times_[i] - times_[0] << "\t" << histogram_[i] << "\n";
      }
      tcf.close();
    }

#if defined(HAVE_FFTW3_H)
    
    // -----------------------------------------------------------------------
    // Determine the number of valid (non-zero count) TCF points.
    // The base class fills histogram_ up to nTimeBins_ = nFrames_, but
    // the longer lag times have fewer samples and count_[i] may be zero
    // for i near nTimeBins_-1.  Use all bins with count > 0.
    // -----------------------------------------------------------------------
    int N = 0;
    for (unsigned int i = 0; i < nTimeBins_; ++i)
      if (count_[i] > 0) N++;

    if (N < 2) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SFG::writeCorrelate: too few TCF points (%d) for FFT.\n", N);
      painCave.isFatal = 1;
      simError();
    }

    // -----------------------------------------------------------------------
    // Allocate FFTW buffers and copy TCF.
    // -----------------------------------------------------------------------
    double*       in  = (double*)      fftw_malloc(sizeof(double)       * N);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

    if (!in || !out) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SFG::writeCorrelate: FFTW memory allocation failed.\n");
      painCave.isFatal = 1;
      simError();
    }

    int idx = 0;
    for (unsigned int i = 0; i < nTimeBins_; ++i)
      if (count_[i] > 0) in[idx++] = static_cast<double>(histogram_[i]);

    // -----------------------------------------------------------------------
    // Execute r2c FFT.
    // FFTW_ESTIMATE avoids expensive measurement; use FFTW_MEASURE if
    // this analysis is run repeatedly on the same N.
    // -----------------------------------------------------------------------
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // -----------------------------------------------------------------------
    // Frequency axis.
    // Use dtMean_ (fs) — the measured mean inter-frame interval — to be
    // consistent with how the base class assigns time bins.
    // -----------------------------------------------------------------------
    const double ps_to_invcm = 33.3564;              // cm⁻¹ per ps⁻¹
    const double dt_ps       = dtMean_ * 1.0e-3;     // fs → ps
    const double dOmega      = ps_to_invcm / (static_cast<double>(N) * dt_ps);

    // -----------------------------------------------------------------------
    // Temperature and quantum correction prefactor.
    // -----------------------------------------------------------------------
    RealType T = 300.0;
    if (info_->getSimParams()->haveTargetTemp())
      T = info_->getSimParams()->getTargetTemp();
    const double kT_invcm = 0.6950356 * static_cast<double>(T);

    // -----------------------------------------------------------------------
    // Open spectrum output file and write header.
    // -----------------------------------------------------------------------
    std::ofstream ofs(outputFilename_.c_str());
    if (!ofs.is_open()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SFG::writeCorrelate: cannot open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs << "# " << getCorrFuncType() << " spectrum\n";
    ofs << "# OpenMD " << r.getFullRevision() << "\n";
    ofs << "# " << r.getBuildDate() << "\n";
    ofs << "# selection script1: \"" << selectionScript1_
        << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
    ofs << "# polarization: " << polarization_ << "\n";
    ofs << "# T = " << T << " K,  kT = " << kT_invcm << " cm-1\n";
    ofs << "# N (TCF points) = " << N << "\n";
    ofs << "# dtMean = " << dtMean_ << " fs,  dOmega = " << dOmega
        << " cm-1\n";
    if (!paramString_.empty())
      ofs << "# parameters: " << paramString_ << "\n";
    ofs << "#omega(cm-1)\tRe(chi2)\tIm(chi2)\t|chi2|^2\n";

    // -----------------------------------------------------------------------
    // Loop over positive-frequency bins, apply quantum correction, write.
    // -----------------------------------------------------------------------
    int nOut = N / 2 + 1;
    for (int k = 0; k < nOut; ++k) {
      double omega = k * dOmega;

      // Normalize the FFT output
      double re = out[k][0] / static_cast<double>(N);
      double im = out[k][1] / static_cast<double>(N);

      // Quantum correction Q(ω) = x / (1 − exp(−x)),  x = ω / kT
      double qCorr = 1.0;
      if (omega > 1.0e-6) {
        double x = omega / kT_invcm;
        qCorr    = (x < 50.0) ? x / (1.0 - std::exp(-x)) : x;
      }

      double im_corr = im * qCorr;
      double mod2    = re * re + im_corr * im_corr;

      ofs << omega    << "\t"
          << re       << "\t"
          << im_corr  << "\t"
          << mod2     << "\n";
    }

    ofs.close();
    fftw_free(in);
    fftw_free(out);
#else
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "SFG::writeCorrelate: FFTW3 is required for SFG spectra.\n");
    painCave.isFatal = 1;
    simError();
#endif
    
  }
} // namespace OpenMD
