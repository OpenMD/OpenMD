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

#include "applications/staticProps/SFGTimeAvg.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

#include "brains/ForceManager.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "io/DumpReader.hpp"
#include "math/Eigenvalue.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

#if defined(HAVE_LAPACK)
extern "C" {
  void dsyevd_(const char* jobz, const char* uplo, const int* n,
               double* a, const int* lda, double* w,
               double* work, const int* lwork,
               int* iwork, const int* liwork, int* info);
}
#endif

namespace OpenMD {

  SFGTimeAvg::SFGTimeAvg(SimInfo* info, const std::string& filename,
                         const std::string& sele1, int nbins,
                         const std::string& polarization, int privilegedAxis,
                         RealType gamma, RealType fc) :
    StaticAnalyser(info, filename, nbins),
    selectionScript1_(sele1), seleMan1_(info_), evaluator1_(info_),
    polarization_(polarization),
    privilegedAxis_(privilegedAxis),
    s1_((privilegedAxis + 1) % 3),
    s2_((privilegedAxis + 2) % 3),
    gamma_(gamma), fc_(fc) {

    setOutputName(getPrefix(filename) + ".sfgta");

    dumpHasElectricFields_ = info_->getSimParams()->getOutputElectricField();

    if (!dumpHasElectricFields_) {
      int atomStorageLayout        = info_->getAtomStorageLayout();
      int rigidBodyStorageLayout   = info->getRigidBodyStorageLayout();
      int cutoffGroupStorageLayout = info->getCutoffGroupStorageLayout();

      atomStorageLayout |= DataStorage::dslElectricField;
      rigidBodyStorageLayout |= DataStorage::dslElectricField;
      info_->setAtomStorageLayout(atomStorageLayout);
      info_->setRigidBodyStorageLayout(rigidBodyStorageLayout);
      info_->setCutoffGroupStorageLayout(cutoffGroupStorageLayout);

      info_->setSnapshotManager(new SimSnapshotManager(info_, atomStorageLayout,
                                                       rigidBodyStorageLayout,
                                                       cutoffGroupStorageLayout));
    }
    info_->getSimParams()->setOutputElectricField(true);

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    chi2_.resize(nBins_);

    // ---------------------------------------------------------------------
    // Spectroscopic maps (identical to the dynamical SFG module).
    // ---------------------------------------------------------------------
    // SPC/E (Auer & Skinner 2008)
    w10_["H_SPCE"]      = std::make_tuple(3761.6, -5060.4, -86225.0);
    muPrime_["H_SPCE"]  = std::make_tuple(0.1333, 14.17, 0.0);
    x10_["H_SPCE"]      = std::make_pair(0.1934,  -1.75e-5);
    p10_["H_SPCE"]      = std::make_pair(1.611,    5.893e-4);
    wintra_["H_SPCE"]   = std::make_tuple(-1789.0, 23852.0, -1.966);
    alphaMap_["H_SPCE"] = std::make_pair(0.185, 0.033);
    tdcLoc_["H_SPCE"]   = 0.58;

    // TIP4P (Gruenbaum et al. 2013)
    w10_["H_TIP4P"]      = std::make_tuple(3760.2,  -3541.7, -152677.0);
    muPrime_["H_TIP4P"]  = std::make_tuple(0.1646,  11.39,    63.41);
    x10_["H_TIP4P"]      = std::make_pair(0.19285, -1.7261e-5);
    p10_["H_TIP4P"]      = std::make_pair(1.6466,   5.7692e-4);
    wintra_["H_TIP4P"]   = std::make_tuple(-1361.0,  27165.0,  -1.887);
    alphaMap_["H_TIP4P"] = std::make_pair(0.185, 0.033);
    tdcLoc_["H_TIP4P"]   = 0.67;

    // TIP4P-Ice (transferability from TIP4P; Takayama 2023)
    w10_["H_TIP4P-Ice"]      = w10_["H_TIP4P"];
    muPrime_["H_TIP4P-Ice"]  = muPrime_["H_TIP4P"];
    x10_["H_TIP4P-Ice"]      = x10_["H_TIP4P"];
    p10_["H_TIP4P-Ice"]      = p10_["H_TIP4P"];
    wintra_["H_TIP4P-Ice"]   = wintra_["H_TIP4P"];
    alphaMap_["H_TIP4P-Ice"] = alphaMap_["H_TIP4P"];
    tdcLoc_["H_TIP4P-Ice"]   = tdcLoc_["H_TIP4P"];

    // HOH bend overtone maps (Ni & Skinner 2015)
    wb01_["TIP4P"]     = std::make_pair(1581.46, 2938.51);
    wb12_["TIP4P"]     = std::make_pair(1551.32, 3147.80);
    wb01_["TIP4P-Ice"] = wb01_["TIP4P"];
    wb12_["TIP4P-Ice"] = wb12_["TIP4P"];
    wb01_["SPCE"]      = wb01_["TIP4P"];
    wb12_["SPCE"]      = wb12_["TIP4P"];

    // Applied external field
    EF_ = V3Zero;
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
                 "SFGTimeAvg: uniformField needs 3 parameters, %zu given.\n",
                 ef.size());
        painCave.isFatal = 1; simError();
      }
      EF_.x() = ef[0]; EF_.y() = ef[1]; EF_.z() = ef[2];
    }
  }

  // =========================================================================
  // process: single pass over the dump file
  // =========================================================================
  void SFGTimeAvg::process() {
    ForceManager* forceMan = nullptr;

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = 0;

    if (!dumpHasElectricFields_) {
      forceMan = new ForceManager(info_);
      forceMan->setDoElectricField(true);
      forceMan->initialize();
    }

    std::fill(chi2_.begin(), chi2_.end(),
              std::complex<double>(0.0, 0.0));

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (!dumpHasElectricFields_) {
        forceMan->setDoElectricField(true);
        forceMan->calcForces();
      }

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }

      std::vector<Vector3d> ohPos;
      std::vector<int>      molIndex;
      std::vector<RealType> intraJ;
      FrameData fd = extractFrame(ohPos, molIndex, intraJ);
      if (fd.N == 0) continue;

      buildHamiltonian(fd, ohPos, molIndex, intraJ);
      accumulateFrame(fd);
      nProcessed_++;
    }

    if (forceMan) delete forceMan;

    // Orientation diagnostic summary
    double avgFree = (diagNFree_ > 0)
      ? diagMuzFree_ / static_cast<double>(diagNFree_) : 0.0;
    double avgHB = (diagNHB_ > 0)
      ? diagMuzHB_ / static_cast<double>(diagNHB_) : 0.0;
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "SFGTimeAvg orientation diagnostic:\n"
             "\t<mu_z> dangling (freq>3650): %.5f  [n=%ld]\n"
             "\t<mu_z> H-bonded (freq<3450): %.5f  [n=%ld]\n"
             "\t(these should have OPPOSITE signs for a real interface)\n",
             avgFree, diagNFree_, avgHB, diagNHB_);
    painCave.isFatal = 0; painCave.severity = OPENMD_INFO; simError();

    writeSpectrum();
  }

  // =========================================================================
  // extractFrame: build chromophore list, transition dipoles/polarizabilities,
  // and Hamiltonian diagonal for the current snapshot.  (Duplicated from the
  // dynamical SFG module.)
  // =========================================================================
  SFGTimeAvg::FrameData
  SFGTimeAvg::extractFrame(std::vector<Vector3d>& ohPos,
                           std::vector<int>&      molIndex,
                           std::vector<RealType>& intraJ) {
    const RealType kcalToAU = 0.0008432975573;
    const RealType chrgToKcal = 23.060548;

    FrameData fd;
    ohPos.clear(); molIndex.clear(); intraJ.clear();

    std::vector<RealType> freqs;
    std::vector<RealType> x10vals, p10vals, mxvals, intraJ_k0, intraJ_kp;

    struct MolData {
      Vector3d Opos;
      Vector3d ohUnit1, ohUnit2;
      Vector3d Efield1, Efield2;
      RealType rOH1 = 0.0, rOH2 = 0.0;
      int ohCount = 0;
      std::string modelName;
    };
    std::map<int, MolData> molMap;

    int ii;
    Molecule* mol;
    for (mol = seleMan1_.beginSelectedMolecule(ii); mol != nullptr;
         mol = seleMan1_.nextSelectedMolecule(ii)) {

      int molID = mol->getGlobalIndex();

      // Locate the oxygen and the (up to two) hydrogens of this water.
      Vector3d Opos;
      bool foundO = false;
      std::vector<Atom*> hAtoms;
      Molecule::AtomIterator ai;
      for (Atom* atom = mol->beginAtom(ai); atom != nullptr;
           atom = mol->nextAtom(ai)) {
        const std::string& aname = atom->getAtomType()->getName();
        if (aname[0] == 'O' && !foundO) {
          Opos = atom->getPos();
          foundO = true;
        } else if (w10_.find(aname) != w10_.end()) {
          // An H atom type for which we have a stretch map
          hAtoms.push_back(atom);
        }
      }
      if (!foundO) continue;
      if (hAtoms.empty()) continue;

      // Build a stretch chromophore for EACH OH bond of this molecule.
      // Iterating over molecules (selected by COM z) and including both
      // OH bonds removes the orientational bias that arises when selecting
      // individual hydrogens within a thin z-slab (which preferentially
      // admits up-pointing OH bonds).
      for (Atom* hAtom : hAtoms) {
        std::string name = hAtom->getAtomType()->getName();

        auto wIt  = w10_.find(name);
        auto mpIt = muPrime_.find(name);
        if (wIt == w10_.end() || mpIt == muPrime_.end()) continue;

        Vector3d Hpos    = hAtom->getPos();
        Vector3d rOH_mic = Hpos - Opos;
        currentSnapshot_->wrapVector(rOH_mic);
        RealType rOH_len = rOH_mic.length();
        Vector3d ohUnit  = rOH_mic; ohUnit.normalize();

        Vector3d sE = hAtom->getElectricField();
        sE += EF_ * chrgToKcal;
        sE *= kcalToAU;
        RealType E = dot(ohUnit, sE);

        auto [a0, a1, a2] = wIt->second;
        RealType freq = a0 + a1*E + a2*E*E;

        auto [m0, m1, m2] = mpIt->second;
        RealType muPrime = m0 + m1*E + m2*E*E;

        auto [x0, x1] = x10_.at(name);
        RealType x10  = x0 + x1*freq;

        auto [pv0, pv1] = p10_.at(name);
        RealType p10  = pv0 + pv1*freq;

        const RealType eaBohr_to_Debye = 2.5417464;
        RealType muMag = muPrime * x10 * eaBohr_to_Debye;

        RealType tdcDist = 0.5 * rOH_len;
        auto tdcIt = tdcLoc_.find(name);
        if (tdcIt != tdcLoc_.end()) tdcDist = tdcIt->second;
        Vector3d dipolePos = Opos + tdcDist * ohUnit;

        auto [apar, aperp] = alphaMap_.at(name);
        RealType x10_gas = x10_.at(name).first;

        fd.N++;
        fd.mu.push_back(muMag * ohUnit);
        fd.alpha.push_back(bondPolarizability(ohUnit, apar, aperp,
                                              x10, x10_gas));

        // Orientation diagnostic: record mu_z by frequency class
        {
          double muz = (muMag * ohUnit)[privilegedAxis_];
          if (freq > 3650.0)      { diagMuzFree_ += muz; diagNFree_++; }
          else if (freq < 3450.0) { diagMuzHB_   += muz; diagNHB_++;   }
        }

        ohPos.push_back(dipolePos);
        molIndex.push_back(molID);
        freqs.push_back(freq);
        x10vals.push_back(x10);
        p10vals.push_back(p10);
        mxvals.push_back(muPrime * x10);

        auto [k0, k1, kp] = wintra_.at(name);
        intraJ.push_back(k0 + k1*E);
        intraJ_k0.push_back(k0);
        intraJ_kp.push_back(kp);

        MolData& md = molMap[molID];
        md.Opos = Opos;
        md.modelName = mol->getType();
        if (md.ohCount == 0) {
          md.ohUnit1 = ohUnit; md.Efield1 = sE; md.rOH1 = rOH_len;
        } else if (md.ohCount == 1) {
          md.ohUnit2 = ohUnit; md.Efield2 = sE; md.rOH2 = rOH_len;
        }
        md.ohCount++;
      }
    }

    int nStretch = fd.N;
    fd.nStretch = nStretch;

    // Bend overtone chromophores
    int nBend = 0;
    if (fc_ != 0.0) {
      for (auto& [molID, md] : molMap) {
        if (md.ohCount != 2) continue;
        auto bIt = wb01_.find(md.modelName);
        if (bIt == wb01_.end()) continue;
        auto cIt = wb12_.find(md.modelName);
        if (cIt == wb12_.end()) continue;

        Vector3d vp = cross(md.ohUnit1, md.ohUnit2);
        RealType vpLen = vp.length();
        if (vpLen < 1.0e-9) continue;
        vp /= vpLen;

        Vector3d ePerp1 = cross(vp, md.ohUnit1);
        Vector3d ePerp2 = cross(md.ohUnit2, vp);

        RealType Eb = dot(ePerp1, md.Efield1) / md.rOH1
	  + dot(ePerp2, md.Efield2) / md.rOH2;

        auto [b0_10, b1_10] = bIt->second;
        auto [b0_21, b1_21] = cIt->second;
        RealType w_2d = (b0_10 + b1_10*Eb) + (b0_21 + b1_21*Eb);

        fd.N++;
        fd.mu.push_back(V3Zero);
        fd.alpha.push_back(Mat3x3d(0.0));

        ohPos.push_back(md.Opos);
        molIndex.push_back(molID);
        freqs.push_back(w_2d);
        x10vals.push_back(0.0);
        p10vals.push_back(0.0);
        mxvals.push_back(0.0);
        intraJ.push_back(0.0);
        intraJ_k0.push_back(0.0);
        intraJ_kp.push_back(0.0);

        nBend++;
      }
    }
    fd.nBend = nBend;
    useFermi_ = (nBend > 0);

    int N = fd.N;
    fd.H = DynamicRectMatrix<RealType>(N, N, 0.0);
    for (int i = 0; i < N; ++i)
      fd.H(i, i) = freqs[i];

    // Pack 6 per-site entries for buildHamiltonian
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
  // buildHamiltonian: fill off-diagonals.  (Duplicated from dynamical SFG.)
  // =========================================================================
  void SFGTimeAvg::buildHamiltonian(FrameData&                   fd,
                                    const std::vector<Vector3d>& ohPos,
                                    const std::vector<int>&      molIndex,
                                    const std::vector<RealType>& intraJ) {
    int N = fd.N;
    if (N == 0) return;
    const int nStr = fd.nStretch;

    for (int i = 0; i < N; ++i) {
      bool i_is_bend = (i >= nStr);
      for (int j = i+1; j < N; ++j) {
        bool j_is_bend = (j >= nStr);
        RealType J = 0.0;
        if (molIndex[i] == molIndex[j]) {
          if (!i_is_bend && !j_is_bend) {
            RealType Ki     = intraJ[6*i + 0];
            RealType xi     = intraJ[6*i + 1];
            RealType pi     = intraJ[6*i + 2];
            RealType k0_val = intraJ[6*i + 3];
            RealType kp_val = intraJ[6*i + 4];
            RealType Kj     = intraJ[6*j + 0];
            RealType xj     = intraJ[6*j + 1];
            RealType pj     = intraJ[6*j + 2];
            J = (Ki + Kj - k0_val) * xi * xj + kp_val * pi * pj;
          } else if (i_is_bend != j_is_bend) {
            J = fc_;     // Fermi coupling stretch <-> bend overtone
          }
        } else {
          Vector3d r_ij = ohPos[j] - ohPos[i];
          currentSnapshot_->wrapVector(r_ij);
          Vector3d e_i = fd.mu[i]; RealType m_i = e_i.length();
          if (m_i > 1.0e-12) e_i /= m_i;
          Vector3d e_j = fd.mu[j]; RealType m_j = e_j.length();
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
  // diagonalize: real symmetric eigensolve.  evecs row-major: evecs[i*N+a].
  // =========================================================================
  void SFGTimeAvg::diagonalize(const DynamicRectMatrix<RealType>& H, int N,
                               std::vector<double>& evals,
                               std::vector<double>& evecs) {
    evals.assign(N, 0.0);
    evecs.assign(N*N, 0.0);

#if defined(HAVE_LAPACK)
    std::vector<double> A(N*N);
    {
      const RealType* src = H.getArrayPointer();
      for (int i = 0; i < N*N; ++i) A[i] = static_cast<double>(src[i]);
    }
    // Workspace query
    int lwork = -1, liwork = -1, info = 0;
    double wq = 0.0; int wiq = 0;
    dsyevd_("V", "U", &N, A.data(), &N, evals.data(),
            &wq, &lwork, &wiq, &liwork, &info);
    lwork  = static_cast<int>(wq);
    liwork = wiq;
    std::vector<double> work(std::max(1, lwork));
    std::vector<int>    iwork(std::max(1, liwork));
    info = 0;
    dsyevd_("V", "U", &N, A.data(), &N, evals.data(),
            work.data(), &lwork, iwork.data(), &liwork, &info);
    if (info == 0) {
      // A holds eigenvectors column-major: A[i + a*N] = component i of state a.
      // Store row-major: evecs[i*N + a].
      for (int i = 0; i < N; ++i)
        for (int a = 0; a < N; ++a)
          evecs[i*N + a] = A[i + a*N];
      return;
    }
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "SFGTimeAvg::diagonalize: dsyevd info=%d; using JAMA.\n", info);
    painCave.isFatal = 0; painCave.severity = OPENMD_WARNING; simError();
#endif

    // JAMA fallback
    JAMA::Eigenvalue<RealType> eig(H);
    DynamicVector<RealType> ev(N);
    DynamicRectMatrix<RealType> V(N, N);
    eig.getRealEigenvalues(ev);
    eig.getV(V);
    for (int a = 0; a < N; ++a) {
      evals[a] = static_cast<double>(ev[a]);
      for (int i = 0; i < N; ++i)
        evecs[i*N + a] = static_cast<double>(V(i, a));
    }
  }

  // =========================================================================
  // accumulateFrame: diagonalize H, project transition moments onto
  // eigenstates, and add each eigenstate's complex Lorentzian to the spectrum.
  // =========================================================================
  void SFGTimeAvg::accumulateFrame(const FrameData& fd) {
    int N = fd.N;
    std::vector<double> evals, evecs;
    diagonalize(fd.H, N, evals, evecs);

    const double dOmega = (maxFreq_ - minFreq_) / static_cast<double>(nBins_);
    // Only evaluate Lorentzian within +/- window bins of each eigenfrequency
    const double cutoff = 20.0 * gamma_;
    const int    halfWin = static_cast<int>(std::ceil(cutoff / dOmega));

    for (int a = 0; a < N; ++a) {
      double wa = evals[a];
      if (wa < minFreq_ - cutoff || wa > maxFreq_ + cutoff) continue;

      // Eigenstate transition dipole and polarizability components.
      // mu_r^a = sum_i V(i,a) mu_{r,i}; alpha_pq^a = sum_i V(i,a) alpha_{pq,i}
      double mu_n = 0.0;          // dipole along privileged (normal) axis
      double a_s1s1 = 0.0, a_s2s2 = 0.0, a_nn = 0.0;
      double a_s1n = 0.0, a_s2n = 0.0;
      for (int i = 0; i < N; ++i) {
        double v = evecs[i*N + a];
        if (v == 0.0) continue;
        const Vector3d& mu = fd.mu[i];
        const Mat3x3d&  al = fd.alpha[i];
        mu_n   += v * mu[privilegedAxis_];
        a_s1s1 += v * al(s1_, s1_);
        a_s2s2 += v * al(s2_, s2_);
        a_nn   += v * al(privilegedAxis_, privilegedAxis_);
        a_s1n  += v * al(s1_, privilegedAxis_);
        a_s2n  += v * al(s2_, privilegedAxis_);
      }

      // Polarization combination amplitude A_a = alpha_pq^a * mu_r^a
      double amp = 0.0;
      if (polarization_ == "ssp") {
        // SSP: 0.5 (alpha_s1s1 + alpha_s2s2) mu_n
        amp = 0.5 * (a_s1s1 + a_s2s2) * mu_n;
      } else if (polarization_ == "ppp") {
        amp = (-0.5*a_s1s1 - 0.5*a_s2s2 + a_nn) * mu_n;
      } else if (polarization_ == "sps") {
        // SPS: alpha_s1n mu_s1 component; use s1-normal cross term
        amp = 0.5 * (a_s1n + a_s2n) * mu_n;
      } else {
        amp = 0.5 * (a_s1s1 + a_s2s2) * mu_n;  // default ssp
      }

      // Add complex Lorentzian to nearby bins.
      // The SFG response carries a leading factor of i:
      //   chi2(w) ∝ i ∫ dt e^{-iwt} <...>,
      // so the absorptive (imaginary) part tracks +A_a.  A bare
      // Lorentzian A/(w - w_a + i*gamma) has Im = -A*gamma/(...),
      // i.e. the wrong sign; we therefore use the conjugate denominator
      // A/(w - w_a - i*gamma), whose imaginary part is +A*gamma/(...).
      int centerBin = static_cast<int>((wa - minFreq_) / dOmega);
      int lo = std::max(0, centerBin - halfWin);
      int hi = std::min(static_cast<int>(nBins_) - 1, centerBin + halfWin);
      for (int b = lo; b <= hi; ++b) {
        double w = minFreq_ + (b + 0.5) * dOmega;
        std::complex<double> denom(w - wa, -gamma_);
        chi2_[b] += amp / denom;
      }
    }
  }

  // =========================================================================
  // writeSpectrum: normalize, apply Bose-Einstein quantum correction, write.
  // =========================================================================
  void SFGTimeAvg::writeSpectrum() {
    const double dOmega = (maxFreq_ - minFreq_) / static_cast<double>(nBins_);

    // Temperature for quantum correction
    RealType T = 300.0;
    if (info_->getSimParams()->haveTargetTemp())
      T = info_->getSimParams()->getTargetTemp();
    const double kT_invcm = 0.6950356 * static_cast<double>(T);

    double invN = (nProcessed_ > 0)
      ? 1.0 / static_cast<double>(nProcessed_) : 1.0;

    std::ofstream ofs(outputFilename_.c_str());
    if (!ofs.is_open()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SFGTimeAvg: cannot open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1; simError();
    }

    ofs << "# SFG spectrum (time-averaging approximation, "
        << polarization_ << ")\n";
    ofs << "# selection: (" << selectionScript1_ << ")\n";
    ofs << "# nFrames = " << nProcessed_ << "\n";
    ofs << "# T = " << T << " K, kT = " << kT_invcm << " cm-1\n";
    ofs << "# Lorentzian half-width gamma = " << gamma_ << " cm-1\n";
    ofs << "# Fermi coupling fc = " << fc_ << " cm-1"
        << (useFermi_ ? "  (bend overtone included)" : "") << "\n";
    ofs << "#omega(cm-1)\tRe(chi2)\tIm(chi2)\t|chi2|^2\n";

    for (int b = 0; b < static_cast<int>(nBins_); ++b) {
      double w = minFreq_ + (b + 0.5) * dOmega;
      std::complex<double> c = chi2_[b] * invN;

      double x = w / kT_invcm;
      double qCorr = (x > 1.0e-6)
        ? ((x < 50.0) ? x / (1.0 - std::exp(-x)) : x) : 1.0;
      double re = c.real() * qCorr;
      double im = c.imag() * qCorr;

      ofs << w << "\t" << re << "\t" << im << "\t"
          << re*re + im*im << "\n";
    }
    ofs.close();
  }

  // =========================================================================
  // bondPolarizability: transition polarizability tensor for one OH bond,
  // scaled by x10/x10_gas (MultiSpec trPol convention).
  // =========================================================================
  Mat3x3d SFGTimeAvg::bondPolarizability(const Vector3d& ohUnit,
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
  // tdcCoupling: transition-dipole coupling in MultiSpec convention.
  //   J [cm-1] = AU_TO_WN * geom * (mu'.x10)_i * (mu'.x10)_j / r_bohr^3
  // =========================================================================
  RealType SFGTimeAvg::tdcCoupling(const Vector3d& r_ij,
                                   const Vector3d& e_i, const Vector3d& e_j,
                                   RealType mx_i, RealType mx_j) {
    RealType r_A = r_ij.length();
    if (r_A < 1.0e-6) return 0.0;
    Vector3d rhat = r_ij / r_A;
    const RealType A_to_bohr = 1.0 / 0.52917721067;
    RealType r_bohr = r_A * A_to_bohr;
    RealType r3 = r_bohr * r_bohr * r_bohr;
    const RealType AU_TO_WN = 219474.6313705;
    RealType geom = dot(e_i, e_j) - 3.0 * dot(e_i, rhat) * dot(e_j, rhat);
    return AU_TO_WN * geom * mx_i * mx_j / r3;
  }

}  // namespace OpenMD
