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

#ifndef APPLICATIONS_STATICPROPS_SFGTIMEAVG_HPP
#define APPLICATIONS_STATICPROPS_SFGTIMEAVG_HPP

#include <complex>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "applications/staticProps/StaticAnalyser.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  /**
   * @class SFGTimeAvg
   *
   * Computes vibrational SFG spectra in the OH-stretch region using the
   * TIME-AVERAGING APPROXIMATION (TAA) of Auer & Skinner, rather than the
   * dynamical adiabatic-propagator method implemented in the SFG dynamic
   * property.  The TAA is a single pass through the dump file: at each
   * frame the instantaneous exciton Hamiltonian is built and diagonalized,
   * and each eigenstate contributes a complex Lorentzian to the spectrum
   * at its eigenfrequency, weighted by the eigenstate transition
   * polarizability and transition dipole.
   *
   *   chi2_pqr(w) = < sum_a (alpha_pq^a mu_r^a) / (w - w_a + i*Gamma) >
   *
   * where for eigenstate a:
   *   mu_r^a    = sum_i V(i,a) mu_{r,i}
   *   alpha_pq^a = sum_i V(i,a) alpha_{pq,i}
   *
   * The TAA neglects motional narrowing (it is the inhomogeneous limit),
   * but is far cheaper than the dynamical method (one diagonalization per
   * frame, no correlation windows, no propagation).  It shares the same
   * spectroscopic maps and Hamiltonian construction as the dynamical SFG
   * module so that the two are directly comparable.
   *
   * Fermi resonance with the HOH bend overtone is included by default
   * (fc = 50 cm-1); set fc = 0 to disable.
   */
  class SFGTimeAvg : public StaticAnalyser {
  public:
    SFGTimeAvg(SimInfo* info, const std::string& filename,
               const std::string& sele1, int nbins,
               const std::string& polarization = "ssp",
               int privilegedAxis = 2,
               RealType gamma = 20.0, RealType fc = 50.0);

    virtual void process();

  private:
    // -----------------------------------------------------------------------
    // Per-frame chromophore data in the local (site) basis.
    // Bend overtones (if Fermi resonance enabled) are appended after the
    // OH stretches: indices [0..nStretch-1] stretches,
    //                       [nStretch..N-1]  bend overtones.
    // -----------------------------------------------------------------------
    struct FrameData {
      int N {0};
      int nStretch {0};
      int nBend {0};
      DynamicRectMatrix<RealType> H;
      std::vector<Vector3d>  mu;
      std::vector<Mat3x3d>   alpha;
    };

    // Build the chromophore list and Hamiltonian diagonal for the current
    // snapshot.  Fills ohPos (TDC dipole positions), molIndex, and the
    // packed per-site array used by buildHamiltonian.
    FrameData extractFrame(std::vector<Vector3d>& ohPos,
                           std::vector<int>&      molIndex,
                           std::vector<RealType>& intraJ);

    void buildHamiltonian(FrameData&                   fd,
                          const std::vector<Vector3d>& ohPos,
                          const std::vector<int>&      molIndex,
                          const std::vector<RealType>& intraJ);

    // Diagonalize a real symmetric matrix; fill evals and row-major
    // eigenvectors evecs (evecs[i*N+a] = component i of eigenstate a).
    void diagonalize(const DynamicRectMatrix<RealType>& H, int N,
                     std::vector<double>& evals,
                     std::vector<double>& evecs);

    Mat3x3d bondPolarizability(const Vector3d& ohUnit,
                               RealType alpha_par, RealType alpha_perp,
                               RealType x10, RealType x10_gas);

    RealType tdcCoupling(const Vector3d& r_ij,
                         const Vector3d& e_i, const Vector3d& e_j,
                         RealType mx_i, RealType mx_j);

    void accumulateFrame(const FrameData& fd);
    void writeSpectrum();

    Snapshot* currentSnapshot_ {};

    std::string selectionScript1_;
    SelectionManager seleMan1_;
    SelectionEvaluator evaluator1_;

    // Spectroscopic maps (keyed by H atom-type name)
    std::map<std::string, std::tuple<RealType,RealType,RealType>> w10_;
    std::map<std::string, std::tuple<RealType,RealType,RealType>> muPrime_;
    std::map<std::string, std::pair<RealType,RealType>>           x10_;
    std::map<std::string, std::pair<RealType,RealType>>           p10_;
    std::map<std::string, std::tuple<RealType,RealType,RealType>> wintra_;
    std::map<std::string, std::pair<RealType,RealType>>           alphaMap_;
    std::map<std::string, RealType>                              tdcLoc_;

    // Bend maps (keyed by water-model name)
    std::map<std::string, std::pair<RealType,RealType>>           wb01_;
    std::map<std::string, std::pair<RealType,RealType>>           wb12_;

    std::string polarization_;
    int  privilegedAxis_ {2};
    int  s1_ {0};
    int  s2_ {1};
    RealType gamma_ {20.0};   // Lorentzian half-width [cm-1]
    RealType fc_    {50.0};   // Fermi coupling [cm-1]
    bool useFermi_  {false};

    int      nProcessed_ {0};
    RealType minFreq_ {2700.0};
    RealType maxFreq_ {4000.0};

    // Orientation diagnostics: accumulate <mu_z> by frequency class to
    // verify dangling vs donor OHs point in opposite directions.
    double   diagMuzFree_ {0.0};   // freq > 3650 (dangling)
    double   diagMuzHB_   {0.0};   // freq < 3450 (H-bonded donor)
    long     diagNFree_   {0};
    long     diagNHB_     {0};

    // Complex spectrum accumulators on the output grid
    std::vector<std::complex<double>> chi2_;

    Vector3d EF_;
    bool dumpHasElectricFields_;
  };
}  // namespace OpenMD

#endif
