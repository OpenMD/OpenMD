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
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifndef APPLICATIONS_STATICPROPS_BONDORDERPARAMETER_HPP
#define APPLICATIONS_STATICPROPS_BONDORDERPARAMETER_HPP

#include "applications/staticProps/StaticAnalyser.hpp"
#include "math/SphericalHarmonic.hpp"
#include "math/Vector3.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  /**
   * @class BondOrderParameter
   * @brief Bond Order Parameter
   *
   * Computes orientational bond order parameters as outlined in:
   *
   *   "Bond-orientaional order in liquids and glasses," by
   *    P. J. Steinhart, D. R. Nelson, and M. Ronchetti,
   *    Phys. Rev. B, 28, 784 (1983).
   *
   * A somewhat more useful reference which has formulae for these order
   * parameters for individual atoms is:
   *
   *   "Numerical calculation of the rate of crystal nucleation in a
   *    Lennard-Jones system at moderate undercooling," by
   *    Pieter Rein ten Wolde, Maria J. Ruiz-Montero, and Daan Frenkel,
   *    J. Chem. Phys. 104, pp. 9932-9947 (1996).
   *
   * Note that this version uses a single cutoff radius to decide
   * membership in the list of neighbors, and does not have use a
   * distance-dependent weighting as used in the second reference above.
   *
   * The selection script can be utilized to look at specific types of
   * central atoms.  A dynamic selector can also be utilized.  By
   * default, this class computes the \f[ Q_{l} \f] and
   * \f[ \hat{W}_{l} \f] parameters up to \f[ l = 12 \f].  The completed
   * configurational averages of these values as well as the
   * distributions of atomic \f[ q_{l} \f] and \f[ \hat{w}_{l} \f]
   * values are then placed in .boq and .bow files.
   */
  class BondOrderParameter : public StaticAnalyser {
  public:
    BondOrderParameter(SimInfo* info, const std::string& filename,
                       const std::string& sele, double rCut, int nbins);

    virtual void process();

  private:
    virtual void initializeHistogram();
    virtual void collectHistogram(std::vector<RealType> q,
                                  std::vector<ComplexType> what);
    void writeOrderParameter(std::vector<RealType> Q,
                             std::vector<ComplexType> What);

    Snapshot* currentSnapshot_;
    std::string selectionScript_;
    SelectionManager seleMan_;
    SelectionEvaluator evaluator_;

    RealType rCut_;
    static const int lMax_ = 12;
    int frameCounter_;
    int nBins_;

    std::map<std::pair<int, int>, int> m2Min;
    std::map<std::pair<int, int>, int> m2Max;
    std::map<std::pair<int, int>, std::vector<RealType>> w3j;

    RealType MinQ_;
    RealType MaxQ_;
    RealType deltaQ_;
    std::vector<int> Qcount_;
    std::map<std::pair<int, int>, int> Q_histogram_;

    RealType MinW_;
    RealType MaxW_;
    RealType deltaW_;
    std::vector<int> Wcount_;
    std::map<std::pair<int, int>, int> W_histogram_;
  };
}  // namespace OpenMD

#endif
