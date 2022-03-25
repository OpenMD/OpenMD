/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#ifndef TYPES_EAMINTERACTIONTYPE_HPP
#define TYPES_EAMINTERACTIONTYPE_HPP

#include "types/NonBondedInteractionType.hpp"

namespace OpenMD {
  /**
   * @class EAMInteractionType
   *
   * EAMInteractionType is one of the basic metallic interactions for
   * representing the bonding in metallic solids. The basic functional
   * form has a non-pairwise density functional and a pair potential
   *
   * \f[ V = \sum_{i} F_i \left[\rho_i\right] + \sum_{i,j} \phi_{ij}(r_{ij}) \f]
   *
   * where the functional depends on a radially-decaying electron density,
   *
   * \f[ \rho_i = \sum_{j \neq i} \rho_{j}(r_{ij})  \f]
   */

  enum EAMiType { eamitTabulated, eamitZhou, eamitOxides };

  class EAMInteractionType : public NonBondedInteractionType {
  public:
    EAMInteractionType(RealType re, RealType alpha, RealType beta, RealType A,
                       RealType B, RealType kappa, RealType lambda) {
      interactionType_ = eamitZhou;
      setEAMZhou();
      re_     = re;
      alpha_  = alpha;
      beta_   = beta;
      A_      = A;
      B_      = B;
      kappa_  = kappa;
      lambda_ = lambda;
    }

    EAMInteractionType(RealType re, RealType alpha, RealType A, RealType Ci,
                       RealType Cj) {
      interactionType_ = eamitOxides;
      setEAMOxides();
      re_    = re;
      alpha_ = alpha;
      A_     = A;
      Ci_    = Ci;
      Cj_    = Cj;
    }

    int getNr() { return nr_; }
    RealType getDr() { return dr_; }
    RealType getRcut() { return rcut_; }
    std::vector<RealType> getPhi() { return phi_; }
    RealType getRe() { return re_; }
    RealType getA() { return A_; }
    RealType getB() { return B_; }
    RealType getAlpha() { return alpha_; }
    RealType getBeta() { return beta_; }
    RealType getKappa() { return kappa_; }
    RealType getLambda() { return lambda_; }
    EAMiType getInteractionType() { return interactionType_; }
    RealType getCi() { return Ci_; }
    RealType getCj() { return Cj_; }

  private:
    // This first set is for parameters read from DYNAMO 86 setfl files:
    int nr_;
    RealType dr_;
    RealType rcut_;
    std::vector<RealType> phi_;  // phi(r)
    // This set is for parameters for the parameterization of EAM described in:
    // Acta mater 49, 4005 (2001), and X. W. Zhou, R. A. Johnson, and
    // H. N. G. Wadley, Phys. Rev. B, 69, 144113 (2004).
    RealType re_;
    RealType alpha_;
    RealType beta_;
    RealType A_;
    RealType B_;
    RealType kappa_;
    RealType lambda_;
    // Extra parameters to support oxide potentials:
    RealType Ci_;
    RealType Cj_;
    EAMiType interactionType_;
  };
}  // namespace OpenMD

#endif
