/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
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

  enum EAMiType {
    eamitTabulated,
    eamitZhou
  };

  class EAMInteractionType : public NonBondedInteractionType {

  public:

    EAMInteractionType(RealType re, RealType alpha,
                       RealType beta, RealType A,
                       RealType B, RealType kappa,
                       RealType lambda) {

      interactionType_ = eamitZhou;
      setEAMZhou();
      re_ = re;
      alpha_ = alpha;
      beta_ = beta;
      A_ = A;
      B_ = B;
      kappa_ = kappa;
      lambda_ = lambda;
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

  private:
    // This first set is for parameters read from DYNAMO 86 setfl files:
    int nr_;
    RealType dr_;
    RealType rcut_;
    std::vector<RealType> phi_;   // phi(r)
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
    EAMiType interactionType_;
  };
}
#endif
