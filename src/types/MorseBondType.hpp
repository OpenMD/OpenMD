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

#ifndef TYPES_MORSEBONDTYPE_HPP
#define TYPES_MORSEBONDTYPE_HPP

#include <cmath>

#include "types/BondType.hpp"

namespace OpenMD {

  /**
 * @class MorseBondType
 *
 * @brief MorseBondType is a more realistic bond potential.
 *
 * The functional form is given by:
   \f[
      V(r) = D_e (1 - e^{-\beta (r - r_0)})^2
    \f]
 * where \f$D_e\f$ is the bond dissociation energy (in
 * kcal / mol), \f$\beta\f$ is an inverse distance parameter related
 * to the force constant. \f$\beta = \sqrt{\frac{k}{2 D_e}}\f$, and
 * \f$r_0\f$ is the equilibrium bond length.
 */
  class MorseBondType : public BondType {
  public:
    MorseBondType(RealType myR0, RealType myD, RealType myBeta) :
        BondType(myR0), De(myD), beta(myBeta) {}

    void setWellDepth(RealType myD) { De = myD; }

    void setBeta(RealType myBeta) { beta = myBeta; }

    void setWellDepthAndForceConstant(RealType myD, RealType myK) {
      De   = myD;
      beta = sqrt(myK / (2.0 * De));
    }

    RealType getWellDepth() { return De; }

    RealType getBeta() { return beta; }

    RealType getForceConstant() { return 2.0 * De * beta * beta; }

    virtual void calcForce(RealType r, RealType& V, RealType& dVdr) {
      RealType dr, eterm, eterm2;

      dr     = r - r0;
      eterm  = exp(-beta * dr);
      eterm2 = eterm * eterm;

      V    = De * (1 - 2.0 * eterm + eterm2);
      dVdr = 2.0 * De * beta * (eterm - eterm2);
    }

  private:
    RealType De;
    RealType beta;
  };
}  // namespace OpenMD

#endif
