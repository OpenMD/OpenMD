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
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#ifndef TYPES_MAWINTERACTIONTYPE_HPP
#define TYPES_MAWINTERACTIONTYPE_HPP

#include "types/NonBondedInteractionType.hpp"

namespace OpenMD {
  /**
   * @class MAWInteractionType
   *

   * MAWInteractionType (Metal-Angular-Water) is one of the basic
   * Metal-to-NonMetal interaction types.
   *
   * \f[ V =  D_e * \exp(-a(r-r_e)) * (\exp(-a(r-r_e)) - 2) *
                         (1 + ca1*(1-\sqrt(3)*\cos(\theta))^2 +
                              cb1*3*(\sin(\theta)*\cos(\phi))^2) \f]

   * The spherical coordinates are defined in the body-fixed frame
   * of a rigid-body water molecule (HO bonds are on the Y-Z plane)
   * and the dipole vector of the water molecule points along the
   * Z-axis.  A metal atom's position is uniquely defined by a set
   * of spherical polar coordinates \f$(r, \theta, \phi)\f$ in the
   * body-fixed frame of each water molecule.
   */

  class MAWInteractionType : public NonBondedInteractionType {
  public:
    MAWInteractionType(RealType myD0, RealType myBeta0, RealType myR0,
                       RealType myCa1, RealType myCb1) {
      D_e  = myD0;
      beta = myBeta0;
      r_e  = myR0;
      ca1  = myCa1;
      cb1  = myCb1;
      setMAW();
    }
    RealType getD() { return D_e; }

    RealType getBeta() { return beta; }

    RealType getR() { return r_e; }

    RealType getCA1() { return ca1; }

    RealType getCB1() { return cb1; }

  private:
    RealType D_e;
    RealType beta;
    RealType r_e;
    RealType ca1;
    RealType cb1;
  };
}  // namespace OpenMD

#endif
