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

/**
 * @file OplsTorsionType.hpp
 * @author    teng lin
 * @date  11/16/2004
 * @version 1.0
 */

#ifndef TYPES_OPLSTORSIONTYPE_HPP
#define TYPES_OPLSTORSIONTYPE_HPP

#include "types/PolynomialTorsionType.hpp"

namespace OpenMD {

  /**
   * @class OplsTorsionType OplsTorsionType.hpp "types/OplsTorsionType.hpp"
   * These torsion types are defined identically with functional form given
   * in the following paper:
   *
   * "Development and Testing of the OPLS All-Atom Force Field on
   * Conformational Energetics and Properties of Organic Liquids,"  by
   * William L. Jorgensen, David S. Maxwell, and Julian Tirado-Rives,
   * J. Am. Chem. Soc.; 1996; 118(45) pp 11225 - 11236;
   * DOI: 10.1021/ja9621760
   *
   * This torsion potential has the form:
   *
   *  Vtors = 0.5* (v1*(1+cos(phi)) + v2*(1-cos(2*phi)) + v3*(1+cos(3*phi)))
   *
   * Notes:
   *
   * 1) OpenMD converts internally to a Polynomial torsion type because
   *    all of the phase angles are zero in the OPLS paper.
   * 2) Coefficients are assumed to be in kcal / mol, and be careful about
   *    that factor of 1/2 when importing the coefficients!
   */
  class OplsTorsionType : public PolynomialTorsionType {
  public:
    OplsTorsionType(RealType v1, RealType v2, RealType v3, bool trans180) :
        PolynomialTorsionType(), v1_(v1), v2_(v2), v3_(v3) {
      // convert OPLS Torsion Type to Polynomial Torsion type
      RealType c0 = v2 + 0.5 * (v1 + v3);
      RealType c1 = 0.5 * (v1 - 3.0 * v3);
      RealType c2 = -v2;
      RealType c3 = 2.0 * v3;

      if (!trans180) {
        c1 = -c1;
        c3 = -c3;
      }

      setCoefficient(0, c0);
      setCoefficient(1, c1);
      setCoefficient(2, c2);
      setCoefficient(3, c3);
    }

    friend std::ostream& operator<<(std::ostream& os, OplsTorsionType& ott);

  private:
    RealType v1_;
    RealType v2_;
    RealType v3_;
  };

  std::ostream& operator<<(std::ostream& os, OplsTorsionType& ott) {
    os << "This OplsTorsionType has below form:" << std::endl;
    os << ott.v1_ << "/2*(1+cos(phi))"
       << " + " << ott.v2_ << "/2*(1-cos(2*phi))"
       << " + " << ott.v3_ << "/2*(1+cos(3*phi))" << std::endl;
    return os;
  }

}  // namespace OpenMD

#endif  // TYPES_OPLSTORSIONTYPE_HPP
