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

#ifndef TYPES_HARMONICTORSIONTYPE_HPP
#define TYPES_HARMONICTORSIONTYPE_HPP

#if defined(_MSC_VER)
#define copysign _copysign
#endif

namespace OpenMD {

  /**
   * @class HarmonicTorsionType
   * These torsion types are defined identically with functional form given
   * in equation 5 in the following paper:
   *
   * "Transferable Potentials for Phase Equilibria. 4. United-Atom
   * Description of Linear and Branched Alkenes and Alkylbenzenes" by
   * Collin D. Wick, Marcus G. Martin and J. Ilja Siepmann,
   * J. Phys. Chem. B; 2000; 104(33) pp 8008 - 8016;
   *
   *    http://pubs.acs.org/doi/abs/10.1021/jp001044x
   *
   * This torsion potential has the form:
   *
   *  \f[
         V_{tors} = \frac{d_0}{2} \left(\phi - \phi_0\right)^2
       \f]
   *
   */
  class HarmonicTorsionType : public TorsionType {
  public:
    HarmonicTorsionType(RealType d0, RealType phi0) :
        TorsionType(), d0_(d0), phi0_(phi0) {}

    virtual void calcForce(RealType cosPhi, RealType& V, RealType& dVdCosPhi) {
      // check roundoff
      if (cosPhi > 1.0) {
        cosPhi = 1.0;
      } else if (cosPhi < -1.0) {
        cosPhi = -1.0;
      }

      RealType phi    = acos(cosPhi);
      RealType sinPhi = sqrt(1.0 - cosPhi * cosPhi);

      // trick to avoid divergence in angles near 0 and pi:

      if (fabs(sinPhi) < 1.0E-6) { sinPhi = copysign(1.0E-6, sinPhi); }

      V = 0.5 * d0_ * pow((phi - phi0_), 2);

      dVdCosPhi = -d0_ * (phi - phi0_) / sinPhi;
    }

    friend std::ostream& operator<<(std::ostream& os, HarmonicTorsionType& ttt);

  private:
    RealType d0_;
    RealType phi0_;
  };

  std::ostream& operator<<(std::ostream& os, HarmonicTorsionType& htt) {
    os << "This HarmonicTorsionType has below form:" << std::endl;
    os << htt.d0_ << "*(phi - " << htt.phi0_ << ")/2" << std::endl;
    return os;
  }

}  // namespace OpenMD

#endif  // TYPES_HARMONICTORSIONTYPE_HPP
