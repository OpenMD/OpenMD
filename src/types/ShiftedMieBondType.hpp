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

#ifndef TYPES_SHIFTEDMIEBONDTYPE_HPP
#define TYPES_SHIFTEDMIEBONDTYPE_HPP

#include "types/BondType.hpp"

namespace OpenMD {

  /**
   * @class ShiftedMieBondType
   *
   * @brief ShiftedMieBondType is used to correct 1-3 bend interactions in
   * the SDK force field.
   *
   * The functional form is given by:
   \f[
   V(r) = V_{Mie}(n,m,sigma,epsilon,r_{ij}) - V_{Mie}(n,m,sigma,epsilon,rs) for
   r_{ij} < r_{s} \f]
   * where
   \f[
   V_{Mie}(n,m,\sigma,\epsilon,r) =  \left(\frac{n}{n-m}\right)
   \left(\frac{n}{m}\right)^{m/(n-m)} \epsilon \left[\left( \frac{\sigma}{r}
   \right)^{n} - \left( \frac{\sigma}{r} \right)^{m}\right] \f]
   * and
   \f[
   r_{s} = e^{\frac{m \log (\sigma )+\log (m)-n \log (\sigma )-\log (n)}{m-n}}
   \f]
   */
  class ShiftedMieBondType : public BondType {
  public:
    ShiftedMieBondType(RealType mySigma, RealType myEpsilon, int myNrep,
                       int myMatt) :
        BondType(0.0) {
      sigma   = mySigma;
      epsilon = myEpsilon;
      n       = myNrep;
      m       = myMatt;

      rS_ = exp(log(pow(sigma, m - n) * RealType(m) / n));

      nmScale_ = n * pow(RealType(n) / m, RealType(m) / (n - m)) / (n - m);

      RealType rss = rS_ / sigma;
      RealType rsi = 1.0 / rss;
      RealType rsn = pow(rsi, n);
      RealType rsm = pow(rsi, m);

      potS_ = nmScale_ * epsilon * (rsn - rsm);

      setEquilibriumBondLength(rS_);
    }

    virtual void calcForce(RealType r, RealType& V, RealType& dVdr) {
      RealType ros, ri, rin, rim, rin1, rim1;

      ros  = r / sigma;
      ri   = 1.0 / ros;
      rin  = pow(ri, n);
      rim  = pow(ri, m);
      rin1 = rin * ri;
      rim1 = rim * ri;

      V    = nmScale_ * epsilon * (rin - rim) - potS_;
      dVdr = nmScale_ * epsilon * (-n * rin1 + m * rim1) / sigma;
    }

    RealType getSigma() { return sigma; }

    RealType getEpsilon() { return epsilon; }

    int getNrep() { return n; }
    int getMatt() { return m; }

  private:
    RealType sigma;
    RealType epsilon;
    int n;
    int m;
    RealType nmScale_;
    RealType rS_;
    RealType potS_;
  };
}  // namespace OpenMD

#endif
