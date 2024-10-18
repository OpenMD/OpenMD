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

#ifndef MATH_INTEGRATION_LEGENDREGAUSS1D_HPP
#define MATH_INTEGRATION_LEGENDREGAUSS1D_HPP

namespace OpenMD {

  template<class returnType, class argumentType>
  class LegendreGauss1d {
  public:
    LegendreGauss1d(size_t points);
    virtual ~LegendreGauss1d();
    virtual returnType integrand(const argumentType& r) const = 0;
    virtual returnType integral(const argumentType& rmin,
                                const argumentType& rmax,
                                RealType length) const;

    std::vector<RealType> xi;
    std::vector<RealType> w;
    size_t points;

  protected:
  };

  struct LegendreGauss1dParams : public LegendreGauss1d<RealType, RealType> {
  public:
    LegendreGauss1dParams(size_t points) :
        LegendreGauss1dParams<RealType, RealType>(points) {};

    virtual RealType integrand(const RealType& r) const { return 0; };

  private:
  };

  template<class returnType, class argumentType>
  LegendreGauss1d<returnType, argumentType>::LegendreGauss1d(size_t p) :
      points(p) {
    xi.resize(points);
    w.resize(points);
    switch (points) {
    case 1:
      xi[0] = 0.0;
      w[0]  = 2.0;
      break;
    case 2:
      xi[0] = 1.0 / sqrt(3.);
      xi[1] = -xi[0];
      w[0]  = 1.0;
      w[1]  = w[0];
      break;
    case 3:
      xi[0] = 0.0;
      xi[1] = sqrt(0.6);
      xi[2] = -xi[1];
      w[0]  = 8.0 / 9.0;
      w[1]  = 5.0 / 9.0;
      w[2]  = w[1];
      break;
    case 4:
      xi[0] = 0.861136;
      xi[1] = -xi[0];
      xi[2] = 0.339981;
      xi[3] = -xi[2];
      w[0]  = 0.347855;
      w[1]  = w[0];
      w[2]  = 0.652145;
      w[3]  = w[2];
      break;
    case 5:
      xi[0] = 0;
      xi[1] = 0.906180;
      xi[2] = -xi[1];
      xi[3] = 0.538469;
      xi[4] = -xi[3];
      w[0]  = 0.568889;
      w[1]  = 0.236927;
      w[2]  = w[1];
      w[3]  = 0.478629;
      w[4]  = w[3];
      break;
    case 6:
      xi[0] = 0.932470;
      xi[1] = -xi[0];
      xi[2] = 0.661209;
      xi[3] = -xi[2];
      xi[4] = 0.238619;
      xi[5] = -xi[4];
      w[0]  = 0.171324;
      w[1]  = w[0];
      w[2]  = 0.360702;
      w[3]  = w[2];
      w[4]  = 0.467914;
      w[5]  = w[4];
      break;
    default:
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "LegendreGauss1d does not implement a %d point algorithm\n",
               points);
      painCave.isFatal = 1;
      ;
      simError();
    }
  }

  template<class returnType, class argumentType>
  LegendreGauss1d<returnType, argumentType>::~LegendreGauss1d() {
    xi.clear();
    xi.shrink_to_fit();
    w.clear();
    w.shrink_to_fit();
  }

  template<class returnType, class argumentType>
  returnType LegendreGauss1d<returnType, argumentType>::integral(
      const argumentType& rmin, const argumentType& rmax,
      RealType length) const {
    returnType s(0);
    for (size_t i = 0; i < points; i++) {
      argumentType r = 0.5 * ((xi[i] + 1.0) * rmax + (1.0 - xi[i]) * rmin);
      s +            = integrand(r) * w[i];
    }
    return 0.5 * length * s;
  }
}  // namespace OpenMD
#endif  // MATH_INTEGRATION_LEGENDREGAUSS1D_HPP
