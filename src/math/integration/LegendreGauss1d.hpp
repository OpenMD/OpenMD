/*
 * Copyright (c) 2004-2022 The University of Notre Dame. All Rights Reserved.
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

#ifndef MATH_INTEGRATION_LEGENDREGAUSS1D_HPP
#define MATH_INTEGRATION_LEGENDREGAUSS1D_HPP

namespace OpenMD {

  template<class returnType, class argumentType>
  class LegendreGauss1d{
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

  struct LegendreGauss1dParams: public LegendreGauss1d<RealType,RealType>{
  public:
    LegendreGauss1dParams(size_t points): LegendreGauss1dParams<RealType,RealType>(points){};

    virtual RealType integrand(const RealType& r) const{return 0;};
  private:
  };
  
  template<class returnType,class argumentType>
  LegendreGauss1d<returnType, argumentType>::LegendreGauss1d(size_t p):points(p){
    xi.resize(points);
    w.resize(points);
    switch(points) {
    case 1:
      xi[0] = 0.0;
      w[0] = 2.0;
      break;
    case 2:
      xi[0] = 1.0/sqrt(3.);
      xi[1] = -xi[0];
      w[0] = 1.0;
      w[1] = w[0];
      break;
    case 3:
      xi[0] = 0.0;
      xi[1] = sqrt(0.6);
      xi[2] = -xi[1];
      w[0] = 8.0/9.0;
      w[1] = 5.0/9.0;
      w[2] = w[1];
      break;
    case 4:
      xi[0] = 0.861136;
      xi[1] = -xi[0];
      xi[2] = 0.339981;
      xi[3] = -xi[2];
      w[0] = 0.347855;
      w[1] = w[0];
      w[2] = 0.652145;
      w[3] = w[2];
      break;
    case 5:
      xi[0] = 0;
      xi[1] = 0.906180;
      xi[2] = -xi[1];
      xi[3] = 0.538469;
      xi[4] = -xi[3];
      w[0] = 0.568889;
      w[1] = 0.236927;
      w[2] = w[1];
      w[3] = 0.478629;
      w[4] = w[3];
      break;
    case 6:
      xi[0] = 0.932470;
      xi[1] = -xi[0];
      xi[2] = 0.661209;
      xi[3] = -xi[2];
      xi[4] = 0.238619;
      xi[5] = -xi[4];
      w[0] = 0.171324;
      w[1] = w[0];
      w[2] = 0.360702;
      w[3] = w[2];
      w[4] = 0.467914;
      w[5] = w[4];
      break;
    default:
      sprintf(painCave.errMsg,
              "LegendreGauss1d does not implement a %d point algorithm\n",
              points);
      painCave.isFatal = 1;;
      simError();
    }
  }
  
  template<class returnType,class argumentType>
  LegendreGauss1d<returnType,argumentType>::~LegendreGauss1d(){
    xi.clear();
    xi.shrink_to_fit();
    w.clear();
    w.shrink_to_fit();
  }
  
  template<class returnType, class argumentType>
  returnType LegendreGauss1d<returnType, argumentType>::
  integral(const argumentType& rmin, const argumentType& rmax,
           RealType length) const{
    returnType s(0);
    for(size_t i = 0; i < points; i++) {
      argumentType r = 0.5*( (xi[i]+1.0) * rmax + (1.0-xi[i]) * rmin);
      s + = integrand(r) * w[i];
    } 
    return 0.5 * length * s;
  }
}
#endif // MATH_INTEGRATION_LEGENDREGAUSS1D_HPP
