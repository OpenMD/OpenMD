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

#ifndef MATH_INTEGRATION_LEGENDREGAUSSTRIANGLE_HPP
#define MATH_INTEGRATION_LEGENDREGAUSSTRIANGLE_HPP

#include "math/Vector3d.hpp"

namespace OpenMD {

  template<class returnType, class argumentType>
  class LegendreGaussTriangle{
  public:
    LegendreGaussTriangle(size_t points);
    virtual ~LegendreGaussTriangle();
    virtual returnType integrand(const argumentType& r) const = 0;
    virtual returnType integral(const argumentType& r0, const argumentType& r1,
                                const argumentType& r2, RealType area) const;

    std::vector<Vector3d> xi;
    std::vector<RealType> w;
    size_t points;
  protected:
  };

  struct LegendreGaussTriangleParams: public LegendreGaussTriangle<RealType,RealType>{
  public:
    LegendreGaussTriangleParams(size_t points): LegendreGaussTriangleParams<RealType,RealType>(points){};

    virtual RealType integrand(const RealType& r) const{return 0;};
  private:
  };
  
  template<class returnType,class argumentType>
  LegendreGaussTriangle<returnType, argumentType>::LegendreGaussTriangle(size_t p):points(p){
    xi.resize(points);
    w.resize(points);
    switch(points) {
    case 1:
      xi[0] = Vector3d(1/3, 1/3, 1/3);
      w[0] = 1.0;
      break;
    case 3:
      xi[0] = Vector3d(0.5, 0.5, 0.0);
      xi[1] = Vector3d(0.0, 0.5, 0.5);
      xi[2] = Vector3d(0.5, 0.0, 0.5);
      w[0] = 1/3;
      w[1] = 1/3;
      w[2] = 1/3;
      break;
    case 4:
      xi[0] = Vector3d(1./3,1./3,1./3);
      xi[1] = Vector3d(0.6,0.2,0.2);
      xi[2] = Vector3d(0.2,0.6,0.2);
      xi[3] = Vector3d(0.2,0.2,0.6);
      w[0] = -9/16;
      w[1] = 25/48;
      w[2] = w[1];
      w[3] = w[1];
      break;
    case 7:
      xi[0] = Vector3d(0.33333333, 0.33333333, 0.33333333);
      xi[1] = Vector3d(0.79742699, 0.10128651, 0.10128651);
      xi[2] = Vector3d(0.10128651, 0.79742699, 0.10128651);
      xi[3] = Vector3d(0.10128651, 0.10128651, 0.79742699);
      xi[4] = Vector3d(0.05971587, 0.47014206, 0.47014206);
      xi[5] = Vector3d(0.47014206, 0.05971587, 0.47014206);
      xi[6] = Vector3d(0.47014206, 0.47014206, 0.05971587);
      w[0] = 0.225;
      w[1] = 0.12593918;
      w[2] = w[1];
      w[3] = w[1];
      w[4] = 0.13239415;
      w[5] = w[4];
      w[6] = w[4];
      break;
    default:
      sprintf(painCave.errMsg,
              "LegendreGaussTriangle does not implement a %d point algorithm\n",
              points);
      painCave.isFatal = 1;;
      simError();
    }
  }
  
  template<class returnType,class argumentType>
  LegendreGaussTriangle<returnType,argumentType>::~LegendreGaussTriangle(){
    xi.clear();
    xi.shrink_to_fit();
    w.clear();
    w.shrink_to_fit();
  }
  
  template<class returnType, class argumentType>
  returnType LegendreGaussTriangle<returnType, argumentType>::
  integral(const argumentType& r0, const argumentType& r1,
           const argumentType& r2, RealType area) const{
    argumentType r01 = r1 - r0;
    argumentType r02 = r2 - r0;
    returnType s(0);
    for(size_t i = 0; i < points; i++) {
      argumentType r = xi[i].x * r01 + xi[i].y * r02 + r0;
      s + = integrand(r) * w[i];
    } 
    return area * s;
  }
}
#endif // MATH_INTEGRATION_LEGENDREGAUSSTRIANGLE_HPP
