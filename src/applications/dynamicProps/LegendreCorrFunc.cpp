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

#include "applications/dynamicProps/LegendreCorrFunc.hpp"
#include "math/LegendrePolynomial.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include <sstream>

namespace OpenMD {
  LegendreCorrFunc::LegendreCorrFunc(SimInfo* info,
                                     const std::string& filename,
                                     const std::string& sele1,
                                     const std::string& sele2, int order)
    : ObjectACF<Vector3d>(info, filename, sele1, sele2,
                          DataStorage::dslAmat), order_(order) {

    setCorrFuncType("Legendre Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".lcorr");

    std::stringstream params;
    params << " order = " << order_;
    const std::string paramString = params.str();
    setParameterString( paramString );
    
    setLabelString( "Pn(costheta_x)\tPn(costheta_y)\tPn(costheta_z)" );
    
    LegendrePolynomial polynomial(order);
    legendre_ = polynomial.getLegendrePolynomial(order);

    rotMats_.resize(nFrames_);
  }
  
  int LegendreCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    rotMats_[frame].push_back( sd->getA() );
    return rotMats_[frame].size() - 1;
  } 

  Vector3d LegendreCorrFunc::calcCorrVal(int frame1, int frame2,
                                         int id1, int id2) {
    
    // The lab frame vector corresponding to the body-fixed 
    // z-axis is simply the second column of A.transpose()
    // or, identically, the second row of A itself.
    // Similar identites give the 0th and 1st rows of A for
    // the lab vectors associated with body-fixed x- and y- axes.

    Vector3d v1x = rotMats_[frame1][id1].getRow(0);
    Vector3d v1y = rotMats_[frame1][id1].getRow(1);    
    Vector3d v1z = rotMats_[frame1][id1].getRow(2);

    Vector3d v2x = rotMats_[frame2][id2].getRow(0);
    Vector3d v2y = rotMats_[frame2][id2].getRow(1);    
    Vector3d v2z = rotMats_[frame2][id2].getRow(2);

    RealType ux = legendre_.evaluate(dot(v1x, v2x)/(v1x.length()*v2x.length()));
    RealType uy = legendre_.evaluate(dot(v1y, v2y)/(v1y.length()*v2y.length()));
    RealType uz = legendre_.evaluate(dot(v1z, v2z)/(v1z.length()*v2z.length()));

    return Vector3d(ux, uy, uz);
  }

  void LegendreCorrFunc::validateSelection(SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;
    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      if (!sd->isDirectional()) {
	sprintf(painCave.errMsg,
                "LegendreCorrFunc::validateSelection Error: "
                "at least one of the selected objects is not Directional\n");
	painCave.isFatal = 1;
	simError();        
      }
    }    
  }
}
