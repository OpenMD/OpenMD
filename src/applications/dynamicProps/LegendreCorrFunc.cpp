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

#include "applications/dynamicProps/LegendreCorrFunc.hpp"

#include <sstream>

#include "math/LegendrePolynomial.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {
  LegendreCorrFunc::LegendreCorrFunc(SimInfo* info, const std::string& filename,
                                     const std::string& sele1,
                                     const std::string& sele2, int order) :
      ObjectACF<Vector3d>(info, filename, sele1, sele2),
    order_(order), doOffset_(false) {
    setCorrFuncType("Legendre Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".lcorr");

    std::stringstream params;
    params << " order = " << order_;
    const std::string paramString = params.str();
    setParameterString(paramString);

    setLabelString("Pn(costheta_x)\tPn(costheta_y)\tPn(costheta_z)");

    LegendrePolynomial polynomial(order);
    legendre_ = polynomial.getLegendrePolynomial(order);

    rotMats_.resize(nFrames_);
  }

  LegendreCorrFunc::LegendreCorrFunc(SimInfo* info, const std::string& filename,
                                     const std::string& sele1,
                                     const std::string& sele2,
				     int seleOffset,
				     int order) :
    ObjectACF<Vector3d>(info, filename, sele1, sele2), doOffset_(true),
    seleOffset_(seleOffset), order_(order) {
    setCorrFuncType("Legendre Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".lcorr");

    std::stringstream params;
    params << " order = " << order_;
    params << " seleoffset = " << seleOffset_;
    const std::string paramString = params.str();
    setParameterString(paramString);

    setLabelString("Pn(costheta_x)\tPn(costheta_y)\tPn(costheta_z)");

    LegendrePolynomial polynomial(order);
    legendre_ = polynomial.getLegendrePolynomial(order);

    vectors_.resize(nFrames_);
  }
  
  int LegendreCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    if (doOffset_) {
      int sd2Index = sd->getGlobalIndex() + seleOffset_;
      StuntDouble* sd2= info_->getIOIndexToIntegrableObject(sd2Index);
      vectors_[frame].push_back(sd->getPos() - sd2->getPos());
      return vectors_[frame].size() - 1;
    } else {
      rotMats_[frame].push_back(sd->getA());
      return rotMats_[frame].size() - 1;
    }
  }

  Vector3d LegendreCorrFunc::calcCorrVal(int frame1, int frame2, int id1,
                                         int id2) {

    if (doOffset_) {
      Vector3d v1z = vectors_[frame1][id1];
      Vector3d v2z = vectors_[frame2][id2];
      RealType uz =
        legendre_.evaluate(dot(v1z, v2z) / (v1z.length() * v2z.length()));
      return Vector3d(0, 0, uz);
      
    } else {

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
      
      RealType ux =
        legendre_.evaluate(dot(v1x, v2x) / (v1x.length() * v2x.length()));
      RealType uy =
        legendre_.evaluate(dot(v1y, v2y) / (v1y.length() * v2y.length()));
      RealType uz =
        legendre_.evaluate(dot(v1z, v2z) / (v1z.length() * v2z.length()));
      
      return Vector3d(ux, uy, uz);
    }
  }

  void LegendreCorrFunc::validateSelection(SelectionManager&) {
    StuntDouble* sd;
    int i;
    if (!doOffset_) {
      for (sd = seleMan1_.beginSelected(i); sd != NULL;
	   sd = seleMan1_.nextSelected(i)) {
	if (!sd->isDirectional()) {
	  snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		   "LegendreCorrFunc::validateSelection Error: "
		   "at least one of the selected objects is not Directional\n");
	  painCave.isFatal = 1;
	  simError();
	}
      }
    }
  }
}  // namespace OpenMD
