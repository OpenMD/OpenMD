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

#include "applications/dynamicProps/VelocityAutoOutProductCorrFunc.hpp"
#include "utils/Revision.hpp"
#include "math/SquareMatrix3.hpp"

namespace OpenMD {
  VelocityAutoOutProductCorrFunc::VelocityAutoOutProductCorrFunc(SimInfo* info,
                                       const std::string& filename,
                                       const std::string& sele1,
                                       const std::string& sele2)
    : ObjectACF<Mat3x3d>(info, filename, sele1, sele2,
                         DataStorage::dslVelocity){

    setCorrFuncType("Velocity - Velocity Auto Outer Product Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".vaOutProdcorr");

    velocity_.resize(nFrames_);
    sumVelocity_ = V3Zero;
    velocityCount_ = 0;
    propertyTemp = V3Zero;
  }

  int VelocityAutoOutProductCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    Vector3d vel = sd->getVel();
    propertyTemp = vel;

    velocity_[frame].push_back( propertyTemp );
    sumVelocity_ += propertyTemp;
    velocityCount_++;
    return velocity_[frame].size() - 1;
  }

  Mat3x3d VelocityAutoOutProductCorrFunc::calcCorrVal(int frame1, int frame2,
                                         int id1, int id2) {
    Mat3x3d tmpMat_1;
    tmpMat_1 = outProduct(velocity_[frame1][id1], velocity_[frame2][id2]);
    Mat3x3d tmpMat_2;
    tmpMat_2 = outProduct(velocity_[frame2][id2], velocity_[frame1][id1]);
    Mat3x3d tmpMat_3;
    tmpMat_3 = 0.5*(tmpMat_1 + tmpMat_2);
    return tmpMat_3;
  }

  void VelocityAutoOutProductCorrFunc::postCorrelate() {
    // Gets the average of the velocities
    sumVelocity_ /= RealType(velocityCount_);

    Mat3x3d correlationOfAverages_ = outProduct(sumVelocity_, sumVelocity_);

    for (unsigned int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= RealType(count_[i]);

        // The outerProduct correlation of the averages is subtracted
        // from the correlation value:
        histogram_[i] -= correlationOfAverages_;
      } else {
        histogram_[i] = M3Zero;
      }
    }
  }

}
