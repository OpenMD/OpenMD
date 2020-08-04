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

#include <string>

#include "applications/dynamicProps/StressCorrFunc.hpp"
#include "applications/dynamicProps/TimeCorrFunc.hpp"
#include "brains/DataStorage.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "utils/Constants.hpp"
#include "utils/Revision.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  // We need all of the positions, velocities, etc. so that we can
  // recalculate pressures and actions on the fly:
  StressCorrFunc::StressCorrFunc(SimInfo* info, const std::string& filename,
				 const std::string& sele1,
                                 const std::string& sele2)
    : SystemACF<Mat3x3d>( info, filename, sele1, sele2,
                          DataStorage::dslPosition |
                          DataStorage::dslVelocity |
                          DataStorage::dslForce ) {

    setCorrFuncType("StressCorrFunc");
    setOutputName(getPrefix(dumpFilename_) + ".action");
    setLabelString( "Txx\tTxy\tTxz\tTyx\tTyy\tTyz\tTzx\tTzy\tTzz" );

    // We'll need the force manager to compute forces for the average pressure
    forceMan_ = new ForceManager(info_);

    // We'll need thermo to compute the pressures from the virial
    thermo_ =  new Thermo(info_);

    action_.resize(nTimeBins_);
    time_.resize(nTimeBins_);
  }

  void StressCorrFunc::computeProperty1(int frame) {

    forceMan_->calcForces();
    RealType vol = thermo_->getVolume();
    RealType pressure = thermo_->getPressure() / Constants::pressureConvert;

    int i;
    StuntDouble* sd;

    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {

      Vector3d r = sd->getPos();
      Vector3d v = sd->getVel();
      RealType m = sd->getMass();

      action_[frame] += m*outProduct(r, v);
    }

    action_[frame] /= vol;
    time_[frame] = info_->getSnapshotManager()->getCurrentSnapshot()->getTime();

    pressure_.add(pressure);
  }

  Mat3x3d StressCorrFunc::calcCorrVal(int frame1, int frame2) {

    Mat3x3d corrTensor(0.0);
    RealType thisTerm;

    RealType pAve = pressure_.getAverage();

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        if (i == j) {
          thisTerm = (action_[frame2](i, j) - action_[frame1](i, j)
                      - pAve * (time_[frame2] - time_[frame1]));
        } else {
          thisTerm = (action_[frame2](i, j) - action_[frame1](i, j));
        }
        corrTensor(i, j) += thisTerm * thisTerm;
      }
    }
    return corrTensor;
  }
}
