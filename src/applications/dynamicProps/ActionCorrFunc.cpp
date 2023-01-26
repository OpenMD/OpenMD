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

#include "applications/dynamicProps/ActionCorrFunc.hpp"

#include <string>

#include "applications/dynamicProps/TimeCorrFunc.hpp"
#include "brains/ForceManager.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/Constants.hpp"
#include "utils/Revision.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  // We need all of the positions, velocities, etc. so that we can
  // recalculate pressures and actions on the fly:
  ActionCorrFunc::ActionCorrFunc(SimInfo* info, const std::string& filename,
                                 const std::string& sele1,
                                 const std::string& sele2) :
      SystemACF<Mat3x3d>(info, filename, sele1, sele2) {
    setCorrFuncType("ActionCorrFunc");
    setOutputName(getPrefix(dumpFilename_) + ".action");
    setLabelString("Txx\tTxy\tTxz\tTyx\tTyy\tTyz\tTzx\tTzy\tTzz");

    // We'll need the force manager to compute forces for the average pressure
    forceMan_ = new ForceManager(info);

    // We'll need thermo to compute the pressures from the virial
    thermo_ = new Thermo(info);

    action_.resize(nTimeBins_);
    time_.resize(nTimeBins_);
  }

  void ActionCorrFunc::computeProperty1(int frame) {
    forceMan_->calcForces();
    RealType vol      = thermo_->getVolume();
    RealType pressure = thermo_->getPressure() / Constants::pressureConvert;

    int i;
    StuntDouble* sd;

    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      Vector3d r = sd->getPos(frame);
      Vector3d v = sd->getVel(frame);
      RealType m = sd->getMass();

      action_[frame] += m * outProduct(r, v);
    }
    action_[frame] /= vol;
    time_[frame] = info_->getSnapshotManager()->getCurrentSnapshot()->getTime();
    pressure_.add(pressure);
  }

  Mat3x3d ActionCorrFunc::calcCorrVal(int frame1, int frame2) {
    Mat3x3d corrTensor(0.0);
    RealType thisTerm;

    RealType pAve = pressure_.getAverage();

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        if (i == j) {
          thisTerm = (action_[frame2](i, j) - action_[frame1](i, j) -
                      pAve * (time_[frame2] - time_[frame1]));
        } else {
          thisTerm = (action_[frame2](i, j) - action_[frame1](i, j));
        }
        corrTensor(i, j) += thisTerm * thisTerm;
      }
    }
    return corrTensor;
  }
}  // namespace OpenMD
