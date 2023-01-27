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

#include "applications/staticProps/DipoleOrientation.hpp"

#include <string>
#include <vector>

#include "applications/staticProps/SpatialStatistics.hpp"
#include "brains/SimInfo.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/Accumulator.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  DipoleOrientation::DipoleOrientation(
      SimInfo* info, const std::string& filename, const std::string& sele,
      const RealType dipoleX, const RealType dipoleY, const RealType dipoleZ,
      int nzbins, int axis) :
      SlabStatistics(info, filename, sele, nzbins, axis),
      axis_(axis) {
    switch (axis_) {
    case 0:
      axisLabel_ = "x";
      refAxis_   = Vector3d(1, 0, 0);
      break;
    case 1:
      axisLabel_ = "y";
      refAxis_   = Vector3d(0, 1, 0);
      break;
    case 2:
    default:
      axisLabel_ = "z";
      refAxis_   = Vector3d(0, 0, 1);
      break;
    }
    setOutputName(getPrefix(filename) + ".Sz");

    dipoleVector_ = Vector3d(dipoleX, dipoleY, dipoleZ);
    dipoleVector_.normalize();

    orderS_               = new OutputData;
    orderS_->units        = "";
    orderS_->title        = "Orientational Order parameter";
    orderS_->dataType     = odtReal;
    orderS_->dataHandling = odhAverage;
    orderS_->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      orderS_->accumulator.push_back(new Accumulator());
    data_.push_back(orderS_);

    orderSCos_               = new OutputData;
    orderSCos_->units        = "";
    orderSCos_->title        = "Orientational Order parameter cosine Theta";
    orderSCos_->dataType     = odtReal;
    orderSCos_->dataHandling = odhAverage;
    orderSCos_->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      orderSCos_->accumulator.push_back(new Accumulator());
    data_.push_back(orderSCos_);
  }

  void DipoleOrientation::processFrame(int) {
    RealType z;

    hmat_ = currentSnapshot_->getHmat();

    for (unsigned int i = 0; i < nBins_; i++) {
      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(axis_, axis_);
      dynamic_cast<Accumulator*>(z_->accumulator[i])->add(z);
    }

    volume_ = currentSnapshot_->getVolume();

    StuntDouble* sd;
    int i;

    std::vector<RealType> binS(nBins_, 0.0);
    std::vector<RealType> binSCos(nBins_, 0.0);

    std::vector<int> count(nBins_, 0.0);

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // loop over the selected atoms:
    SquareMatrix3<RealType> rotMat;
    Vector3d rotatedDipoleVector;
    RealType ctheta;
    RealType orderParameter;

    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {
      // figure out where that object is:
      Vector3d pos = sd->getPos();

      int bin = getBin(pos);

      if (sd->isDirectional() || sd->isRigidBody()) {
        rotMat              = sd->getA();
        rotatedDipoleVector = rotMat * dipoleVector_;
        rotatedDipoleVector.normalize();
        ctheta = dot(rotatedDipoleVector, refAxis_);

        orderParameter = (3 * (ctheta * ctheta) - 1) / 2;

        binS[bin] += orderParameter;
        binSCos[bin] += ctheta;

        count[bin] += 1;
      }
    }

    for (unsigned int i = 0; i < nBins_; i++) {
      count[i] != 0 ?
          dynamic_cast<Accumulator*>(orderS_->accumulator[i])
              ->add(binS[i] / count[i]) :
          dynamic_cast<Accumulator*>(orderS_->accumulator[i])->add(binS[i]);
      count[i] != 0 ? dynamic_cast<Accumulator*>(orderSCos_->accumulator[i])
                          ->add(binSCos[i] / count[i]) :
                      dynamic_cast<Accumulator*>(orderSCos_->accumulator[i])
                          ->add(binSCos[i]);
    }
  }
}  // namespace OpenMD
