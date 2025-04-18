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

#include "applications/staticProps/DipoleOrientation.hpp"

#include <string>
#include <vector>

#include "applications/staticProps/SpatialStatistics.hpp"
#include "brains/SimInfo.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/Accumulator.hpp"
#include "utils/AccumulatorView.hpp"
#include "utils/BaseAccumulator.hpp"
#include "utils/StringUtils.hpp"

using namespace OpenMD::Utils;

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

    // Pre-load the OutputData
    data_.resize(DipoleOrientation::ENDINDEX);

    OutputData z;
    z.units        = "Angstroms";
    z.title        = axisLabel_;
    z.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      z.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[Z] = std::move(z);

    OutputData orderS;
    orderS.units        = "";
    orderS.title        = "Orientational Order parameter";
    orderS.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      orderS.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[ORDERS] = std::move(orderS);

    OutputData orderSCos;
    orderSCos.units        = "";
    orderSCos.title        = "Orientational Order parameter cosine Theta";
    orderSCos.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      orderSCos.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[ORDERSCOS] = std::move(orderSCos);
  }

  void DipoleOrientation::processFrame(int) {
    RealType z;

    hmat_ = currentSnapshot_->getHmat();

    for (unsigned int i = 0; i < nBins_; i++) {
      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(axis_, axis_);
      data_[Z].accumulator[i]->add(z);
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
      count[i] != 0 ? data_[ORDERS].accumulator[i]->add(binS[i] / count[i]) :
                      data_[ORDERS].accumulator[i]->add(binS[i]);
      count[i] != 0 ?
          data_[ORDERSCOS].accumulator[i]->add(binSCos[i] / count[i]) :
          data_[ORDERSCOS].accumulator[i]->add(binSCos[i]);
    }
  }
}  // namespace OpenMD
