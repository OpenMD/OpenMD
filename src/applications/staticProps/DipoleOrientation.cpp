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
#include <vector>

#include "applications/staticProps/DipoleOrientation.hpp"
#include "applications/staticProps/SpatialStatistics.hpp"
#include "brains/SimInfo.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/Accumulator.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  DipoleOrientation::DipoleOrientation(SimInfo* info, const std::string& filename, const std::string& sele,
				       const RealType dipoleX, const RealType dipoleY, const RealType dipoleZ,
				       int nzbins, int axis)
    : SlabStatistics(info, filename, sele, nzbins, axis), axis_(axis) {

    switch(axis_) {
    case 0:
      axisLabel_ = "x";
      refAxis_ = Vector3d(1,0,0);
      break;
    case 1:
      axisLabel_ = "y";
      refAxis_ = Vector3d(0,1,0);
      break;
    case 2:
    default:
      axisLabel_ = "z";
      refAxis_ = Vector3d(0,0,1);
      break;
    }
    setOutputName(getPrefix(filename) + ".Sz");

    dipoleVector_ = Vector3d(dipoleX, dipoleY, dipoleZ);
    dipoleVector_.normalize();

    orderS_ = new OutputData;
    orderS_->units =  "";
    orderS_->title =  "Orientational Order parameter";
    orderS_->dataType = odtReal;
    orderS_->dataHandling = odhAverage;
    orderS_->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      orderS_->accumulator.push_back( new Accumulator() );
    data_.push_back(orderS_);


    orderSCos_ = new OutputData;
    orderSCos_->units =  "";
    orderSCos_->title =  "Orientational Order parameter cosine Theta";
    orderSCos_->dataType = odtReal;
    orderSCos_->dataHandling = odhAverage;
    orderSCos_->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      orderSCos_->accumulator.push_back( new Accumulator() );
    data_.push_back(orderSCos_);
  }

  void DipoleOrientation::processFrame(int istep) {

    RealType z;

    hmat_ = currentSnapshot_->getHmat();

    for (unsigned int i = 0; i < nBins_; i++) {
      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(axis_,axis_);
      dynamic_cast<Accumulator*>(z_->accumulator[i])->add(z);
    }

    volume_ = currentSnapshot_->getVolume();

    StuntDouble* sd;
    int i;

    std::vector<RealType> binS(nBins_, 0.0);
    std::vector<RealType> binSCos(nBins_, 0.0);

    std::vector<int> count(nBins_,0.0);

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
        rotMat = sd->getA();
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
      count[i] !=0 ? dynamic_cast<Accumulator *>(orderS_->accumulator[i])->add(binS[i]/count[i]) : dynamic_cast<Accumulator *>(orderS_->accumulator[i])->add(binS[i]);
      count[i] !=0 ? dynamic_cast<Accumulator *>(orderSCos_->accumulator[i])->add(binSCos[i]/count[i]) : dynamic_cast<Accumulator *>(orderSCos_->accumulator[i])->add(binSCos[i]);
    }
  }
}
