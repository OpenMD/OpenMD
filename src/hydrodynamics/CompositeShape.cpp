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

#include "hydrodynamics/CompositeShape.hpp"

#include "hydrodynamics/HydrodynamicsModel.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {

  CompositeShape::CompositeShape() { origin_ = V3Zero; }

  CompositeShape::~CompositeShape() { Utils::deletePointers(shapes_); }

  bool CompositeShape::isInterior(Vector3d pos) {
    bool result = false;
    std::vector<Shape*>::iterator iter;
    for (iter = shapes_.begin(); iter != shapes_.end(); ++iter) {
      if ((*iter)->isInterior(pos)) {
        result = true;
        break;
      }
    }

    return result;
  }

  template<class Cont, class Predict>
  void swap_if(Cont& b1, Cont& b2, Predict predict) {
    unsigned int size = b1.size();
    assert(size == b2.size());
    for (unsigned int i = 0; i < size; ++i) {
      if (predict(b1[i], b2[i])) std::swap(b1[i], b2[i]);
    }
  }

  std::pair<Vector3d, Vector3d> CompositeShape::getBoundingBox() {
    std::vector<Shape*>::iterator iter     = shapes_.begin();
    std::pair<Vector3d, Vector3d> boundary = (*iter)->getBoundingBox();
    for (++iter; iter != shapes_.end(); ++iter) {
      std::pair<Vector3d, Vector3d> currBoundary = (*iter)->getBoundingBox();
      swap_if(boundary.first, currBoundary.first, std::greater<RealType>());
      swap_if(boundary.second, currBoundary.second, std::less<RealType>());
    }

    return boundary;
  }

  HydroProp* CompositeShape::getHydroProp(RealType viscosity) {
    HydroProp* props = new HydroProp();
    props->setCenterOfResistance(V3Zero);
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "CompositeShape was asked to return an analytic HydroProps.\n");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal  = 1;
    simError();
    return props;
  }
}  // namespace OpenMD
