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

#include "visitors/LipidTransVisitor.hpp"

#include <memory>

#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {
  LipidTransVisitor::LipidTransVisitor(SimInfo* info,
                                       const std::string& originSeleScript,
                                       const std::string& refSeleScript) :
      BaseVisitor(),
      info_(info), originEvaluator_(info), originSeleMan_(info),
      refEvaluator_(info), refSeleMan_(info), refSd_(NULL) {
    visitorName = "LipidTransVisitor";

    originEvaluator_.loadScriptString(originSeleScript);
    if (!originEvaluator_.isDynamic()) {
      originSeleMan_.setSelectionSet(originEvaluator_.evaluate());
      if (originSeleMan_.getSelectionCount() == 1) {
        int i;
        originDatom_ =
            dynamic_cast<DirectionalAtom*>(originSeleMan_.beginSelected(i));
        if (originDatom_ == NULL) {
          sprintf(painCave.errMsg,
                  "LipidTransVisitor: origin selection must select an "
                  "directional atom");
          painCave.isFatal = 1;
          simError();
        }
      } else {
        sprintf(
            painCave.errMsg,
            "LipidTransVisitor: origin selection must select an directional "
            "atom");
        painCave.isFatal = 1;
        simError();
      }
    }

    refEvaluator_.loadScriptString(refSeleScript);
    if (!refEvaluator_.isDynamic()) {
      refSeleMan_.setSelectionSet(refEvaluator_.evaluate());
      if (refSeleMan_.getSelectionCount() == 1) {
        int i;
        refSd_ = refSeleMan_.beginSelected(i);

      } else {
        // error
      }
    }
  }

  void LipidTransVisitor::update() {
    Vector3d ref = refSd_->getPos();
    origin_      = originDatom_->getPos();
    Vector3d v1  = ref - origin_;
    info_->getSnapshotManager()->getCurrentSnapshot()->wrapVector(v1);

    MultipoleAdapter ma = MultipoleAdapter(originDatom_->getAtomType());
    Vector3d zaxis;
    if (ma.isDipole()) {
      zaxis = originDatom_->getDipole();
    } else {
      zaxis = originDatom_->getA().transpose() * V3Z;
    }

    Vector3d xaxis = cross(v1, zaxis);
    Vector3d yaxis = cross(zaxis, xaxis);

    xaxis.normalize();
    yaxis.normalize();
    zaxis.normalize();

    rotMat_.setRow(0, xaxis);
    rotMat_.setRow(1, yaxis);
    rotMat_.setRow(2, zaxis);
  }

  void LipidTransVisitor::internalVisit(StuntDouble* sd) {
    std::shared_ptr<GenericData> data;
    std::shared_ptr<AtomData> atomData;
    std::shared_ptr<AtomInfo> atomInfo;
    std::vector<std::shared_ptr<AtomInfo>>::iterator i;

    data = sd->getPropertyByName("ATOMDATA");

    if (data != nullptr) {
      atomData = std::dynamic_pointer_cast<AtomData>(data);

      if (atomData == nullptr) return;
    } else
      return;

    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();

    for (atomInfo = atomData->beginAtomInfo(i); atomInfo;
         atomInfo = atomData->nextAtomInfo(i)) {
      Vector3d tmp = atomInfo->pos - origin_;
      currSnapshot->wrapVector(tmp);
      atomInfo->pos = rotMat_ * tmp;
      ;
      atomInfo->vec = rotMat_ * atomInfo->vec;
    }
  }

  const std::string LipidTransVisitor::toString() {
    char buffer[65535];
    std::string result;

    sprintf(
        buffer,
        "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer, "Visitor Description: rotate the whole system\n");
    result += buffer;

    sprintf(
        buffer,
        "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

}  // namespace OpenMD
