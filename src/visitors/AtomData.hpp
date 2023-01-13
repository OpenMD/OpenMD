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

#ifndef VISITOR_ATOMDATA_HPP
#define VISITOR_ATOMDATA_HPP

#include <memory>
#include <string>
#include <vector>

#include "math/Vector3.hpp"
#include "utils/GenericData.hpp"

namespace OpenMD {

  struct AtomInfo {
    AtomInfo() :
        pos(V3Zero), vec(V3Zero), vel(V3Zero), frc(V3Zero), eField(V3Zero),
        charge(0.0), hasCharge(false), hasVector(false), hasVelocity(false),
        hasForce(false), hasElectricField(false), hasGlobalID(false) {}

    std::string atomTypeName;
    int globalID;
    Vector3d pos;
    Vector3d vec;
    Vector3d vel;
    Vector3d frc;
    Vector3d eField;
    RealType charge;
    bool hasCharge;
    bool hasVector;
    bool hasVelocity;
    bool hasForce;
    bool hasElectricField;
    bool hasGlobalID;
  };

  class AtomData : public GenericData {
  public:
    AtomData(const std::string& id = "ATOMDATA") : GenericData(id) {}

    void addAtomInfo(std::shared_ptr<AtomInfo> info) { data.push_back(info); }

    void clearAllAtomInfo();

    std::shared_ptr<AtomInfo> beginAtomInfo(
        std::vector<std::shared_ptr<AtomInfo>>::iterator& i) {
      i = data.begin();
      return i != data.end() ? *i : nullptr;
    }

    std::shared_ptr<AtomInfo> nextAtomInfo(
        std::vector<std::shared_ptr<AtomInfo>>::iterator& i) {
      ++i;
      return i != data.end() ? *i : nullptr;
    }

    std::vector<std::shared_ptr<AtomInfo>> getData() { return data; }

    size_t getSize() { return data.size(); }

  protected:
    std::vector<std::shared_ptr<AtomInfo>> data;
  };
}  // namespace OpenMD

#endif  // VISITOR_ATOMDATA_HPP
