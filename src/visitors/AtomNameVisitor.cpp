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
#include "visitors/AtomNameVisitor.hpp"

#include <fstream>
#include <memory>
#include <sstream>

#include "brains/SimInfo.hpp"
#include "utils/StringTokenizer.hpp"

namespace OpenMD {
  AtomNameVisitor::AtomNameVisitor(SimInfo* info) : BaseVisitor(), info_(info) {
    visitorName = "AtomNameVisitor";
    ff_         = info_->getForceField();
  }

  void AtomNameVisitor::visitAtom(Atom* atom) {
    std::shared_ptr<AtomData> atomData;
    std::shared_ptr<GenericData> data = atom->getPropertyByName("ATOMDATA");

    if (data != nullptr) {
      atomData = std::dynamic_pointer_cast<AtomData>(data);
      if (atomData == nullptr) {
        std::cerr << "can not get Atom Data from " << atom->getType()
                  << std::endl;
        atomData = std::make_shared<AtomData>();
      }
    } else {
      atomData = std::make_shared<AtomData>();
    }

    std::vector<std::shared_ptr<AtomInfo>>::iterator i;
    for (std::shared_ptr<AtomInfo> atomInfo = atomData->beginAtomInfo(i);
         atomInfo != nullptr; atomInfo      = atomData->nextAtomInfo(i)) {
      // query the force field for the AtomType associated with this
      // atomTypeName:
      AtomType* at = ff_->getAtomType(atomInfo->atomTypeName);
      // Sometimes the atomInfo is set to a fictitious type, so we'll
      // only do base replacement if this AtomType is known in the forceField:
      if (at != NULL) {
        // get the chain of base types for this atom type:
        std::vector<AtomType*> ayb = at->allYourBase();
        // use the last type in the chain of base types for the name:
        std::string bn = ayb[ayb.size() - 1]->getName();

        atomInfo->atomTypeName = bn;
      }
    }
  }

  void AtomNameVisitor::visit(RigidBody* rb) {
    std::vector<Atom*>::iterator i;

    for (Atom* atom = rb->beginAtom(i); atom != NULL; atom = rb->nextAtom(i)) {
      visit(atom);
    }
  }

  const std::string AtomNameVisitor::toString() {
    char buffer[65535];
    std::string result;

    snprintf(
        buffer, 65535,
        "------------------------------------------------------------------\n");
    result += buffer;

    snprintf(buffer, 65535, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    snprintf(buffer, 65535, "Visitor Description: print base atom types\n");
    result += buffer;

    snprintf(
        buffer, 65535,
        "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

}  // namespace OpenMD
