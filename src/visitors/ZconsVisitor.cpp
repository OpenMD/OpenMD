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

#include "visitors/ZconsVisitor.hpp"

#include <cmath>
#include <memory>

#include "primitives/Molecule.hpp"
#include "types/ZconsStamp.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  ZConsVisitor::ZConsVisitor(SimInfo* info) :
      BaseVisitor(), zconsReader_(NULL), info_(info) {
    visitorName       = "ZConsVisitor";
    currSnapshot_     = info_->getSnapshotManager()->getCurrentSnapshot();
    Globals* simParam = info_->getSimParams();

    if (simParam->haveZconsTime()) {
      zconsTime_ = simParam->getZconsTime();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ZConstraint error: If you use a ZConstraint,\n"
               "\tyou must set zconsTime.\n");
      painCave.isFatal = 1;
      simError();
    }

    if (simParam->haveZconsTol()) {
      zconsTol_ = simParam->getZconsTol();
    } else {
      zconsTol_ = 0.01;
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ZConstraint Warning: Tolerance for z-constraint method is not "
               "specified.\n"
               "\tOpenMD will use a default value of %f.\n"
               "\tTo set the tolerance, use the zconsTol variable.\n",
               zconsTol_);
      painCave.isFatal = 0;
      simError();
    }

    int nZconstraints              = simParam->getNZconsStamps();
    std::vector<ZConsStamp*> stamp = simParam->getZconsStamps();
    for (int i = 0; i < nZconstraints; i++) {
      int zmolIndex = stamp[i]->getMolIndex();
      zmolStates_.insert(std::make_pair(zmolIndex, zsMoving));
    }

    // fill zatomToZmol_ array
    /** @todo only works for single version now*/
    std::map<int, ZConsState>::iterator j;
    for (j = zmolStates_.begin(); j != zmolStates_.end(); ++j) {
      Molecule* mol = info_->getMoleculeByGlobalIndex(j->first);
      assert(mol != NULL);
      Molecule::AtomIterator ai;
      Atom* at;
      for (at = mol->beginAtom(ai); at != NULL; at = mol->nextAtom(ai)) {
        zatomToZmol_.insert(
            std::make_pair(at->getGlobalIndex(), mol->getGlobalIndex()));
      }
    }

    zconsFilename_ = getPrefix(info_->getFinalConfigFileName()) + ".fz";

    zconsReader_ = new ZConsReader(info);

    if (zconsReader_->hasNextFrame()) zconsReader_->readNextFrame();
  }

  ZConsVisitor::~ZConsVisitor() { delete zconsReader_; }

  void ZConsVisitor::visit(Atom* atom) {
    std::string prefix;
    if (isZconstraint(atom->getGlobalIndex(), prefix))
      internalVisit(atom, prefix);
  }

  void ZConsVisitor::visit(DirectionalAtom* datom) {
    std::string prefix;

    if (isZconstraint(datom->getGlobalIndex(), prefix))
      internalVisit(datom, prefix);
  }

  void ZConsVisitor::visit(RigidBody* rb) {
    std::string prefix;
    std::vector<Atom*> atoms;

    atoms = rb->getAtoms();

    if (isZconstraint(atoms[0]->getGlobalIndex(), prefix))
      internalVisit(rb, prefix);
  }

  void ZConsVisitor::update() {
    Vector3d com;
    std::map<int, ZConsState>::iterator i;
    for (i = zmolStates_.begin(); i != zmolStates_.end(); ++i) {
      i->second = zsMoving;
    }

    readZconsFile(currSnapshot_->getTime());

    const std::vector<ZconsData>& fixedZmolData =
        zconsReader_->getFixedZMolData();
    std::vector<ZconsData>::const_iterator j;
    for (j = fixedZmolData.begin(); j != fixedZmolData.end(); ++j) {
      std::map<int, ZConsState>::iterator k = zmolStates_.find(j->zmolIndex);
      assert(k != zmolStates_.end());
      k->second = zsFixed;
    }
  }

  void ZConsVisitor::readZconsFile(RealType time) {
    RealType tempTime;
    while (zconsReader_->hasNextFrame()) {
      tempTime = zconsReader_->getCurTime();
      if (tempTime >= time) { return; }

      zconsReader_->readNextFrame();
    }
  }

  void ZConsVisitor::internalVisit(StuntDouble* sd, const std::string& prefix) {
    std::shared_ptr<GenericData> data;
    std::shared_ptr<AtomData> atomData;
    std::shared_ptr<AtomInfo> atomInfo;
    std::vector<std::shared_ptr<AtomInfo>>::iterator iter;

    // if there is not atom data, just skip it
    data = sd->getPropertyByName("ATOMDATA");
    if (data != nullptr) {
      atomData = std::dynamic_pointer_cast<AtomData>(data);
      if (atomData == nullptr) return;
    } else
      return;

    for (atomInfo = atomData->beginAtomInfo(iter); atomInfo;
         atomInfo = atomData->nextAtomInfo(iter))
      (atomInfo->atomTypeName).insert(0, prefix);
  }

  bool ZConsVisitor::isZconstraint(int atomIndex, std::string& prefix) {
    std::string prefixString[]     = {"ZF", "ZM"};
    std::map<int, int>::iterator i = zatomToZmol_.find(atomIndex);
    if (i == zatomToZmol_.end()) {
      prefix = "";
      return false;
    } else {
      std::map<int, ZConsState>::iterator j = zmolStates_.find(i->second);
      assert(j != zmolStates_.end());
      prefix = prefixString[j->second];
      return true;
    }
  }

  const std::string ZConsVisitor::toString() {
    char buffer[65535];
    std::string result;

    snprintf(
        buffer, 65535,
        "------------------------------------------------------------------\n");
    result += buffer;

    snprintf(buffer, 65535, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    snprintf(buffer, 65535, "number of zconstraint molecule: %d\n",
             (int)zmolStates_.size());
    result += buffer;

    snprintf(buffer, 65535, "zconstraint tolerance = %lf\n", zconsTol_);
    result += buffer;

    snprintf(buffer, 65535, "zconstraint sample time = %lf\n", zconsTime_);
    result += buffer;

    snprintf(buffer, 65535, "zconstraint output filename = %s\n",
             zconsFilename_.c_str());
    result += buffer;

    std::map<int, ZConsState>::iterator i;
    int j = 0;
    for (i = zmolStates_.begin(); i != zmolStates_.end(); ++i) {
      snprintf(buffer, 65535, "zconstraint molecule[%d] = %d\n", j++, i->first);
      result += buffer;
    }

    snprintf(
        buffer, 65535,
        "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

}  // namespace OpenMD
