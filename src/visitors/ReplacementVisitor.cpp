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

#include "visitors/ReplacementVisitor.hpp"

#include <cstring>
#include <memory>

#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"

namespace OpenMD {

  void ReplacementVisitor::addReplacedAtomName(const std::string& repName) {
    myTypes_.insert(repName);
  }

  bool ReplacementVisitor::isReplacedAtom(const std::string& atomType) {
    std::set<std::string>::iterator strIter;
    strIter = myTypes_.find(atomType);
    return strIter != myTypes_.end() ? true : false;
  }

  void ReplacementVisitor::addSite(const std::string& name,
                                   const Vector3d& refPos) {
    std::shared_ptr<AtomInfo> atomInfo = std::make_shared<AtomInfo>();
    atomInfo->atomTypeName             = name;
    atomInfo->pos                      = refPos;
    sites_->addAtomInfo(atomInfo);
  }
  void ReplacementVisitor::addSite(const std::string& name,
                                   const Vector3d& refPos,
                                   const Vector3d& refVec) {
    std::shared_ptr<AtomInfo> atomInfo = std::make_shared<AtomInfo>();
    atomInfo->atomTypeName             = name;
    atomInfo->pos                      = refPos;
    atomInfo->vec                      = refVec;
    atomInfo->hasVector                = true;
    sites_->addAtomInfo(atomInfo);
  }

  void ReplacementVisitor::visit(DirectionalAtom* datom) {
    RotMat3x3d A;
    RotMat3x3d Atrans;
    Mat3x3d I;
    Vector3d pos;
    Vector3d vel;
    Vector3d frc;
    Vector3d trq;
    Vector3d j;
    Mat3x3d skewMat;

    Vector3d newVec;
    std::shared_ptr<AtomInfo> atomInfo;
    std::shared_ptr<AtomData> atomData;
    std::shared_ptr<GenericData> data;
    bool haveAtomData;

    // if atom is not one of our recognized atom types, just skip it
    if (!isReplacedAtom(datom->getType())) return;

    data = datom->getPropertyByName("ATOMDATA");

    if (data != nullptr) {
      atomData = std::dynamic_pointer_cast<AtomData>(data);

      if (atomData == nullptr) {
        std::cerr << "can not get Atom Data from " << datom->getType()
                  << std::endl;
        atomData     = std::make_shared<AtomData>();
        haveAtomData = false;
      } else
        haveAtomData = true;
    } else {
      atomData     = std::make_shared<AtomData>();
      haveAtomData = false;
    }

    pos = datom->getPos();
    vel = datom->getVel();

    j = datom->getJ();
    I = datom->getI();
    A = datom->getA();

    skewMat(0, 0) = 0;
    skewMat(0, 1) = j[2] / I(2, 2);
    skewMat(0, 2) = -j[1] / I(1, 1);
    skewMat(1, 0) = -j[2] / I(2, 2);
    skewMat(1, 1) = 0;
    skewMat(1, 2) = j[0] / I(0, 0);
    skewMat(2, 0) = j[1] / I(1, 1);
    skewMat(2, 1) = -j[0] / I(0, 0);
    skewMat(2, 2) = 0;
    Mat3x3d mat   = (A * skewMat).transpose();

    // We need A^T to convert from body-fixed to space-fixed:
    Atrans = A.transpose();

    std::shared_ptr<AtomInfo> siteInfo;
    std::vector<std::shared_ptr<AtomInfo>>::iterator iter;

    for (siteInfo = sites_->beginAtomInfo(iter); siteInfo;
         siteInfo = sites_->nextAtomInfo(iter)) {
      newVec = Atrans * siteInfo->pos;

      atomInfo               = std::make_shared<AtomInfo>();
      atomInfo->atomTypeName = siteInfo->atomTypeName;
      atomInfo->pos          = pos + newVec;

      if (siteInfo->hasVector) {
        newVec        = Atrans * siteInfo->vec;
        atomInfo->vec = newVec;
      } else {
        atomInfo->vec = V3Zero;
      }

      atomInfo->vel         = vel + mat * siteInfo->pos;
      atomInfo->hasVelocity = true;

      atomData->addAtomInfo(atomInfo);
    }
    if (!haveAtomData) {
      atomData->setID("ATOMDATA");
      datom->addProperty(atomData);
    }

    setVisited(datom);
  }

  const std::string ReplacementVisitor::toString() {
    char buffer[65535];
    std::string result;

    snprintf(
        buffer, 65535,
        "------------------------------------------------------------------\n");
    result += buffer;

    snprintf(buffer, 65535, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    snprintf(buffer, 65535,
             "Visitor Description: replace atom with other sites\n");
    result += buffer;

    snprintf(
        buffer, 65535,
        "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }
}  // namespace OpenMD
