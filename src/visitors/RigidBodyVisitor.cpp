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

#include "visitors/RigidBodyVisitor.hpp"

#include <memory>

#include "primitives/RigidBody.hpp"

namespace OpenMD {

  void LipidHeadVisitor::visit(RigidBody* rb) {
    int globalID;
    Vector3d pos;
    Vector3d u(0, 0, 1);
    Vector3d newVec;
    std::shared_ptr<GenericData> data;
    std::shared_ptr<AtomData> atomData;
    std::shared_ptr<AtomInfo> atomInfo;
    bool haveAtomData;
    RotMat3x3d rotMatrix;

    if (!canVisit(rb->getType())) return;

    globalID  = rb->getGlobalIndex();
    pos       = rb->getPos();
    rotMatrix = rb->getA();
    // matVecMul3(rotMatrix, u, newVec);
    newVec = rotMatrix * u;

    data = rb->getPropertyByName("ATOMDATA");

    if (data != nullptr) {
      atomData = std::dynamic_pointer_cast<AtomData>(data);

      if (atomData == nullptr) {
        std::cerr << "can not get Atom Data from " << rb->getType()
                  << std::endl;

        atomData     = std::make_shared<AtomData>();
        haveAtomData = false;

      } else
        haveAtomData = true;

    } else {
      atomData     = std::make_shared<AtomData>();
      haveAtomData = false;
    }

    atomInfo               = std::make_shared<AtomInfo>();
    atomInfo->atomTypeName = "X";
    atomInfo->globalID     = globalID;
    atomInfo->pos[0]       = pos[0];
    atomInfo->pos[1]       = pos[1];
    atomInfo->pos[2]       = pos[2];
    atomInfo->vec[0]       = newVec[0];
    atomInfo->vec[1]       = newVec[1];
    atomInfo->vec[2]       = newVec[2];

    atomData->addAtomInfo(atomInfo);

    if (!haveAtomData) {
      atomData->setID("ATOMDATA");
      rb->addProperty(atomData);
    }
  }

  void LipidHeadVisitor::addLipidHeadName(const std::string& name) {
    lipidHeadName.insert(name);
  }

  bool LipidHeadVisitor::canVisit(const std::string& name) {
    return lipidHeadName.find(name) != lipidHeadName.end() ? true : false;
  }

  const std::string LipidHeadVisitor::toString() {
    char buffer[65535];
    std::string result;
    std::set<std::string>::iterator i;

    sprintf(
        buffer,
        "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    // print the ignore type list
    sprintf(buffer, "lipidHeadName list contains below types:\n");
    result += buffer;

    for (i = lipidHeadName.begin(); i != lipidHeadName.end(); ++i) {
      sprintf(buffer, "%s\t", i->c_str());
      result += buffer;
    }

    sprintf(buffer, "\n");
    result += buffer;

    sprintf(
        buffer,
        "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

  void RBCOMVisitor::visit(RigidBody* rb) {
    std::shared_ptr<AtomData> atomData;
    std::shared_ptr<AtomInfo> atomInfo;
    Vector3d pos;
    pos = rb->getPos();

    atomInfo               = std::make_shared<AtomInfo>();
    atomInfo->atomTypeName = "X";
    atomInfo->pos[0]       = pos[0];
    atomInfo->pos[1]       = pos[1];
    atomInfo->pos[2]       = pos[2];
    atomInfo->vec[0]       = 0;
    atomInfo->vec[1]       = 0;
    atomInfo->vec[2]       = 0;

    atomData = std::make_shared<AtomData>();
    atomData->setID("ATOMDATA");
    atomData->addAtomInfo(atomInfo);

    rb->addProperty(atomData);
  }

  const std::string RBCOMVisitor::toString() {
    char buffer[65535];
    std::string result;

    sprintf(
        buffer,
        "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    // print the ignore type list
    sprintf(buffer, "Visitor Description: add a pseudo atom at the center of "
                    "the mass of the "
                    "rigidbody\n");
    result += buffer;

    sprintf(
        buffer,
        "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

}  // namespace OpenMD
