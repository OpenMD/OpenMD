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

#include "visitors/CompositeVisitor.hpp"

#include <cstring>

#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"

namespace OpenMD {

  CompositeVisitor::~CompositeVisitor() {
    VisitorIterator i;
    BaseVisitor* curVisitor;

    for (curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i))
      delete curVisitor;

    visitorList.clear();
  }
  void CompositeVisitor::addVisitor(BaseVisitor* newVisitor, int priority) {
    VisitorIterator i;
    int curPriority;

    for (i = visitorList.begin(); i != visitorList.end(); ++i) {
      curPriority = (*i).second;
      // if new visitor has higher priority, just insert it before current
      // visitor
      if (priority > curPriority) {
        visitorList.insert(i, std::make_pair(newVisitor, priority));
      }
    }

    // if new visitor has lowest priority, insert it at the end of the list
    visitorList.insert(visitorList.end(), std::make_pair(newVisitor, priority));
  }

  BaseVisitor* CompositeVisitor::beginVisitor(VisitorIterator& i) {
    i = visitorList.begin();
    return i != visitorList.end() ? (*i).first : NULL;
  }

  BaseVisitor* CompositeVisitor::nextVisitor(VisitorIterator& i) {
    ++i;
    return i != visitorList.end() ? (*i).first : NULL;
  }

  void CompositeVisitor::visit(Atom* atom) {
    VisitorIterator i;
    BaseVisitor* curVisitor;

    for (curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i))
      atom->accept(curVisitor);
  }

  void CompositeVisitor::visit(DirectionalAtom* datom) {
    VisitorIterator i;
    BaseVisitor* curVisitor;

    for (curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i))
      datom->accept(curVisitor);
  }
  void CompositeVisitor::visit(RigidBody* rb) {
    VisitorIterator i;
    BaseVisitor* curVisitor;
    std::vector<Atom*> myAtoms;
    std::vector<Atom*>::iterator atomIter;

    myAtoms = rb->getAtoms();

    for (curVisitor = beginVisitor(i); curVisitor;
         curVisitor = nextVisitor(i)) {
      rb->accept(curVisitor);

      for (atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
        (*atomIter)->accept(curVisitor);
    }
  }

  const std::string CompositeVisitor::toString() {
    VisitorIterator i;
    std::string result;
    char buffer[65535];

    sprintf(
        buffer,
        "******************************************************************\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(
        buffer,
        "Visitor Description: visitor manager  maintaining a priority  list\n");
    result += buffer;

    sprintf(buffer, "visitors in current priority list:\n");
    result += buffer;

    for (i = visitorList.begin(); i != visitorList.end(); ++i) {
      sprintf(buffer, "Priority = %d\tvisitor = %s\n", (*i).second,
              ((*i).first->getVisitorName()).c_str());
      result += buffer;
    }

    sprintf(buffer, "Detail information about every visitor:\n");

    for (i = visitorList.begin(); i != visitorList.end(); ++i)
      result += ((*i).first)->toString();

    sprintf(
        buffer,
        "******************************************************************\n");
    result += buffer;

    return result;
  }

  void CompositeVisitor::update() {
    VisitorIterator i;
    BaseVisitor* curVisitor;

    for (curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i))
      curVisitor->update();
  }

}  // namespace OpenMD
