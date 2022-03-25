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

#ifndef VISITORS_REPLACEMENTVISITOR_HPP
#define VISITORS_REPLACEMENTVISITOR_HPP

#include <memory>
#include <set>

#include "visitors/AtomVisitor.hpp"

namespace OpenMD {

  /**
   * @class ReplacementVisitor
   *
   * Replaces an atomic object with a collection atomic sites.  These
   * sites are specified with reference location to the object, as well as
   * a name.
   */
  class ReplacementVisitor : public BaseAtomVisitor {
  public:
    using BaseVisitor::visit;
    ReplacementVisitor(SimInfo* info) : BaseAtomVisitor(info) {
      visitorName = "ReplacementVisitor";
      sites_      = std::make_shared<AtomData>();
    }

    void visit(Atom*) {}
    void visit(DirectionalAtom* datom);
    void visit(RigidBody*) {}

    const std::string toString();

    void addReplacedAtomName(const std::string& repName);
    void addSite(const std::string& name, const Vector3d& refPos);
    void addSite(const std::string& name, const Vector3d& refPos,
                 const Vector3d& refVec);

  private:
    inline bool isReplacedAtom(const std::string& atomType);
    std::set<std::string> myTypes_;
    std::shared_ptr<AtomData> sites_;
  };

  class SSDAtomVisitor : public ReplacementVisitor {
  public:
    using BaseVisitor::visit;
    SSDAtomVisitor(SimInfo* info) : ReplacementVisitor(info) {
      visitorName = "SSDAtomVisitor";

      /// these are the atom names we can replace with a fixed structure
      addReplacedAtomName("SSD");
      addReplacedAtomName("SSD_E");
      addReplacedAtomName("SSD_RF");
      addReplacedAtomName("SSD1");
      addReplacedAtomName("SSDQ");
      addReplacedAtomName("SSDQO");
      addReplacedAtomName("TAP");
      addReplacedAtomName("TRED");

      // this is the reference structure we'll use for the replacement:
      addSite("H", Vector3d(0.0, -0.75695, 0.5206));
      addSite("H", Vector3d(0.0, 0.75695, 0.5206));
      addSite("O", Vector3d(0.0, 0.0, -0.0654));
      addSite("X", Vector3d(0.0, 0.0, 0.0), Vector3d(0, 0, 1));
    }
  };

  class GBtailVisitor : public ReplacementVisitor {
  public:
    using BaseVisitor::visit;
    GBtailVisitor(SimInfo* info) : ReplacementVisitor(info) {
      visitorName = "GBtailVisitor";

      /// these are the atom names we can replace with a fixed structure
      addReplacedAtomName("GBtail");

      // this is the reference structure we'll use for the replacement:
      addSite("C", Vector3d(0.0, 0.0, 9.0));
      addSite("C", Vector3d(0.0, 0.0, 0.0));
      addSite("C", Vector3d(0.0, 0.0, -9.0));
    }
  };

  class GBheadVisitor : public ReplacementVisitor {
  public:
    using BaseVisitor::visit;
    GBheadVisitor(SimInfo* info) : ReplacementVisitor(info) {
      visitorName = "GBheadVisitor";

      /// these are the atom names we can replace with a fixed structure
      addReplacedAtomName("GBhead");

      // this is the reference structure we'll use for the replacement:
      addSite("N", Vector3d(0.0, 0.0, 3.5));
      addSite("C", Vector3d(0.0, 0.0, 0.0));
      addSite("P", Vector3d(0.0, 0.0, -3.5));
    }
  };
}  // namespace OpenMD

#endif
