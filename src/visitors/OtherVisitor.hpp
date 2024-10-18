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

#ifndef VISITORS_OTHERVISITOR_HPP
#define VISITORS_OTHERVISITOR_HPP

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "primitives/StuntDouble.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "visitors/AtomData.hpp"
#include "visitors/BaseVisitor.hpp"

namespace OpenMD {

  class SimInfo;

  class WrappingVisitor : public BaseVisitor {
  public:
    using BaseVisitor::visit;
    WrappingVisitor(SimInfo* info, bool useCom = true) :
        BaseVisitor(), useCom_(useCom) {
      this->info  = info;
      visitorName = "WrappingVisitor";
    }
    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);

    virtual const std::string toString();

    virtual void update();

  private:
    void internalVisit(StuntDouble* sd);
    SimInfo* info {nullptr};
    Vector3d origin_;
    bool useCom_ {false};
  };

  class ReplicateVisitor : public BaseVisitor {
  public:
    using BaseVisitor::visit;
    ReplicateVisitor(SimInfo* info, Vector3i opt);
    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);

    virtual const std::string toString();

  protected:
    void internalVisit(StuntDouble* sd);
    void replicate(std::vector<std::shared_ptr<AtomInfo>>& infoList,
                   std::shared_ptr<AtomData> data, const Mat3x3d& box);

  private:
    std::vector<Vector3d> dir;
    SimInfo* info {nullptr};
    Vector3i replicateOpt;
  };

  class XYZVisitor : public BaseVisitor {
  public:
    using BaseVisitor::visit;

    XYZVisitor(SimInfo* info);

    XYZVisitor(SimInfo* info, const std::string& script);

    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);

    virtual void update();

    virtual const std::string toString();

    void writeFrame(std::ostream& outStream);
    void clear() { frame.clear(); }
    void doPositions(bool pos) { doPositions_ = pos; }
    void doVelocities(bool vel) { doVelocities_ = vel; }
    void doForces(bool frc) { doForces_ = frc; }
    void doVectors(bool vec) { doVectors_ = vec; }
    void doCharges(bool chg) { doCharges_ = chg; }
    void doElectricFields(bool efl) { doElectricFields_ = efl; }
    void doGlobalIDs(bool gid) { doGlobalIDs_ = gid; }

  protected:
    void internalVisit(StuntDouble* sd);
    bool isSelected(StuntDouble* sd);

  private:
    std::string trimmedName(const std::string& atomType);

    SimInfo* info {nullptr};
    SelectionManager seleMan;
    SelectionEvaluator evaluator;
    std::vector<std::string> frame;
    bool doPositions_ {false};
    bool doVelocities_ {false};
    bool doForces_ {false};
    bool doVectors_ {false};
    bool doCharges_ {false};
    bool doElectricFields_ {false};
    bool doGlobalIDs_ {false};
  };

  class PrepareVisitor : public BaseVisitor {
  public:
    using BaseVisitor::visit;
    PrepareVisitor() : BaseVisitor() { visitorName = "prepareVisitor"; }

    virtual void visit(Atom* atom) { internalVisit(atom); }
    virtual void visit(DirectionalAtom* datom) {
      internalVisit(reinterpret_cast<Atom*>(datom));
    }
    virtual void visit(RigidBody* rb) { internalVisit(rb); }

    virtual const std::string toString();

  protected:
    void internalVisit(Atom* atom);
    void internalVisit(RigidBody* rb);
  };

  class WaterTypeVisitor : public BaseVisitor {
  public:
    using BaseVisitor::visit;
    WaterTypeVisitor();
    virtual void visit(Atom*) {}
    virtual void visit(DirectionalAtom*) {}
    virtual void visit(RigidBody* rb);

    virtual const std::string toString();

  private:
    std::string trimmedName(const std::string& atomType);

    std::set<std::string> waterTypeList;
  };
}  // namespace OpenMD

#endif  // _OTHERVISITOR_H_
