/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
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
 
#ifndef VISITORS_OTHERVISITOR_HPP
#define VISITORS_OTHERVISITOR_HPP

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "visitors/BaseVisitor.hpp"
#include "primitives/StuntDouble.hpp"
#include "visitors/AtomData.hpp"
#include "selection/SelectionManager.hpp"
#include "selection/SelectionEvaluator.hpp"

namespace OpenMD {

  class SimInfo;

  class WrappingVisitor : public BaseVisitor{
  public:
    using BaseVisitor::visit; 
    WrappingVisitor(SimInfo* info, bool useCom = true) : BaseVisitor(), useCom_(useCom) {
      this->info = info;
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

  class ReplicateVisitor : public BaseVisitor{
  public:
    using BaseVisitor::visit; 
    ReplicateVisitor(SimInfo* info, Vector3i opt);
    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);

    virtual const std::string toString();
  protected:
    void internalVisit(StuntDouble* sd);
    void replicate(std::vector<AtomInfo*>& infoList,
		   std::shared_ptr<AtomData> data,
                   const Mat3x3d& box);
    
  private:
    std::vector<Vector3d> dir;
    SimInfo* info {nullptr};
    Vector3i replicateOpt;
  };

  class XYZVisitor : public BaseVisitor{
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
    void clear() {frame.clear();}
    void doPositions(bool pos) {doPositions_ = pos;}
    void doVelocities(bool vel) {doVelocities_ = vel;}
    void doForces(bool frc) {doForces_ = frc;}
    void doVectors(bool vec) {doVectors_ = vec;}
    void doCharges(bool chg) {doCharges_ = chg;}
    void doElectricFields(bool efl) {doElectricFields_ = efl;}
    void doGlobalIDs(bool gid) {doGlobalIDs_ = gid;}

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

  class PrepareVisitor : public BaseVisitor{
  public:
    using BaseVisitor::visit; 
    PrepareVisitor() : BaseVisitor() {visitorName = "prepareVisitor";}

    virtual void visit(Atom* atom) {internalVisit(atom);}
    virtual void visit(DirectionalAtom* datom) {internalVisit(reinterpret_cast<Atom*>(datom));}
    virtual void visit(RigidBody* rb) {internalVisit(rb);}

    virtual const std::string toString();

  protected:
    void internalVisit(Atom* atom);
    void internalVisit(RigidBody* rb);
  };

  class WaterTypeVisitor : public BaseVisitor{
  public:
    using BaseVisitor::visit; 
    WaterTypeVisitor() ;
    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom) {}
    virtual void visit(RigidBody* rb);
    
    virtual const std::string toString();
    
  private:
    std::string trimmedName(const std::string& atomType);
    
    std::set<std::string> waterTypeList;
  };
} // namespace OpenMD

#endif // _OTHERVISITOR_H_
