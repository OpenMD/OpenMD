/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef VISITORS_REPLACEMENTVISITOR_HPP
#define VISITORS_REPLACEMENTVISITOR_HPP

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
  class ReplacementVisitor : public BaseAtomVisitor{
  public:
    using BaseVisitor::visit;
    ReplacementVisitor(SimInfo* info) : BaseAtomVisitor(info) {
      visitorName = "ReplacementVisitor";   
      sites_ = new AtomData; 
    }
    ~ReplacementVisitor();
    
    void visit(Atom* atom) {}
    void visit(DirectionalAtom* datom);       
    void visit(RigidBody* rb) {}
    
    const std::string toString();

    void addReplacedAtomName(const std::string &repName);
    void addSite(const std::string &name, const Vector3d &refPos);
    void addSite(const std::string &name, const Vector3d &refPos, const Vector3d &refVec);
  private:
    inline bool isReplacedAtom(const std::string& atomType);
    std::set<std::string> myTypes_;
    AtomData* sites_;
  };
  
  class SSDAtomVisitor : public ReplacementVisitor{
  public:
    using BaseVisitor::visit;
    SSDAtomVisitor(SimInfo* info) : ReplacementVisitor(info) {
      visitorName = "SSDAtomVisitor";
      
      /// these are the atom names we can replace with a fixed structure
      addReplacedAtomName("SSD");
      addReplacedAtomName("SSD_E");
      addReplacedAtomName("SSD_RF");
      addReplacedAtomName("SSD1");
      addReplacedAtomName("TAP");
      addReplacedAtomName("TRED");
      
      // this is the reference structure we'll use for the replacement:
      addSite("H", Vector3d(0.0, -0.75695, 0.5206));
      addSite("H", Vector3d(0.0,  0.75695, 0.5206));
      addSite("O", Vector3d(0.0,  0.0,    -0.0654));
      addSite("X", Vector3d(0.0,  0.0,     0.0   ), Vector3d(0,0,1));
    }
  };
  
  class GBtailVisitor : public ReplacementVisitor{
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
  
  class GBheadVisitor : public ReplacementVisitor{
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
}//namespace OpenMD
#endif
