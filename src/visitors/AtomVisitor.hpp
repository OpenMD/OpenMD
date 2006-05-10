/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#ifndef VISITORS_BASEATOMVISITOR_HPP
#define VISITORS_BASEATOMVISITOR_HPP

#include <set>

#include "visitors/BaseVisitor.hpp"
#include "visitors/AtomData.hpp"

namespace oopse {

  /**
   * @class BaseAtomVisitor
   * @todo document
   */
  class BaseAtomVisitor : public BaseVisitor{
  public:
    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom) {}
    virtual void visit(RigidBody* rb);
    void setVisited(Atom* atom);
    bool isVisited(Atom* atom);

  protected:
    BaseAtomVisitor(SimInfo* info) : BaseVisitor() {}    
    SimInfo* info;
  };

  /**
   * @class SSDAtomVisitor
   * @todo document
   */
  class SSDAtomVisitor : public BaseAtomVisitor{
  public:
    SSDAtomVisitor(SimInfo* info) : BaseAtomVisitor(info) {
      visitorName = "SSDAtomVisitor";
      ssdAtomType.insert("SSD");
      ssdAtomType.insert("SSD_E");
      ssdAtomType.insert("SSD_RF");
      ssdAtomType.insert("SSD1");
      ssdAtomType.insert("TAP");
    }

    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom);       
    virtual void visit(RigidBody* rb) {}

    virtual const std::string toString();
  private:
    inline bool isSSDAtom(const std::string& atomType);
    std::set<std::string> ssdAtomType;   
  };

  class LinearAtomVisitor : public BaseAtomVisitor{
  public:
    LinearAtomVisitor(SimInfo* info) : BaseAtomVisitor(info) {
      visitorName = "LinearAtomVisitor";
      linearAtomType.insert("linear");
    }

    void addGayBerneAtomType(const std::string& atomType); 
    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom);       
    virtual void visit(RigidBody* rb) {}

    virtual const std::string toString();
  private:
    inline bool isLinearAtom(const std::string& atomType);
    std::set<std::string> linearAtomType;   
  };

  class GBLipidAtomVisitor : public BaseAtomVisitor{
  public:
    GBLipidAtomVisitor(SimInfo* info) : BaseAtomVisitor(info) {
      visitorName = "GBLipidAtomVisitor";
      GBLipidAtomType.insert("GBlipid");
    }

    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom);       
    virtual void visit(RigidBody* rb) {}

    virtual const std::string toString();
  private:
    inline bool isGBLipidAtom(const std::string& atomType);
    std::set<std::string> GBLipidAtomType;   
  };


  class DefaultAtomVisitor : public BaseAtomVisitor{
  public:
    DefaultAtomVisitor(SimInfo* info) : BaseAtomVisitor(info) { visitorName = "DefaultAtomVisitor";}

    virtual void visit(Atom* atom);    
    virtual void visit(DirectionalAtom* datom);    
    virtual void visit(RigidBody* rb) {}

    virtual const std::string toString();

  };

}//namespace oopse
#endif
