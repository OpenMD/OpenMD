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
 
#ifndef VISITORS_ZCONSVISITOR_HPP
#define VISITORS_ZCONSVISITOR_HPP

#include "io/ZConsReader.hpp"
#include "primitives/RigidBody.hpp"
#include "visitors/BaseVisitor.hpp"
#include "visitors/AtomData.hpp"



namespace oopse {

  /**
   * @class ZConsVisitor
   * @note 
   */
  class ZConsVisitor : public BaseVisitor{
  public:
    enum ZConsState{zsFixed = 0, zsMoving};

    ZConsVisitor(SimInfo* info);
    ~ZConsVisitor();

    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);

    virtual void update();

    bool haveZconsMol() {return !zmolStates_.empty(); }

    virtual const std::string toString();
  protected:
    void internalVisit(StuntDouble* sd, const std::string& prefix);
    bool isZconstraint(int index, std::string& prefix);
    void readZconsFile(double time);

  private:  
    std::vector<double> zconsPos;
    std::map<int, ZConsState> zmolStates_;
    double zconsTol_;
    double zconsTime_;
    std::string zconsFilename_;
    ZConsReader* zconsReader_;
    SimInfo* info_;
    Snapshot* currSnapshot_;
    std::map<int, int> zatomToZmol_;
  };

}//namespace oopse
#endif // _ZCONS_VISITOR_H_


