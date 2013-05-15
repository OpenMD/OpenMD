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
 
#include <cstring>
#include "visitors/CompositeVisitor.hpp"
#include "primitives/RigidBody.hpp"
#include "primitives/DirectionalAtom.hpp"

namespace OpenMD {

  CompositeVisitor::~CompositeVisitor(){
    VisitorIterator i;
    BaseVisitor* curVisitor;
  
    for(curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i))
      delete curVisitor;   

    visitorList.clear();
  
  }
  void CompositeVisitor::addVisitor(BaseVisitor* newVisitor, int priority){
    VisitorIterator i;
    int curPriority;
  
    for(i = visitorList.begin(); i != visitorList.end(); ++i){
      curPriority = (*i).second;
      //if new visitor has higher priority, just insert it before current visitor
      if(priority > curPriority){      
	visitorList.insert(i, std::make_pair(newVisitor, priority));
      }
    }

    //if new visitor has lowest priority, insert it at the end of the list
    visitorList.insert(visitorList.end(), std::make_pair(newVisitor, priority));
  }

  BaseVisitor* CompositeVisitor::beginVisitor(VisitorIterator& i){
    i = visitorList.begin();
    return i != visitorList.end() ? (*i).first : NULL;
  }

  BaseVisitor* CompositeVisitor::nextVisitor(VisitorIterator& i){
    ++i;
    return i != visitorList.end() ? (*i).first : NULL;

  }

  void CompositeVisitor::visit(Atom* atom){
    VisitorIterator i;
    BaseVisitor* curVisitor;

    for(curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i))
      atom->accept(curVisitor);
  }

  void CompositeVisitor::visit(DirectionalAtom* datom){
    VisitorIterator i;
    BaseVisitor* curVisitor;

    for(curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i))
      datom->accept(curVisitor);
  }
  void CompositeVisitor::visit(RigidBody* rb){
    VisitorIterator i;
    BaseVisitor* curVisitor;
    std::vector<Atom*> myAtoms;
    std::vector<Atom*>::iterator atomIter;

    myAtoms = rb->getAtoms();
  
    for(curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i)){
      rb->accept(curVisitor);
    
      for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
	(*atomIter)->accept(curVisitor);
    }

  
  
  }

  const  std::string CompositeVisitor::toString(){
    VisitorIterator i;
    std::string result;
    char buffer[65535];

    sprintf(buffer ,"******************************************************************\n");
    result += buffer;
  
    sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer , "Visitor Description: visitor manager  maintaining a priority  list\n");
    result += buffer;

    sprintf(buffer , "visitors in current priority list:\n");
    result += buffer;
  
    for(i = visitorList.begin(); i != visitorList.end(); ++i){
      sprintf(buffer, "Priority = %d\tvisitor = %s\n",  (*i).second, ((*i).first->getVisitorName()).c_str());
      result += buffer;
    }

    sprintf(buffer, "Detail information about every visitor:\n");


    for(i = visitorList.begin(); i != visitorList.end(); ++i)
      result += ((*i).first)->toString();

    sprintf(buffer ,"******************************************************************\n");
    result += buffer;
  
    return result;
  }

  void CompositeVisitor::update(){
    VisitorIterator i;
    BaseVisitor* curVisitor;  
  
    for(curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i))
      curVisitor->update();
  }

}//namespace OpenMD
