#include <cstring>
#include "CompositeVisitor.hpp"
#include "RigidBody.hpp"
#include "DirectionalAtom.hpp"
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
  
  for(i = visitorList.begin(); i != visitorList.end(); i++){
    curPriority = (*i).second;
    //if new visitor has higher priority, just insert it before current visitor
    if(priority > curPriority){      
      visitorList.insert(i, make_pair(newVisitor, priority));
    }
  }

  //if new visitor has lowest priority, insert it at the end of the list
  visitorList.insert(visitorList.end(), make_pair(newVisitor, priority));
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
  vector<Atom*> myAtoms;
  vector<Atom*>::iterator atomIter;

  myAtoms = rb->getAtoms();
  
  for(curVisitor = beginVisitor(i); curVisitor; curVisitor = nextVisitor(i)){
    rb->accept(curVisitor);
    
    for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
      (*atomIter)->accept(curVisitor);
  }

  
  
}

const string CompositeVisitor::toString(){
  VisitorIterator i;
  string result;
  char buffer[65535];

  sprintf(buffer ,"******************************************************************\n");
  result += buffer;
  
  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer , "Visitor Description: visitor manager  maintaining a priority  list\n");
  result += buffer;

  sprintf(buffer , "visitors in current priority list:\n");
  result += buffer;
  
  for(i = visitorList.begin(); i != visitorList.end(); i++){
    sprintf(buffer, "Priority = %d\tvisitor = %s\n",  (*i).second, ((*i).first->getVisitorName()).c_str());
    result += buffer;
  }

  sprintf(buffer, "Detail information about every visitor:\n");


  for(i = visitorList.begin(); i != visitorList.end(); i++)
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
