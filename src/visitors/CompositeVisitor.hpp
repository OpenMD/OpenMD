#ifndef _COMPOSITEVISITOR_H_
#define _COMPOSITEVISITOR_H_

#include <list>
#include "visitors/BaseVisitor.hpp"
#include "primitives/StuntDouble.hpp"
using namespace std;

typedef list<pair<BaseVisitor*, int> >::iterator VisitorIterator;

class CompositeVisitor: public BaseVisitor{
  public:
    CompositeVisitor() : BaseVisitor() { visitorName = "CompositeVisitor";}
    ~CompositeVisitor();
    
    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom); 
    virtual void visit(RigidBody* rb); 
     virtual void update();
    
    void addVisitor(BaseVisitor* v, int priority = 0);
    BaseVisitor* beginVisitor(VisitorIterator& i);
    BaseVisitor* nextVisitor(VisitorIterator& i);

    const string toString();
  protected:
    list<pair<BaseVisitor*, int> > visitorList;
};

#endif //_COMPOSITEVISITOR_H_
