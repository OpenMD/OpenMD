#ifndef _RIGIDBODYVISITOR_H_
#define _RIGIDBODYVISITOR_H_
#include <iostream>
#include <set>
#include <string>
#include "BaseVisitor.hpp"
#include "GenericData.hpp"

using namespace std;

class BaseRigidBodyVisitor: public BaseVisitor{
  public:

  protected:
     BaseRigidBodyVisitor(SimInfo* info) : BaseVisitor(){ this->info = info;}

      SimInfo* info;
};


//LipidHeadVisitor  adds a pesudo atom into rigidbody which holds a dipole moment
class LipidHeadVisitor : public BaseRigidBodyVisitor{
  public:
     LipidHeadVisitor(SimInfo* info) : BaseRigidBodyVisitor(info){ visitorName = "LipidHeadVisitor";}

    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom) {}
    virtual void visit(RigidBody* rb);
    
    void addLipidHeadName(const string& name);
    virtual const string toString();
    
  protected:
     bool canVisit(const string& name);

     set<string> lipidHeadName;
};

class RBCOMVisitor : public BaseRigidBodyVisitor{
  public:
    RBCOMVisitor(SimInfo* info) : BaseRigidBodyVisitor(info){ visitorName = "RBCOMVisitor";}

    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom) {}
    virtual void visit(RigidBody* rb);

    virtual const string toString();
};
#endif
