#ifndef _BASEATOMVISITOR_H_
#define _BASEATOMVISITOR_H_


#include <vector>
#include "visitors/BaseVisitor.hpp"
#include "visitors/AtomData.hpp"

using namespace std;

namespace oopse {
    
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


class SSDAtomVisitor : public BaseAtomVisitor{
  public:
    SSDAtomVisitor(SimInfo* info) : BaseAtomVisitor(info) {
      visitorName = "SSDAtomVisitor";
      ssdAtomType.push_back("SSD");
      ssdAtomType.push_back("SSD_E");
      ssdAtomType.push_back("SSD_RF");
      ssdAtomType.push_back("SSD1");
    }

    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom);       
    virtual void visit(RigidBody* rb) {}
    
    virtual const string toString();
  private:
    inline bool isSSDAtom(const string& atomType);
    vector<string> ssdAtomType;   
};

class LinearAtomVisitor : public BaseAtomVisitor{
  public:
    LinearAtomVisitor(SimInfo* info) : BaseAtomVisitor(info) {
      visitorName = "LinearAtomVisitor";
      linearAtomType.push_back("linear");
    }

    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom);       
    virtual void visit(RigidBody* rb) {}
    
    virtual const string toString();
  private:
    inline bool isLinearAtom(const string& atomType);
    vector<string> linearAtomType;   
};




class DefaultAtomVisitor : public BaseAtomVisitor{
  public:
    DefaultAtomVisitor(SimInfo* info) : BaseAtomVisitor(info) { visitorName = "DefaultAtomVisitor";}

    virtual void visit(Atom* atom);    
    virtual void visit(DirectionalAtom* datom);    
    virtual void visit(RigidBody* rb) {}
    
    virtual const string toString();
    
};

}//namespace oopse
#endif
