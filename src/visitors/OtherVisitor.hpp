#ifndef _OTHERVISITOR_H_
#define _OTHERVISITOR_H_
#include <set>
#include <string>
#include <vector>

#include "visitors/BaseVisitor.hpp"
#include "primitives/StuntDouble.hpp"
#include "visitors/AtomData.hpp"

using namespace std;

namespace oopse {

//IgnoreVisitor will turn on the ignoring flag of the stuntdouble
class IgnoreVisitor : public BaseVisitor{
  public:
    IgnoreVisitor() : BaseVisitor() {visitorName = "IgnoreVisitor";}

    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);
    
    virtual const string toString();

    void addIgnoreType(const string& type) {itList.insert(type);}
    
  protected:
    bool isIgnoreType(const string& name);
    void internalVisit(StuntDouble* sd);
    set<string> itList; //ignore type list;
};


class WrappingVisitor : public BaseVisitor{
  public:
    WrappingVisitor(SimInfo* info) : BaseVisitor() {
      this->info = info;
      visitorName = "WrappingVisitor";
    }
    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);

    virtual const string toString();

  protected:
    void internalVisit(StuntDouble* sd);
    SimInfo* info;
};


class IntVec3 {
  public:
    IntVec3(){}
    IntVec3(int i, int j, int k){
      vec[0] = i;
      vec[1] = j;
      vec[2] = k;
    }
    
    int vec[3];
    int& operator[](int index) {return vec[index];}  
};

class ReplicateVisitor : public BaseVisitor{
  public:
    ReplicateVisitor(SimInfo* info, IntVec3 opt);
    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);

    virtual const string toString();
  protected:
    void internalVisit(StuntDouble* sd);
    void replicate(vector<AtomInfo*>& infoList,  AtomData* data, double boxM[3][3]);
    
  private:
    vector<IntVec3> dir;
    SimInfo* info;
    IntVec3 replicateOpt;
};

class XYZVisitor : public BaseVisitor{
  public:
    XYZVisitor(SimInfo* info, bool printDipole = true);
    
    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);

    virtual const string toString();
    
    void writeFrame(ostream& outStream);    
    void clear() {frame.clear();}
    
  protected:
    void internalVisit(StuntDouble* sd);
    bool isIgnore(StuntDouble* sd);

  private:  
    SimInfo* info;
    vector<string> frame;
    bool printDipole;
};


class PrepareVisitor : public BaseVisitor{
  public:
    PrepareVisitor() : BaseVisitor() {visitorName = "prepareVisitor";}

    virtual void visit(Atom* atom) {internalVisit(atom);}
    virtual void visit(DirectionalAtom* datom) {internalVisit((Atom*)datom);}
    virtual void visit(RigidBody* rb) {internalVisit(rb);}

    virtual const string toString();

  protected:
    void internalVisit(Atom* atom);
    void internalVisit(RigidBody* rb);
};

class WaterTypeVisitor : public BaseVisitor{
  public:
    WaterTypeVisitor() ;
    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom) {}
    virtual void visit(RigidBody* rb);

    virtual const string toString();
    
  private:
    void replaceType(string& atomType);
      
    set<string> waterTypeList;
};


}//namespace oopse
#endif //_OTHERVISITOR_H_
