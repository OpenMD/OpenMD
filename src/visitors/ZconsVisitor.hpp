#ifndef _ZCONS_VISITOR_H_
#define _ZCONS_VISITOR_H_

#include "visitors/BaseVisitor.hpp"
#include "io/ZConsReader.hpp"
#include "visitors/AtomData.hpp"
#include "constraints/ZconsData.hpp"
enum ZConsState{zsFixed, zsMoving};
namespace oopse {
    
class ZConsVisitor : public BaseVisitor{
  public:

    
    ZConsVisitor(SimInfo* info);
    ~ZConsVisitor();
    
    virtual void visit(Atom* atom);
    virtual void visit(DirectionalAtom* datom);
    virtual void visit(RigidBody* rb);
    
    virtual void update();
    
    bool haveZconsMol() {return haveZcons;}

    virtual const string toString();
  protected:
    void internalVisit(StuntDouble* sd, const string& prefix);
    bool isZconstraint(int index, string& prefix);
    Molecule* findZconsMol(int index);
    void getZconsPos(double time);
    
  private:  
    vector<Molecule*> zconsMol;
    vector<double> zconsPos;
    map<int, ZConsState> zconsState;
    bool haveZcons;
    double zconsTol;
    double zconsTime;
    string zconsFilename;
    ZConsReader* zconsReader;
    SimInfo* info;
};

}//namespace oopse
#endif // _ZCONS_VISITOR_H_


