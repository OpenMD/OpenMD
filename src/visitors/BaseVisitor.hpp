#ifndef _BASEVISITOR_H_
#define _BASEVISITOR_H_
#include <iostream>
#include <string>
using namespace std;

class SimInfo;
class Atom;
class DirectionalAtom;
class RigidBody;

class BaseVisitor{
  public:
    virtual ~BaseVisitor() {}
    virtual void visit(Atom* atom) {}
    virtual void visit(DirectionalAtom* datom) {}
    virtual void visit(RigidBody* rb) {}

    virtual void update() {}

    const string& getVisitorName() {return visitorName;}
    virtual const string toString() {
      string result;
      char buffer[65535];

      sprintf(buffer,"------------------------------------------------------------------\n");
      result += buffer;      
      
      sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
      result += buffer;      

      sprintf(buffer, "------------------------------------------------------------------\n");
      result += buffer;      

      return result;
    }

  protected:
    BaseVisitor() {}
    string visitorName;
};

#endif
