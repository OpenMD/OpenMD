#ifndef __STUNTDOUBLE_HPP__
#define __STUNTDOUBLE_HPP__

#include <map>
#include "GenericData.hpp"

#define OT_ATOM 0
#define OT_DATOM 1
#define OT_RIGIDBODY 2

using namespace std;
class BaseVisitor;

class StuntDouble {
 public:
  virtual ~StuntDouble();
  
  int getObjType();
  
  bool isAtom(){
    return objType == OT_ATOM || objType == OT_DATOM;
  }

  bool isDirectionalAtom(){
    return objType == OT_DATOM;
  }

  bool isRigidBody(){
    return objType == OT_RIGIDBODY;
  }

  bool isDirectional(){
    return isDirectionalAtom() || isRigidBody();
  }

  virtual double getMass(void);
  
  virtual void   getPos(double pos[3]);
  virtual void   setPos(double pos[3]);
  
  virtual void   getVel(double vel[3]);
  virtual void   setVel(double vel[3]);
  
  virtual void   getFrc(double frc[3]);
  virtual void   addFrc(double frc[3]);

  virtual void   getA(double A[3][3]);
  virtual void   setA(double A[3][3]);

  virtual void   getJ(double j[3]);
  virtual void   setJ(double j[3]);

  virtual void getQ( double q[4] ); // get the quanternions
  virtual void setQ( double q[4] );

  virtual void setType(char* type) = 0;
  virtual char* getType() = 0;
  

  virtual void   getTrq(double trq[3]);
  virtual void   addTrq(double trq[3]);

  virtual void   getI(double I[3][3]);
  virtual void   lab2Body(double vec[3]);

  virtual void   getGrad(double grad[6]);
  virtual void   setEuler(double phi, double theta, double psi);
  virtual void   getEulerAngles(double eulers[3]);

  virtual bool isLinear() {return false;}
  virtual int linearAxis() {return -1;}

  virtual double   getZangle();
  virtual void   setZangle(double zAngle);
  virtual void   addZangle(double zAngle);

  virtual void accept(BaseVisitor* v) = 0;

  void addProperty(GenericData* data);
  void removeProperty(const string& propName);
  GenericData* getProperty(const string& propName);
  
 protected:
  StuntDouble(){}

  //prevent default copy constructor copy information from properties which will cause problem
  StuntDouble(const StuntDouble& sd){
    objType = sd.objType;
  }
  
  int objType;

  map<string, GenericData*> properties;
};

#endif
