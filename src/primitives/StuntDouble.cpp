#include "primitives/StuntDouble.hpp"
#include "primitives/Atom.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"
#include "utils/simError.h"

/*   
        "Don't move, or you're dead! Stand up! Captain, we've got them!"
     
        "Spectacular stunt, my friends, but all for naught. Turn around
         please. Ha. What a pity. What a pity. So, Princess, you thought
         you could outwit the imperious forces of...."

        "You idiots! These are not them. You've captured their stunt
         doubles! Search the area. Find them! Find them!"

      StuntDouble is a very strange idea.  A StuntDouble stands in for
      some object that can be manipulated by the Integrators or
      Minimizers.  Some of the manipulable objects are Atoms, some are
      DirectionalAtoms, and some are RigidBodies.  StuntDouble
      provides an interface for the Integrators and Minimizers to use,
      and does some preliminary sanity checking so that the program
      doesn't try to do something stupid like torque an Atom.  (The
      quotes above are from Spaceballs...)

*/

StuntDouble::~StuntDouble(){
  map<string, GenericData*>::iterator iter;

  for(iter = properties.begin(); iter != properties.end(); ++iter ){
    delete iter->second;
    properties.erase(iter);
  }
    
}

int StuntDouble::getObjType(){
  return objType;
}

double StuntDouble::getMass(){
  switch (objType)
    {
    case OT_ATOM :
    case OT_DATOM:      
      return ((Atom*)this)->getMass();
      break;
    case OT_RIGIDBODY:
      return ((RigidBody*)this)->getMass();
      break;      
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getMass.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
      return 0.0;
    }
}


void StuntDouble::getPos(double pos[3]){
  switch (objType)
    {
    case OT_ATOM :
    case OT_DATOM:      
      ((Atom*)this)->getPos(pos);
      break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getPos(pos);
      break;      
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getPos.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }  
}

void StuntDouble::getVel(double vel[3]){
  switch (objType)
    {
    case OT_ATOM :
    case OT_DATOM:      
      ((Atom*)this)->getVel(vel);
      break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getVel(vel);
      break;      
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getVel.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }  
}

void StuntDouble::getFrc(double frc[3]){
  switch (objType)
    {
    case OT_ATOM :
    case OT_DATOM:      
      ((Atom*)this)->getFrc(frc);
      break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getFrc(frc);
      break;      
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getFrc.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }  
}


void StuntDouble::setPos(double pos[3]){
  switch (objType)
    {
    case OT_ATOM :
    case OT_DATOM:      
      ((Atom*)this)->setPos(pos);
      break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->setPos(pos);
      break;      
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::setPos.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }  
}

void StuntDouble::setVel(double vel[3]){

  switch (objType)
    {
    case OT_ATOM :
    case OT_DATOM:      
      ((Atom*)this)->setVel(vel);
      break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->setVel(vel);
      break;      
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::setVel.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }  
}

void StuntDouble::addFrc(double frc[3]){
  switch (objType)
    {
    case OT_ATOM :
    case OT_DATOM:      
      ((Atom*)this)->addFrc(frc);
      break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->addFrc(frc);
      break;      
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::addFrc.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }  
}

void StuntDouble::getA(double A[3][3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::getA was called for a regular atom.\n"
               "\tRegular Atoms don't have rotation matrices.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      // Return unit matrix
      A[0][0] = 1;  A[0][1] = 0;  A[0][2] = 0;
      A[1][0] = 0;  A[1][1] = 1;  A[1][2] = 0;
      A[2][0] = 0;  A[2][1] = 0;  A[2][2] = 1;
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->getA(A);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getA(A);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getA.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
  }
}

void StuntDouble::setA(double A[3][3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::setA was called for a regular atom.\n"
               "\tRegular Atoms don't have rotation matrices.  Be smarter.\n");
      painCave.isFatal = 1;
      simError();
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->setA(A);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->setA(A);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::setA.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }
}

void StuntDouble::getJ(double j[3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::getJ was called for a regular atom.\n"
               "\tRegular Atoms don't have angular momentum.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      // Return zeros.
      j[0] = 0;
      j[1] = 0;
      j[2] = 0;
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->getJ(j);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getJ(j);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getJ.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
  }
}

void StuntDouble::setJ(double j[3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::setJ was called for a regular atom.\n"
               "\tRegular Atoms don't have angular momentum.  Be smarter.\n");
      painCave.isFatal = 1;
      simError();
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->setJ(j);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->setJ(j);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::setJ.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }
}

void StuntDouble::getQ(double q[4] ){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::getJ was called for a regular atom.\n"
               "\tRegular Atoms don't have angular momentum.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      // Return zeros.
      q[0] = 0;
      q[1] = 0;
      q[2] = 0;
      q[3] = 0;
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->getQ(q);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getQ(q);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getJ.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
  }
}

void StuntDouble::setQ(double q[4] ){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::setJ was called for a regular atom.\n"
               "\tRegular Atoms don't have angular momentum.  Be smarter.\n");
      painCave.isFatal = 1;
      simError();
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->setJ(q);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->setJ(q);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::setJ.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }
}
void StuntDouble::getTrq(double trq[3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::getTrq was called for a regular atom.\n"
               "\tRegular Atoms don't have torques.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      // Return zeros.
      trq[0] = 0;
      trq[1] = 0;
      trq[2] = 0;
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->getTrq(trq);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getTrq(trq);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getTrq.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
  }
}

void StuntDouble::addTrq(double trq[3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::addTrq was called for a regular atom.\n"
               "\tRegular Atoms don't have torques.  Be smarter.\n");
      painCave.isFatal = 1;
      simError();
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->addTrq(trq);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->addTrq(trq);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::addTrq.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }
}

void StuntDouble::getI(double I[3][3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::getI was called for a regular atom.\n"
               "\tRegular Atoms don't have moments of inertia.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      // Return zero matrix
      I[0][0] = 0;  I[0][1] = 0;  I[0][2] = 0;
      I[1][0] = 0;  I[1][1] = 0;  I[1][2] = 0;
      I[2][0] = 0;  I[2][1] = 0;  I[2][2] = 0;
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->getI(I);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getI(I);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getI.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
  }
}

void StuntDouble::lab2Body(double vec[3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::lab2Body was called for a regular atom.\n"
               "\tRegular Atoms don't have reference frames.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->lab2Body(vec);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->lab2Body(vec);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::lab2Body.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
  }
}

void StuntDouble::getGrad(double grad[6]){
  double frc[3];

  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::getGrad was called for a regular atom.\n"
               "\tRegular Atoms don't have 6 coordinates.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      ((Atom*)this)->getFrc(frc);
      grad[0] = -frc[0];
      grad[1] = -frc[1];
      grad[2] = -frc[2];
      grad[3] = 0.0;
      grad[4] = 0.0;
      grad[5] = 0.0;
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->getGrad(grad);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getGrad(grad);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getGrad.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
  }
}

void StuntDouble::setEuler(double phi, double theta, double psi){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::setEuler was called for a regular atom.\n"
               "\tRegular Atoms don't have Euler Angles.  Be smarter.\n");
      painCave.isFatal = 1;
      simError();
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->setEuler(phi, theta, psi);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->setEuler(phi, theta, psi);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::setA.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }
}

void StuntDouble::getEulerAngles(double eulers[3]){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::getEulerAngles was called for a regular atom.\n"
               "\tRegular Atoms don't have Euler angles.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      // Return zeros.
      eulers[0] = 0;
      eulers[1] = 0;
      eulers[2] = 0;
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->getEulerAngles(eulers);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->getEulerAngles(eulers);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getEulerAngles.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
  }
}

bool StuntDouble::isLinear() {
  int i;
  double momI[3][3];
  bool linearTest = false;
  double tolerance = 0.001;

  getI(momI);
  
  for (i=0; i<3; i++){
    if (momI[i][i]<tolerance){
      linearTest = true;
      zeroAxis = i;
    }
  }

  return linearTest;
}

double StuntDouble::getZangle(){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::getZangle was called for a regular atom.\n"
               "\tRegular Atoms don't have zAngles.  Be smarter.\n");
      painCave.isFatal = 0;
      simError();
      // Return zeros.
      return 0;
      break;
    case OT_DATOM:
      return ((DirectionalAtom*)this)->getZangle();
    break;
    case OT_RIGIDBODY:
      return ((RigidBody*)this)->getZangle();
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::getZangle.\n",
               objType );
      painCave.isFatal = 1;
      simError();    
      return 0;
  }
}

void StuntDouble::setZangle(double zAngle){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::setZangle was called for a regular atom.\n"
               "\tRegular Atoms don't have zAngles.  Be smarter.\n");
      painCave.isFatal = 1;
      simError();
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->setZangle(zAngle);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->setZangle(zAngle);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::setZangle.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }
}

void StuntDouble::addZangle(double zAngle){
  switch (objType) 
    {
    case OT_ATOM:
      sprintf( painCave.errMsg,
               "StuntDouble::addZangle was called for a regular atom.\n"
               "\tRegular Atoms don't have zAngles.  Be smarter.\n");
      painCave.isFatal = 1;
      simError();
      break;
    case OT_DATOM:
      ((DirectionalAtom*)this)->addZangle(zAngle);
    break;
    case OT_RIGIDBODY:
      ((RigidBody*)this)->addZangle(zAngle);
    break;
    default:
      sprintf( painCave.errMsg,
               "Unrecognized ObjType (%d) in StuntDouble::addZangle.\n",
               objType );
      painCave.isFatal = 1;
      simError();     
    }
}

void StuntDouble::addProperty(GenericData* data){
  map<string, GenericData*>::iterator result;
  result = properties.find(data->getID());
  
  //we can't simply use  properties[prop->getID()] = prop,
  //it will cause memory leak if we already contain a propery which has the same name of prop
  
  if(result != properties.end()){
    delete (*result).second;
    (*result).second = data;      
  }
  else
    properties[data->getID()] = data;

  
}
void StuntDouble::removeProperty(const string& propName){
  map<string, GenericData*>::iterator result;
    
  result = properties.find(propName);
  
  if(result != properties.end()){
    delete result->second;
    properties.erase(result);
    
  }
  
}
GenericData* StuntDouble::getProperty(const string& propName){
  map<string, GenericData*>::iterator result;
  
  
  result = properties.find(propName);
  
  if(result != properties.end()) 
    return (*result).second;  
  else   
    return NULL;    
}
