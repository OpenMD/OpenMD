#include "RigidBodyVisitor.hpp"
#include "RigidBody.hpp"
#include "MatVec3.h"

void LipidHeadVisitor::visit(RigidBody* rb){
  double pos[3];
  double u[3] = {0, 0, 1};
  double newVec[3];
  GenericData* data;
  AtomData* atomData;
  AtomInfo* atomInfo;
  bool haveAtomData;
  double rotMatrix[3][3];

  if(!canVisit(rb->getType()))
    return;

  rb->getPos(pos);
  rb->getA(rotMatrix);
  matVecMul3(rotMatrix, u, newVec);

  data = rb->getProperty("ATOMDATA");
  if(data != NULL){

    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData == NULL){
      cerr << "can not get Atom Data from " << rb->getType() << endl;
      atomData = new AtomData; 
      haveAtomData = false;      
    }
    else
      haveAtomData = true;
  }
  else{
    atomData = new AtomData;
    haveAtomData = false;
  }

  atomInfo = new AtomInfo;
  atomInfo->AtomType = "X";
  atomInfo->pos[0] = pos[0];
  atomInfo->pos[1] = pos[1];
  atomInfo->pos[2] = pos[2];
  atomInfo->dipole[0] = newVec[0];
  atomInfo->dipole[1] = newVec[1];
  atomInfo->dipole[2] = newVec[2];

  atomData->addAtomInfo(atomInfo);

  if(!haveAtomData){
    atomData->setID("ATOMDATA");
    rb->addProperty(atomData);
  }
    
}

void LipidHeadVisitor::addLipidHeadName(const string& name){
  lipidHeadName.insert(name);

}

bool LipidHeadVisitor::canVisit(const string& name){
  return lipidHeadName.find(name) != lipidHeadName.end() ? true : false;

}

const string LipidHeadVisitor::toString(){
  char buffer[65535];
  string result;
  set<string>::iterator i;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;
  
  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  //print the ignore type list
  sprintf(buffer , "lipidHeadName list contains below types:\n");
  result += buffer;

  for(i = lipidHeadName.begin(); i != lipidHeadName.end(); ++i){
    sprintf(buffer ,"%s\t", i->c_str());
    result += buffer;
  }

  sprintf(buffer ,"\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;

}

void RBCOMVisitor::visit(RigidBody* rb){
  AtomData* atomData;
  AtomInfo* atomInfo;
  double pos[3];
  
  rb->getPos(pos);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "X";
  atomInfo->pos[0] = pos[0];
  atomInfo->pos[1] = pos[1];
  atomInfo->pos[2] = pos[2];
  atomInfo->dipole[0] = 0;
  atomInfo->dipole[1] = 0;
  atomInfo->dipole[2] = 0;

  atomData = new AtomData; 
  atomData->setID("ATOMDATA");
  atomData->addAtomInfo(atomInfo);

  rb->addProperty(atomData);
}

const string RBCOMVisitor::toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;
  
  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  //print the ignore type list
  sprintf(buffer , "Visitor Description: add a pseudo atom at the center of the mass of the rigidbody\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;

}

