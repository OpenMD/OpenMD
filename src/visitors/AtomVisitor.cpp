#include <cstring>
#include "visitors/AtomVisitor.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "math/MatVec3.h"
#include "primitives/RigidBody.hpp"

namespace oopse {

void BaseAtomVisitor::visit(RigidBody* rb){
  //vector<Atom*> myAtoms;
  //vector<Atom*>::iterator atomIter;

  //myAtoms = rb->getAtoms();
  
  //for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
  //  (*atomIter)->accept(this);
}

void BaseAtomVisitor::setVisited(Atom* atom){
  GenericData* data;
  data = atom->getProperty("VISITED");

  //if visited property is not existed, add it as new property
  if(data == NULL){
    data = new GenericData();
    data->setID("VISITED");
    atom->addProperty(data);  
  }
}

bool BaseAtomVisitor::isVisited(Atom* atom){
  GenericData* data;
  data = atom->getProperty("VISITED");
  return data == NULL ?  false : true;
}

bool SSDAtomVisitor::isSSDAtom(const string& atomType){
  vector<string>::iterator strIter;
  
  for(strIter = ssdAtomType.begin(); strIter != ssdAtomType.end(); ++strIter)
   if(*strIter == atomType)
    return true;
  
  return false;  
}

void SSDAtomVisitor::visit(DirectionalAtom* datom){

  vector<AtomInfo*> atoms;

  //we need to convert SSD into 4 differnet atoms
  //one oxygen atom, two hydrogen atoms and one pseudo atom which is the center of the mass
  //of the water with a dipole moment
  double h1[3] = {0.0, -0.75695, 0.5206};
  double h2[3] = {0.0, 0.75695, 0.5206};
  double ox[3] = {0.0, 0.0, -0.0654};
  double u[3] = {0, 0, 1};
  double rotMatrix[3][3];
  double rotTrans[3][3];
  AtomInfo* atomInfo;
  double pos[3];
  double newVec[3];
  double q[4];
  AtomData* atomData;
  GenericData* data;
  bool haveAtomData;
  
  //if atom is not SSD atom, just skip it
  if(!isSSDAtom(datom->getType()))
    return;

  data = datom->getProperty("ATOMDATA");
  if(data != NULL){

    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData == NULL){
      cerr << "can not get Atom Data from " << datom->getType() << endl;
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
   
  
  datom->getPos(pos);
  datom->getQ(q);
  datom->getA(rotMatrix);

  // We need A^T to convert from body-fixed to space-fixed:
  transposeMat3(rotMatrix, rotTrans);
  
  //center of mass of the water molecule
  matVecMul3(rotTrans, u, newVec);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "X";
  atomInfo->pos[0] = pos[0];
  atomInfo->pos[1] = pos[1];
  atomInfo->pos[2] = pos[2];
  atomInfo->dipole[0] = newVec[0];
  atomInfo->dipole[1] = newVec[1];
  atomInfo->dipole[2] = newVec[2];

  atomData->addAtomInfo(atomInfo);

  //oxygen
  matVecMul3(rotTrans, ox, newVec);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "O";
  atomInfo->pos[0] = pos[0] + newVec[0];
  atomInfo->pos[1] = pos[1] + newVec[1];
  atomInfo->pos[2] = pos[2] + newVec[2];
  atomInfo->dipole[0] = 0.0;
  atomInfo->dipole[1] = 0.0;
  atomInfo->dipole[2] = 0.0;
  atomData->addAtomInfo(atomInfo);


  //hydrogen1
    matVecMul3(rotTrans, h1, newVec);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "H";
  atomInfo->pos[0] = pos[0] + newVec[0];
  atomInfo->pos[1] = pos[1] + newVec[1];
  atomInfo->pos[2] = pos[2] + newVec[2];
  atomInfo->dipole[0] = 0.0;
  atomInfo->dipole[1] = 0.0;
  atomInfo->dipole[2] = 0.0;
  atomData->addAtomInfo(atomInfo);

  //hydrogen2
  matVecMul3(rotTrans, h2, newVec);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "H";
  atomInfo->pos[0] = pos[0] + newVec[0];
  atomInfo->pos[1] = pos[1] + newVec[1];
  atomInfo->pos[2] = pos[2] + newVec[2];
  atomInfo->dipole[0] = 0.0;
  atomInfo->dipole[1] = 0.0;
  atomInfo->dipole[2] = 0.0;
  atomData->addAtomInfo(atomInfo);

  //add atom data into atom's property

  if(!haveAtomData){
    atomData->setID("ATOMDATA");
    datom->addProperty(atomData);
  }

  setVisited(datom);

}

const string SSDAtomVisitor::toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer , "Visitor Description: Convert SSD into 4 different atoms\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

bool LinearAtomVisitor::isLinearAtom(const string& atomType){
  vector<string>::iterator strIter;
  
  for(strIter = linearAtomType.begin(); strIter != linearAtomType.end(); 
      ++strIter)
    if(*strIter == atomType)
      return true;
  
  return false;  
}

void LinearAtomVisitor::visit(DirectionalAtom* datom){

  vector<AtomInfo*> atoms;

  //we need to convert linear into 4 different atoms
  double c1[3] = {0.0, 0.0, -1.8};
  double c2[3] = {0.0, 0.0, -0.6};
  double c3[3] = {0.0, 0.0,  0.6};
  double c4[3] = {0.0, 0.0,  1.8};
  double rotMatrix[3][3];
  double rotTrans[3][3];
  AtomInfo* atomInfo;
  double pos[3];
  double newVec[3];
  double q[4];
  AtomData* atomData;
  GenericData* data;
  bool haveAtomData;
  
  //if atom is not SSD atom, just skip it
  if(!isLinearAtom(datom->getType()))
    return;
  
  data = datom->getProperty("ATOMDATA");
  if(data != NULL){

    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData == NULL){
      cerr << "can not get Atom Data from " << datom->getType() << endl;
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
   
  
  datom->getPos(pos);
  datom->getQ(q);
  datom->getA(rotMatrix);

  // We need A^T to convert from body-fixed to space-fixed:
  transposeMat3(rotMatrix, rotTrans);
  
  matVecMul3(rotTrans, c1, newVec);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "C";
  atomInfo->pos[0] = pos[0] + newVec[0];
  atomInfo->pos[1] = pos[1] + newVec[1];
  atomInfo->pos[2] = pos[2] + newVec[2];
  atomInfo->dipole[0] = 0.0;
  atomInfo->dipole[1] = 0.0;
  atomInfo->dipole[2] = 0.0;
  atomData->addAtomInfo(atomInfo);

  matVecMul3(rotTrans, c2, newVec);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "C";
  atomInfo->pos[0] = pos[0] + newVec[0];
  atomInfo->pos[1] = pos[1] + newVec[1];
  atomInfo->pos[2] = pos[2] + newVec[2];
  atomInfo->dipole[0] = 0.0;
  atomInfo->dipole[1] = 0.0;
  atomInfo->dipole[2] = 0.0;
  atomData->addAtomInfo(atomInfo);

  matVecMul3(rotTrans, c3, newVec);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "C";
  atomInfo->pos[0] = pos[0] + newVec[0];
  atomInfo->pos[1] = pos[1] + newVec[1];
  atomInfo->pos[2] = pos[2] + newVec[2];
  atomInfo->dipole[0] = 0.0;
  atomInfo->dipole[1] = 0.0;
  atomInfo->dipole[2] = 0.0;
  atomData->addAtomInfo(atomInfo);

  matVecMul3(rotTrans, c4, newVec);
  atomInfo = new AtomInfo;
  atomInfo->AtomType = "C";
  atomInfo->pos[0] = pos[0] + newVec[0];
  atomInfo->pos[1] = pos[1] + newVec[1];
  atomInfo->pos[2] = pos[2] + newVec[2];
  atomInfo->dipole[0] = 0.0;
  atomInfo->dipole[1] = 0.0;
  atomInfo->dipole[2] = 0.0;
  atomData->addAtomInfo(atomInfo);

  //add atom data into atom's property

  if(!haveAtomData){
    atomData->setID("ATOMDATA");
    datom->addProperty(atomData);
  }

  setVisited(datom);

}

const string LinearAtomVisitor::toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer , "Visitor Description: Convert linear into 4 different atoms\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

//----------------------------------------------------------------------------//

void DefaultAtomVisitor::visit(Atom* atom){
  AtomData* atomData;
  AtomInfo* atomInfo;
  double pos[3];

  if(isVisited(atom))
    return;

 atomInfo =new AtomInfo;

  atomData = new AtomData; 
  atomData->setID("ATOMDATA");
 
  atom->getPos(pos);
  atomInfo->AtomType = atom->getType();
  atomInfo->pos[0] = pos[0];
  atomInfo->pos[1] = pos[1];
  atomInfo->pos[2] = pos[2];
  atomInfo->dipole[0] = 0.0;
  atomInfo->dipole[1] = 0.0;
  atomInfo->dipole[2] = 0.0;


  atomData->addAtomInfo(atomInfo);
  
  atom->addProperty(atomData);

  setVisited(atom);
}
void DefaultAtomVisitor::visit(DirectionalAtom* datom){
  AtomData* atomData;
  AtomInfo* atomInfo;
  double pos[3];
  double u[3];

  if(isVisited(datom))
    return;
  
  datom->getPos(pos);
  datom->getU(u);

  atomData = new AtomData; 
  atomData->setID("ATOMDATA");
  atomInfo =new AtomInfo;
  
  atomInfo->AtomType = datom->getType();
  atomInfo->pos[0] = pos[0];
  atomInfo->pos[1] = pos[1];
  atomInfo->pos[2] = pos[2];
  atomInfo->dipole[0] = u[0];
  atomInfo->dipole[1] = u[1];
  atomInfo->dipole[2] = u[2];  

  atomData->addAtomInfo(atomInfo);

  datom->addProperty(atomData);

  setVisited(datom);
}


const string DefaultAtomVisitor::toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer , "Visitor Description: copy atom infomation into atom data\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}    

}//namespace oopse
