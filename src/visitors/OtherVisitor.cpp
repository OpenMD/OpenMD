#include "OtherVisitor.hpp"
#include "DirectionalAtom.hpp"
#include "RigidBody.hpp"
#include "Molecule.hpp"
#include "SimInfo.hpp"
//----------------------------------------------------------------------------//
void IgnoreVisitor::visit(Atom* atom){
  if(isIgnoreType(atom->getType()))
    internalVisit(atom);
}

void IgnoreVisitor::visit(DirectionalAtom* datom){
  if(isIgnoreType(datom->getType()))
    internalVisit(datom);
}

void IgnoreVisitor::visit(RigidBody* rb){
  vector<Atom*> myAtoms;
  vector<Atom*>::iterator atomIter;
  AtomInfo* atomInfo;
  
  if(isIgnoreType(rb->getType())){
    
    internalVisit(rb);

    myAtoms = rb->getAtoms();    
    
    for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
      internalVisit(*atomIter);

  }
  
}

bool IgnoreVisitor::isIgnoreType(const string& name){
  return itList.find(name) != itList.end() ? true : false;
}

void IgnoreVisitor::internalVisit(StuntDouble* sd){
  GenericData* data;
  data = sd->getProperty("IGNORE");

  //if this stuntdoulbe is already marked as ignore just skip it
  if (data == NULL){
    data = new GenericData;
    data->setID("IGNORE");
    sd->addProperty(data);
  }
    
}

const string IgnoreVisitor::toString(){
  char buffer[65535];
  string result;
  set<string>::iterator i;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer ,"Visitor Description: ignore  stuntdoubles\n");
  result += buffer;

  //print the ignore type list
  sprintf(buffer , "Ignore type list contains below types:\n");
  result += buffer;

  for(i = itList.begin(); i != itList.end(); ++i){
    sprintf(buffer ,"%s\t", i->c_str());
    result += buffer;

  }
  sprintf(buffer ,"\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

//----------------------------------------------------------------------------//

void WrappingVisitor::visit(Atom* atom){
  internalVisit(atom);
}
void WrappingVisitor::visit(DirectionalAtom* datom){
  internalVisit(datom);
}

void WrappingVisitor::visit(RigidBody* rb){
  internalVisit(rb);
}

void WrappingVisitor::internalVisit(StuntDouble* sd){
  GenericData* data;
  AtomData* atomData;
  AtomInfo* atomInfo;
  vector<AtomInfo*>::iterator i;

  data = sd->getProperty("ATOMDATA");
  if(data != NULL){
    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData == NULL)
      return;
  }
  else
    return;

  for(atomInfo = atomData->beginAtomInfo(i); atomInfo; atomInfo = atomData->nextAtomInfo(i))
    info->wrapVector(atomInfo->pos);

   
}

const string WrappingVisitor::toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer ,"Visitor Description: wrapping atoms back to periodic box\n");
  result += buffer;
 
  sprintf(buffer,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

//----------------------------------------------------------------------------//

ReplicateVisitor::ReplicateVisitor(SimInfo* info, IntVec3 opt) : BaseVisitor(){
  this->info = info;
  visitorName = "ReplicateVisitor";
  this->replicateOpt = opt;
  //generate the replicate directions
  for(int i = 0; i <= replicateOpt[0]; i ++)
    for(int j = 0; j <= replicateOpt[1]; j ++)
      for(int k = 0; k <= replicateOpt[2]; k ++)
        //skip original frame
        if(i == 0 && j ==0 && k ==0)
          continue;
        else
          dir.push_back(IntVec3(i, j, k));
  
}
void ReplicateVisitor::visit(Atom* atom){
  internalVisit(atom);
}
void ReplicateVisitor::visit(DirectionalAtom* datom){
  internalVisit(datom);
}

void ReplicateVisitor::visit(RigidBody* rb){
  internalVisit(rb);
}

void ReplicateVisitor::internalVisit(StuntDouble* sd){
  GenericData* data;
  AtomData* atomData;
  AtomInfo* atomInfo;
  double box[3][3];
  vector<AtomInfo*> atomInfoList;
  IntVec3 dir;
  
  //if there is not atom data, just skip it
  data = sd->getProperty("ATOMDATA");
  if(data != NULL){
    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData == NULL)
      return;
  }
  else
    return;

  
  info->getBoxM(box);

  atomInfoList = atomData->getData();
  
  replicate(atomInfoList, atomData, box);

}

void ReplicateVisitor::replicate(vector<AtomInfo*>& infoList, AtomData* data, double boxM[3][3]){
  AtomInfo * newAtomInfo;
  vector<IntVec3>::iterator dirIter;
  vector<AtomInfo*>::iterator i;
  
  for(dirIter = dir.begin(); dirIter != dir.end(); ++dirIter){
    for(i = infoList.begin(); i != infoList.end(); i++){
      newAtomInfo = new AtomInfo;
      *newAtomInfo = *(*i);    

      for(int j = 0; j < 3; j++)
        newAtomInfo->pos[j] +=   (*dirIter)[0] * boxM[j][0] + (*dirIter)[1]* boxM[j][1] + (*dirIter)[2] * boxM[j][2];

      data->addAtomInfo(newAtomInfo);
    }
  }// end for(dirIter)  
}

const string ReplicateVisitor::toString(){
  char buffer[65535];
  string result;
  set<string>::iterator i;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer ,"Visitor Description: replicate the atoms in different direction\n");
  result += buffer;
  
  //print the replicate direction
  sprintf(buffer , "repeatX = %d:\n", replicateOpt[0]);
  result += buffer;

  sprintf(buffer , "repeatY = %d:\n", replicateOpt[1]);
  result += buffer;

  sprintf(buffer , "repeatZ = %d:\n", replicateOpt[2]);
  result += buffer;


  sprintf(buffer,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

//----------------------------------------------------------------------------//

XYZVisitor::XYZVisitor(SimInfo* info,  bool printDipole) : BaseVisitor(){
  this->info = info;
  visitorName = "XYZVisitor";
  this->printDipole = printDipole;
}

void XYZVisitor::visit(Atom* atom){
  if(!isIgnore(atom))
    internalVisit(atom);
}

void XYZVisitor::visit(DirectionalAtom* datom){
  if(!isIgnore(datom))
    internalVisit(datom);
}

void XYZVisitor::visit(RigidBody* rb){
  if(!isIgnore(rb))
    internalVisit(rb);

}

void XYZVisitor::internalVisit(StuntDouble* sd){
  GenericData* data;
  AtomData* atomData;
  AtomInfo* atomInfo;
  vector<AtomInfo*>::iterator i;
  char buffer[1024];
  
  //if there is not atom data, just skip it
  data = sd->getProperty("ATOMDATA");
  if(data != NULL){
    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData == NULL)
      return;
  }
  else
    return;

  for(atomInfo = atomData->beginAtomInfo(i); atomInfo; atomInfo = atomData->nextAtomInfo(i)){

    if(printDipole)
      sprintf(buffer, "%s%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f",
                  atomInfo->AtomType.c_str(),
                  atomInfo->pos[0],
                  atomInfo->pos[1],
                  atomInfo->pos[2],
                  atomInfo->dipole[0],
                  atomInfo->dipole[1],
                  atomInfo->dipole[2]);
    else
      sprintf(buffer, "%s%15.8f%15.8f%15.8f",
                  atomInfo->AtomType.c_str(),
                  atomInfo->pos[0],
                  atomInfo->pos[1],
                  atomInfo->pos[2]);   

    frame.push_back(buffer);
              
  }

}

bool XYZVisitor::isIgnore(StuntDouble* sd){
  GenericData* data;
  
  data = sd->getProperty("IGNORE");
  return data ==NULL ? false : true;
}

void XYZVisitor::writeFrame(ostream& outStream){
  vector<string>::iterator i;
  double box[3][3];
  char buffer[1024];

  if(frame.size() == 0)
    cerr << "Current Frame does not contain any atoms" << endl;

  //total number of atoms  
  outStream << frame.size() << endl;

  //write comment line
  info->getBoxM(box);
  sprintf(buffer,"%15.8f;%15.8f%15.8f%15.8f;%15.8f%15.8f%15.8f;%15.8f%15.8f%15.8f",
              info->getTime(),
              box[0][0], box[0][1], box[0][2],
              box[1][0], box[1][1], box[1][2],
              box[2][0], box[2][1], box[2][2]);

  outStream << buffer << endl;  

  for(i = frame.begin(); i != frame.end(); ++i)
    outStream << *i << endl;
}

const string XYZVisitor::toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer ,"Visitor Description: assemble the atom data and output xyz file\n");
  result += buffer;
 
  sprintf(buffer,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

//----------------------------------------------------------------------------//

void PrepareVisitor::internalVisit(Atom * atom){
  GenericData* data;
  AtomData* atomData;

  //if visited property is  existed, remove it
  data = atom->getProperty("VISITED");
  if(data != NULL){
    atom->removeProperty("VISITED");  
  }

  //remove atomdata
  data = atom->getProperty("ATOMDATA");
  if(data != NULL){
    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData != NULL)
      atom->removeProperty("ATOMDATA");
  }
  
}

void PrepareVisitor::internalVisit(RigidBody * rb){
  GenericData* data;
  AtomData* atomData;
  vector<Atom*> myAtoms;
  vector<Atom*>::iterator atomIter;
  
  //if visited property is  existed, remove it
  data = rb->getProperty("VISITED");
  if(data != NULL){
    rb->removeProperty("VISITED");  
  }

  //remove atomdata
  data = rb->getProperty("ATOMDATA");
  if(data != NULL){
    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData != NULL)
      rb->removeProperty("ATOMDATA");
  }

  myAtoms = rb->getAtoms();
    
  for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
   internalVisit (*atomIter);  
}

const string PrepareVisitor::toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s", visitorName.c_str());
  result += buffer;

  sprintf(buffer ,"Visitor Description: prepare for operation of other vistors\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

//----------------------------------------------------------------------------//

WaterTypeVisitor:: WaterTypeVisitor(){
  visitorName = "WaterTypeVisitor";
  waterTypeList.insert("TIP3P_RB_0");
  waterTypeList.insert("TIP4P_RB_0");
  waterTypeList.insert("TIP5P_RB_0");
  waterTypeList.insert("SPCE_RB_0");  
}


void WaterTypeVisitor:: visit(RigidBody* rb){
  string rbName;
  vector<Atom*> myAtoms;
  vector<Atom*>::iterator atomIter;
  GenericData* data;
  AtomData* atomData;
  AtomInfo* atomInfo;
  vector<AtomInfo*>::iterator i;
  
  rbName = rb->getType();

  if(waterTypeList.find(rbName) != waterTypeList.end()){

    myAtoms = rb->getAtoms();    
    for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter){

      data = (*atomIter)->getProperty("ATOMDATA");
      if(data != NULL){
        atomData = dynamic_cast<AtomData*>(data);  
        if(atomData == NULL)
          continue;
      }
      else
        continue;
      
      for(atomInfo = atomData->beginAtomInfo(i); atomInfo; atomInfo = atomData->nextAtomInfo(i)){
        replaceType(atomInfo->AtomType);
      }//end for(atomInfo)

    }//end for(atomIter)
      
  }//end if (waterTypeList.find(rbName) != waterTypeList.end())
  
}

void WaterTypeVisitor:: replaceType(string& atomType){
  atomType = atomType.substr(0, atomType.find('_'));
}

const string WaterTypeVisitor:: toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer ,"Visitor Description: Replace the atom type in water model\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

