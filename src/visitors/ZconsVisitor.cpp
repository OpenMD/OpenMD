#include "ZconsVisitor.hpp"
#include <cmath>

ZConsVisitor::ZConsVisitor(SimInfo* info) : BaseVisitor(), zconsReader(NULL){
  GenericData* data;
  DoubleData* tolerance;  
  ZConsParaData* zConsParaData;
  StringData* filename; 
  DoubleData* sampleTime;
  vector<ZConsParaItem>* parameters;

  this->info = info;
  visitorName = "ZConsVisitor";
  
  //retrieve tolerance for z-constraint molecuels
  data = info->getProperty(ZCONSTOL_ID);

  if (!data){
    cerr << "Can not get zconstraint tolerance from SimInfo" << endl;
    haveZcons = false;  
    return;
  }
  else{
    tolerance = dynamic_cast<DoubleData*>(data);

    if (!tolerance){
      cerr << "Can not get zconstraint tolerance  from SimInfo" << endl;
      haveZcons = false;    
      return;
    }
    else{
      zconsTol = tolerance->getData();
    }
  }

  //retrieve sample time of z-contraint 
  data = info->getProperty(ZCONSTIME_ID);

  if (!data){
    cerr << "Can not get zcons time  from SimInfo" << endl;
    haveZcons = false;    
    return;
  }
  else{
    sampleTime = dynamic_cast<DoubleData*>(data);

    if (!sampleTime){
      cerr << "Can not get zcons time  from SimInfo" << endl;
      haveZcons = false;    
      return;
    }
    else{
     zconsTime = sampleTime->getData();
    }
  }

  //retrieve index of z-constraint molecules
  data = info->getProperty(ZCONSPARADATA_ID);
  if (!data){
    cerr << "Can not get index of z-constraint molecules from SimInfo" << endl;
    haveZcons = false;
    return;
  }
  else{
    zConsParaData = dynamic_cast<ZConsParaData*>(data);

    if (!zConsParaData){
      cerr << "Can not get index of z-constraint molecules from SimInfo" << endl;
      haveZcons = false;
      return;
    }
    else{
      vector<ZConsParaItem>::iterator i;
      Molecule* mol;

      parameters = zConsParaData->getData();
      for(i = parameters->begin(); i != parameters->end(); ++i){

        mol = findZconsMol(i->zconsIndex);
        if(mol != NULL)
          zconsMol.push_back(mol);

      }

      if(zconsMol.size() < 1){
        cerr << "number of zconstraint molecules is less than one" << endl;
        haveZcons = false;
        return;
      }
      
    }//end if (!zConsParaData)

  }//end  if (!data  
  
  //retrieve output filename of z force
  data = info->getProperty(ZCONSFILENAME_ID);
  if (!data){
    cerr << "Can not get filename of z-constraint output from SimInfo" << endl;
    haveZcons = false;
    return;
  }
  else{
    filename = dynamic_cast<StringData*>(data);

    if (!filename){
      cerr << "Can not get filename of z-constraint output from SimInfo" << endl;
      haveZcons = false;
      return;
    }
    else{
      zconsFilename = filename->getData();
    }
  }

  zconsReader = new ZConsReader(info);

  if (zconsReader->hasNextFrame())
  zconsReader->readNextFrame();
  
  haveZcons = true;
}

ZConsVisitor::~ZConsVisitor(){
  if(!zconsReader) 
    delete zconsReader;
  
}

void ZConsVisitor::visit(Atom* atom){
  string prefix;
  if(isZconstraint(atom->getIndex(), prefix))
    internalVisit(atom, prefix);
}

void ZConsVisitor::visit(DirectionalAtom* datom){
  string prefix;
  
  if(isZconstraint(datom->getIndex(), prefix))
    internalVisit(datom, prefix);
}

void ZConsVisitor::visit(RigidBody* rb){
  string prefix;
  vector<Atom*> atoms;
  
  atoms = rb->getAtoms();
    
  if(isZconstraint(atoms[0]->getIndex(), prefix))
    internalVisit(rb, prefix);
}

Molecule* ZConsVisitor::findZconsMol(int index){
  Molecule* mols;
  mols = info->molecules;
  for(int i =0; i < info->n_mol; i++)
    if(index == mols[i].getMyIndex())
      return &mols[i];

  return NULL;
}

void ZConsVisitor::update(){
  Molecule* mol;
  double com[3];
  ZConsState state;
  Atom** atoms;
  int natoms;

  getZconsPos(info->currentTime);
  
  for(size_t i = 0; i < zconsMol.size(); i++){
    zconsMol[i]->getCOM(com);

    if(fabs(com[2] - zconsPos[i]) < zconsTol)
      state = zsFixed;
    else
      state = zsMoving;

    //update the state of zconstraint atom
    natoms = zconsMol[i]->getNAtoms();
    atoms = zconsMol[i]->getMyAtoms();    
    for(int j =0; j < natoms; j++)
      zconsState[atoms[j]->getIndex()] = state;
  }
  
}

void ZConsVisitor::getZconsPos(double time){

  vector<double> tempPos;
  vector<double> prevPos;
  double tempTime;

  while(true){
    tempTime = zconsReader->getCurTime();

    zconsPos = zconsReader->getZConsPos();
    
    if(tempTime >= time)
      return;
    
    zconsReader->readNextFrame();

  } 

}

void ZConsVisitor::internalVisit(StuntDouble* sd, const string& prefix){
  GenericData* data;
  AtomData* atomData;
  AtomInfo* atomInfo;
  vector<AtomInfo*>::iterator iter;

  
  //if there is not atom data, just skip it
  data = sd->getProperty("ATOMDATA");
  if(data != NULL){
    atomData = dynamic_cast<AtomData*>(data);  
    if(atomData == NULL)
      return;
  }
  else
    return;

  for(atomInfo  = atomData->beginAtomInfo(iter); atomInfo; atomInfo = atomData->nextAtomInfo(iter))
    (atomInfo->AtomType).insert(0, prefix);
}


bool ZConsVisitor::isZconstraint(int index, string& prefix){
  map<int, ZConsState>::iterator i;
  string prefixString[] = {"ZF", "ZM"};
  
  i = zconsState.find(index);
  if(i !=zconsState.end()){
    prefix = prefixString[(*i).second];
    return true;
  }
  else{
    prefix = "";
    return false;
  }
}

const string ZConsVisitor::toString(){
  char buffer[65535];
  string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;
  
  sprintf(buffer , "number of zconstraint molecule: %d\n", zconsMol.size());
  result += buffer;

  sprintf(buffer , "zconstraint tolerance = %lf\n", zconsTol);
  result += buffer;

  sprintf(buffer , "zconstraint sample time = %lf\n", zconsTime);
  result += buffer;

  sprintf(buffer , "zconstraint output filename = %s\n", zconsFilename.c_str());
  result += buffer;
  
  for(size_t i = 0; i < zconsMol.size(); ++i){
    sprintf(buffer ,"zconstraint molecule[%d] = %d\n", i, zconsMol[i]->getMyIndex());
    result += buffer;

  }

  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}
