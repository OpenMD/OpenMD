#include "io/ZConsReader.hpp"
#include "utils/simError.h"

ZConsReader::ZConsReader(SimInfo* info)
                   :istream(NULL){

  GenericData* data;
  StringData* filename; 
  this->info = info;

 //retrieve output filename of z force
  data = info->getProperty(ZCONSFILENAME_ID);
  if(!data) {

      
    sprintf( painCave.errMsg,
               "ZConsReader error: If you use an ZConstraint\n"
               " , you must set output filename of z-force.\n");
    painCave.isFatal = 1;
    simError();  

  }
  else{

    filename = dynamic_cast<StringData*>(data);
    
    if(!filename){

      sprintf( painCave.errMsg,
                 "ZConsReader error: Can not get property from SimInfo\n");
      painCave.isFatal = 1;
      simError();  
        
    }
    else{
      zconsFileName = filename->getData();
    }
    
  }

  istream = new ifstream(zconsFileName.c_str());

  if (!istream){
    cerr << "open " << filename << "error" << endl;
    exit(1);
  }
  
  readHeader();  
}

ZConsReader::ZConsReader(const string& filename){
  istream = new ifstream(zconsFileName.c_str());

  if (!istream){
    cerr << "open " << filename << "error" << endl;
    exit(1);
  }
  
  readHeader();  
}

ZConsReader::~ZConsReader(){
  istream->close();
}

int ZConsReader::getNumZMol(){
  return nZMol;
}

vector<int> ZConsReader::getZConsIndex(){
  return index;
}

vector<double> ZConsReader::getZConsPos(){
  return zconsPos;
}

//double ZConsReader::getKRatio(){
//  return kRatio;
//}
    
double ZConsReader::getCurTime(){
  return curTime;
}

vector<double> ZConsReader::getCurZPos(){
  return curZPos;
}

vector<double> ZConsReader::getCurFZ(){
  return curFZ;
}

void ZConsReader::readHeader(){
  const int MAXBUFFERSIZE = 2000;
  char line[MAXBUFFERSIZE];
  int zmolIndex;
  float zmolPos;
  int sscanfCount;
  
  istream->getline(line, MAXBUFFERSIZE);

  cout << line << endl;
  //skip the comment lines
  while(line[0] == '#')
      istream->getline(line, MAXBUFFERSIZE);

  sscanfCount = sscanf(line, "%d", &nZMol);

  if (sscanfCount != 1){
    cerr << "ZConsReader Error : reading file error" << endl;
    exit(1);
  }
  
  for(int i = 0 ; i < nZMol; i++){

    istream->getline(line, MAXBUFFERSIZE);
      
    sscanfCount = sscanf(line, "%d\t%f", &zmolIndex, &zmolPos);
    if (sscanfCount != 2){
      cerr << "ZConsReader Error : reading file error" << endl;
      exit(1);    
    }
    
    index.push_back(zmolIndex);
    zconsPos.push_back(zmolPos);
  }

  curZPos.resize(nZMol);
  curFZ.resize(nZMol);
}

void ZConsReader::readNextFrame(){
  const int MAXBUFFERSIZE = 2000;
  char line[MAXBUFFERSIZE];  
  int tempNZMol;
  int sscanfCount;
  int tempIndex;
  float tempCurTime;
  float tempFZ;
  float tempCurZPos;
  float tempZconsPos;
  
  istream->getline(line, MAXBUFFERSIZE);
  sscanfCount = sscanf(line, "%f", &tempCurTime);
  if (sscanfCount != 1){
    cerr << "ZConsReader Error : reading file error" << endl;
    exit(1);
  }
  curTime = tempCurTime;
  
  istream->getline(line, MAXBUFFERSIZE);
  sscanfCount = sscanf(line, "%d", &tempNZMol);
  if (sscanfCount != 1){
    cerr << "ZConsReader Error : reading file error" << endl;
    exit(1);
  }
  
  if (tempNZMol != nZMol){
    cerr << "ZConsReader Error: reading file error" << endl;
    exit(1);
  }

  for(int i = 0; i < nZMol; i++){
    istream->getline(line, MAXBUFFERSIZE);
    sscanfCount = sscanf(line, "%d\t%f\t%f\t%f", &tempIndex, &tempFZ, &tempCurZPos,&tempZconsPos);
    if (sscanfCount != 4){
      cerr << "ZConsReader Error : reading file error" << endl;
      exit(1);
    }

    index[i] = tempIndex;
    curFZ[i] = tempFZ;
    curZPos[i]= tempCurZPos;
    zconsPos[i] = tempZconsPos;
  }

}

bool ZConsReader::hasNextFrame(){
  return istream->peek() != EOF ? true : false;
}

