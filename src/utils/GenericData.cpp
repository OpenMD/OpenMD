#include "utils/GenericData.hpp"
////////////////////////////////////////////////////////////////////////////////
//Implementation of GenericData
////////////////////////////////////////////////////////////////////////////////
GenericData::GenericData(){

  id = "undefined";

}

GenericData& GenericData::operator = (const GenericData& rhs){
  
  if(this == &rhs)
    return (*this);
  
  id = rhs.id;
  
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
//Implementation of ZConsParaData
////////////////////////////////////////////////////////////////////////////////
ZConsParaData::ZConsParaData(){
  id = ZCONSPARADATA_ID;
}

void ZConsParaData::sortByIndex(){
  sort(data.begin(), data.end(), ZConsParaSortCriterion());
}
bool ZConsParaData::isIndexUnique(){
  
  for(int i = 0; i < (int)(data.size() - 1); i++)
    for(int j = i + 1; j < (int)(data.size()); j++)
      if(data[i].zconsIndex == data[j].zconsIndex)
        return false;  

  return true;
}

////////////////////////////////////////////////////////////////////////////////
//Implementation of AtomData
////////////////////////////////////////////////////////////////////////////////
AtomData::~AtomData(){
  vector<AtomInfo*>::iterator i;
  AtomInfo* atomInfo;

  for(atomInfo = beginAtomInfo(i); atomInfo; atomInfo  = nextAtomInfo(i))
    delete atomInfo;

  data.clear();
}