#ifndef VISITOR_ATOMDATA_HPP
#define VISITOR_ATOMDATA_HPP
#include "utils/GenericData.hpp"

using namespace std;

namespace oopse {

struct AtomInfo {
  public:
    string AtomType;
    double pos[3];
    double dipole[3];  
};

class AtomData : public GenericData{
  public:
    AtomData(const string& id = "ATOMDATA") : GenericData(id) {}
    ~AtomData() {
        vector<AtomInfo*>::iterator i;
        AtomInfo* atomInfo;

        for(atomInfo = beginAtomInfo(i); atomInfo; atomInfo  = nextAtomInfo(i))
            delete atomInfo;

        data.clear();
    }
    void addAtomInfo(AtomInfo* info) {data.push_back(info);}
    void clearAllAtomInfo();
    AtomInfo* beginAtomInfo(vector<AtomInfo*>::iterator& i){
      i = data.begin();
      return i != data.end()? *i : NULL;
    }
    AtomInfo* nextAtomInfo(vector<AtomInfo*>::iterator& i){
      ++i;
      return i != data.end()? *i: NULL;
    }
    vector<AtomInfo*> getData() {return data;}
    int getSize() {return data.size();}
  protected:
    vector<AtomInfo*> data;
};


}
#endif //VISITOR_ATOMDATA_HPP
