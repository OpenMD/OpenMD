#ifndef USETHEFORCE_FORCEFIELD_HPP
#define USETHEFORCE_FORCEFIELD_HPP

#define MK_STR(s) # s
#define STR_DEFINE(t, s) t = MK_STR(s)

#include <utilities>

#include "utils/Tuple.hpp"
#include "types/ShapeAtomType.hpp"
#include "io/basic_ifstrstream.hpp"

using namespace std;
using namespace oopse;

class ForceField{

public:

  ForceField(){ 
    hasVariant=false;
    ffPath = getenv("FORCE_PARAM_PATH");
    if( ffPath.empty() ) {
      STR_DEFINE(ffPath, FRC_PATH );
    }   
  }

  virtual ~ForceFields(){}

  void setVariant(const string &variant) { hasVariant = true; theVariant = variant; }
  virtual void readParams( void ) = 0;  
  
  AtomType* getMatchingAtomType(const string &at);
  BondType* getMatchingBondType(const string &at1, const string &at2);
  BendType* getMatchingBendType(const string &at1, const string &at2,
                                const string &at3);
  TorsionType* getMatchingTorsionType(const string &at1, const string &at2,
                                      const string &at3, const string &at4);

  double getRcutForAtomType(AtomType* at);
  
protected:
  
  string ffPath;
  ifstrstream forceFile;
  bool hasVariant;
  string variant;
  map<string, AtomType*> atomTypeMap;
  map<pair<string,string>, BondType*> bondTypeMap;
  map<tuple3<string,string,string>, BendType*> bendTypeMap;
  map<tuple4<string,string,string,string>, TorsionType*> torsionTypeMap;
  string wildCardAtomTypeName;

};


#endif

