#ifndef _LATTICECREATOR_H_
#define _LATTICECREATOR_H_
#include <string>
#include "applications/simpleBuilder/BaseLattice.hpp"

using namespace std;

class BaseLatticeCreator{
  public:
    virtual BaseLattice* createLattice() = 0;
    const string& getType() {return latticeType;}
    
  protected:
    BaseLatticeCreator(const string& latType);
  private:
    string latticeType;
};

template<class LatticeClass>
class LatticeCreator : public BaseLatticeCreator
{
public:
	LatticeCreator(const string& latticeType): BaseLatticeCreator(latticeType) {}
	virtual BaseLattice* createLattice()  { return new LatticeClass();}
};

#endif
