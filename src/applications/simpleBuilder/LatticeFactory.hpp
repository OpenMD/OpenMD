#ifndef _LATTICEFACTORY_H_
#define _LATTICEFACTORY_H_
#include <map>
#include <string>

using namespace std;

class BaseLatticeCreator;
class BaseLattice;

class LatticeFactory{
public:
	~LatticeFactory();

	static LatticeFactory* getInstance();

	bool registerCreator( BaseLatticeCreator*  latCreator );


    bool hasLatticeCreator( const string& latticeType );

    const string toString();

    BaseLattice* createLattice( const string& latticeType );
    
private:
	LatticeFactory(){}
    static LatticeFactory* instance;
    map<string,  BaseLatticeCreator*> creatorMap;
};
#endif 
