#include "applications/simpleBuilder/LatticeFactory.hpp"
#include "applications/simpleBuilder/BaseLattice.hpp"
#include "applications/simpleBuilder/LatticeCreator.hpp"

LatticeFactory* LatticeFactory::instance = NULL;
LatticeFactory::~LatticeFactory(){
  map<string, BaseLatticeCreator*>::iterator mapIter;
	for (mapIter = creatorMap.begin(); mapIter == creatorMap.end(); ++mapIter) {
		delete mapIter->second;
	}  
}

 LatticeFactory* LatticeFactory::getInstance(){
	if (instance == NULL)
      instance = new LatticeFactory();
	
  	return instance;
}

bool LatticeFactory::registerCreator( BaseLatticeCreator*  latCreator ){  
  string latticeType = latCreator->getType();
  map<string, BaseLatticeCreator*>::iterator mapIter;

  mapIter = creatorMap.find(latticeType);

  if (mapIter == creatorMap.end()) {
    creatorMap[ latticeType ] = latCreator;
    return true;
  }
  else{
    delete mapIter->second;
    mapIter->second = latCreator;
    return false;
  }
}

BaseLattice* LatticeFactory::createLattice( const string& latticeType ){ 
  map<string, BaseLatticeCreator*>::iterator mapIter; 
  
  mapIter = creatorMap.find(latticeType);

  if (mapIter != creatorMap.end()) {
     return (mapIter->second)->createLattice();
  }
  else
    return NULL;
}

bool LatticeFactory::hasLatticeCreator( const string& latticeType ){
  map<string, BaseLatticeCreator*>::iterator mapIter;

  mapIter = creatorMap.find(latticeType);

  if (mapIter != creatorMap.end())
    return true;
  else
    return false;
}

const string LatticeFactory::toString(){
  string result;
  map<string, BaseLatticeCreator*>::iterator mapIter;

  result = "Avaliable lattice creators in LatticeFactory are:\n";
  
  for(mapIter = creatorMap.begin(); mapIter != creatorMap.end(); ++mapIter){
    result += mapIter->first + " ";
  }
  
  result += "\n";

  return result;
}

