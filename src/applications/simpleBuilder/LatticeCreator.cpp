#include "applications/simpleBuilder/LatticeCreator.hpp"
#include "applications/simpleBuilder/LatticeFactory.hpp"

BaseLatticeCreator::BaseLatticeCreator(const string& latType){
  latticeType = latType;
  LatticeFactory::getInstance()->registerCreator(this);
}

