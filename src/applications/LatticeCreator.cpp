#include "applications/LatticeCreator.hpp"
#include "applications/LatticeFactory.hpp"

BaseLatticeCreator::BaseLatticeCreator(const string& latType){
  latticeType = latType;
  LatticeFactory::getInstance()->registerCreator(this);
}

