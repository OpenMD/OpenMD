#include "LatticeCreator.hpp"
#include "LatticeFactory.hpp"

BaseLatticeCreator::BaseLatticeCreator(const string& latType){
  latticeType = latType;
  LatticeFactory::getInstance()->registerCreator(this);
}

