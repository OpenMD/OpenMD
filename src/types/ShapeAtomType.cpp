#include "types/ShapeAtomType.hpp"

using namespace oopse;

ShapeAtomType::~ShapeAtomType() {
  vector<RealSphericalHarmonic*>::iterator iter;
  for (iter = contactFuncs.begin(); iter != contactFuncs.end(); ++iter) 
    delete (*iter);
  for (iter = rangeFuncs.begin(); iter != rangeFuncs.end(); ++iter) 
    delete (*iter);
  for (iter = strengthFuncs.begin(); iter != strengthFuncs.end(); ++iter) 
    delete (*iter);
  contactFuncs.clear();
  rangeFuncs.clear();
  strengthFuncs.clear();
}

double ShapeAtomType::getContactValueAt(double costheta, double phi) {
  
  vector<RealSphericalHarmonic*>::iterator contactIter;
  double contactVal;
  
  contactVal = 0.0;
  
  for(contactIter = contactFuncs.begin();  contactIter != contactFuncs.end(); 
      ++contactIter) 
    contactVal += (*contactIter)->getValueAt(costheta, phi);

  return contactVal;
}

double ShapeAtomType::getRangeValueAt(double costheta, double phi) {
  
  vector<RealSphericalHarmonic*>::iterator rangeIter;
  double rangeVal;
  
  rangeVal = 0.0;
  
  for(rangeIter = rangeFuncs.begin();  rangeIter != rangeFuncs.end(); 
      ++rangeIter)     
    rangeVal += (*rangeIter)->getValueAt(costheta, phi);
  
  return rangeVal;
}

double ShapeAtomType::getStrengthValueAt(double costheta, double phi) {
  
  vector<RealSphericalHarmonic*>::iterator strengthIter;
  double strengthVal;
  
  strengthVal = 0.0;
  
  for(strengthIter = strengthFuncs.begin();  
      strengthIter != strengthFuncs.end(); 
      ++strengthIter)     
    strengthVal += (*strengthIter)->getValueAt(costheta, phi);
  
  return strengthVal;
}
