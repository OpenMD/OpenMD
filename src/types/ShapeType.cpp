#include "types/ShapeType.hpp"

using namespace oopse;

ShapeType::ShapeType(void) {
  mass = 0.0;
  I = new Matrix3x3d();
  
}

ShapeType::~ShapeType(void) {
}

double ShapeType::getContactValueAt(double costheta, double phi) {

  vector<RealSphericalHarmonic*>::iterator contactIter;
  double contactVal;

  contactVal = 0.0;

  for(contactIter = contactFuncs.begin();  contactIter != contactFuncs.end(); 
      ++contactIter)     
    contactVal += (*contactIter)->getValueAt(costheta, phi);
  
  return contactVal;
}

double ShapeType::getRangeValueAt(double costheta, double phi) {
  
  vector<RealSphericalHarmonic*>::iterator rangeIter;
  double rangeVal;
  
  rangeVal = 0.0;
  
  for(rangeIter = rangeFuncs.begin();  rangeIter != rangeFuncs.end(); 
      ++rangeIter)     
    rangeVal += (*rangeIter)->getValueAt(costheta, phi);
  
  return rangeVal;
}

double ShapeType::getStrengthValueAt(double costheta, double phi) {

  vector<RealSphericalHarmonic*>::iterator strengthIter;
  double strengthVal;

  strengthVal = 0.0;

  for(strengthIter = strengthFuncs.begin();  
      strengthIter != strengthFuncs.end(); 
      ++strengthIter)     
    strengthVal += (*strengthIter)->getValueAt(costheta, phi);

  return strengthVal;
}
