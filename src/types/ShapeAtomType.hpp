/**
 * @file ShapeAtomType.hpp
 * @author Dan Gezelter
 * @date 10/18/2004
 * @version 1.0
 */

#ifndef TYPES_SHAPEATOMTYPE_HPP
#define TYPES_SHAPEATOMTYPE_HPP

#include <vector>
#include "math/RealSphericalHarmonic.hpp"
#include "math/SquareMatrix3.hpp"
#include "types/DirectionalAtomType.hpp"

namespace oopse {
  using namespace std;
  
  class ShapeAtomType : public DirectionalAtomType {
    
  public: 
    
    ShapeAtomType() : DirectionalAtomType() { atp.is_Shape = 1; }
    ~ShapeAtomType();
        
    vector<RealSphericalHarmonic*> getContactFuncs(void) {return contactFuncs;}
    vector<RealSphericalHarmonic*> getRangeFuncs(void) {return rangeFuncs;}
    vector<RealSphericalHarmonic*> getStrengthFuncs(void) {return strengthFuncs;}
    
    void setContactFuncs(vector<RealSphericalHarmonic*> cf) {
      contactFuncs = cf;
    }
    void setRangeFuncs(vector<RealSphericalHarmonic*> rf) {
      rangeFuncs = rf;
    }
    void setStrengthFuncs(vector<RealSphericalHarmonic*> sf) {
      strengthFuncs = sf;
    }
    
    /**
     * Gets the value of the contact function at a particular orientation
     * @param costheta
     * @param phi
     */
    double getContactValueAt(double costheta, double phi);
    
    /**
     * Gets the value of the range function at a particular orientation
     * @param costheta
     * @param phi
     */
    double getRangeValueAt(double costheta, double phi);
    
    /**
     * Gets the value of the strength function at a particular orientation
     * @param costheta
     * @param phi
     */
    double getStrengthValueAt(double costheta, double phi);
    
  private:
    
    vector<RealSphericalHarmonic*> contactFuncs;  // The contact functions
    vector<RealSphericalHarmonic*> rangeFuncs;    // The range functions
    vector<RealSphericalHarmonic*> strengthFuncs; // The strength functions
    
  }; 
}
#endif

