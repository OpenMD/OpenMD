/**
 * @file ShapeType.hpp
 * @author Dan Gezelter
 * @date 10/18/2004
 * @version 1.0
 */

#ifndef TYPES_SHAPETYPE_HPP
#define TYPES_SHAPETYPE_HPP

#include <fstream>
#include <vector>
#include "math/RealSphericalHarmonic.hpp"
#include "math/SquareMatrix3.hpp"

namespace oopse {
  using namespace std;
  
  class ShapeType {
    
  public: 
    
    ShapeType(void);
    ~ShapeType(void);  
    
    char *getName(void) {return shapeName;}
    void setName(char * name) {shapeName = strdup(name);}
    
    double getMass(void) {return mass;}
    void setMass(double m) {mass = m;}
    
    Mat3x3d getI(void) {return I;}
    void setI(Mat3x3d theI) {I = theI;}
    
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
    
    char *shapeName;             // The name of the shape
    double mass;                 // The mass
    Mat3x3d I;                // The moment of inertia tensor
    vector<RealSphericalHarmonic*> contactFuncs;  // The contact functions
    vector<RealSphericalHarmonic*> rangeFuncs;    // The range functions
    vector<RealSphericalHarmonic*> strengthFuncs; // The strength functions
    
  }; 
}
#endif

