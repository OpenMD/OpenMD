/**
 * @file RealSphericalHarmonic.hpp
 * @author Dan Gezelter
 * @date 10/18/2004
 * @version 1.0
 */

#ifndef MATH_REALSPHERICALHARMONIC_HPP
#define MATH_REALSPHERICALHARMONIC_HPP

#include <string.h>

#define SH_SIN  0
#define SH_COS  1

namespace oopse {
  class RealSphericalHarmonic {
  public:
    
    RealSphericalHarmonic();
    virtual ~RealSphericalHarmonic() {} 
    
    void setL(int theL) { L = theL; };
    int getL() { return L; }
    
    void setM(int theM) { M = theM; };
    int getM() { return M; }
    
    void setCoefficient(double co) {coefficient = co;}
    double getCoefficient() {return coefficient;}

    void setFunctionType(short int theType) {functionType = theType;}
    short int getFunctionType() { return functionType; }

    void makeSinFunction() {functionType = SH_SIN;}
    void makeCosFunction() {functionType = SH_COS;}

    bool isSinFunction() { return functionType == SH_SIN ? true : false;}
    bool isCosFunction() { return functionType == SH_COS ? true : false;}
    
    double getValueAt(double costheta, double phi);
    
  protected:
    
    double LegendreP (int l, int m, double x);
    
    int L;
    int M;
    short int functionType;
    double coefficient;
    
  };
}

#endif
