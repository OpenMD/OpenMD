/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */
 
/**
 * @file RealSphericalHarmonic.hpp
 * @author Dan Gezelter
 * @date 10/18/2004
 * @version 1.0
 */

#ifndef MATH_REALSPHERICALHARMONIC_HPP
#define MATH_REALSPHERICALHARMONIC_HPP

#include <string.h>

#define RSH_SIN  0
#define RSH_COS  1

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

    void makeSinFunction() {functionType = RSH_SIN;}
    void makeCosFunction() {functionType = RSH_COS;}

    bool isSinFunction() { return functionType == RSH_SIN ? true : false;}
    bool isCosFunction() { return functionType == RSH_COS ? true : false;}
    
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
