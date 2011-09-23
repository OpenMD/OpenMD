/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
#ifndef BRAINS_THERMO_HPP
#define BRAINS_THERMO_HPP

#include "primitives/Atom.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {

  class Thermo{

  public:

    Thermo( SimInfo* info) : info_(info) {}

    // note: all the following energies are in kcal/mol

    RealType getKinetic(); // the total kinetic energy 
    RealType getPotential(); // the total potential energy
    RealType getTotalE(); // gets the total energy

    RealType getTemperature(); // gives the instant temp. in K

    RealType getPressure(); // gives the instant pressure in atm;
    RealType getPressureX() { return getPressure(0); }
    RealType getPressureY() { return getPressure(1); }
    RealType getPressureZ() { return getPressure(2); }

    // gives the pressure tensor in amu*fs^-2*Ang^-1
    Mat3x3d getPressureTensor(); 
    RealType getVolume();   // gives the volume in Ang^3 

    // accumulate and return the simulation box dipole moment in C*m
    Vector3d getBoxDipole(); 
    // accumulate and return helfand thermal momoment
    Vector3d getThermalHelfand();
    
    void saveStat();
    
  private:
    RealType getPressure(int direction);
    
    SimInfo* info_;
  };
  
} //end namespace OpenMD
#endif
