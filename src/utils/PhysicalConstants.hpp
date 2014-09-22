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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef UTILS_PHYSICALCONSTANTS_HPP
#define UTILS_PHYSICALCONSTANTS_HPP
#include "config.h"
#include <cmath>

namespace OpenMD {

  namespace PhysicalConstants {
    const RealType kb = 1.9872156E-3; // boltzman's constant in kcal/(mol K)
    const RealType kB = 8.31451e-7;   // boltzmann constant amu*Ang^2*fs^-2/K
    const RealType energyConvert = 4.184E-4; // convert kcal/mol -> (amu A^2)/fs^2
    const RealType rotationalEnergyConvert = energyConvert*2.0*M_PI;
    
    const RealType pressureConvert = 1.63882576e8; // converts amu*fs^-2*Ang^-1 -> atm

    //! \name chargeFieldConvert Converts electron-volts to kcal/mol
    const RealType chargeFieldConvert = 23.0609; 
    //! \name dipoleFieldConvert  Converts Debye*Volts/Angstroms to kcal/mol
    const RealType dipoleFieldConvert = 4.8018969509; 

    /* 
     *  surfaceTensionConvert   
     *    multiplies standard input file units of 
     *      surfaceTension (Newton / meter)
     *    returns values of
     *      kcal mol^-1 Angstrom^-2
     */
    const RealType surfaceTensionConvert = 1.439326479; // convert N/m to kcal/mol*Ang^-2
    
    /* 
     *  viscoConvert   
     *    used for products of:
     *      viscosity (Poise) * distance (Angstroms) * velocity (Angstrom / fs)
     *    returns values of:
     *      force in (kcal mol^-1 Angstrom^-1) 
     */
    const RealType viscoConvert = 1.439326479e4; 

    /*
     *  densityConvert
     *    used for converting amu / Angstroms^3 into  g / cm^3
     */
    const RealType densityConvert = 1.66053886;

    /* 
     *  thermalConductivityConvert   
     *    multiplies standard input file units of 
     *      themalConductivity (watts meter^-1 Kelvin^-1)  
     *    returns values of:
     *      kcal mol^-1 Angstrom^-1 fs^-1 Kelvin^-1
     */
    const RealType thermalConductivityConvert = 1.439326479e-5; 

    /* Atomic Units are used in the Slater overlap code, and we need
     * to get distances back and forth to angstroms and energies back 
     * and forth to kcal / mol 
     */
    const RealType angstromToBohr = 1.88972612;
    const RealType bohrToAngstrom = 0.52917721092;
    const RealType hartreeToKcal   = 627.509469;
  }
}
#endif 
