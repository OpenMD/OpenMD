/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifndef UTILS_CONSTANTS_HPP
#define UTILS_CONSTANTS_HPP

#include <config.h>

#ifdef _MSC_VER
#include <math.h>
#else
#include <cmath>
#endif

namespace OpenMD::Constants {

  static constexpr RealType PI     = M_PI;
  static constexpr RealType TWO_PI = 2.0 * PI;

  // Internal units in OpenMD:
  // distance = Angstroms (1.0e-10 m)
  // time = fs (1.0e-15 s)
  // mass = amu
  // charge = electron charge (1.602176634e-19 Coulombs)
  // energy = kcal/mol
  // temperature = K
  // velocity = Angstroms / fs
  // force = kcal/mol / Angstroms

  inline constexpr RealType kb =
      1.9872156E-3;  //!< boltzman's constant in kcal/(mol K)
  inline constexpr RealType kB =
      8.31451e-7;  //!< boltzmann constant amu*Ang^2*fs^-2/K
  inline constexpr RealType c =
      299792458. * 1e-5;  //! speed of light in Ang fs^-1
  inline constexpr RealType epsilon0 =
      2.396451e-4;  //! vacuum permittivity in e^2 Ang^-1 (kcal/mol)^-1
  inline constexpr RealType mu0 =
      1.25663706212e-6;  //! vacuum permeability in N Amp^-2

  inline constexpr RealType energyConvert =
      4.184E-4;  //!< convert kcal/mol -> (amu A^2)/fs^2
  inline constexpr RealType rotationalEnergyConvert = energyConvert * TWO_PI;

  inline constexpr RealType pressureConvert =
      1.63882576e8;  //!< converts amu*fs^-2*Ang^-1 -> atm
  inline constexpr RealType elasticConvert =
      1.66053386e4;  //!< converts amu*fs^-2*Ang^-1 -> GPa
  inline constexpr RealType energyElasticConvert =
      6.947695345;  //!< converts kcal*mol^-1*Ang^-3 -> GPa

  //! \name chargeFieldConvert Converts electron-volts to kcal/mol
  inline constexpr RealType chargeFieldConvert = 23.0609;
  //! \name dipoleFieldConvert  Converts Debye*Volts/Angstroms to kcal/mol
  inline constexpr RealType dipoleFieldConvert = 4.8018969509;
  //! \name dipoleConvert Converts Debye to electron Angstroms
  inline constexpr RealType dipoleConvert = 0.2081943;

  //!\name magneticFieldConvert Converts Tesla to Volts fs/Ang^2
  inline constexpr RealType magneticFieldConvert = 1.0e-5;

  /**
   *  surfaceTensionConvert
   *    multiplies standard input file units of
   *      surfaceTension (Newton / meter)
   *    returns values of
   *      kcal mol^-1 Angstrom^-2
   */
  inline constexpr RealType surfaceTensionConvert =
      1.439326479;  //!< converts N/m to kcal/mol*Ang^-2

  /**
   *  viscoConvert
   *    used for products of:
   *      viscosity (Poise) * distance (Angstroms) * velocity (Angstrom / fs)
   *    returns values of:
   *      force in (kcal mol^-1 Angstrom^-1)
   */
  inline constexpr RealType viscoConvert = 1.439326479e4;

  /**
   *  densityConvert
   *    used for converting amu / Angstroms^3 into  g / cm^3
   */
  inline constexpr RealType densityConvert = 1.66053886;

  /**
   *  thermalConductivityConvert
   *    multiplies standard input file units of
   *      themalConductivity (watts meter^-1 Kelvin^-1)
   *    returns values of:
   *      kcal mol^-1 Angstrom^-1 fs^-1 Kelvin^-1
   */
  inline constexpr RealType thermalConductivityConvert = 1.439326479e-5;

  /**
   * currentConvert
   *   multiplies standard input file units of
   *     electricalCurrent (Amperes)
   *   returns values of:
   *     electrons fs^-1
   */
  inline constexpr RealType currentConvert = 6241.573027317;

  /**
   * currentDensityConvert
   *   multiplies standard input file units of
   *     currentDensity (Amperes m^-2)
   *   returns values of:
   *     electrons fs^-1 Angstrom^-2
   */
  inline constexpr RealType currentDensityConvert = 6.241573027317e-17;

  /**
   * chargeDensityConvert
   *   multiplies standard input file units of
   *     chargeDensity (Coulombs m^-2)
   *   returns values of:
   *     electrons Angstrom^-2
   */
  inline constexpr RealType chargeDensityConvert = 6.241573027317e-2;

  /**
   * concentrationConvert
   *   multiplies standard number density units (Angstrom^-3)
   *   returns values of molarity (1 M = 1 mole / Liter)
   */
  inline constexpr RealType concentrationConvert = 1660.5390404272;

  /**
   * Atomic Units are used in the Slater overlap code, and we need
   * to get distances back and forth to angstroms and energies back
   * and forth to kcal / mol
   */
  inline constexpr RealType angstromToBohr = 1.88972612;
  inline constexpr RealType bohrToAngstrom = 0.52917721092;
  inline constexpr RealType hartreeToKcal  = 627.509469;
}  // namespace OpenMD::Constants

#endif
