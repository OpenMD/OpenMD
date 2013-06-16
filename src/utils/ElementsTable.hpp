/**********************************************************************

This basic Periodic Table class was originally taken from the data.h 
file in OpenBabel. The code has been modified to match the OpenMD coding style.

We have retained the OpenBabel copyright and GPL license on this class:  

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/**
 * @file ElementsTable.hpp
 * @author gezelter
 * @date 12/21/2007
 * @version 1.0
 */

#ifndef UTILS_ELEMENTSTABLE_HPP
#define UTILS_ELEMENTSTABLE_HPP

#include "config.h"
#include <vector>
#include "primitives/Element.hpp"

namespace OpenMD {

  /**
   * @class ElementsTable
   * @brief Periodic Table of the Elements 
   * Using element data is a place holder when we lack information about
   * a specific atom type.  In particular, the Langevin algorithms must
   * assume specific atomic radii to predict drag and random forces on those
   * atoms.  For force fields which do not specify Lennard-Jones radii,
   * the element's van der Waals radius is used instead.
   * The ElementsTable class (etab) is declared as external in
   * ElementsTable.cpp. Source files that include the header file 
   * ElementsTable.hpp automatically have an extern definition to etab. 
   * The following code sample demonstrates the use of the ElementsTable class:
   * @code
   * cout << "The symbol for element 6 is " << etab.GetSymbol(6) << endl;
   * cout << "The atomic number for Sulfur is " << etab.GetAtomicNum(16) << endl;
   * cout << "The van der Waal radius for Nitrogen is " << etab.GetVdwRad(7);
   * @endcode
   * Stored information in the OBElementTable includes elemental:
   *   - symbols
   *   - covalent radii
   *   - van der Waal radii
   *   - expected maximum bonding valence
   *   - molar mass (by IUPAC recommended atomic masses)
   *   - electronegativity
   *   - ionization potential
   *   - electron affinity
   *   - RGB colors for visualization programs
   *   - names (by IUPAC recommendation)
   */
  class ElementsTable {
  public:
    /** Constructor */
    ElementsTable(); 
    /** Destructor */
    ~ElementsTable();
    
    /**
     * Read in the data file.
     */
    void  Init();
    /**
     * Set the directory before calling Init()
     */
    void  SetReadDirectory(char *dir) { dir_ = dir; }
    /**
     * Set the environment variable to use before calling Init()
     */
    void  SetEnvironmentVariable(char *var) { envvar_ = var; }
    /**
     * Specified by particular table classes (parses an individual data line)
     * @param line the data line to parse
     */
    void ParseLine(const char *line);
    /**
     * @return the number of elements in the periodic table
     */
    unsigned int GetNumberOfElements();
    unsigned int GetSize() { return GetNumberOfElements(); }
    /**
     * @return the atomic number matching the element symbol passed
     * or 0 if not defined. 
     * @param str the element symbol
     */
    int GetAtomicNum(const char *str);
    /**
     * @return the atomic number matching the element symbol passed
     * or 0 if not defined. For 'D' or 'T' hydrogen isotopes, will return
     * a value in the second argument
     * @param str the element symbol
     * @param iso the isotope index for Deuterium or Tritium
     */
    int GetAtomicNum(const char *str, int &iso);
    /**
     * @return the element symbol matching the atomic number passed
     * @param atomicnum the atomic number of the element
     */
    const char *GetSymbol(int atomicnum);
    /**
     * @return the van der Waals radius for this atomic number
     * @param atomicnum the atomic number of the element
     */
    RealType GetVdwRad(int atomicnum);
    /**
     * @return the covalent radius for this atomic number
     * @param atomicnum the atomic number of the element
     */
    RealType GetCovalentRad(int atomicnum);
    /**
     * @return the average atomic mass for this element.
     * @param atomicnum the atomic number of the element
     */
    RealType GetMass(int atomicnum);
    /**
     * @return a "corrected" bonding radius based on the hybridization.
     * Scales the covalent radius by 0.95 for sp2 and 0.90 for sp hybrids
     * @param atomicnum the atomic number of the element
     * @param hyb the hybridization of the element
     */
    RealType CorrectedBondRad(int atomicnum, int hyb = 3); 
    /**
     * @return a "corrected" vdW radius based on the hybridization.
     * Scales the van der Waals radius by 0.95 for sp2 and 0.90 for sp hybrids
     * @param atomicnum the atomic number of the element
     * @param hyb the hybridization of the element
     */
    RealType CorrectedVdwRad(int atomicnum, int hyb = 3);
    /**
     * @return the maximum expected number of bonds to this element
     * @param atomicnum the atomic number of the element
     */
    int	GetMaxBonds(int atomicnum);
    /**
     * @return the Pauling electronegativity for this element
     * @param atomicnum the atomic number of the element
     */
    RealType GetElectroNeg(int atomicnum);
    /**
     * @return the ionization potential (in eV) for this element
     * @param atomicnum the atomic number of the element
     */
    RealType GetIonization(int atomicnum);
    /**
     * @return the electron affinity (in eV) for this element
     * @param atomicnum the atomic number of the element
     */
    RealType GetElectronAffinity(int atomicnum);
    /**
     * @return a vector with red, green, blue color values for this element
     * @param atomicnum the atomic number of the element
     */
    std::vector<RealType> GetRGB(int atomicnum);

    /**
     * @return the name of this element
     * @param atomicnum the atomic number of the element
     */
    std::string GetName(int atomicnum);
    
  protected:
    bool         init_;		//!< whether the data been read already
    std::string  filename_;	//!< file to search for
    std::string  dir_;		//!< data directory for file if _envvar fails
    std::string  subdir_;	//!< subdirectory (if using environment variable)
    std::string  envvar_;	//!< environment variable to check first
    std::vector<Element*> elements_;
    const char  *dataptr_;      //!< default data table if file is unreadable

  };

  extern ElementsTable etab;
}

#endif
