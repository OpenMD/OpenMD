/**********************************************************************

This basic Element data-holding class was originally taken from the data.h 
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
 * @file Element.hpp
 * @author gezelter
 * @date  12/21/2007
 * @version 1.0
 */ 

#ifndef PRIMITIVES_ELEMENT_HPP
#define PRIMITIVES_ELEMENT_HPP
#include <string>
#include <cstring>
#include "config.h"

namespace OpenMD{
  class Element {
  public:
    Element()    {}
    /** Constructor
        @param num     Atomic number
        @param sym     Elemental symbol (maximum 3 characters)
        @param ARENeg  Allred-Rochow electronegativity
        @param rcov    Covalent radius (in Angstrom)
        @param rvdw    van der Waals radius (in Angstrom)
        @param maxbo   Maximum bonding valence
        @param mass    Atomic mass (in amu)
        @param elNeg   Electronegativity (in Pauling units)
        @param ionize  Ionization potential (in eV)
        @param elAffin Electron affinity (in eV)
        @param red     RGB value for a suggest visualization color (0 .. 1)
        @param green   RGB value for a suggest visualization color (0 .. 1)
        @param blue    RGB value for a suggest visualization color (0 .. 1)
        @param name Full IUPAC name
    **/
    Element(int num, const char *sym, RealType ARENeg, RealType rcov,
            RealType rvdw, int maxbo, RealType mass, RealType elNeg,
            RealType ionize, RealType elAffin, RealType red, RealType green,
            RealType blue, std::string name) :
      num_(num), name_(name), Rcov_(rcov), Rvdw_(rvdw), mass_(mass), 
      elNeg_(elNeg), ARENeg_(ARENeg), ionize_(ionize), elAffinity_(elAffin), 
      red_(red), green_(green), blue_(blue),
      maxbonds_(maxbo)
    {
      symbol_.assign(sym, 3);
    }
        
    /**
     * Returns the atomic number of this element
     * @return the atomic number of this element
     */
    int GetAtomicNum() {
      return(num_);    
    }
    
    /**
     * Returns the atomic symbol for this element
     * @return the atomic symbol for this element
     */
    const char *GetSymbol() {
      return(symbol_.c_str()); 
    }
 
    /**
     * Returns the covalent radius of this element
     * @return the covalent radius of this element
     */
    RealType GetCovalentRad() {       
      return(Rcov_);
    }

    /**
     * Returns the van der Waals radius of this element
     * @return the van der Waals radius of this element
     */
    RealType GetVdwRad() {
      return(Rvdw_);
    }

    /**
     * Returns the standard atomic mass for this element (in amu)
     * @return the standard atomic mass for this element (in amu)
     */
    RealType GetMass() {
      return(mass_);   
    }

    /**
     * Returns the maximum expected number of bonds to this element
     * @return the maximum expected number of bonds to this element
     */
    int GetMaxBonds() {
      return(maxbonds_);
    }

    /**
     * Returns the Pauling electronegativity for this element
     * @return the Pauling electronegativity for this element
     */
    RealType GetElectroNeg() {
      return(elNeg_);  
    }

    /**
     * Returns the Allred-Rochow electronegativity for this element
     * @return the Allred-Rochow electronegativity for this element
     */
    RealType GetAllredRochowElectroNeg() {
      return(ARENeg_);  
    }

    /**
     * Returns the ionization potential (in eV) of this element
     * @return the ionization potential (in eV) of this element
     */
    RealType GetIonization() {
      return(ionize_);  
    }

    /**
     * Returns the electron affinity (in eV) of this element
     * @return the electron affinity (in eV) of this element
     */
    RealType GetElectronAffinity() {      
      return(elAffinity_);  
    }

    /**
     * Returns the name of this element (in English)
     * @return the name of this element (in English)
     */
    std::string GetName() {
      return(name_);    
    }

    /**
     * Returns the red component of this element's default visualization color
     * @return the red component of this element's default visualization color
     */
    RealType GetRed() {
      return(red_);
    }

    /**
     * Returns the green component of this element's default color
     * @return the green component of this element's default color
     */
    RealType GetGreen() {
      return(green_);   
    }

    /**
     * Returns the blue component of this element's default color
     * @return the blue component of this element's default color
     */
    RealType GetBlue() {
      return(blue_);
    }
    
  protected: 
    int num_;
    std::string symbol_;
    std::string name_;
    RealType Rcov_, Rvdw_, mass_, elNeg_, ARENeg_, ionize_, elAffinity_;
    RealType red_, green_, blue_;
    int maxbonds_;
    
  };
}
#endif
