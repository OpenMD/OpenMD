/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#ifndef OB_AMBERFORMAT_HPP
#define OB_AMBERFORMAT_HPP

#include "mol.hpp"
#include "obconversion.hpp"
#include "obmolecformat.hpp"

using namespace std;
namespace OpenBabel
{

class AmberPrepFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    AmberPrepFormat()
    {
        OBConversion::RegisterFormat("prep",this);
    }

  virtual const char* Description() //required
  {
    return
      "Amber Prep format\n \
       Read Options e.g. -as\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n\n";
  };

  virtual const char* SpecificationURL()
  {return "http://www.amber.ucsf.edu/amber/formats.html";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return NOTWRITABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
};

}
#endif
