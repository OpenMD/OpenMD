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
#ifndef OB_XYZFORMAT_HPP
#define OB_XYZFORMAT_HPP

#include "mol.hpp"
#include "obconversion.hpp"
#include "obmolecformat.hpp"

#ifdef HAVE_SSTREAM
#include <sstream>
#else
#include <strstream>
#endif

using namespace std;
namespace OpenBabel
{

class XYZFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID
  XYZFormat()
  {
    OBConversion::RegisterFormat("xyz", this, "chemical/x-xyz");
  }

  virtual const char* Description() //required
  {
    return
      "XYZ cartesian coordinates format\n \
       Read Options e.g. -as\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n\n";
  };

  virtual const char* SpecificationURL()
  {return "http://openbabel.sourceforge.net/formats/xyz.shtml";}; //optional

  virtual const char* GetMIMEType() 
  { return "chemical/x-xyz"; };

  //*** This section identical for most OBMol conversions ***
  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
};

}

#endif
