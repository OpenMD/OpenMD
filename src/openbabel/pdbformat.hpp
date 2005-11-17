/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2003-2005 Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#ifndef OB_PDBFORMAT_HPP
#define OB_PDBFORMAT_HPP

#include "config.h"
#include "mol.hpp"
#include "obconversion.hpp"
#include "obmolecformat.hpp"

#if !HAVE_SNPRINTF
extern "C" int snprintf( char *, size_t, const char *, /* args */ ...);
#endif

#include <vector>
#include <map>

#ifdef HAVE_SSTREAM
#include <sstream>
#else
#include <strstream>
#endif

using namespace std;
namespace OpenBabel
{

class PDBFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    PDBFormat()
    {
        OBConversion::RegisterFormat("pdb",this, "chemical/x-pdb");
        OBConversion::RegisterFormat("ent",this, "chemical/x-pdb");
    }

  virtual const char* Description() //required
  {
    return
      "Protein Data Bank format\n \
       Read Options e.g. -as\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n\n";
  };

  virtual const char* SpecificationURL()
  { return "http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html";};

  virtual const char* GetMIMEType() 
  { return "chemical/x-pdb"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return READONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

};

}
#endif
