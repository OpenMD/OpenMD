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
#ifndef OB_TINKERFORMAT_HPP
#define OB_TINKERFORMAT_HPP

#include "mol.hpp"
#include "obconversion.hpp"

using namespace std;
namespace OpenBabel
{

class TinkerFormat : public OBFormat
{
public:
    //Register this format type ID
    TinkerFormat()
    {
        OBConversion::RegisterFormat("txyz",this);
    }

    virtual const char* Description() //required
    {
        return
            "Tinker MM2 format\n \
            No comments yet\n";
    };

  virtual const char* SpecificationURL()
  {return "http://dasher.wustl.edu/tinker/";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return NOTREADABLE | WRITEONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ////////////////////////////////////////////////////
    /// The "Convert" interface functions
    virtual bool WriteChemObject(OBConversion* pConv)
    {
        //Retrieve the target OBMol
        OBBase* pOb = pConv->GetChemObject();
        OBMol* pmol = dynamic_cast<OBMol*> (pOb);
        bool ret=false;
        if(pmol)
            ret=WriteMolecule(pmol,pConv);

	std::string auditMsg = "OpenBabel::Write molecule ";
	std::string description(Description());
        auditMsg += description.substr( 0, description.find('\n') );
        obErrorLog.ThrowError(__func__,
                              auditMsg,
                              obAuditMsg);

        delete pOb;
        return ret;
    };
};
}
#endif

