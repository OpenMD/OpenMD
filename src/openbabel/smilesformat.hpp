/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
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
//Contains SMIFormat and FIXFormat classes

#ifndef OB_SMILESFORMAT_HPP
#define OB_SMILESFORMAT_HPP

#include "mol.hpp"
#include "obconversion.hpp"
#include "obmolecformat.hpp"

using namespace std;

namespace OpenBabel
{
class SMIFormat : public OBMoleculeFormat
{
public:
	//Register this format type ID
	SMIFormat()
	{
		OBConversion::RegisterFormat("smi",this, "chemical/x-daylight-smiles");
		OBConversion::RegisterOptionParam("n", this);
		OBConversion::RegisterOptionParam("t", this);
	}

  virtual const char* GetMIMEType() 
  { return "chemical/x-daylight-smiles"; };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ///////////////////////////////////////////////////////

    virtual const char* Description()
    {
        return
            "SMILES format\n \
            A linear text format which can describe the connectivity\n \
            and chirality of a molecule\n \
            Write Options e.g. -xt\n \
            -n no molecule name\n \
            -t molecule name only\n \
	    -r radicals lower case eg ethyl is Cc\n\n";
    };

  virtual const char* SpecificationURL()
  {return "http://www.daylight.com/dayhtml/smiles/";}; //optional

		virtual int SkipObjects(int n, OBConversion* pConv)
		{
			if(n==0) return 1; //already points after current line
			string temp;
			istream& ifs = *pConv->GetInStream();
			int i;
			for(i=0;i<n && ifs.good();i++)
				getline(ifs, temp);
			return ifs.good() ? 1 : -1;	
		};	
};

}

#endif

