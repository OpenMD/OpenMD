/**********************************************************************
Copyright (C) 2002 by Steffen Reith <streit@streit.cc>
Some portions Copyright (C) 2003-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_POVRAYFORMAT_HPP
#define OB_POVRAYFORMAT_HPP

/* ---- C includes ---- */
#include <math.h>
#include <time.h>
#include <stdlib.h>

/* ---- OpenBabel include ---- */
#include "config.h"
#include "mol.hpp"
#include "obconversion.hpp"

/* ---- C++ includes ---- */
#include <string>
#if defined(HAVE_SSTREAM)
#include <sstream>
#else
#include <strstream>
#endif

/* ---- Max. length of a atom-label ---- */
#define StrLen 32

/* ---- Define max. length of domainname ---- */
#define MAXDOMAINNAMELEN 256

/* ---- Maximal radius of an atom. Needed for bounding box ---- */
#define MAXRADIUS (double) 3.0

/* ---- Define index of first atom if needed ---- */
#ifndef MIN_ATOM
#define MIN_ATOM 1
#endif

/* ---- Size of time-string ---- */
#define TIME_STR_SIZE 64

/* ---- if x < = EPSILON then x = 0.0 ---- */
#define EPSILON (double) 1e-4

/* ---- Define makro for calculating x^2 ---- */
#ifdef SQUARE
#undef SQUARE
#endif
#define SQUARE(x) ((x) * (x))

/* ---- Define PI (if needed) ---- */
#ifndef PI
#define PI ((double) 3.1415926535897932384626433)
#endif

/* ---- Convert RAD to DEG ---- */
#define RAD2DEG(r) (((double) 180.0 * r) / PI)

using namespace std;
namespace OpenBabel
{

class PovrayFormat : public OBFormat
{
public:
    //Register this format type ID
    PovrayFormat()
    {
        OBConversion::RegisterFormat("pov",this);
    }

    virtual const char* Description() //required
    {
        return
            "POV-Ray input format\n \
            No comments yet\n";
    };

  virtual const char* SpecificationURL()
  {return "http://www.povray.org/";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return NOTREADABLE | WRITEONEONLY;
    };

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
        obErrorLog.ThrowError(__FUNCTION__,
                              auditMsg,
                              obAuditMsg);
        delete pOb;
        return ret;
    };
};

}
#endif
