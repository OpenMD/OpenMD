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
#ifndef OB_OOPSEFORMAT_HPP
#define OB_OOPSEFORMAT_HPP

#include <set>

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

class OOPSEFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID
  OOPSEFormat()
  {
    OBConversion::RegisterFormat("in", this, "chemical/x-in");
  }

  virtual const char* Description() //required
  {
    return
      "oopse cartesian coordinates format\n";
  };

  virtual const char* SpecificationURL()
  { return "";}; //optional

  virtual const char* GetMIMEType() 
  { return "chemical/x-in"; };

   virtual unsigned int Flags() { return WRITEONEONLY;}

  //*** This section identical for most OBMol conversions ***
  ////////////////////////////////////////////////////
  /// The "API" interface functions
  //virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
      bool AreSameFragments(OBMol& mol, vector<int>& frag1, vector<int>& frag2);
      void findAngles(OBMol& mol);
      OBMol* createMolFromFragment(OBMol& mol, vector<int>& fragment);
      void WriteMDFile(vector<OBMol*> mols, vector<int> numMols);
      void WriteINFile(OBMol& mol, ostream& ofs, vector<int>& indices);
};
}
#endif
