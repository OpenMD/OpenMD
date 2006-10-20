/**********************************************************************
Copyright (C) 2005 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#ifndef OB_FINGERPRINTFORMAT_HPP
#define OB_FINGERPRINTFORMAT_HPP

#include "mol.hpp"
#include "obconversion.hpp"
#include "obmolecformat.hpp"
#include "fingerprint.hpp"
#include <vector>
#include <string>
#include <iomanip>

using namespace std;
namespace OpenBabel {

/// \brief Constructs and displays fingerprints. For details see OBFingerprint class
class FingerprintFormat : public OBMoleculeFormat
{
public:
	//Register this format type ID
	FingerprintFormat() {OBConversion::RegisterFormat("fpt",this);}

	virtual const char* Description() //required
	{ return
"Fingerprint format\n \
Constructs and displays fingerprints and (for multiple input objects)\n \
the Tanimoto coefficient and whether a superstructure of the first object\n \
Options e.g. -xfFP3 -xn128\n \
 f<id> fingerprint type\n \
 N# fold to specified number of bits, 32, 64, 128, etc.\n \
 h  hex output when multiple molecules\n \
 F  displays the available fingerprint types\n \
";
	};

	virtual unsigned int Flags(){return NOTREADABLE;};
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

private:
	vector<unsigned int> firstfp;
	string firstname;
	bool IsPossibleSubstructure(vector<unsigned int>Mol, vector<unsigned int>Frag);
};

}
#endif

