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

#include "tinkerformat.hpp"

using namespace std;
namespace OpenBabel
{

bool TinkerFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    OBBond *bond;
    vector<OBEdgeBase*>::iterator j;

    sprintf(buffer,"%6d %-20s   MM2 parameters",mol.NumAtoms(),mol.GetTitle());
    ofs << buffer << endl;

    ttab.SetFromType("INT");

    OBAtom *atom;
    string str,str1;
    for(i = 1;i <= mol.NumAtoms(); i++)
    {
        atom = mol.GetAtom(i);
        str = atom->GetType();
        ttab.SetToType("MM2");
        ttab.Translate(str1,str);
        sprintf(buffer,"%6d %2s  %12.6f%12.6f%12.6f %5d",
                i,
                etab.GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ(),
                atoi((char*)str1.c_str()));
        ofs << buffer;

        for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
        {
            sprintf(buffer,"%6d", (bond->GetNbrAtom(atom))->GetIdx());
            ofs << buffer;
        }

        ofs << endl;
    }

    return(true);
}

} //namespace OpenBabel
