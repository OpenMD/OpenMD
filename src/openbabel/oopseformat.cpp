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

#include "oopseformat.hpp"
#include <fstream>
namespace OpenBabel
{

bool OOPSEFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;
    
    vector<vector<int> > fragmentLists;
    pmol->ContigFragList(fragmentLists);
    OBBitVec unused;
    vector<bool> used(fragmentLists.size(), 0);
    vector<vector<int> > molecules;
    vector<int> indices;
    for(int i =0; i < used.size(); ++i) {
        if (used[i])
        {
            continue;
        }
        used[i] = true;
        vector<int> sameMolTypes;
        sameMolTypes.push_back(i);
        indices.insert(indices.end(), fragmentLists[i].begin(), fragmentLists[i].end());
        for (int j = i + 1;j < used.size(); ++j)
        {
            if (used[j])
            {
                continue;
            }
            
            if (AreSameFragments(*pmol, fragmentLists[i], fragmentLists[j]))
            {
                sameMolTypes.push_back(j);
                indices.insert(indices.end(), fragmentLists[j].begin(), fragmentLists[j].end());
                used[j]=true;
            }
        }
        molecules.push_back(sameMolTypes);
        
    }

    //
    vector<OBMol*> mdMols;    
    vector<int> numMols;
    for(vector<vector<int> >::iterator  i = molecules.begin(); i != molecules.end(); ++i) 
    {
        mdMols.push_back(createMolFromFragment(*pmol, fragmentLists[i->front()]));
        numMols.push_back((*i).size());
    }

    string OutputFileName = pConv->GetInFilename();
    unsigned int pos = OutputFileName.rfind(".");
    if(pos==string::npos)
        OutputFileName += ".md";
    else
        OutputFileName = OutputFileName.substr(0, pos) + ".md";       
    ofstream ofs(OutputFileName.c_str());
    if(!ofs)
    {
    	 cerr << "Cannot write to " << OutputFileName <<endl;
   	 return false;
    }
                
    WriteMDFile(mdMols, numMols, ofs);    

    for(vector<OBMol*>::iterator  i = mdMols.begin(); i != mdMols.end(); ++i) 
    {
        delete *i;
    }

    //    
    WriteINFile(*pmol, *pConv->GetOutStream(), indices);



    return(true);
}

bool OOPSEFormat::AreSameFragments(OBMol& mol, vector<int>& frag1, vector<int>& frag2)
{
    if (frag1.size() != frag2.size())
    {
        return false;
    }

    //exact graph matching is a NP complete problem, 
    for (unsigned int i =0 ; i < frag1.size(); ++i)
    {
        if (strcmp(mol.GetAtom(frag1[i])->GetType(), mol.GetAtom(frag2[i])->GetType()) )
        {
            return false;
        }
    }
    return true;
}

struct SameAngle
{
  bool operator()(const triple<OBAtom*,OBAtom*,OBAtom*> t1, const triple<OBAtom*,OBAtom*,OBAtom*> t2) const
  {
    return (t1.second == t2.second) && ( (t1.first == t2.first && t1.third == t2.third) || (t1.first == t2.third && t1.third == t2.first));
  }
};

void OOPSEFormat::findAngles(OBMol& mol)
{
    /*
    //if already has data return
    if(mol.HasData(OBGenericDataType::AngleData))
        return;

    vector<OBEdgeBase*>::iterator bi1,bi2;
    OBBond* bond;
    OBAtom *a,*b,*c;

    set<triple<OBAtom*,OBAtom*,OBAtom*>, SameAngle> uniqueAngles;
    //loop through all bonds generating torsions
    for(bond = mol.BeginBond(bi1);bond;bond = mol.NextBond(bi1))
    {
        b = bond->GetBeginAtom();
        c = bond->GetEndAtom();
        if(b->IsHydrogen())
            continue;

        for(a = b->BeginNbrAtom(bi2);a;a = b->NextNbrAtom(bi2))
        {
            if(a == c)
                continue;          
            
            uniqueAngles.insert(triple<OBAtom*,OBAtom*,OBAtom*>(a, b, c));
        }
    }

    //get new data and attach it to molecule
    OBAngleData *angles = new OBAngleData;
    mol.SetData(angles);
    set<triple<OBAtom*,OBAtom*,OBAtom*>, SameAngle>::iterator i;

    for (i = uniqueAngles.begin(); i != uniqueAngles.end(); ++i) {
        OBAngle angle;
        angle.SetAtoms(i->first, i->second, i->second);
        angles->SetData(angle);
    }
    */
}
OBMol* OOPSEFormat::createMolFromFragment(OBMol& mol, vector<int>& fragment)
{
    OBMol* newMol = new OBMol();
    newMol->ReserveAtoms(fragment.size());
    newMol->BeginModify();
    for(vector<int>::iterator i = fragment.begin(); i != fragment.end(); ++i)
    {
        OBAtom* newAtom = newMol->NewAtom();
        *newAtom = *mol.GetAtom(*i);
    }
    newMol->EndModify();
    newMol->ConnectTheDots();
    findAngles(*newMol);
    newMol->FindTorsions();
    return newMol;
}
void OOPSEFormat::WriteMDFile(vector<OBMol*> mols, vector<int> numMols, ostream& os)
{
    std::string identLevel1("\t");
    std::string identLevel2("\t\t");
    std::string molPrefix("MolName");
    const int BUFFLEN = 1024;
    char buffer[BUFFLEN];
    
    for(unsigned int i = 0; i < mols.size(); ++i)
    {
        OBMol* pmol = mols[i];
        map<OBAtom*, int> atomMap;
        os << "molecule {\n";
        sprintf(buffer, "%d", i);
        os << identLevel1 << "name = " << "\"" << molPrefix << buffer << "\"" << ";\n";

        
        //atom
        int ai = 0;
        FOR_ATOMS_OF_MOL(atom, *pmol ) {
            os << identLevel1 << "atom[" << ai << "] {\n";
            os << identLevel2 << "type = " << "\"" << atom->GetType() << "\"" << ";\n";
            os << identLevel1 << "}\n";
            atomMap[&(*atom)] = ai++;
        }        
        os << "\n";

        //bond
        FOR_BONDS_OF_MOL(bond, *pmol ) {
            os << identLevel1 << "bond {\n";
            os << identLevel2 << "members(" << atomMap[bond->GetBeginAtom()] <<  ", " << atomMap[bond->GetEndAtom()] << ");\n";
            os << identLevel1 << "}\n";
        }  
        /*
        //bend
        OBGenericData* pGenericData = pmol->GetData(OBGenericDataType::AngleData);
        OBAngleData* pAngleData = dynamic_cast<OBAngleData*>(pGenericData);
        vector<OBAngle> angles = pAngleData->GetData();

        os << identLevel1 << "nBends = " << angles.size() << ";\n";        
        int bendIndex = 0;
        for (vector<OBAngle>::iterator ti = angles.begin(); ti != angles.end(); ++ti)
        {
            triple<OBAtom*, OBAtom*, OBAtom*> bendAtoms = ti->getAtoms();
            os << identLevel1 << "bend[" << bendIndex++ << "] {\n";
            os << identLevel2 << "member(" << atomMap[bendAtoms.first] <<  ", " << atomMap[bendAtoms.second] << atomMap[bendAtoms.third] <<");\n";
            os << identLevel1 << "}\n";            
        }
        
        //torsion
        pGenericData = pmol->GetData(OBGenericDataType::TorsionData);
        OBTorsionData* pTorsionData = dynamic_cast<OBTorsionData*>(pGenericData);
        vector<OBTorsion> torsions = pTorsionData->GetData();
        vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > torsionArray;
        for (vector<OBTorsion>::iterator ti = torsions.begin(); ti != torsions.end(); ++ti)
        {
            vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > tmpTorsions = ti->getTorsions();
            torsionArray.insert(torsionArray.end(), tmpTorsions.begin(), tmpTorsions.end());            
        }

        os << identLevel1 << "nTorsions = " << torsionArray.size() << ";\n";
        int torsionIndex = 0;
        for (vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> >::iterator ti = torsionArray.begin(); ti != torsionArray.end(); ++ti)
        {
            os << identLevel1 << "torsion[" << torsionIndex++ << "] {\n";
            os << identLevel2 << "member(" << atomMap[ti->first] <<  ", " << atomMap[ti->second] <<", " << atomMap[ti->third] <<", " << atomMap[ti->forth] << ");\n";
            os << identLevel1 << "}\n";          
        }
        */
        os << "}\n";
        os << "\n";

    }

    os << "\n";
    os << "nComponents = " << mols.size() << ";\n";
    
    for(unsigned int i =0; i < mols.size(); ++i)
    {
        os << "component{\n";
        sprintf(buffer, "%d", i);
        os << "type = " << molPrefix << buffer << ";\n";
        os << "nMol = " << numMols[i]<< ";\n";
        os << "}\n";
    }
}
void OOPSEFormat::WriteINFile(OBMol& mol, ostream& ofs, vector<int>& indices)
{
    unsigned int i;
    char buffer[BUFF_SIZE];
    
    sprintf(buffer,"%d", mol.NumAtoms());
    ofs << buffer << endl;
    //sprintf(buffer,"0;%f, 0, 0; 0, %15.7f, 0; 0, 0, %15.7f;", boxx, boxy, boxz);
    sprintf(buffer, "0;%f, 0, 0; 0, %15.7f, 0; 0, 0, %15.7f;", 100.0, 100.0, 100.0);
    ofs << buffer << endl;

    OBAtom *atom;
    string str,str1;
    
    for(vector<int>::iterator i = indices.begin();i != indices.end(); ++i)
    {
        atom = mol.GetAtom(*i);
        sprintf(buffer,"%15s%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f",
                atom->GetType(),
                atom->GetX(), atom->GetY(), atom->GetZ(),
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
                );
        ofs << buffer << endl;
    }
    
}

} //namespace OpenBabel

