/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2008 by J. Daniel Gezelter
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/obiter.h>
#include <openbabel/mol.h>
#include <openbabel/chains.h>
#include <fstream>

#include "utils/StringUtils.hpp"

using namespace std;
namespace OpenBabel
{
  
  class OOPSEFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    OOPSEFormat() 
    {      
      OBConversion::RegisterFormat("md",this);
    }
    
    virtual const char* Description() //required
    {
      return
        "OOPSE combined meta-data / cartesian coordinates format\nNo comments yet\n";
    };
    
    virtual const char* SpecificationURL()
    {return "http://www.oopse.org";}; //optional
    
    virtual const char* GetMIMEType() 
    {return "chemical/x-md"; };
    
    virtual unsigned int Flags() 
    { 
      return NOTREADABLE | WRITEONEONLY;
    }
    
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    bool AreSameFragments(OBMol& mol, vector<int>& frag1, vector<int>& frag2);
    OBMol* createMolFromFragment(OBMol& mol, vector<int>& fragment);
    void WriteMDFile(vector<OBMol*> mols, vector<int> numMols, ostream& os, 
                     OBMol& mol, vector<int>& indices);
    void CalcBoundingBox(OBMol &mol,
                         double &min_x, double &max_x,
                         double &min_y, double &max_y,
                         double &min_z, double &max_z);
      
  };
  
  //Make an instance of the format class
  OOPSEFormat theOOPSEFormat;
 
  bool OOPSEFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv) {
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
        continue;

      used[i] = true;
      vector<int> sameMolTypes;
      sameMolTypes.push_back(i);
      indices.insert(indices.end(), fragmentLists[i].begin(), 
                     fragmentLists[i].end());
      for (int j = i + 1;j < used.size(); ++j) {
        if (used[j])
          continue;
                
        if (AreSameFragments(*pmol, fragmentLists[i], fragmentLists[j])) {
          sameMolTypes.push_back(j);
          indices.insert(indices.end(), fragmentLists[j].begin(), 
                         fragmentLists[j].end());
          used[j]=true;
        }
      }
      molecules.push_back(sameMolTypes);      
    }
    
    vector<OBMol*> mdMols;    
    vector<int> numMols;
    for(vector<vector<int> >::iterator  i = molecules.begin(); 
        i != molecules.end(); ++i) {
      
      mdMols.push_back(createMolFromFragment(*pmol, 
                                             fragmentLists[i->front()]));
      numMols.push_back((*i).size());
    }
    
    string OutputFileName = pConv->GetInFilename();
    size_t pos = OutputFileName.rfind(".");
    if(pos!=string::npos)
      OutputFileName = OutputFileName.substr(0, pos) + ".md";       
    else
      OutputFileName += ".md";
    
    ofstream ofs(OutputFileName.c_str());
    if(!ofs) {
        cerr << "Cannot write to " << OutputFileName <<endl;
        return false;
    }
    
    WriteMDFile(mdMols, numMols, ofs, *pmol, indices);
    
    for(vector<OBMol*>::iterator  i = mdMols.begin(); i != mdMols.end(); ++i) {
      delete *i;
    }
    
    return(true);
  }
  
  bool OOPSEFormat::AreSameFragments(OBMol& mol, vector<int>& frag1, 
                                     vector<int>& frag2) {
    if (frag1.size() != frag2.size())
      return false;
    
    // Exact graph matching is an NP complete problem.
    // This just matches all of the atom atomic numbers and may falsely 
    // detect identical fragments which aren't really identical.
    // @todo using sparse matrix to store the connectivities 

    for (unsigned int i =0 ; i < frag1.size(); ++i) {
        OBAtom* atom1 = mol.GetAtom(frag1[i]);
        OBAtom* atom2 = mol.GetAtom(frag2[i]);
        
        if (atom1->GetAtomicNum() != atom2->GetAtomicNum()) 
          return false;
        
    }
    return true;
  }
  
  struct SameAngle {
    bool operator()(const triple<OBAtom*,OBAtom*,OBAtom*> t1, 
                    const triple<OBAtom*,OBAtom*,OBAtom*> t2) const {
      return (t1.second == t2.second) && ( (t1.first == t2.first && t1.third == t2.third) || (t1.first == t2.third && t1.third == t2.first));
    }
  };


  OBMol* OOPSEFormat::createMolFromFragment(OBMol& mol, 
                                            vector<int>& fragment) {
    
    OBMol* newMol = new OBMol();
    newMol->ReserveAtoms(fragment.size());
    newMol->BeginModify();
    for(vector<int>::iterator i = fragment.begin(); i != fragment.end(); ++i) {
      OBAtom* newAtom = newMol->NewAtom();
      *newAtom = *mol.GetAtom(*i);
    }

    newMol->EndModify();
    newMol->ConnectTheDots();
    newMol->PerceiveBondOrders();

    return newMol;
  }
  
  void OOPSEFormat::WriteMDFile(vector<OBMol*> mols, vector<int> numMols, 
                                ostream& os, OBMol& mol, 
                                vector<int>& indices) {
    
    std::string molPrefix("MolName");
    std::string resName;
    unsigned int i;
    const int BUFFLEN = 1024;
    char buffer[BUFFLEN];
    string str, str1;
    OBAtom *a, *b, *c, *d;    
    bool molIsWater = false;
    OBResidue *r;
    int resKey;
    char type_name[10];
    char *element_name;
    int res_num;
    OBChainsParser* chainParser = new OBChainsParser();   
    double min_x, max_x, min_y, max_y, min_z, max_z; /* Edges of bounding box */
    
    os << "<OOPSE version=4>" << endl;
    os << "  <MetaData>" << endl << endl;
    
    for(i = 0; i < mols.size(); ++i) {
      OBMol* pmol = mols[i];
      map<OBAtom*, int> atomMap;

      chainParser->PerceiveChains(*pmol, false);
      molIsWater = false;
      FOR_RESIDUES_OF_MOL(residue, *pmol) {
        std::cerr << "residue = " << residue->GetName() << "\n";
        if (residue->GetName().compare("HOH") == 0) {
          molIsWater = true;
        }
      }
      
      if (molIsWater) {
        // water include files define all of the known water types
        os << "#include \"water.md\";\n";
        pmol->SetTitle("HOH");
      } else {

        os << "molecule {\n";
        sprintf(buffer, "%d", i);
        os << "  name = \"" << molPrefix << buffer << "\";\n";
        
        int ai = 0;
        FOR_ATOMS_OF_MOL(atom, *pmol ) {
          str = atom->GetType();
          r = atom->GetResidue();
          
          if (r == NULL) 
            resName = "NULL";
          else 
            resName = r->GetName();
          
          if (resName.compare("NULL") ==0 || 
              resName.compare("LIG") == 0 || 
              resName.compare("UNK") == 0) {
            // Either couldn't find a residue at all or couldn't find a
            // reasonable residue name to use.  We'll punt and use
            // OpenBabel's internal atom typing:
            ttab.SetFromType("INT");
            ttab.SetToType("INT");
            ttab.Translate(str1, str);
          } else {                        
            
            // If we know what residue we've got, the specific atom name can
            // be used to help specify partial charges. 
            
            str = r->GetAtomID(&*atom);
            size_t startpos = str.find_first_not_of(" ");
            size_t endpos = str.find_last_not_of(" ");
            str = str.substr( startpos, endpos-startpos+1 );
            str1 = resName + "-" + str;
          }       
          os << "  atom[" << ai << "] { ";
          os << "type = " << "\"" << str1 << "\"" << "; ";
          os << "}\n";
          atomMap[&(*atom)] = ai++;
        }        
        os << "\n";
        
        //bond
        
        int b1, b2;
        FOR_BONDS_OF_MOL(bond, *pmol ) {
          b1 = atomMap[bond->GetBeginAtom()];
          b2 = atomMap[bond->GetEndAtom()];
          
          os << "  bond { ";
          
          if (b1 < b2) 
            os << "members(" << b1 <<  ", " << b2 << "); ";
          else
            os << "members(" << b2 <<  ", " << b1 << "); ";
          
          os << "}" << endl;
        } 
        
        os << endl;
       
        os << "}" << endl;
        os << endl;
      }
    }
    
    os << endl;
        
    for(i=0; i < mols.size(); ++i) {
      OBMol* pmol = mols[i];      
      os << "component{" << endl;
      if (std::string(pmol->GetTitle()).compare("HOH") == 0) {
        os << "  type = " << "HOH" << ";" << endl;
      } else {
        sprintf(buffer, "%d", i);
        os << "  type = " << molPrefix << buffer << ";" << endl;
      }
      os << "  nMol = " << numMols[i]<< ";" << endl;
      os << "}" << endl;
    }
    
    os << "  </MetaData>" << endl;
    os << "  <Snapshot>" << endl;
    os << "    <FrameData>" << endl;
    
    sprintf(buffer, "        Time: %.10g", 0.0);
    
    os << buffer << endl;

    CalcBoundingBox(mol, min_x, max_x, min_y, max_y, min_z, max_z);

    // still to do: should compute a bounding box here
    sprintf(buffer, "        Hmat: {{ %.10g, %.10g, %.10g }, { %.10g, %.10g, %.10g }, { %.10g, %.10g, %.10g }}", 
            max_x - min_x, 0.0, 0.0, 0.0, max_y - min_y, 0.0, 0.0, 0.0, max_z - min_z);
    
    os << buffer << endl;
    os << "    </FrameData>" << endl;
    os << "    <StuntDoubles>" << endl;
    
    OBAtom *atom;
    
    // still to do: intercept waters and recompute pvqj lines

    for(vector<int>::iterator i = indices.begin();i != indices.end(); ++i) {     
      atom = mol.GetAtom(*i);
      sprintf(buffer, "%10d %7s %18.10g %18.10g %18.10g %13e %13e %13e", *i - 1, 
              "pv", atom->GetX(), atom->GetY(), atom->GetZ(), 0.0, 0.0, 0.0);
      os << buffer << endl;
    }
    os << "    </StuntDoubles>" << endl;
    os << "  </Snapshot>" << endl;
    os << "</OOPSE>" << endl;
  }

  void OOPSEFormat::CalcBoundingBox(OBMol &mol,
                                    double &min_x, double &max_x,
                                    double &min_y, double &max_y,
                                    double &min_z, double &max_z
                                    )
  {
    /* ---- Init bounding-box variables ---- */
    min_x = (double) 0.0;
    max_x = (double) 0.0;
    min_y = (double) 0.0;
    max_y = (double) 0.0;
    min_z = (double) 0.0;
    max_z = (double) 0.0;
    
    /* ---- Check all atoms ---- */
    for(unsigned int i = 1; i <= mol.NumAtoms(); ++i)
      {
        
        /* ---- Get a pointer to ith atom ---- */
        OBAtom *atom = mol.GetAtom(i);
        
        /* ---- Check for minimal/maximal x-position ---- */
        if (atom -> GetX() < min_x)
          min_x = atom -> GetX();
        if (atom -> GetX() > max_x)
          max_x = atom -> GetX();
        
        /* ---- Check for minimal/maximal y-position ---- */
        if (atom -> GetY() < min_y)
          min_y = atom -> GetY();
        if (atom -> GetY() > max_y)
          max_y = atom -> GetY();
        
        /* ---- Check for minimal/maximal z-position ---- */
        if (atom -> GetZ() < min_z)
          min_z = atom -> GetZ();
        if (atom -> GetZ() > max_z)
          max_z = atom -> GetZ();
        
      }    
  }
} //namespace OpenBabel

