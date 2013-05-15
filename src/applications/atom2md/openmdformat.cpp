/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2008-2009 by J. Daniel Gezelter
 
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
#include <openbabel/data.h>
#include <fstream>

#include "utils/StringUtils.hpp"

using namespace std;
namespace OpenBabel
{
  
  class OpenMDFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    OpenMDFormat() 
    {      
      OBConversion::RegisterFormat("md",this);
    }
    
    virtual const char* Description() //required
    {
      return
        "OpenMD combined meta-data / cartesian coordinates format\n\
        No comments yet\n";
    };
    
    virtual const char* SpecificationURL()
    {return "http://www.openmd.net";}; //optional
    
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
  OpenMDFormat theOpenMDFormat;
 
  bool OpenMDFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv) {
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
  
  bool OpenMDFormat::AreSameFragments(OBMol& mol, vector<int>& frag1, 
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


  OBMol* OpenMDFormat::createMolFromFragment(OBMol& mol, 
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
  
  void OpenMDFormat::WriteMDFile(vector<OBMol*> mols, vector<int> numMols, 
                                ostream& os, OBMol& mol, 
                                vector<int>& indices) {
    
    std::string molPrefix("MolName");
    std::string resName;
    unsigned int i;
    const int BUFFLEN = 1024;
    char buffer[BUFFLEN];
    string str, str1, str2, str3;
    bool molIsWater = false;
    OBResidue *r;
    double min_x, max_x, min_y, max_y, min_z, max_z; /* Edges of bounding box */
    
    os << "<OpenMD version=2>" << endl;
    os << "  <MetaData>" << endl << endl;
    
    for(i = 0; i < mols.size(); ++i) {
      OBMol* pmol = mols[i];
      map<OBAtom*, int> atomMap;

      molIsWater = false;
      FOR_RESIDUES_OF_MOL(residue, *pmol) {
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
        sprintf(buffer, "%u", i);
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

            //resdat.SetResName(resName);
            
            // atom type from residue: 
            str = r->GetAtomID(&*atom);
           
	    // arginine has separate indices for chemically-identical
	    // nitrogen atoms:
	    if (resName.compare("ARG") == 0) {
	      if (str.compare("NH1") == 0 || str.compare("NH2") == 0) {
		str = "NH";
	      }
	    }
	    if (resName.compare("VAL") == 0) {
	      if (str.compare("CG1") == 0 || str.compare("CG2") == 0) {
		str = "CG";
	      }
	    }
	    if (resName.compare("LEU") == 0) {
	      if (str.compare("CD1") == 0 || str.compare("CD2") == 0) {
		str = "CD";
	      }
	    }
	    if (resName.compare("ASP") == 0) {
	      if (str.compare("OD1") == 0 || str.compare("OD2") == 0) {
		str = "OD";
	      }
	    }
	    if (resName.compare("GLU") == 0) {
	      if (str.compare("OE1") == 0 || str.compare("OE2") == 0) {
		str = "OE";
	      }
	    }
	    if (resName.compare("TYR") == 0) {
	      if (str.compare("CD1") == 0 || str.compare("CD2") == 0) {
		str = "CD";
	      }
	      if (str.compare("CE1") == 0 || str.compare("CE2") == 0) {
		str = "CE";
	      }
	    }
	    

            if ((&*atom)->IsHydrogen()) {
               FOR_NBORS_OF_ATOM(nbr, *atom) {
                 str2 = r->GetAtomID(&*nbr);
                 size_t startpos = str2.find_first_not_of(" ");
                 size_t endpos = str2.find_last_not_of(" ");
                 if ((endpos - startpos) < 1) {
                   // if the bonded atom type has only one character (i.e. N)
                   // then the hydrogen will be labeled "HN" to show what
                   // kind of proton it is:
                   str3 = str2;
                 } else {
                   if (str2.compare("OH") == 0) {
                      str3 = "O";
                   } else {
                     // When the bonded atom type is more specific, we drop
                     // the first character:  i.e. H bonded to OG1 is HG1 type:
                     str3 = str2.substr(startpos+1, endpos-startpos);
                   }
                 }
                str = "H" + str3;
               }
	       // same problem with arginine NH atoms, but now for connected hydrogens
	       if (resName.compare("ARG") == 0) {
		 if (str.compare("HH1") == 0 || str.compare("HH2") == 0) {
		   str = "HH";
		 }
	       }
	       if (resName.compare("VAL") == 0) {
		 if (str.compare("HG1") == 0 || str.compare("HG2") == 0) {
		   str = "HG";
		 }
	       }
	       if (resName.compare("LEU") == 0) {
		 if (str.compare("HD1") == 0 || str.compare("HD2") == 0) {
		   str = "HD";
		 }
	       }
	       if (resName.compare("TYR") == 0) {
		 if (str.compare("HD1") == 0 || str.compare("HD2") == 0) {
		   str = "HD";
		 }
		 if (str.compare("HE1") == 0 || str.compare("HE2") == 0) {
		   str = "HE";
		 }
	       }

	    }

            // atom type from residue table:
            //resdat.LookupType(str, str2, hyb);
            size_t startpos = str.find_first_not_of(" ");
            size_t endpos = str.find_last_not_of(" ");
            str = str.substr( startpos, endpos-startpos+1 );
            str1 = resName + "-" + str;
          }      
          os << "  atom[" << ai << "] { ";
          os << "type = " << "\"" << str1 << "\"" << "; ";
          os << "position( " << (&*atom)->GetX() << ", " << (&*atom)->GetY() << ", " << (&*atom)->GetZ() << ");";
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
        os << "  type = " << "\"HOH\"" << "; // change to appropriate water model" << endl;
      } else {
        sprintf(buffer, "%u", i);
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
         
    for(vector<int>::iterator i = indices.begin();i != indices.end(); ++i) {     
      
      atom = mol.GetAtom(*i);
      sprintf(buffer, "%10d %7s %18.10g %18.10g %18.10g %13e %13e %13e", *i - 1, 
              "pv", atom->GetX(), atom->GetY(), atom->GetZ(), 0.0, 0.0, 0.0);
      os << buffer << endl;
    }
    os << "    </StuntDoubles>" << endl;
    os << "  </Snapshot>" << endl;
    os << "</OpenMD>" << endl;
  }

  void OpenMDFormat::CalcBoundingBox(OBMol &mol,
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

