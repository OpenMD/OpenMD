#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

namespace OpenMD {
  struct XYZAtom {
    Vector3d pos;
    std::string type;
    int atomicNum;
    RealType charge;
  };

  class XYZFormat {     
  public:
    XYZFormat() {}
    bool ReadMolecule(istream& ifs);
    std::vector <XYZAtom*> mol_;
    std::string title_;
  };
  
  bool XYZFormat::ReadMolecule(istream& ifs) {
    char buffer[BUFF_SIZE]; 
    unsigned int natoms;
    stringstream errorMsg;
    
    if (!ifs)
      return false; // we're attempting to read past the end of the file
    
    if (!ifs.getline(buffer,BUFF_SIZE)) {
      strcpy(painCave.errMsg,
	     "Problems reading an XYZ file: Cannot read the first line.\n");
      painCave.isFatal = 1;
      simError();
    }
    
    if (sscanf(buffer, "%d", &natoms) == 0 || !natoms) {
      strcpy(painCave.errMsg, 
	     "Problems reading an XYZ file: The first line must contain the number of atoms.\n");
      painCave.isFatal = 1;
      simError();
    }
    
    mol_.reserve(natoms);
    
    // The next line contains a title string for the molecule. Use this
    // as the title for the molecule if the line is not
    // empty. Otherwise, use the title given by the calling function.
    if (!ifs.getline(buffer,BUFF_SIZE)) {
      strcpy(painCave.errMsg, "Problems reading an XYZ file: Could not read the second line (title/comments).\n");
      painCave.isFatal = 1;
      simError();
    }
    string readTitle(buffer);
    string::size_type location = readTitle.find("Energy");
    if (location != string::npos)
      readTitle.erase(location);
    trim(readTitle);
    
    location = readTitle.find_first_not_of(" \t\n\r");
    if (location != string::npos)
      title_ = readTitle;
    else
      title_ = "";
    
    // The next lines contain four items each, separated by white
    // spaces: the atom type, and the coordinates of the atom
    for (unsigned int i = 1; i <= natoms; i ++) {
      if (!ifs.getline(buffer,BUFF_SIZE)) {
	errorMsg << "Problems reading an XYZ file: "
		 << "Could not read line #" << i+2 << ", file error." << endl
		 << " According to line one, there should be " << natoms
		 << " atoms, and therefore " << natoms+2 << " lines in the file.";
	
	strcpy(painCave.errMsg, errorMsg.str().c_str() );
	painCave.isFatal  = 0;
	painCave.severity = OPENMD_WARNING;
	simError();
	return(false);
      }
      StringTokenizer tokenizer(buffer, " ;,\t\n\r");
      std::vector<std::string> vs = tokenizer.getAllTokens();
      if (vs.size() < 4) // ignore extra columns which some applications add
	{
	  errorMsg << "Problems reading an XYZ file: "
		   << "Could not read line #" << i+2 << "." << endl
		   << "OpenBabel found the line '" << buffer << "'" << endl
		   << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
		   << "However, OpenBabel found " << vs.size() << " items.";
	
	  strcpy(painCave.errMsg, errorMsg.str().c_str() );
	  painCave.isFatal  = 0;
	  painCave.severity = OPENMD_WARNING;
	  simError();
	  return(false);
	}

      // Atom Type: get the atomic number from the element table, using
      // the first entry in the currently read line. If the entry makes
      // sense, set the atomic number and leave the atomic type open
      // (the type is then later faulted in when atom->GetType() is
      // called). If the entry does not make sense to use, set the atom
      // type manually, assuming that the author of the xyz-file had
      // something "special" in mind.
   
      XYZAtom* atom = new XYZAtom();
      mol_.push_back(atom);
    
      int atomicNum = etab.GetAtomicNum(vs[0].c_str());
      //set atomic number, or '0' if the atom type is not recognized
      if (atomicNum == 0) {
	// Sometimes people call this an XYZ file, but it's actually Unichem
	// i.e., the first column is the atomic number, not a symbol
	// so we'll first check if we can convert this to an element number
	atomicNum = atoi(vs[0].c_str());
      }
    
      atom->atomicNum = atomicNum;
      if (atomicNum == 0) // still strange, try using an atom type
	atom->type = vs[0];
    
      // Read the atom coordinates
      char *endptr;
      double x = strtod((char*)vs[1].c_str(),&endptr);
      if (endptr == (char*)vs[1].c_str()) {
      
	errorMsg << "Problems reading an XYZ file: "
		 << "Could not read line #" << i+2 << "." << endl
		 << "OpenBabel found the line '" << buffer << "'" << endl
		 << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
		 << "OpenBabel could not interpret item #1 as a number.";
      
	strcpy(painCave.errMsg, errorMsg.str().c_str());
	painCave.isFatal = 1;
	simError();
      }
      double y = strtod((char*)vs[2].c_str(),&endptr);
      if (endptr == (char*)vs[2].c_str()) {

	errorMsg << "Problems reading an XYZ file: "
		 << "Could not read line #" << i+2 << "." << endl
		 << "OpenBabel found the line '" << buffer << "'" << endl
		 << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
		 << "OpenBabel could not interpret item #2 as a number.";
      
	strcpy(painCave.errMsg, errorMsg.str().c_str());
	painCave.isFatal = 1;
	simError();
      }
      double z = strtod((char*)vs[3].c_str(),&endptr);
      if (endptr == (char*)vs[3].c_str()) {
      
	errorMsg << "Problems reading an XYZ file: "
		 << "Could not read line #" << i+2 << "." << endl
		 << "OpenBabel found the line '" << buffer << "'" << endl
		 << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
		 << "OpenBabel could not interpret item #3 as a number.";
      
	strcpy(painCave.errMsg, errorMsg.str().c_str() );
	painCave.isFatal  = 0;
	painCave.severity = OPENMD_WARNING;
	simError();
	return(false);
      }
      atom->pos = Vector3d(x,y,z); //set coordinates
    
      // OK, sometimes there's sym x y z charge -- accepted by Jmol
      if (vs.size() > 5) {
	string::size_type decimal = vs[4].find('.');
	if (decimal !=string::npos) { // period found
	  double charge = strtod((char*)vs[4].c_str(),&endptr);
	  if (endptr != (char*)vs[4].c_str())
	    atom->charge = charge;
	}
      } // attempt to parse charges
    }
  
    // clean out any remaining blank lines
    std::streampos ipos;
    do {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);  
    return(true);
  }
}
