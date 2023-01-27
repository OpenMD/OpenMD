/**********************************************************************

This basic Periodic Table class was originally taken from the data.cpp
file in OpenBabel. The code has been modified to match the OpenMD coding style.

We have retained the OpenBabel copyright and GPL license on this class:

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "io/XYZFormat.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "math/Vector3.hpp"
#include "utils/ElementsTable.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/Trim.hpp"
#include "utils/simError.h"

namespace OpenMD {

  bool XYZFormat::ReadMolecule(std::istream& ifs) {
    char buffer[BUFF_SIZE];
    unsigned int natoms;
    std::stringstream errorMsg;

    if (!ifs)
      return false;  // we're attempting to read past the end of the file

    if (!ifs.getline(buffer, BUFF_SIZE)) {
      strcpy(painCave.errMsg,
             "Problems reading an XYZ file: Cannot read the first line.\n");
      painCave.isFatal = 1;
      simError();
    }

    if (sscanf(buffer, "%d", &natoms) == 0 || !natoms) {
      strcpy(painCave.errMsg, "Problems reading an XYZ file: The first line "
                              "must contain the number of atoms.\n");
      painCave.isFatal = 1;
      simError();
    }

    mol_.reserve(natoms);

    // The next line contains a title string for the molecule. Use this
    // as the title for the molecule if the line is not
    // empty. Otherwise, use the title given by the calling function.
    if (!ifs.getline(buffer, BUFF_SIZE)) {
      strcpy(painCave.errMsg, "Problems reading an XYZ file: Could not read "
                              "the second line (title/comments).\n");
      painCave.isFatal = 1;
      simError();
    }
    std::string readTitle(buffer);
    std::string::size_type location = readTitle.find("Energy");
    if (location != std::string::npos) readTitle.erase(location);
    Utils::trim(readTitle);

    location = readTitle.find_first_not_of(" \t\n\r");
    if (location != std::string::npos)
      title_ = readTitle;
    else
      title_ = "";

    // The next lines contain four items each, separated by white
    // spaces: the atom type, and the coordinates of the atom
    for (unsigned int i = 1; i <= natoms; i++) {
      if (!ifs.getline(buffer, BUFF_SIZE)) {
        errorMsg << "Problems reading an XYZ file: "
                 << "Could not read line #" << i + 2 << ", file error.\n"
                 << " According to line one, there should be " << natoms
                 << " atoms, and therefore " << natoms + 2
                 << " lines in the file.";

        strcpy(painCave.errMsg, errorMsg.str().c_str());
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_WARNING;
        simError();
        return (false);
      }
      StringTokenizer tokenizer(buffer, " ;,\t\n\r");
      std::vector<std::string> vs = tokenizer.getAllTokens();
      if (vs.size() < 4)  // ignore extra columns which some applications add
      {
        errorMsg << "Problems reading an XYZ file: "
                 << "Could not read line #" << i + 2 << ".\n"
                 << "OpenBabel found the line '" << buffer << "'\n"
                 << "According to the specifications, this line should contain "
                    "exactly 4 entries, separated by white space.\n"
                 << "However, OpenBabel found " << vs.size() << " items.";

        strcpy(painCave.errMsg, errorMsg.str().c_str());
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_WARNING;
        simError();
        return (false);
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
      // set atomic number, or '0' if the atom type is not recognized
      if (atomicNum == 0) {
        // Sometimes people call this an XYZ file, but it's actually Unichem
        // i.e., the first column is the atomic number, not a symbol
        // so we'll first check if we can convert this to an element number
        atomicNum = atoi(vs[0].c_str());
      }

      atom->atomicNum = atomicNum;
      if (atomicNum == 0)  // still strange, try using an atom type
        atom->type = vs[0];

      // Read the atom coordinates
      char* endptr;
      double x = strtod((char*)vs[1].c_str(), &endptr);
      if (endptr == (char*)vs[1].c_str()) {
        errorMsg << "Problems reading an XYZ file: "
                 << "Could not read line #" << i + 2 << ".\n"
                 << "OpenBabel found the line '" << buffer << "'\n"
                 << "According to the specifications, this line should contain "
                    "exactly 4 entries, separated by white space.\n"
                 << "OpenBabel could not interpret item #1 as a number.";

        strcpy(painCave.errMsg, errorMsg.str().c_str());
        painCave.isFatal = 1;
        simError();
      }
      double y = strtod((char*)vs[2].c_str(), &endptr);
      if (endptr == (char*)vs[2].c_str()) {
        errorMsg << "Problems reading an XYZ file: "
                 << "Could not read line #" << i + 2 << ".\n"
                 << "OpenBabel found the line '" << buffer << "'\n"
                 << "According to the specifications, this line should contain "
                    "exactly 4 entries, separated by white space.\n"
                 << "OpenBabel could not interpret item #2 as a number.";

        strcpy(painCave.errMsg, errorMsg.str().c_str());
        painCave.isFatal = 1;
        simError();
      }
      double z = strtod((char*)vs[3].c_str(), &endptr);
      if (endptr == (char*)vs[3].c_str()) {
        errorMsg << "Problems reading an XYZ file: "
                 << "Could not read line #" << i + 2 << ".\n"
                 << "OpenBabel found the line '" << buffer << "'\n"
                 << "According to the specifications, this line should contain "
                    "exactly 4 entries, separated by white space.\n"
                 << "OpenBabel could not interpret item #3 as a number.";

        strcpy(painCave.errMsg, errorMsg.str().c_str());
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_WARNING;
        simError();
        return (false);
      }
      atom->pos = Vector3d(x, y, z);  // set coordinates

      // OK, sometimes there's sym x y z charge -- accepted by Jmol
      if (vs.size() > 5) {
        std::string::size_type decimal = vs[4].find('.');
        if (decimal != std::string::npos) {  // period found
          double charge = strtod((char*)vs[4].c_str(), &endptr);
          if (endptr != (char*)vs[4].c_str()) atom->charge = charge;
        }
      }  // attempt to parse charges
    }

    // clean out any remaining blank lines
    std::streampos ipos;
    do {
      ipos = ifs.tellg();
      ifs.getline(buffer, BUFF_SIZE);
    } while (strlen(buffer) == 0 && !ifs.eof());
    ifs.seekg(ipos);
    return (true);
  }
}  // namespace OpenMD
