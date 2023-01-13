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

#ifndef IO_XYZFORMAT_HPP
#define IO_XYZFORMAT_HPP

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

#include <iostream>
#include <string>
#include <vector>

#include "math/Vector3.hpp"

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
    bool ReadMolecule(std::istream& ifs);
    std::vector<XYZAtom*> mol_;
    std::string title_;
  };
}  // namespace OpenMD

#endif
