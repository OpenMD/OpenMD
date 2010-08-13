/*******************************************************************\

              Copyright (C) 2004 Joseph Coffland

    This program is free software; you can redistribute it and/or
     modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
        of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
             GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
                           02111-1307, USA.

            For information regarding this software email
                   jcofflan@users.sourceforge.net

\*******************************************************************/
#ifndef PIPE_H
#define PIPE_H

#include <iostream>

namespace ProtoMol {

  class Pipe {
    int pipeFDs[2];
    bool fdOpen[2];
    std::istream *outStream;
    std::ostream *inStream;

  public:
    Pipe();
    ~Pipe();

    void closeOut();
    void closeIn();

    int getOutFD() {return pipeFDs[0];}
    int getInFD() {return pipeFDs[1];}

    void duplicateOutFD(int newFD);
    void duplicateInFD(int newFD);

    void moveOutFD(int newFD);
    void moveInFD(int newFD);

    std::istream *getOutStream();
    std::ostream *getInStream();
  };
}
#endif // PIPE_H
