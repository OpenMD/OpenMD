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

#include <protomol/debug/Pipe.h>

#include <protomol/base/Exception.h>
#include <protomol/base/Zap.h>

#include <stdio.h>
#include <unistd.h>
#include <ext/stdio_filebuf.h> // NOTE: This only works in GCC 3.2 and newer

using namespace std;
using namespace ProtoMol;

Pipe::Pipe() :
  outStream(0), inStream(0) {
  if (pipe(pipeFDs)) THROW("Error creating pipe!");
  fdOpen[0] = fdOpen[1] = true;
}

Pipe::~Pipe() {
  if (fdOpen[0]) closeOut();
  if (fdOpen[1]) closeIn();

  zap(outStream);
  zap(inStream);
}

void Pipe::closeOut() {
  ASSERT_OR_THROW("Pipe output not open!", fdOpen[0]);
  close(pipeFDs[0]);
  fdOpen[0] = false;
}

void Pipe::closeIn() {
  ASSERT_OR_THROW("Pipe input not open!", fdOpen[1]);
  close(pipeFDs[1]);
  fdOpen[1] = false;
}

void Pipe::duplicateOutFD(int newFD) {
  if (dup2(pipeFDs[0], newFD) != newFD)
    THROW("Error duplicating file descriptor!");
}

void Pipe::moveOutFD(int newFD) {
  duplicateOutFD(newFD);

  if (close(pipeFDs[0])) THROW("Error closing old file descriptor!");

  pipeFDs[0] = newFD;
}

void Pipe::duplicateInFD(int newFD) {
  if (dup2(pipeFDs[1], newFD) != newFD)
    THROW("Error duplicating file descriptor!");
}

void Pipe::moveInFD(int newFD) {
  duplicateInFD(newFD);

  if (close(pipeFDs[1])) THROW("Error closing old file descriptor!");

  pipeFDs[1] = newFD;
}

istream *Pipe::getOutStream() {
  ASSERT_OR_THROW("Pipe output not open!", fdOpen[0]);

  if (!outStream) {
    // C file descriptor to C++ streams magic
    // NOTE: This only works in GCC 3.2 and newer
    //       Hopefully they will leave the API alone now!
    __c_file *pipeFile = fdopen(pipeFDs[0], "r");
    __gnu_cxx::stdio_filebuf<char> *pipeBuf =
      new __gnu_cxx::stdio_filebuf<char>(pipeFile, ios::in);

    // FIXME There are two leaks above.  fdopen must be closed
    //       and the stdio_filebuf must be freed.
    //       Also in the getInStream()

    outStream = new istream(pipeBuf);
  }

  return outStream;
}

ostream *Pipe::getInStream() {
  ASSERT_OR_THROW("Pipe input not open!", fdOpen[1]);

  if (!inStream) {
    // C file descriptor to C++ streams magic
    // NOTE: This only works in GCC 3.2 and newer
    //       Hopefully they will leave the API alone now!
    __c_file *pipeFile = fdopen(pipeFDs[1], "w");
    __gnu_cxx::stdio_filebuf<char> *pipeBuf =
      new __gnu_cxx::stdio_filebuf<char>(pipeFile, ios::out);

    inStream = new ostream(pipeBuf);
  }

  return inStream;
}

