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

#include <protomol/base/Exception.h>

#ifdef HAVE_DEBUGGER
#include <protomol/debug/Debugger.h>
#include <protomol/debug/Process.h>
#include <protomol/type/String.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>

#include <sstream>

using namespace std;
using namespace ProtoMol;

string Debugger::executableName;
int Debugger::numTraces = 0;
bool Debugger::traceFiltering = true;
int Debugger::maxTraces = 10;

void Debugger::initStackTrace(string executableName) {
  Debugger::executableName = executableName;
  Exception::enableStackTraces = true;
}

bool Debugger::printStackTrace(ostream &stream) {
  trace_t trace;
  bool retVal = getStackTrace(trace);

  trace_t::iterator it;
  for (it = trace.begin(); it != trace.end(); it++)
    stream << *it << endl;

  return retVal;
}

#define BUF_SIZE 2048

bool Debugger::_getStackTrace(trace_t &trace) {
  if (executableName == "") {
    trace.push_back("Stack Trace Error: Stack dumper not initialized!");
    return false;
  }

  numTraces++;
  if (maxTraces && numTraces > maxTraces) {
    trace.push_back("Stack Trace Error: Exceeded maxTraces of " +
      String(maxTraces));
    return false;
  }

  // Spawn gdb process
  int argc = 0;
  const char *argv[5];

  argv[argc++] = "gdb";
  argv[argc++] = (char *)executableName.c_str();
  String pid(getpid());
  argv[argc++] = (char *)pid.c_str();
  argv[argc] = 0;

  try {
    Process debugProc;
    Pipe *inPipe = debugProc.getChildPipe(Process::TO_CHILD, 0);
    Pipe *outPipe = debugProc.getChildPipe(Process::FROM_CHILD, 1);
    Pipe *errPipe = debugProc.getChildPipe(Process::FROM_CHILD, 2);

    // Run gdb commands
    string debugCmd =
      string("set width ") + String(BUF_SIZE - 1) + "\nwhere\nquit\n";
    ssize_t size =
      write(inPipe->getInFD(), debugCmd.c_str(), debugCmd.length());
    size = size;

    // Execute debugger process
    debugProc.exec((char **)argv);

    // Read output
    FILE *out = fdopen(outPipe->getOutFD(), "r");
    FILE *err = fdopen(errPipe->getOutFD(), "r");
    if (!out || !err) {
      trace.push_back("Stack Trace Error: Opening debugger output streams!");

      return false;
    }

    char buf[BUF_SIZE + 1];
    int offset = 0;
    int count = 0;
    while (fgets(buf, BUF_SIZE, out))
      if (buf[0] == '#') {
        if (traceFiltering) {
          count++;

          if (strstr(buf, "Debugger::") ||
              strstr(buf, "Exception::init") ||
              strstr(buf, "Exception (")) {
            offset = count;
            trace.clear();
            continue;
          }
        }

        int line = atoi(&buf[1]) - offset;
        char *start = strchr(buf, ' ');
        int len = strlen(buf);

        if (buf[len - 1] == '\n' || buf[len - 1] == '\r') buf[len - 1] = 0;
        trace.push_back(string("#") + String(line) + start);
      }

#ifdef DEBUGGER_PRINT_ERROR_STREAM
    while (fgets(buf, BUF_SIZE, err)) {
      int len = strlen(buf);
      if (buf[len - 1] == '\n' || buf[len - 1] == '\r') buf[len - 1] = 0;
      if (buf[0] != 0) trace.push_back(buf);
    }
#endif


    // Clean up
    fclose(out);
    fclose(err);

    debugProc.wait();
    if (debugProc.getReturnCode()) {
      trace.push_back("Stack Trace Error: gdb returned an error.");

      return false;
    }

    return true;
  } catch (const Exception &e) {
    trace.push_back(string("Stack Trace Error: ") + e.getMessage());
  }

  return false;
}

bool Debugger::getStackTrace(trace_t &trace) {
  static bool inStackTrace = false;
  bool ret;

  if (inStackTrace) {
    trace.push_back("Stack Trace Error: Already in stack trace!");
    return false;
  }

  inStackTrace = true;

  try {
    ret = _getStackTrace(trace);
  } catch (...) {
    inStackTrace = false;
    throw;
  }
  inStackTrace = false;

  return ret;
}
#endif // HAVE_DEBUGGER

#ifdef DEBUGGER_TEST
#include <iostream>

void b(int x) {
  THROW("Test cause!");
}

void a(char *x) {
  try {
    b(10);
  } catch (const Exception &e) {
    THROWC("Test exception!", e);
  }
}

int main(int argc, char *argv[]) {
  Debugger::initStackTrace(argv[0]);

  try {
    a("test");
  } catch (const Exception &e) {
    cerr << "Exception: " << e << endl;
  }

  return 0;
}

#endif // DEBUGGER_TEST
