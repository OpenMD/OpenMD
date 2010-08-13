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
#ifndef PROCESS_H
#define PROCESS_H

#include <iostream>
#include <string>
#include <vector>
#include <list>

#include <protomol/debug/Pipe.h>

namespace ProtoMol {

  class ProcessFunctor {
  public:
	virtual ~ProcessFunctor() {}
    virtual void child() = 0;
    virtual void parent() = 0;
  };

  class PipeProcessFunctor : public ProcessFunctor {
    Pipe pipe;
    int direction;
    int fd;
  public:

    PipeProcessFunctor(int direction, int fd) : direction(direction), fd(fd) {}
    virtual ~PipeProcessFunctor() {}

    virtual void child();
    virtual void parent();

    Pipe *getPipe() {return &pipe;}
  };

  class FDReplaceProcessFunctor : public ProcessFunctor {
    int fd;
    int replacement;

  public:
    FDReplaceProcessFunctor(int fd, int replacement) :
      fd(fd), replacement(replacement) {}
    virtual ~FDReplaceProcessFunctor() {}

    virtual void child();
    virtual void parent() {}
  };

  class Process {
    typedef std::vector<ProcessFunctor *> functors_t;
    functors_t functors;

    int pid;
    bool running;
    int returnCode;

  public:
    Process();
    ~Process();

    typedef enum {TO_CHILD, FROM_CHILD} dir_t;

    void exec(std::list<std::string> &args);
    void exec(const char *args);
    void exec(char *args);
    void exec(char *argv[]);

    static void parseArgs(char *args, int &argc, char *argv[], int n);

    Pipe *getChildPipe(dir_t dir, int childFD = -1);
    void replaceChildFD(int fd, int replacement);

    int getPID() {return pid;}
    void kill(int sig);
    int wait(int options = 0);

    bool isRunning();
    int getReturnCode() {return returnCode;}
  };
}
#endif // PROCESS_H
