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
#include <protomol/debug/Process.h>

#include <protomol/base/Exception.h>
#include <protomol/debug/Debugger.h>
#include <protomol/debug/Pipe.h>
#include <protomol/type/String.h>
#include <protomol/base/Zap.h>

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>
#include <string.h>
#include <errno.h>
#include <string>

using namespace std;
using namespace ProtoMol;

void PipeProcessFunctor::child() {
  if (direction == Process::TO_CHILD) {
    pipe.closeIn();
    if (fd >= 0) pipe.moveOutFD(fd);
  } else {
    pipe.closeOut();
    if (fd >= 0) pipe.moveInFD(fd);
  }
}

void PipeProcessFunctor::parent() {
  if (direction == Process::TO_CHILD) pipe.closeOut();
  else pipe.closeIn();
}

void FDReplaceProcessFunctor::child() {
  if (dup2(replacement, fd) != fd)
    THROW("Error replacing file descriptor in child process!");
}

Process::Process() :
  pid(0), running(false), returnCode(0) {}

Process::~Process() {
  if (isRunning()) {
    kill(SIGHUP);
    wait();
  }

  for (unsigned int i = 0; i < functors.size(); i++)
    delete functors[i];
}

void Process::exec(list<string> &args) {
  SmartPointer<char *>::Array argv = new char *[args.size() + 1];
  list<string>::iterator it;
  unsigned int i = 0;

  for (i = 0, it = args.begin(); it != args.end(); i++, it++)
    argv[i] = (char *)it->c_str();

  argv[i] = 0;

  exec(argv.get());
}

void Process::exec(const char *args) {
  SmartPointer<char, SP_MALLOC> tmp = strdup(args);
  exec(tmp.get());
}

void Process::exec(char *args) {
  char *argv[256];
  int argc = 0;

  parseArgs(args, argc, argv, 256);
  exec(argv);
}

void Process::exec(char *argv[]) {
  ASSERT_OR_THROW("Process all ready running!", !isRunning());

  pid = fork();

  if (pid == 0) {
    try {
      for (unsigned int i = 0; i < functors.size(); i++)
        functors[i]->child();
    } catch (const Exception &e) {
      perror(e.getMessage().c_str());
    }

    execvp(argv[0], argv);

    // Execution failed!
    string errorStr = "Process() executing '";
    for (unsigned i = 0; argv[i]; i++) {
      if (i) errorStr += " ";
      errorStr += argv[i];
    }

    errorStr += "'";

    perror(errorStr.c_str());
    exit(1);
  } else if (pid == -1)
    THROW("Failed to spawn child!");

  running = true;

  for (unsigned int i = 0; i < functors.size(); i++)
    functors[i]->parent();
}

void Process::parseArgs(char *args, int &argc, char *argv[], int n) {
  if (args) {
    bool inArg = false;
    bool inSQuote = false;
    bool inDQuote = false;
    bool addChar;
    int i = 0;

    for (char *ptr = args; *ptr; ptr++) {
      addChar = false;

      switch (*ptr) {
      case '\\':
        if (ptr[1] != '\0') {
          ptr++;
          switch (*ptr) {
          case 'n': *ptr = '\n';
            break;
          case 't': *ptr = '\t';
            break;
          case 'r': *ptr = '\r';
            break;
          default: break;
          }

          addChar = true;
        }
        break;

      case '\'':
        if (inDQuote) addChar = true;
        else
        if (inSQuote) inSQuote = false;
        else inSQuote = true;
        break;

      case '"':
        if (inSQuote) addChar = true;
        else
        if (inDQuote) inDQuote = false;
        else inDQuote = true;
        break;

      case '\t':
      case ' ':
      case '\n':
      case '\r':
        if (inArg) {
          if (inSQuote || inDQuote) addChar = true;
          else {
            args[i++] = 0;
            inArg = false;
          }
        }
        break;

      default: addChar = true;
        break;
      }

      if (addChar) {
        if (!inArg) {
          ASSERT_OR_THROW("Too many arguments!", argc < n - 1);
          argv[argc++] = &args[i];
          inArg = true;
        }
        args[i++] = *ptr;
      }
    }

    if (inArg) args[i] = 0;
  }

  argv[argc] = 0;
}

Pipe *Process::getChildPipe(dir_t dir, int childFD) {
  ASSERT_OR_THROW("Process all ready running!", !running);

  PipeProcessFunctor *functor = new PipeProcessFunctor(dir, childFD);
  functors.push_back(functor);

  return functor->getPipe();
}

void Process::replaceChildFD(int fd, int replacement) {
  functors.push_back(new FDReplaceProcessFunctor(fd, replacement));
}

void Process::kill(int sig) {
  ASSERT_OR_THROW("Process not running!", running);

  if (::kill(pid, sig) != 0)
    THROW(string("Failed to kill process ") + String(pid) + ":" +
      strerror(errno));
}

int Process::wait(int options) {
  ASSERT_OR_THROW("Process not running!", running);

  int status = 0;
  int retVal = waitpid(pid, &status, options);

  if (retVal == -1) {
    running = false;
    THROW(string("Failed to wait on process ") + String(pid) + ":" +
      strerror(errno));
  }

  if (retVal) {
    if (WIFEXITED(status)) {
      returnCode = WEXITSTATUS(status);
      running = false;
    } else if (WIFSIGNALED(status))
      // TODO report on signal
      running = false;
  }

  return status;
}

bool Process::isRunning() {
  if (running) wait(WNOHANG);
  return running;
}

