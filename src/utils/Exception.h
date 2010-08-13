/*******************************************************************\

              Copyright (C) 2003 Joseph Coffland

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


#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <protomol/config.h>
#include <protomol/base/FileLocation.h>
#include <protomol/base/SmartPointer.h>
#include <protomol/base/Zap.h>

#ifdef HAVE_DEBUGGER
#include <protomol/debug/Debugger.h>
#endif
#ifdef HAVE_BACKTRACE
#include <protomol/debug/Backtrace.h>
#endif

#include <vector>
#include <string>
#include <iostream>

#ifdef HAVE_NO_SSTREAM
#include <protomol/base/sstream_local.h>
#else
#include <sstream>
#endif

#if defined(HAVE_DEBUGGER) || defined(HAVE_BACKTRACE)
#define HAVE_STACK_TRACE
#endif


namespace ProtoMol {
  // Forward Declarations
  template <class T, class DEALLOC_T>
  class SmartPointer;

  /**
   * Exception is a general purpose exception class.  It is similar to
   * the java Exception class.  A Exception can carry a message,
   * a FileLocation and/or a pointer to an exception which was the
   * original cause.
   *
   * There are preprocessor macros that can be used to as a convient way
   * to add the current file, line and column where the exception occured.
   * These are:
   *
   * THROW(const string message)
   * 
   * and
   * 
   * ASSERT_OR_THROW(const string message, const bool condition)
   *
   * The latter can be used to in place of assert(const bool condition).  
   * Throwing an exception instead of aborting overcomes some of the 
   * limitations of the standard assert.
   */
  class Exception {
#ifdef HAVE_STACK_TRACE
  public:
    typedef std::vector<std::string> trace_t;
#endif

  private:
    std::string message;
    FileLocation location;
    SmartPointer<Exception> cause;
#ifdef HAVE_STACK_TRACE
    SmartPointer<trace_t> trace;
#endif

  public:
    static unsigned int causePrintLevel;
#ifdef HAVE_STACK_TRACE
    static bool enableStackTraces;
#endif

    Exception() {init();}

    Exception(const std::string message) : message(message) {
      init();
    }

    Exception(const std::string message, const FileLocation &location) 
      : message(message), location(location) {
      init();
    }
  
    Exception(const std::string message, Exception &cause) :
      message(message) {
      this->cause = new Exception(cause);
      init();
    }

    Exception(const std::string message, const FileLocation &location, 
              const Exception &cause) :
      message(message), location(location) {
      this->cause = new Exception(cause);
      init();
    }

#ifdef HAVE_STACK_TRACE
    Exception(const Exception &e, const Exception &cause) :
      message(e.message), location(e.location), cause(new Exception(cause)),
      trace(e.trace) {}

    /// Copy constructor
    Exception(const Exception &e) :
      message(e.message), location(e.location), cause(e.cause),
      trace(e.trace) {}
#else
    Exception(const Exception &e, const Exception &cause) :
      message(e.message), location(e.location), cause(new Exception(cause)) {}

    /// Copy constructor
    Exception(const Exception &e) :
      message(e.message), location(e.location), cause(e.cause) {}
#endif

    virtual ~Exception() {}

    const std::string &getMessage() const {return message;}
    const FileLocation &getLocation() const {return location;}

    /**
     * @return A SmartPointer to the Exception that caused this 
     *         exception or NULL.
     */  
    SmartPointer<Exception> getCause() const {return cause;}

#ifdef HAVE_STACK_TRACE
    SmartPointer<trace_t> getTrace() const {return trace;}
#endif

    /** 
     * Prints the complete exception recuring down to the cause exception if
     * not null.  WARNING: If there are many layers of causes this function
     * could print a very large amount of data.  This can be limited by
     * setting the causePrintLevel variable.
     * 
     * @param stream The output stream.
     * @param printLocations Print file locations.
     * @param printLevel The current cause print level.
     * 
     * @return A reference to the passed stream.
     */  
    std::ostream &print(std::ostream &stream,
                        bool printLocations = true,
                        unsigned int printLevel = 0) const {

      stream << message << std::endl;

      if (printLocations && !location.isEmpty())
        stream << "@ " << location << std::endl;

#ifdef HAVE_STACK_TRACE
      if (enableStackTraces && !trace.isNull()) {
        trace_t::iterator it;
        for (it = trace->begin(); it != trace->end(); it++)
          stream << std::endl << "  " << *it;
      }
#endif

      if (!cause.isNull()) {
        stream << std::endl << " ";

        if (printLevel > causePrintLevel) {
          stream << "Aborting exception dump due to causePrintLevel limit! "
                 << "Increase Exception::causePrintLevel to see more.";

        } else {
          stream << "caused by: ";
          cause->print(stream, printLocations, printLevel);
        }
      }

      return stream;
    }

  protected:
    void init() {
#ifdef HAVE_STACK_TRACE
      if (enableStackTraces) {
        trace = new trace_t();

        // When Optimization is turned on functions such as this
        // one are often inlined and not visable to the debugger.
        // This means stack traces for optimized code will often
        // be incomplete. 
#ifdef __OPTIMIZE__
        trace->push_back("Warning: Optimization can cause incomplete traces.");
#endif

#if defined(HAVE_DEBUGGER) && defined(HAVE_BACKTRACE)
        if (!Debugger::getStackTrace(*trace)) {
          trace->clear();
          Backtrace::getStackTrace(*trace);
        }

#else
#if defined(HAVE_DEBUGGER)
        Debugger::getStackTrace(*trace);
#else
        Backtrace::getStackTrace(*trace);
#endif
#endif
      }
#endif // HAVE_STACK_TRACE
    }

    friend std::ostream &operator<<(std::ostream &, const Exception &);
  };

  /** 
   * An stream output operator for Exception.  This allows you to print the
   * text of an exception to a stream like so:
   *
   * . . .
   * } catch (Exception &e) {
   *   cout << e << endl;
   *   return 0;
   * }
   */
  inline std::ostream &operator<<(std::ostream &stream, const Exception &e) {
    e.print(stream);
    return stream;
  }
}

#ifndef THROW
#define THROW(msg) throw ProtoMol::Exception((msg), FILE_LOCATION)
#define THROWC(msg, cause) \
  throw ProtoMol::Exception((msg), FILE_LOCATION, (cause))
#define ASSERT_OR_THROW(msg, condition) {if (!(condition)) THROW(msg);}

#define THROWS(msgs) {                          \
    std::ostringstream errMsg;                  \
    errMsg << msgs;                             \
    THROW(errMsg.str());                        \
  }

#define THROWSC(msgs, cause) {                  \
    std::ostringstream errMsg;                  \
    errMsg << msgs;                             \
    THROWC(errMsg.str(), cause);                \
  }
#endif // THROW

#endif
