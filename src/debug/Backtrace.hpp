#ifndef DEBUG_BACKTRACE_H
#define DEBUG_BACKTRACE_H

#include <vector>
#include <string>

namespace OpenMD {
  class Backtrace {
  public:
    typedef std::vector<std::string> trace_t;
  
    static void getStackTrace(trace_t &trace);
  };
}
#endif // BACKTRACE_H
