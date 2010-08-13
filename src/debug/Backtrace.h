#ifndef BACKTRACE_H
#define BACKTRACE_H

#include <vector>
#include <string>

namespace ProtoMol {
  class Backtrace {
  public:
    typedef std::vector<std::string> trace_t;
  
    static void getStackTrace(trace_t &trace);
  };
}
#endif // BACKTRACE_H
