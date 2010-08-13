#include <protomol/debug/Backtrace.h>
#include <protomol/config.h>

#ifdef HAVE_BACKTRACE
#include <execinfo.h>

using namespace std;
using namespace ProtoMol;

#define MAX_STACK 256

void Backtrace::getStackTrace(trace_t &trace) {
  void *stack[MAX_STACK];
  char **strings;

  int n = backtrace(stack, MAX_STACK);
  strings = backtrace_symbols(stack, n);
  
  for (int i = 0; i < n; i++)
    trace.push_back(strings[i]);
}
#endif // HAVE_BACKTRACE
