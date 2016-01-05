#ifndef UTILS_REVISION_HPP
#define UTILS_REVISION_HPP
#include <string>

extern "C" {
  extern const char version[];
  extern const char revision[];
  extern const char build_date[];
}

namespace OpenMD {

  class Revision {
  public:
    std::string getVersion();
    std::string getRevision();
    std::string getFullRevision();
    std::string getBuildDate();
  };
}

#endif
