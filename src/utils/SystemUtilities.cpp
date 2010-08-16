#include "utils/SystemUtilities.hpp"

#include "utils/StringUtils.hpp"
#include "utils/Exception.hpp"

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h> /* for GetFullPathName */

//____ Define the missing symbols from <unistd.h> for M$ ....
#include <direct.h>
#define CHDIR _chdir
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#define access _access
#include <fcntl.h>
#include <io.h>
#define F_OK 0
#define W_OK 2
#define R_OK 4

#else // _WIN32

#include <libgen.h>
#include <unistd.h>
#define CHDIR chdir
#define PATHSEP '/'
#define PATHSEPSTR "/"
#endif

#include <fstream>
#include <sys/stat.h>

#ifdef BUILD_FOR_FAH
#include <fah/core/chksum/overrides.h>
#endif

using namespace std;

namespace OpenMD {
//____ changeDirectory
  bool changeDirectory(const string &fileName) {
    char *confFile = (char *)fileName.c_str();
    char *currentdir = confFile;
    char *tmp = NULL;

#ifdef WIN32
    // Replace all '/' by '\'
    for (tmp = confFile; *tmp; ++tmp)
      if (*tmp == '/')
        *tmp = '\\';
#endif


    for (tmp = confFile; *tmp; ++tmp) ;

    // find final null
    for (; tmp != confFile && *tmp != PATHSEP; --tmp) ;

    // find last '/'
    if (tmp != confFile) {
      *tmp = 0;
      confFile = tmp + 1;
      if (CHDIR(currentdir))
        return false;
    } else if (*tmp == PATHSEP) // config file in / is odd, but it might happen
      if (CHDIR(PATHSEPSTR))
        return false;

    return true;
  }

//____ isAccessible
  bool isAccessible(const string &fileName) {
    return ::access(fileName.c_str(), F_OK) == 0;
  }
  
  void splitFileName(const string &filename, string &dirname,
                     string &basename, string &extension) {
    string::size_type pos = filename.rfind(PATHSEP);

    if (pos == string::npos) basename = filename;
    else {
      dirname = filename.substr(0, pos);
      basename = filename.substr(pos + 1);
    }

    pos = basename.rfind('.');

    if (pos != string::npos) {
      extension = basename.substr(pos + 1);
      basename = basename.substr(0, pos);
    }
  }

  unsigned int getFileSize(const string &filename) {
    ifstream f(filename.c_str());
    if (!f.is_open()) THROWS("Failed to open '" << filename << "'");

    f.seekg(0, ios::end);
    return (unsigned int)f.tellg();
  }

//____ openMDAbort

  static void (*myAbortFunction)() = NULL;

  void openMDAbort() {
    if (myAbortFunction) (*myAbortFunction)();

    THROW("EXIT");
  }

//____ setOpenMDAbort
  void setOpenMDAbort(void (*abortFunction)()) {
    myAbortFunction = abortFunction;
  }

//____ openMDExit

  static void (*myExitFunction)() = NULL;

  void openMDExit() {
    if (myExitFunction) (*myExitFunction)();

    THROW("EXIT");
  }

//____ setOpenMDExit
  void setOpenMDExit(void (*exitFunction)()) {
    myExitFunction = exitFunction;
  }

//____ openMDStartSerial
  static void (*myStartSerial)(bool) = NULL;

  void openMDStartSerial(bool exludeMaster) {
    if (myStartSerial != NULL)
      (*myStartSerial)(exludeMaster);
  }

//____ setOpenMDStartSerial
  void setOpenMDStartSerial(void (*startSerialFunction)(bool)) {
    myStartSerial = startSerialFunction;
  }

//____ openMDEndSerial
  static void (*myEndSerial)(bool) = NULL;

  void openMDEndSerial(bool exludeMaster) {
    if (myEndSerial != NULL)
      (*myEndSerial)(exludeMaster);
  }

//____ setOpenMDExit
  void setOpenMDEndSerial(void (*endSerialFunction)(bool)) {
    myEndSerial = endSerialFunction;
  }

//____ ISLITTLEENDIAN
  struct Endian {
    // Helper class to make sure that we get endianess correct ... M$
    static bool isLittleEndian() {
      unsigned int tmp = 1;
      return 0 != *(reinterpret_cast<const char *> ( &tmp));
    }
  };
  const bool ISLITTLEENDIAN = Endian::isLittleEndian();

  string getCanonicalPath(const string &path) {
    char buf[4096];

#ifdef _WIN32
    char *finalpart;
    DWORD len = GetFullPathName(path.c_str(), 4096, buf, &finalpart);
    if (len == 0 || len > 4095)
      THROW(string("GetFullPathName '") + path + "' failed.");

    return buf;

#else
    char tmp[path.length() + 3];

    // The file might not exist yet but its directory must.
    strcpy(tmp, path.c_str());
    char *dir = dirname(tmp);

    if (!realpath(dir, buf))
      THROW(string("realpath '") + path + "' failed.");

    strcpy(tmp, path.c_str());

    return string(buf) + "/" + basename(tmp);
#endif
  }

  bool SystemUtilities::unlink(const string &path) {
#ifdef BUILD_FOR_FAH
    return fah_unlink(path.c_str()) == 0;
#else
    return ::unlink(path.c_str()) == 0;
#endif
  }

  void SystemUtilities::rename(const string &src, const string &dst) {
    unlink(dst);
#ifdef BUILD_FOR_FAH
    fah_rename(src.c_str(), dst.c_str());
#else
    ::rename(src.c_str(), dst.c_str());
#endif
  }
}









