/*  -*- c++ -*-  */
#ifndef UTILS_SYSTEMUTILITIES_HPP
#define UTILS_SYSTEMUTILITIES_HPP

#include "config.h"
#include <string>
#include <algorithm>

#include "utils/Report.hpp"


namespace OpenMD {
  /// Changes to the actual directory of the file name
  bool changeDirectory(const std::string &fileName);

  /// Test if the file is accessible
  bool isAccessible(const std::string &fileName);

  void splitFileName(const std::string &filename, std::string &dirname,
                     std::string &basename, std::string &extension);

  unsigned int getFileSize(const std::string &filename);

  /// Does an abort, calling the adequate abort system function
  void openMDAbort();

  /// Sets function to be called when calling OpenMDAbort()
  void setopenMDAbort(void (*abortFunction)());

  /// Does an exit, calling the adequate exit system function
  void openMDExit();

  /// Sets function to be called when calling OpenMDxit()
  void setOpenMDExit(void (*exitFunction)());

  /// Initiates the serialization block, e.g., used by Report
  void openMDStartSerial(bool exludeMaster);

  /// Sets function to be called to start serialization
  void setOpenMDStartSerial(void (*startSerialFunction)(bool));

  /// Finializes the serialization block
  void openMDEndSerial(bool exludeMaster);

  /// Sets function to be called to end serialization
  void setopenMDEndSerial(void (*endSerialFunction)(bool));

  /// Return an resolved absolute path
  std::string getCanonicalPath(const std::string &path);

  /// Swap function to change endianess
  template<typename T> inline
  void swapBytes(T &t) {
    if (sizeof(T) % 2 != 0)
      Report::report << Report::error << "Cannot swap types of uneven size." <<
      Report::endr;
    char *res = reinterpret_cast<char *> (&t);
    std::reverse(res, res + sizeof(T));
  }

  /// Shift left of four Real's
  inline void shift(RealType &a, RealType &b, RealType &c, const RealType d) {
    a = b; b = c; c = d;
  }

  /// bool constant if the machine is littleEndian or not
  extern const bool ISLITTLEENDIAN;

  /// Clears a container explicitly
  template<typename T> inline
  void realclear(T &t) {T tmp; t.swap(tmp);}

  /// Shrinks the capacity of a container explicitly
  template<typename T> inline
  void shrink(T &t) {T tmp(t); t.swap(tmp);}

  namespace SystemUtilities {
    bool unlink(const std::string &path);
    void rename(const std::string &src, const std::string &dst);
  };
}
#endif /* SYSTEMUTILITIES_H */
