/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

/**
 * @file ifstrstream.cpp
 * @author Teng Lin
 * @date 10/14/2004
 * @version 1.0
 */

#include "io/ifstrstream.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {

/**  Constructs an object of class ifstream.  */
#ifdef IS_MPI
  ifstrstream::ifstrstream() :
      std::basic_istream<char, std::char_traits<char>>(0), internalStringBuf_(),
      isRead(false) {
    this->init(&internalStringBuf_);
  }
#else
  ifstrstream::ifstrstream() :
      std::basic_istream<char, std::char_traits<char>>(0), internalFileBuf_(),
      isRead(false) {
    this->init(&internalFileBuf_);
  }
#endif

/**
 * Explicit constructor
 * \param filename String containing the name of the file to be opened
 * \param mode Flags describing the requested i/o mode for the file,
 * default value is ios_base::in
 * \param checkFilename Flags indicating checking the file name in parallel
 */
#ifdef IS_MPI
  ifstrstream::ifstrstream(const char* filename, std::ios_base::openmode mode,
                           bool checkFilename) :
      std::basic_istream<char, std::char_traits<char>>(0),
      internalStringBuf_(), isRead(false) {
    this->init(&internalStringBuf_);
#if defined(_MSC_VER)
    isRead =
        internalOpen(filename, mode | std::ios_base::in | std::ios_base::binary,
                     checkFilename);
#else
    isRead = internalOpen(filename, mode | std::ios_base::in, checkFilename);
#endif
  }
#else
  ifstrstream::ifstrstream(const char* filename, std::ios_base::openmode mode,
                           bool checkFilename) :
      std::basic_istream<char, std::char_traits<char>>(0),
      internalFileBuf_(), isRead(false) {
    this->init(&internalFileBuf_);
#if defined(_MSC_VER)
    isRead =
        internalOpen(filename, mode | std::ios_base::in | std::ios_base::binary,
                     checkFilename);
#else
    isRead = internalOpen(filename, mode | std::ios_base::in, checkFilename);
#endif
  }
#endif
  /**
   * virtual destructor will close the file (in single mode) and clear
   * the stream buffer
   */
  ifstrstream::~ifstrstream() { close(); }

  /**
   * Opens a file and associates a buffer with the specified file to
   * perform the i/o operations (single mode). The primary node reads a
   * file and broadcasts its content to the secondary nodes. After
   * broadcasting, all nodes fall back to stringstream (parallel
   * mode).
   * \param filename String containing the name of the file to be opened
   * \param mode Flags describing the requested i/o mode for the file
   * \param checkFilename Flags indicating checking the file name in parallel
   */
  void ifstrstream::open(const char* filename, std::ios_base::openmode mode,
                         bool checkFilename) {
    if (!isRead) {
#if defined(_MSC_VER)
      isRead =
          internalOpen(filename, mode | std::ios_base::binary, checkFilename);
#else
      isRead = internalOpen(filename, mode, checkFilename);
#endif
    }
  }

  /**
   * Tests if the stream is currently associated with a valid buffer.
   * \return true if a file has successfully been opened (single mode)
   * or the whole file has been read and spread among all of the
   * processors (parallel mode), otherwise false is returned
   */
  bool ifstrstream::is_open() {
#ifdef IS_MPI
    return isRead;
#else
    // single version fall back to ifstream
    return internalFileBuf_.is_open();
#endif
  }

  /**
   * In single mode, closes a file. The stream's file buffer is
   * released from its association with the currently open file. In
   * parallel mode, cleans up.
   */
  void ifstrstream::close() {
#ifndef IS_MPI
    // single version fall back to ifstream
    if (!internalFileBuf_.close()) this->setstate(std::ios_base::failbit);
#endif

    isRead = false;
  }

  /**
   * Gets the stream buffer object associated with the stream
   * \return A pointer to the stream buffer object(filebuf in single
   * mode, stringbuf in parallel mode) associated with the stream.
   */
  std::basic_streambuf<char, std::char_traits<char>>* ifstrstream::rdbuf() {
#ifdef IS_MPI
    return static_cast<_StringBuf*>(&internalStringBuf_);
#else
    return static_cast<_FileBuf*>(&internalFileBuf_);
#endif
  }

  /**
   * Internal function used to open the file
   * \return true if succesfully opens a file (single mode) or gets
   * the file content (parallel mode) otherwise return false
   * \param filename String containing the name of the file to be opened
   * \param mode Flags describing the requested i/o mode for the file
   * \param checkFilename Flags indicating checking the file name in parallel
   * \todo use try - catch syntax to make the program more readable
   */
#ifdef IS_MPI
  bool ifstrstream::internalOpen(const char* filename,
                                 std::ios_base::openmode mode,
                                 bool checkFilename) {
    // int commStatus;
    long fileSize;
    char* fbuf;
    int filenameLen;
    int diffFilename;
    int error;
    int myRank;
    int primaryNode;

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    primaryNode = 0;

    if (myRank == primaryNode) {
      if (checkFilename) {
        // check the filename is the same
        filenameLen = strlen(filename);
        MPI_Bcast(&filenameLen, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
        MPI_Bcast((void*)filename, filenameLen, MPI_CHAR, primaryNode,
                  MPI_COMM_WORLD);

        diffFilename = 0;
        MPI_Allreduce(&diffFilename, &error, 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);

        // if file names are different just return false
        if (error > 0) return false;
      }

      std::ifstream fin(filename, mode);

      if (fin.is_open()) {
        fin.seekg(0, std::ios::end);
        fileSize = fin.tellg();
        fin.seekg(0, std::ios::beg);

        // '\0' need one more char
        fbuf = new char[fileSize + 1];

        assert(fbuf != 0);

        fin.read(fbuf, fileSize);

        if (fin.fail()) fileSize = FileIOError;

        // broadcast the file size
        MPI_Bcast(&fileSize, 1, MPI_LONG, primaryNode, MPI_COMM_WORLD);

        if (fileSize < 0) {
          fin.close();
          delete[] fbuf;

          return false;
        }

        // make a c-style  std::string and broadcast it
        fbuf[fileSize] = '\0';
        MPI_Bcast(fbuf, fileSize + 1, MPI_CHAR, primaryNode, MPI_COMM_WORLD);

        // close the file and delete the buffer
        fin.close();
        internalStringBuf_.str(fbuf);
        delete[] fbuf;
      } else {
        fileSize = FileNotExists;
        MPI_Bcast(&fileSize, 1, MPI_LONG, primaryNode, MPI_COMM_WORLD);
        return false;
      }

    } else {  // secondary nodes

      // check file name
      if (checkFilename) {
        MPI_Bcast(&filenameLen, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);

        char* primaryFilename = new char[filenameLen];
        MPI_Bcast(primaryFilename, filenameLen, MPI_CHAR, primaryNode,
                  MPI_COMM_WORLD);

        if (strcmp(primaryFilename, filename) == 0)
          diffFilename = 0;
        else
          diffFilename = 1;

        delete[] primaryFilename;

        MPI_Allreduce(&diffFilename, &error, 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);

        if (error > 0) return false;
      }
      // get file size
      MPI_Bcast(&fileSize, 1, MPI_LONG, primaryNode, MPI_COMM_WORLD);

      if (fileSize >= 0) {
        fbuf = new char[fileSize + 1];
        assert(fbuf);

        // receive file content
        MPI_Bcast(fbuf, fileSize + 1, MPI_CHAR, primaryNode, MPI_COMM_WORLD);

        internalStringBuf_.str(fbuf);
        delete[] fbuf;

      } else if (fileSize == FileNotExists) {
        return false;

      } else if (fileSize == FileIOError) {
        return false;
      }
    }

    this->clear();
    return true;
  }
#else
  bool ifstrstream::internalOpen(const char* filename,
                                 std::ios_base::openmode mode, bool) {
    // in single version, fall back to ifstream
    if (!internalFileBuf_.open(filename, mode)) {
      this->setstate(std::ios_base::failbit);
      return false;
    }

    this->clear();
    return true;
  }
#endif
}  // namespace OpenMD
